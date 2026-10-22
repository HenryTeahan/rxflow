"""Concurrent workers must upgrade a historical queue exactly once."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from pathlib import Path
import sqlite3
import tempfile
from threading import Barrier
import unittest
from unittest.mock import patch

from hpc_workflow import connect_database


class DatabaseMigrationTests(unittest.TestCase):
    def test_workers_reread_schema_after_acquiring_migration_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "workflow.sqlite"
            db = connect_database(path)
            db.execute("ALTER TABLE jobs DROP COLUMN ts_parent_id")
            db.execute("ALTER TABLE jobs DROP COLUMN ts_sampling_json")
            db.execute("""INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                input_xyz,smiles,charge,settings_json,status,gibbs_free_energy_eh)
                VALUES ('1',1,1,'TS',1,'existing.xyz','graph',0,'{}','succeeded',-10.)""")
            original = dict(db.execute("SELECT * FROM jobs").fetchone())
            db.close()

            workers = 16
            barrier = Barrier(workers)
            real_connect = sqlite3.connect

            class SimultaneousSchemaReaders(sqlite3.Connection):
                first_schema_read = True

                def execute(self, sql, *args, **kwargs):
                    cursor = super().execute(sql, *args, **kwargs)
                    if sql == "PRAGMA table_info(jobs)" and self.first_schema_read:
                        self.first_schema_read = False
                        # Every worker has seen the same missing-column schema.
                        barrier.wait(timeout=15)
                    return cursor

            def open_database(_):
                connection = connect_database(path)
                try:
                    return dict(connection.execute("SELECT * FROM jobs").fetchone())
                finally:
                    connection.close()

            def simultaneous_connect(*args, **kwargs):
                return real_connect(*args, factory=SimultaneousSchemaReaders, **kwargs)

            with patch("hpc_workflow.sqlite3.connect", side_effect=simultaneous_connect):
                with ThreadPoolExecutor(max_workers=workers) as pool:
                    rows = list(pool.map(open_database, range(workers)))
            for row in rows:
                self.assertEqual({key: row[key] for key in original}, original)
                self.assertIsNone(row["ts_parent_id"])
                self.assertIsNone(row["ts_sampling_json"])

    def test_failed_migration_rolls_back_all_added_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "workflow.sqlite"
            db = connect_database(path)
            db.execute("ALTER TABLE jobs DROP COLUMN ts_parent_id")
            db.execute("ALTER TABLE jobs DROP COLUMN ts_sampling_json")
            db.close()
            real_connect = sqlite3.connect

            class FailedMigration(sqlite3.Connection):
                def execute(self, sql, *args, **kwargs):
                    if sql == "ALTER TABLE jobs ADD COLUMN ts_sampling_json TEXT":
                        raise sqlite3.OperationalError("simulated failure")
                    return super().execute(sql, *args, **kwargs)

            with patch("hpc_workflow.sqlite3.connect", side_effect=lambda *args, **kwargs:
                       real_connect(*args, factory=FailedMigration, **kwargs)):
                with self.assertRaisesRegex(sqlite3.OperationalError, "simulated failure"):
                    connect_database(path)
            with closing(real_connect(path)) as db:
                columns = {row[1] for row in db.execute("PRAGMA table_info(jobs)")}
                self.assertNotIn("ts_parent_id", columns)
                self.assertNotIn("ts_sampling_json", columns)
            repaired = connect_database(path)
            repaired.close()


if __name__ == "__main__":
    unittest.main()
