"""Refresh notebook outputs once its selected DFT single points are ready."""
import os
import argparse
from pathlib import Path
import sys
import time

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from hpc_precursors import load_diagram_inputs


def main():
    import nbformat
    from nbclient import NotebookClient

    parser = argparse.ArgumentParser()
    parser.add_argument('--wait-for-dft', action='store_true')
    args = parser.parse_args()
    cache = Path('/data/tmp/rxflow-notebook-visuals')
    for variable, folder in [('MPLCONFIGDIR', 'matplotlib'), ('XDG_CACHE_HOME', 'cache'),
                             ('IPYTHONDIR', 'ipython'), ('JUPYTER_RUNTIME_DIR', 'jupyter')]:
        directory = cache / folder
        directory.mkdir(parents=True, exist_ok=True)
        os.environ[variable] = str(directory)
    path = BASE / 'reaction_diagrams.ipynb'
    notebook = nbformat.read(path, as_version=4)
    setup = next(cell.source for cell in notebook.cells if cell.id == 'cfa56857')
    database_line = next(line for line in setup.splitlines() if line.startswith('DB_PATH = '))
    database = BASE / database_line.split('Path("', 1)[1].split('"', 1)[0]
    for _ in range(60):
        data = load_diagram_inputs(database, database.with_name('precursors.sqlite'),
                                   BASE / 'solution_references/references.sqlite', energy_mode='auto')
        progress = data['singlepoint_summary']
        print(f"DFT single points: {progress['ready']}/{progress['selected_states']}", flush=True)
        if progress['active_mode'] == 'dft' or not args.wait_for_dft:
            break
        if progress['counts'].get('failed'):
            raise RuntimeError(f"Selected single points failed: {progress['counts']}")
        time.sleep(30)
    else:
        raise TimeoutError('Selected DFT corrections are still incomplete after 30 minutes')

    def starting(cell, cell_index, **kwargs):
        if cell.cell_type == 'code':
            print(f"Executing cell {cell_index}: {cell.id}", flush=True)

    client = NotebookClient(notebook, timeout=900, kernel_name='python3',
                            resources={'metadata': {'path': str(BASE)}}, on_cell_start=starting)
    manager = client.create_kernel_manager()
    manager.kernel_spec.argv = [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}']
    client.execute()
    temporary = path.with_name('reaction_diagrams.refreshed.ipynb')
    nbformat.write(notebook, temporary)
    temporary.replace(path)
    print('Saved notebook with refreshed profiles and active energy-model labels.', flush=True)


if __name__ == '__main__':
    main()
