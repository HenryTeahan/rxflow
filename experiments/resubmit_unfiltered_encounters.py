"""Submit the four existing runs after removing encounter geometry rejection."""
import json
import os
from pathlib import Path
import sys
import tomllib

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from hpc_workflow import submit_workflow


def main():
    executables = dict(
        ORCA_EXE='/lustre/hpc/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca',
        ORCA_EXTERNAL_TOOLS='/lustre/hpc/kemi/thf862/bin/orca-external-tools/bin/oet_gxtb_alpb',
        XTB_EXE='/lustre/hpc/kemi/hteahan/opt/xtb-6.7.1/xtb-6.7.1/bin/xtb')
    for name, value in executables.items():
        if not os.access(value, os.X_OK):
            raise ValueError(f'{name} is not executable: {value}')
    os.environ.update(executables)
    directory = str(Path(executables['ORCA_EXE']).parent)
    os.environ['LD_LIBRARY_PATH'] = directory + (os.pathsep + os.environ['LD_LIBRARY_PATH']
                                               if os.environ.get('LD_LIBRARY_PATH') else '')
    submitted = {}
    report = BASE / 'experiments/unfiltered_encounter_submission.json'
    for name in ('hpc.toml', 'hpc_ED_dF.toml', 'hpc_ED_pyr.toml', 'hpc_ni.toml'):
        path = BASE / name
        config = tomllib.loads(path.read_text())
        submitted[name] = submit_workflow(path, config, resume=True, workers=100, partition='kemi1')
        report.write_text(json.dumps(submitted, indent=2) + '\n')
        print(name, submitted[name], flush=True)


if __name__ == '__main__':
    main()
