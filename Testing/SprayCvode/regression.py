#!/usr/bin/env python3
"""Run unmodified spray regression inputs in isolated output directories.

No third-party Python packages are required. Preserve baseline executables before
building the change; use --compare to compare every plotfile FAB and ASCII particle
history and checkpoint state, excluding build metadata.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

LMEX = Path(os.environ.get(
    'PELELMEX_HOME', Path(__file__).resolve().parents[3] / 'PeleLMeX'
)).expanduser().resolve()
SINGLE = LMEX / 'Exec/RegTests/SingleDropEvap'
SPRAY = LMEX / 'Exec/RegTests/SprayTest'
GROUPS = {
    'gcm2': ['WongLin', 'Nomura', 'Daif', 'RungeHep', 'RungeDec', 'RungeMix'],
    'gcm67': ['RungeJP8', 'Burger1Bar', 'Burger10Bar', 'Burger50Bar'],
    'gcm1': ['RungeJP8-H', 'Burger1Bar-H', 'Burger10Bar-H', 'Burger50Bar-H'],
    'detailed67': ['RungeJP8-D'],
    'mp2': ['WongLin', 'Nomura', 'Daif', 'RungeHep', 'RungeDec', 'RungeMix'],
    'manifold2': ['Daif', 'WongLin'],
    'input-gcm': [],
}


def run(exe, work, args, np):
    cmd = ([os.environ.get('MPIEXEC', 'mpiexec'), '-n', str(np)] if np > 1 else [])
    cmd += [str(exe), *map(str, args)]
    (work / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
    with (work / 'run.log').open('w') as log:
        code = subprocess.call(cmd, cwd=work, stdout=log, stderr=subprocess.STDOUT)
    if code:
        raise RuntimeError(f'{work}: exit {code}; see run.log')
    print(f'PASS {work.name}', flush=True)


def payloads(path):
    files = sorted(p for p in path.rglob('*') if p.is_file() and
                   (p.suffix == '.p3d' or
                    (any(part.startswith('chk') for part in
                         p.relative_to(path).parts) and
                     p.name != 'PeleLMeX_job_info') or
                    (p.name == 'Header' and p.parent.name.startswith('plt')) or
                    ((p.name.startswith('Cell_D_') or p.name == 'Cell_H') and
                     p.parent.name.startswith('Level_'))))
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files}


def compare(before, after):
    pending = [str(p) for root in (before, after) for p in root.rglob('run.log')
               if not any(line.startswith('AMReX (') and 'finalized' in line
                          for line in p.read_text().splitlines())]
    if pending:
        raise SystemExit('Runs have not completed: ' + ', '.join(pending))
    a, b = payloads(before), payloads(after)
    missing = sorted(a.keys() ^ b.keys())
    different = sorted(k for k in a.keys() & b.keys() if a[k] != b[k])
    result = {'baseline_payloads': len(a), 'changed_payloads': len(b),
              'missing': missing, 'different': different,
              'bitwise_equal': bool(a) and not missing and not different}
    (after / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    if not result['bitwise_equal']:
        raise SystemExit(1)


def check_grid(output, executable):
    variables = {'density', 'rho.Y(NC10H22)', 'x_velocity', 'temp', 'rhoh',
                 'spray_vol'}
    results = []
    for first, second in [('32_continuous', '64_continuous'),
                          ('32_continuous', '32_restart'),
                          ('64_continuous', '64_restart')]:
        command = [str(executable.resolve()), '-n', '2', '-a',
                   str(output / first / 'plt00009'),
                   str(output / second / 'plt00009')]
        process = subprocess.run(command, capture_output=True, text=True)
        text = process.stdout + process.stderr
        (output / f'fcompare_{first}_{second}.txt').write_text(text)
        errors = {name: [] for name in variables}
        for line in process.stdout.splitlines():
            words = line.split()
            if len(words) >= 3 and words[0] in variables:
                errors[words[0]].append(float(words[2]))
        # fcompare exits nonzero for harmless roundoff too. Apply the original
        # SprayTest's six-variable, 1e-12 relative-error acceptance criterion.
        passed = all(values and all(math.isfinite(x) and x <= 1.e-12
                                    for x in values)
                     for values in errors.values())
        results.append({'first': first, 'second': second, 'passed': passed,
                        'relative_errors': errors})
    (output / 'grid_checks.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))
    if not all(result['passed'] for result in results):
        raise SystemExit(1)


def check_backend(output, executable):
    common = [str(executable.resolve()), str(SPRAY / 'input.2d'),
              'particles.integrator=cvode_bdf',
              'particles.cvode.atol_velocity=1.e-12',
              'particles.cvode.atol_temperature=1.e-8',
              'particles.cvode.rtol=1.e-10', 'particles.cvode.max_order=2']
    cases = [('balances', [], None),
             ('invalid_method', ['particles.integrator=unknown'],
              'particles.integrator must be'),
             ('invalid_tolerance', ['particles.cvode.rtol=-1'],
              'Invalid particles.cvode'),
             ('step_budget', ['particles.cvode.max_steps=1'],
              'failed with status -1')]
    results = []
    for name, extra, failure in cases:
        command = common + extra
        process = subprocess.run(command, cwd=output, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output / f'{name}.log').write_text(process.stdout)
        passed = ((process.returncode == 0) if failure is None else
                  (process.returncode != 0 and failure in process.stdout))
        results.append({'name': name, 'command': command,
                        'returncode': process.returncode, 'passed': passed})
    (output / 'backend_checks.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))
    if not all(result['passed'] for result in results):
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--group', choices=['spray', *GROUPS])
    parser.add_argument('--np', type=int, default=2)
    parser.add_argument('--max-steps', type=int)
    parser.add_argument('--integrator', choices=['explicit', 'cvode_bdf'])
    parser.add_argument('--runtime', action='append', default=[],
                        help='Additional runtime override; may be repeated')
    parser.add_argument('--cmlm', type=Path,
                        help='CMLM source directory for manifold table generation')
    parser.add_argument('--compare', type=Path)
    parser.add_argument('--check-grid', type=Path,
                        help='Run the original grid/restart checks with fcompare')
    parser.add_argument('--check-backend', type=Path,
                        help='Run analytical balances and expected failure checks')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.compare:
        return compare(args.compare.resolve(), out)
    if args.check_grid:
        return check_grid(out, args.check_grid)
    if args.check_backend:
        return check_backend(out, args.check_backend)
    exe = args.exe.resolve()
    flags = list(args.runtime)
    if args.integrator:
        flags += [f'particles.integrator={args.integrator}']
    if args.max_steps:
        flags += [f'amr.max_step={args.max_steps}']
    if args.group == 'spray':
        for size in [32, 64]:
            common = [SPRAY / 'input.2d', f'amr.max_grid_size={size}',
                      f'amr.blocking_factor={size}',
                      f'amr.initial_grid_file={SPRAY}/two_d_gridfiles/gridfile_{size}_1',
                      f'amr.regrid_file={SPRAY}/two_d_gridfiles/gridfile_{size}_2',
                      'amr.max_step=9', 'amr.plot_int=100',
                      'particles.write_ascii_files=1', *flags]
            first = out / f'{size}_continuous'; first.mkdir()
            run(exe, first, [*common, 'amr.check_int=5'], args.np)
            restart = out / f'{size}_restart'; restart.mkdir()
            run(exe, restart, [*common, f'amr.restart={first}/chk00005',
                               'amr.check_int=-1'], args.np)
    elif args.group == 'input-gcm':
        work = out / 'EmbeddedGCM'
        work.mkdir()
        run(exe, work, [SINGLE / 'single-drop-evap-gcm-test.inp', *flags], args.np)
    else:
        sys.path.insert(0, str(SINGLE))
        from CaseInfo import SpecifyCase, CreateInputFile, CreateManifoldFiles
        for name in GROUPS[args.group]:
            props = 'mp' if args.group == 'mp2' else 'gcm'
            case = SpecifyCase(name, props)
            case.gen_input_file = 'single-drop-evap.inp'
            if 'Burger' in name:
                fuel = 'mixture_posf10325' if name.endswith('-H') else 'posf10325'
            elif 'JP8' in name:
                fuel = 'mixture_jp8' if name.endswith('-H') else 'jp8'
            else:
                fuel = 'heptane-decane'
            case.spray_input_file = f'sprayProps{props.upper()}_{fuel}.inp'
            if name.endswith('-D'):
                case.spray_input_file = ('../../../Submodules/PelePhysics/Mechanisms/'
                    'fuellib_posf_nonreacting/spray_input_files/'
                    f'sprayProps{props.upper()}_posf10264.inp')
            ref = SINGLE / 'ref_files' / case.ref_name / 'refdvals.csv'
            with ref.open() as f:
                times = [float(row[0]) for row in list(csv.reader(f))[1:]
                         if row and row[0]]
            case.set_end_time(max(times) / case.xconv)
            work = out / name; work.mkdir()
            case.case_path = case.case_dir = str(work)
            case.input_file = str(work / 'input.inp')
            case.input_spray = str(work / f'input_{case.name}_spray.inp')
            CreateInputFile(case)
            extra = list(flags)
            if args.group == 'manifold2':
                if not args.cmlm:
                    parser.error('--cmlm is required for manifold table generation')
                # Deterministic component order for independently generated tables.
                case.droplet.unique_dep_fuel_names = list(dict.fromkeys(
                    case.droplet.dep_fuel_names))
                previous = Path.cwd()
                try:
                    os.chdir(SINGLE)
                    CreateManifoldFiles(case, str(args.cmlm.resolve()))
                finally:
                    os.chdir(previous)
            run(exe, work, [case.input_file, *extra], args.np)
    manifest = payloads(out)
    (out / 'payloads.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
