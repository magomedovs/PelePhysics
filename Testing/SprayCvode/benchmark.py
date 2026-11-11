#!/usr/bin/env python3
"""Compare particle-integrator settings serially on identical SprayTest steps.

Use --steps 1 for an identical physical interval in every variant. Longer runs
can select slightly different flow timesteps as the particle state changes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import statistics

from regression import run, SPRAY


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True, type=Path)
    parser.add_argument('--after', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--steps', type=int, default=1)
    parser.add_argument('--box-size', choices=[32, 64], type=int, default=32)
    parser.add_argument('--np', type=int, default=2)
    args = parser.parse_args()
    if min(args.repeats, args.steps, args.np) < 1:
        parser.error('repeats, steps, and np must be positive')
    before, after = args.before.resolve(), args.after.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    variants = [('explicit', after, ['particles.integrator=explicit']),
                ('old_order1', before, ['particles.cvode.max_order=1']),
                ('old_order2', before, ['particles.cvode.max_order=2']),
                ('scaled_order1', after, ['particles.cvode.max_order=1']),
                ('scaled_order2', after, ['particles.cvode.max_order=2'])]
    size = args.box_size
    common = [SPRAY / 'input.2d', f'amr.max_grid_size={size}',
              f'amr.blocking_factor={size}',
              f'amr.initial_grid_file={SPRAY}/two_d_gridfiles/gridfile_{size}_1',
              f'amr.regrid_file={SPRAY}/two_d_gridfiles/gridfile_{size}_2',
              f'amr.max_step={args.steps}', 'amr.check_int=-1',
              'particles.write_ascii_files=1']
    records = []
    for repeat in range(args.repeats):
        # Rotate the order to reduce systematic warmup and background-load bias.
        order = variants[repeat % len(variants):] + variants[:repeat % len(variants)]
        for name, exe, flags in order:
            work = output / f'{repeat + 1}_{name}'
            work.mkdir()
            method = [] if name == 'explicit' else ['particles.integrator=cvode_bdf']
            run(exe, work, [*common, *method, *flags], args.np)
            log = (work / 'run.log').read_text()
            advances = [float(x) for x in re.findall(
                r'PeleLMeX::Advance\(\) --> Time:\s*([^\n]+)', log)]
            if len(advances) != args.steps or 'finalized' not in log:
                raise RuntimeError(f'{work}: incomplete run')
            records.append({'variant': name, 'repeat': repeat + 1,
                            'advance_seconds': sum(advances),
                            'total_seconds': float(re.findall(
                                r'Total Time:\s*([^\n]+)', log)[-1]),
                            'physical_end_time': float(re.findall(
                                r'Final simulation time:\s*([^\n]+)', log)[-1])})
    medians = {name: {key: statistics.median(row[key] for row in records
                                           if row['variant'] == name)
                      for key in ('advance_seconds', 'total_seconds')}
               for name, _, _ in variants}
    result = {'platform': platform.platform(), 'mpi_ranks': args.np,
              'box_size': size, 'steps': args.steps, 'repeats': args.repeats,
              'executable_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in (before, after)},
              'runs': records, 'medians': medians}
    (output / 'timings.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(medians, indent=2))


if __name__ == '__main__':
    main()
