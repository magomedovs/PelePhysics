#!/usr/bin/env python3
"""Link the particle relaxation test using a completed SprayTest GNU make build.

Usage: python3 build_relaxation.py BUILD_LOG OUTPUT_EXECUTABLE
The log must include a compile command (make -B -n can provide one).
"""
from pathlib import Path
import os
import shlex
import subprocess
import sys

lmex = Path(os.environ.get(
    'PELELMEX_HOME', Path(__file__).resolve().parents[3] / 'PeleLMeX'
)).expanduser().resolve()
case = lmex / 'Exec/RegTests/SprayTest'
lines = Path(sys.argv[1]).read_text().splitlines()
compile_line = next(x for x in lines if x.startswith('mpicxx ') and ' -c ' in x)
link_line = next(x for x in reversed(lines) if x.startswith('mpicxx ') and ' -o ' in x and ' -c ' not in x)
output = Path(sys.argv[2]).resolve()
source = Path(__file__).with_name('relaxation.cpp')
args = shlex.split(compile_line)
i = args.index('-c'); args[i+1] = str(source)
i = args.index('-o'); args[i+1] = str(output.with_suffix('.o'))
subprocess.run(args, cwd=case, check=True)
args = shlex.split(link_line)
args = [str(output.with_suffix('.o')) if x.endswith('/main.o') else x for x in args]
i = args.index('-o'); args[i+1] = str(output)
subprocess.run(args, cwd=case, check=True)
