#!/usr/bin/env python3
"""A38 symbolic benchmark (SA69 T158).

Task: print A38_Symbolic_Benchmark_Spec.txt exactly. That text is a complete symbolic statement of the Standard Model
with general relativity: gauge group, matter fields and charges, the action, the quantum rule, and 28 declared inputs at
stated precision.

Declared language: Python 3 with its standard library (lzma, base64). This script compresses the target with raw LZMA2,
writes the counted program A38_Symbolic_Benchmark_Program.py, runs it, checks its output against the target by SHA-256,
and writes A38_Symbolic_Benchmark_Results.json.

A working program is an upper bound for its declared task in its declared language. It is not the shortest description,
and it computes no predictions: a program that derives observable predictions from these laws is the next, predictive
task. Inputs are rounded values near current central values, declared for the count; they are not a new fit.

usage: python3 A38_Symbolic_Benchmark.py
"""
import base64, hashlib, json, lzma, os, platform, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
spec_path = os.path.join(here, 'A38_Symbolic_Benchmark_Spec.txt')
prog_path = os.path.join(here, 'A38_Symbolic_Benchmark_Program.py')
res_path = os.path.join(here, 'A38_Symbolic_Benchmark_Results.json')

target = open(spec_path, 'rb').read()
payload = lzma.compress(target, format=lzma.FORMAT_RAW, filters=[{'id': lzma.FILTER_LZMA2, 'preset': 9 | lzma.PRESET_EXTREME}])
program = ("import lzma,base64;print(lzma.decompress(base64.b85decode('" + base64.b85encode(payload).decode()
           + "'),lzma.FORMAT_RAW,filters=[{'id':lzma.FILTER_LZMA2}]).decode(),end='')").encode()
open(prog_path, 'wb').write(program)

output = subprocess.run([sys.executable, '-I', prog_path], capture_output=True, check=True).stdout
sha = lambda b: hashlib.sha256(b).hexdigest()
result = {
    'task': 'Print A38_Symbolic_Benchmark_Spec.txt exactly: the Standard Model with general relativity, stated completely.',
    'declared_language': 'Python 3 with its standard library (lzma, base64)',
    'environment': {'python': platform.python_version(), 'implementation': platform.python_implementation()},
    'target': {'bytes': len(target), 'bits': 8 * len(target), 'sha256': sha(target)},
    'compressed_statement': {'format': 'raw LZMA2, preset 9 extreme', 'bytes': len(payload), 'bits': 8 * len(payload)},
    'program': {'file': 'A38_Symbolic_Benchmark_Program.py', 'bytes': len(program), 'bits': 8 * len(program), 'sha256': sha(program)},
    'output_matches_target': sha(output) == sha(target),
    'tolerance': 'exact match',
    'scope': ('An upper bound for the symbolic task in the declared language. The program states the laws and inputs; it '
              'computes no predictions. A shorter program may exist; the predictive task is the next benchmark.'),
}
open(res_path, 'w', encoding='utf-8').write(json.dumps(result, indent=2) + '\n')
print(json.dumps({'program_bits': result['program']['bits'], 'compressed_statement_bits': result['compressed_statement']['bits'],
                  'target_bits': result['target']['bits'], 'output_matches_target': result['output_matches_target']}))
if not result['output_matches_target']:
    sys.exit('program output differs from the target')
