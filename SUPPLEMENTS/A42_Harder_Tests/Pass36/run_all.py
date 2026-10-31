"""Regenerate the experiment, figure and numerical summary."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
with (root / 'run.log').open('w') as log:
    subprocess.run([sys.executable, str(root / 'run.py')], cwd=root,
                   stdout=log, stderr=subprocess.STDOUT, check=True)
for script in ('plot.py', 'report.py'):
    subprocess.run([sys.executable, str(root / script)], cwd=root, check=True)
print('Finished. Read README.txt, results.json and quantum_backreaction.png.')

