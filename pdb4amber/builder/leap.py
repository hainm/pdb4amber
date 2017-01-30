from __future__ import absolute_import
import os
import subprocess
import shutil

def run(command, verbose=False):
    fn = 'jamber_tmp.in'
    if 'quit' not in command:
        command = command + '\nquit'
    with open(fn, 'w') as fh:
        fh.write(command)
    build_command = '{} -f {}'.format(shutil.which('tleap'), fn).split()
    output = subprocess.check_output(build_command).decode()
    os.unlink(fn)
    if verbose:
        print(output)
