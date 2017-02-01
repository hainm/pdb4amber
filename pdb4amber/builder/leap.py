from __future__ import absolute_import
import os
from ..utils import which, easy_call

def run(command, verbose=False):
    fn = 'jamber_tmp.in'
    if 'quit' not in command:
        command = command + '\nquit'
    with open(fn, 'w') as fh:
        fh.write(command)
    build_command = '{} -f {}'.format(which('tleap'), fn).split()
    output = easy_call(build_command)
    os.unlink(fn)
    if verbose:
        print(output)
