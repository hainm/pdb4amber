import os
from ..utils import which, easy_call

def run(command):
    root = 'tmp_jamber'
    nabin = root + '.nab'
    nabout = root + '.out'
    nabc = root + '.c'

    with open(nabin, 'w') as fh:
        fh.write(command)
    nab_bin = which('nab')
    prefix = os.path.abspath(
        os.path.join(os.path.dirname(nab_bin), '..'))
    os.environ['AMBERHOME'] = prefix
    build_command = [
            nab_bin,
            nabin,
            '-o',
            nabout
    ]
    easy_call(build_command)
    easy_call(['./{}'.format(nabout)])
    os.unlink(nabin)
    os.unlink(nabout)
    os.unlink(nabc)
