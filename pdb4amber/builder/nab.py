import os
import subprocess
import shutil

def run(command):
    root = 'tmp_jamber'
    nabin = root + '.nab'
    nabout = root + '.out'
    nabc = root + '.c'

    with open(nabin, 'w') as fh:
        fh.write(command)
    nab_bin = shutil.which('nab')
    prefix = os.path.abspath(
        os.path.join(os.path.dirname(nab_bin), '..'))
    os.environ['AMBERHOME'] = prefix
    build_command = [
            nab_bin,
            nabin,
            '-o',
            nabout
    ]
    subprocess.check_output(build_command)
    try:
        subprocess.check_output(['./{}'.format(nabout)])
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    os.unlink(nabin)
    os.unlink(nabout)
    os.unlink(nabc)
