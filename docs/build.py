import sys

# output : path to amber-md.github.io/pdb4amber folder
output = sys.argv[1]
import subprocess

subprocess.check_call(['make', 'html'])
subprocess.check_call('cp -rf build/html/* {}'.format(output), shell=True)
