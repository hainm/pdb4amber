import os
import subprocess
from contextlib import contextmanager
import tempfile
from shutil import rmtree

@contextmanager
def tempfolder():
  """run everything in temp folder
  """
  my_temp = tempfile.mkdtemp()
  cwd = os.getcwd()
  os.chdir(my_temp)
  yield
  os.chdir(cwd)
  rmtree(my_temp)


def get_fn(basename):
  fn = os.path.join(os.path.dirname(__file__),
                      'files',
                      basename)
  assert os.path.exists(fn), 'File must exists {}'.format(fn)
  return fn

def _has_program(pname):
    try:
        subprocess.check_call(['which', pname])
    except subprocess.CalledProcessError:
        return False
