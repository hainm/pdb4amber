import os
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

def amberbin(program_str):
    amberhome = os.environ.get('AMBERHOME', '')
    program = os.path.join(amberhome, 'bin', program_str)
    if os.path.exists(program):
        return program
    else:
        return ''

def which(program):
    from distutils.spawn import find_executable
    return find_executable(program)
