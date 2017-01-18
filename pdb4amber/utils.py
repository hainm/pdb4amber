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

# copied from $AMBERHOME/updateutils/utils.py
def which(program):
   """ Looks for a program in the PATH """
   def is_exe(prog):
      return os.path.exists(prog) and os.access(prog, os.X_OK)

   # See if program has a full path. If it is, do not search PATH
   fpath, fname = os.path.split(program)
   if fpath:
      return is_exe(program)

   for d in os.environ['PATH'].split(os.pathsep):
      if is_exe(os.path.join(d, program)):
         return os.path.join(d, program)

   return None

def amberbin(program_str):
    amberhome = os.environ.get('AMBERHOME', '')
    program = os.path.join(amberhome, 'bin', program_str)
    if os.path.exists(program):
        return program
    else:
        return ''
