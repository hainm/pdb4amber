import parmed as pmd
from pdb4amber import Leapify
from .utils import get_fn

def test_minimization():
    pdb_fn = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fn)
    leapify = Leapify(parm)
    leapify.minimize()
