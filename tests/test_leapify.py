import parmed as pmd
from pdb4amber import Leapify

from utils import get_fn

def test_minimization():
    pdb_fn = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fn)
    fixer = Leapify(parm)

    def callback(xyz):
        print('hello')

    fixer.leapify()
    fixer.minimize(igb=None, saltcon=None,
            cutoff=8.,
            tol=1E-6,
            maxcyc=2,
            disp=False,
            callback=callback)
