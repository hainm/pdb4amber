import parmed as pmd
from pdb4amber import Leapify
from pdb4amber.leap_runner import run_tleap

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

leap_input = """
source leaprc.protein.ff14SB
source leaprc.water.tip3p
pdb = loadpdb {input_pdb}
solvateBox pdb TIP3PBOX 15
saveamberparm pdb {prmtop} {rst7}
"""
def test_run_tleap():
    pdb_fn = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fn)
    parm.strip(':HOH')
    print(parm)

    new_parm = run_tleap(parm, ns_names=[],
            gaplist=[],
            sslist=[],
            leap_input=leap_input)
    assert len(new_parm.atoms) == 22884
    assert 'WAT' in set(res.name for res in new_parm.residues)
