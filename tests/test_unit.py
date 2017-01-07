from utils import get_fn
from pdb4amber.pdb4amber import assign_his, constph
import parmed as pmd

def test_assign_his():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    his_residues = [res.name for res in parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HIS']

    assign_his(parm)
    his_residues = [res.name for res in parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HID']

def test_constph():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    resnames_before = set(res.name for res in parm.residues)

    constph(parm)
    resnames_after = set(res.name for res in parm.residues)
    assert sorted(resnames_after - resnames_before) == sorted({'HIP', 'AS4', 'GL4'})
