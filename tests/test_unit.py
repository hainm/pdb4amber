from utils import get_fn
from pdb4amber.pdb4amber import assign_his, constph
from pdb4amber import pdb4amber
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

def test_find_disulfide():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    cys_cys_set = pdb4amber.find_disulfide(parm)
    assert sorted(cys_cys_set) == [(5, 126), (29, 114), (63, 79), (75, 93)]
