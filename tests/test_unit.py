
from utils import get_fn
from pdb4amber.pdb4amber import assign_his
import parmed as pmd

def test_assign_his():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    his_residues = [res.name for res in parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HIS']

    assign_his(parm)
    his_residues = [res.name for res in parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HID']
