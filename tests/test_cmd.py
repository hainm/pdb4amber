import subprocess
import parmed as pmd
from parmed.residue import WATER_NAMES

# local
from utils import tempfolder, get_fn

pdb_fn = get_fn('4lzt/4lzt_h.pdb')

def test_dry():
    option = '--dry'
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        resnames = set(res.name for res in orig_parm.residues)
        assert resnames.intersection(WATER_NAMES)

        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        resnames = set(res.name for res in parm.residues)
        assert not resnames.intersection(WATER_NAMES)
