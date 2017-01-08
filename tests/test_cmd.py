import os
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

def test_no_hydrogen():
    option = '--nohyd'
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        atom_names = set(atom.name for atom in orig_parm.atoms if atom.atomic_number == 1)
        assert atom_names

        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert not atom_names

def test_reduce():
    option = '--reduce'
    pdb_fn = get_fn('2igd/2igd.pdb')
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 
    print(' '.join(command))

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        atom_names = set(atom.name for atom in orig_parm.atoms if atom.atomic_number == 1)
        assert not atom_names

        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert atom_names

def test_write_other_formats_like_mol2():
    # mol2
    pdb_out = 'out.mol2'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out] 
    with tempfolder():
        subprocess.check_call(command)
        with open(pdb_out) as fh:
            assert fh.read().startswith('@<TRIPOS>MOLECULE')
