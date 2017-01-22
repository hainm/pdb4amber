import os
import subprocess
import unittest
import parmed as pmd
import numpy as np
from parmed.residue import WATER_NAMES
import pytest
try:
    from io import StringIO
except ImportError:
    from cStringIO import StringIO

from pdb4amber import pdb4amber
# local
from utils import tempfolder, get_fn

pdb_fn = get_fn('4lzt/4lzt_h.pdb')

try:
    # check internet
    pmd.download_PDB('1l2y')
    internet_ok = True
except OSError:
    internet_ok = False

@unittest.skipUnless(internet_ok, 'must have internet connection to rcsb')
def test_write_model():
    orig_parm = pmd.download_PDB('1l2y')
    pdb_out = 'out.pdb'

    # default
    with tempfolder():
        subprocess.check_call(['pdb4amber', '1l2y', '--pdbid', '-o', pdb_out])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (38, 304, 3)

    # model 1
    with tempfolder():
        subprocess.check_call(['pdb4amber', '1l2y', '--pdbid', '-o', pdb_out,
            '--model', '1'])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (1, 304, 3)
        np.testing.assert_almost_equal(parm.coordinates,
                                       orig_parm.get_coordinates()[1])

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

        # water
        water_parm = pmd.load_file('out_water.pdb')
        assert set(res.name for res in water_parm.residues) == {'HOH'}

def test_write_sslink():
    pdb_out = 'out.pdb'
    pdb_fn = get_fn('4lzt/4lzt_h.pdb')
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out]
    sslink_name = 'out_sslink'
    sslink_pair = [(6, 127), (30, 115), (64, 80), (76, 94)]

    with tempfolder():
        subprocess.check_call(command)
        with open(sslink_name) as fh:
            for index, line in enumerate(fh):
                id0, idx1 = [int(i) for i in line.split()]
                assert (id0, idx1)== sslink_pair[index]

def test_constantph():
    option = '--constantph'
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        # just run to increase code coverage
        # we already test in another file
        subprocess.check_call(command)

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

def test_prot_only():
    option = '--pro'
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        res_names = set(res.name for res in orig_parm.residues)
        assert 'NO3' in res_names
        assert 'HOH' in res_names

        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        res_names = set(res.name for res in parm.residues)
        assert 'NO3' not in res_names
        assert 'HOH' not in res_names

def test_writing_renum():
    pdb_fn = get_fn('2igd/2igd_2_residues.pdb')
    pdb_out = 'out.pdb'
    command = ['pdb4amber', pdb_fn] 
    expected_lines = """
MET     1    MET     0
PRO     3    PRO     1
    """.strip().split('\n')

    with tempfolder():
        subprocess.check_call(command)
        with open('stdout_renum.txt') as fh:
            output_lines = [line.strip() for line in fh.readlines()]
            assert expected_lines == output_lines

def test_reduce_with_pdb_input():
    option = '--reduce'
    pdb_fn = get_fn('2igd/2igd.pdb')
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        atom_names = set(atom.name for atom in orig_parm.atoms if atom.atomic_number == 1)
        assert not atom_names

        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert atom_names

def test_strip_atoms():
    pdb_fn = get_fn('2igd/2igd.cif')
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out,
            '--strip', ':3-500'] 
    with tempfolder():
        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        assert len(parm.residues) == 2

def test_reduce_with_cif_input():
    option = '--reduce'
    pdb_fn = get_fn('2igd/2igd.cif')
    pdb_out = 'out.pdb'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out, option] 

    with tempfolder():
        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert atom_names

def test_stdin_stdout():
    ''' e.g: cat my.pdb | pdb4amber '''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['cat', pdb_fn, '|', 'pdb4amber'] 

    with tempfolder():
        # use shell=True since check_output return exit 1 with |
        # not sure why.
        output = subprocess.check_output(' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 574

def test_fetch_pdbid():
    ''' e.g: pdb4amber 1l2y --pdbid '''
    pdb_fn = '1l2y'
    command = ['pdb4amber', pdb_fn, '--pdbid']

    with tempfolder():
        output = subprocess.check_output(command).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 304
        
def test_fetch_pdbid_and_use_reduce():
    ''' e.g: pdb4amber 1tsu --pdbid --reduce'''
    pdb_fn = '1tsu'
    command = ['pdb4amber', pdb_fn, '--pdbid', '--reduce']

    with tempfolder():
        output = subprocess.check_output(command).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 3174

def test_simplest_command_pdb4amber_mypdb():
    # pdb4amber my.pdb
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['pdb4amber', pdb_fn] 

    with tempfolder():
        output = subprocess.check_output(' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 574

def test_simplest_command_pdb4amber():
    command = ['pdb4amber']
    with tempfolder():
        output = subprocess.check_output(command).decode()
        assert 'usage: pdb4amber' in output


def test_stdin_stdout_with_reduce():
    ''' e.g: cat my.pdb | pdb4amber --reduce '''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['cat', pdb_fn, '|', 'pdb4amber', '--reduce'] 

    with tempfolder():
        # use shell=True since check_output return exit 1 with |
        # not sure why.
        output = subprocess.check_output(' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 1033

def test_write_other_formats_like_mol2():
    # mol2
    pdb_out = 'out.mol2'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out] 
    with tempfolder():
        subprocess.check_call(command)
        with open(pdb_out) as fh:
            assert fh.read().startswith('@<TRIPOS>MOLECULE')

def test_keep_altlocs():
    ''' e.g: pdb4amber 2igd.pdb --keep-altlocs --reduce'''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['pdb4amber', pdb_fn, '--keep-altlocs', '--reduce']

    with tempfolder():
        output = subprocess.check_output(command).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        res4 = parm.residues[4]
        for atom in res4.atoms:
            if atom.name.startswith('CB') or atom.name.startswith('CG'):
                assert atom.other_locations

def test_increase_code_coverage_for_small_stuff():
    with pytest.raises(RuntimeError):
        pdb4amber.run('fake.pdb', 'fake.pdb')
