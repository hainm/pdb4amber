import sys
import subprocess
from pdb4amber.pdb4amber import assign_his, constph, StringIO
from pdb4amber import pdb4amber
import parmed as pmd

# local
from utils import get_fn, tempfolder

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

    one_cys_parm = parm[':CYS'][':1']
    assert not pdb4amber.find_disulfide(one_cys_parm)

def test_strip_water():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    assert 'HOH' in set(res.name for res in parm.residues)

    water_mask = ':' + ','.join(pmd.residue.WATER_NAMES)
    parm.strip(water_mask)
    assert 'HOH' not in set(res.name for res in parm.residues)

def test_find_non_starndard_resnames():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    assert pdb4amber.find_non_starndard_resnames(parm) == {'NO3'}

def test_run_with_StringIO_log():
    stringio_file = StringIO()
    with tempfolder():
         pdb4amber.run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'),
            arg_logfile=stringio_file)
    stringio_file.seek(0)
    assert 'Summary of pdb4amber' in stringio_file.read()

def test_run_with_stderr_log():
    # dummy
    stringio_file = StringIO()
    with tempfolder():
         output = subprocess.check_call([
             'pdb4amber',
             get_fn('4lzt/4lzt_h.pdb'),
             '-o',
             'out.pdb',
             '--logfile=stderr',
         ], stderr=subprocess.STDOUT)

def test_run_with_filename_log():
    stringio_file = StringIO()

    # default
    logfile = 'pdb4amber.log'
    with tempfolder():
         pdb4amber.run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'))
         with open(logfile) as fh:
              assert 'Summary of pdb4amber' in fh.read()

    # given name
    logfile = 'hello.log'
    with tempfolder():
         pdb4amber.run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'),
                 arg_logfile=logfile)
         with open(logfile) as fh:
              assert 'Summary of pdb4amber' in fh.read()
