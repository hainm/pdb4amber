import subprocess
from pdb4amber.pdb4amber import StringIO
from pdb4amber import AmberPDBFixer, run
import parmed as pmd

# local
from utils import get_fn, tempfolder

def test_assign_histidine():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    his_residues = [res.name for res in parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HIS']

    pdbfixer = AmberPDBFixer(parm)
    pdbfixer.assign_histidine()
    his_residues = [res.name for res in pdbfixer.parm.residues if res.name in {'HIS', 'HIE', 'HID', 'HIP'}]
    assert his_residues == ['HID']

def test_constph():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    resnames_before = set(res.name for res in parm.residues)
    pdbfixer = AmberPDBFixer(parm)
    pdbfixer.constph()
    resnames_after = set(res.name for res in pdbfixer.parm.residues)
    assert sorted(resnames_after - resnames_before) == sorted({'HIP', 'AS4', 'GL4'})

def test_find_disulfide():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    pdbfixer = AmberPDBFixer(parm)
    cys_cys_set = pdbfixer.find_disulfide()
    assert sorted(cys_cys_set) == [(5, 126), (29, 114), (63, 79), (75, 93)]

    one_cys_parm = parm[':CYS'][':1']
    pdbfixer_2 = AmberPDBFixer(one_cys_parm)
    assert not pdbfixer_2.find_disulfide()

def test_find_missing_heavy_atoms():
    fn = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(fn)
    parm2 = parm[':1-2&!@CB']
    parm2.save("test.pdb", overwrite=True)
    pdbfixer = AmberPDBFixer(parm2)
    assert len(pdbfixer.find_missing_heavy_atoms()) == 2
    assert 'CB' not in set(atom.name for atom in pdbfixer.parm.atoms)

    pdbfixer.add_missing_atoms()
    assert len(pdbfixer.find_missing_heavy_atoms()) == 0
    assert 'CB' in set(atom.name for atom in pdbfixer.parm.atoms)

def test_strip_water():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    assert 'HOH' in set(res.name for res in parm.residues)

    water_mask = ':' + ','.join(pmd.residue.WATER_NAMES)
    parm.strip(water_mask)
    assert 'HOH' not in set(res.name for res in parm.residues)

def test_find_non_standard_resnames():
    fn = get_fn('4lzt/4lzt_h.pdb')
    parm = pmd.load_file(fn)
    pdbfixer = AmberPDBFixer(parm)
    assert pdbfixer.find_non_starndard_resnames() == {'NO3'}

def test_run_with_StringIO_log():
    stringio_file = StringIO()
    with tempfolder():
         run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'),
            arg_logfile=stringio_file)
    stringio_file.seek(0)
    assert 'Summary of pdb4amber' in stringio_file.read()

def test_run_with_stderr_stdout_log():
    # dummy
    with tempfolder():
         subprocess.check_call([
             'pdb4amber',
             get_fn('4lzt/4lzt_h.pdb'),
             '-o',
             'out.pdb',
             '--logfile=stderr',
         ])

    with tempfolder():
         subprocess.check_call([
             'pdb4amber',
             get_fn('4lzt/4lzt_h.pdb'),
             '-o',
             'out.pdb',
             '--logfile=stdout',
         ])

def test_run_with_filename_log():
    # default
    logfile = 'pdb4amber.log'
    with tempfolder():
         run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'))
         with open(logfile) as fh:
              assert 'Summary of pdb4amber' in fh.read()

    # given name
    logfile = 'hello.log'
    with tempfolder():
         run(arg_pdbout='out.pdb', arg_pdbin=get_fn('4lzt/4lzt_h.pdb'),
                 arg_logfile=logfile)
         with open(logfile) as fh:
              assert 'Summary of pdb4amber' in fh.read()

def test_find_gaps():
    pdb_fh = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fh)
    parm_gap = parm[':1,3']
    pdbfixer = AmberPDBFixer(parm_gap)
    assert pdbfixer.find_gaps() == [(4.134579301452567, 'MET', 1, 'PRO', 2)]

def test_mutate():
    pdb_fh = get_fn('ala3_alpha.pdb')
    parm = pmd.load_file(pdb_fh)
    assert set(res.name for res in parm.residues) == {'ALA'}

    pdbfixer = AmberPDBFixer(parm)
    pdbfixer.mutate([(1, 'ARG'),])
    pdbfixer.add_missing_atoms()
    assert [res.name for res in pdbfixer.parm.residues] == ['ALA', 'ARG', 'ALA']
    assert ([atom.name for atom in pdbfixer.parm.residues[1] if atom.atomic_number==6] == 
            ['CA', 'CB', 'CG', 'CD', 'CZ', 'C'])
