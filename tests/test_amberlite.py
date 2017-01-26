import subprocess
import pytest
import parmed as pmd

from utils import tempfolder, get_fn

# copied from AmberTools/test/amberlite
commands = """
pdb4amber -i 3ptb.pdb -o 3ptb_amber_dpy.pdb -dpy 
pdb4amber -i 3ptb.pdb -o 3ptb_amber.pdb
pdb4amber -i 3ptb.pdb -o 3ptb_amber_ph.pdb --constantph 
pdb4amber -i 3ptb_H.pdb -o 3ptb_amber_H.pdb
cat 3ptb.pdb | pdb4amber > 3ptb_amber_pipe.pdb 2>>error_output.txt
pdb4amber -i 4lzt.pdb -o 4lzt_amber_mostpop.pdb --most-populous -dpy 
pdb4amber -i 4lzt.pdb -o 4lzt_amber.pdb -dpy
pdb4amber -i HIStest.pdb -o HIStest_amber.pdb
pdb4amber -i 3ptb.pdb -o 3ptb_amber_reduce.pdb --reduce
pdb4amber --reduce < 3ptb.pdb > 3ptb_amber_reduce_stdin.pdb
pdb4amber -i 2mpi.pdb -o 2mpi_amber_model.pdb --model 2
""".strip().split('\n')

# pdb_pairs : List[[orig_pdb, expected_pdb]]
pdb_pairs = []
for line in commands:
    pdb_pairs.append([fn for fn in line.split() if fn.endswith('.pdb') and line])

func_inputs = [tuple([command, ] + [pdb_in, pdb_out]) for command, (pdb_in, pdb_out)
        in zip(commands, pdb_pairs)]

@pytest.mark.parametrize('command, pdb_in, pdb_out', func_inputs
)
def test_amberlite(command, pdb_in, pdb_out):
    saved_pdb_fn_abspath = get_fn('amberlite/' + pdb_out + '.save')
    pdb_in_abspath = get_fn('amberlite/' + pdb_in)
    with tempfolder():
        command = command.replace(pdb_in, pdb_in_abspath)
        subprocess.check_call(command, shell=True) 
        if 'stdin' in pdb_out:
            # e.g: pdb4amber --reduce < 3ptb.pdb > 3ptb_amber_reduce_stdin.pdb
            # Not sure why getting diff error for spacing, ack
            parm_out = pmd.load_file(pdb_out)
            saved_parmed = pmd.load_file(saved_pdb_fn_abspath)
            for atom0, atom1 in zip(parm_out.atoms, saved_parmed.atoms):
                assert atom0.name == atom1.name
                assert atom0.xx == atom1.xx
                assert atom0.xy == atom1.xy
                assert atom0.xz == atom1.xz
                assert atom0.residue.name == atom1.residue.name
        else:
            with open(pdb_out) as fh, open(saved_pdb_fn_abspath) as fh_saved:
                assert fh.read() == fh_saved.read()
