import os
import pytraj
import parmed
from . import leap
from ..utils import tempfolder
from . import nab

__all__ = [
        'build_protein',
        'build_adna',
        'build_bdna',
        'build_arna'
]

def build_protein(seq, command=None):
    '''

    Parameters
    ----------
    seq : protein sequence
    command : None, str or list of str
        pytraj command for building secondary structure for a given residue range
        if None, build linear sequence

    Requires
    --------
    tleap, pytraj and ParmEd

    Examples
    --------
    >>> # Ala10
    >>> seq = "ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA"
    >>> parm = build_protein(seq, ['alpha:1-10'])
    >>> # visualize in Jupyter notebook
    >>> parm.visualize()
    >>> # save to pdb file
    >>> parm.save("new.pdb")

    Returns
    -------
    parmed.Structure
    '''
    pdb_fn = 'seq.pdb'
    leap_command = """
    source leaprc.protein.ff14SB
    x = sequence{%s}
    savepdb x %s
    """ % (seq, pdb_fn)

    if command is None:
        command = ''

    command_str = ' '.join(command) if isinstance(command, list) else command

    amberhome = os.getenv('AMBERHOME', '')
    if  os.path.exists(amberhome + '/dat/leap/cmd/leaprc.ff14SB'):
        leap_command = leap_command.replace('protein.ff14SB', 'ff14SB')

    with tempfolder():
        pdb_out = 'pdb_out.pdb'
        leap.run(leap_command)
        traj = pytraj.load(pdb_fn)
        pytraj.make_structure(traj, command_str)
        traj.save(pdb_out, overwrite=True)
        return parmed.load_file(pdb_out)

def solvate(parm, buffer=8.):
    pdb_out = 'my.pdb'
    leap_command = """
    source leaprc.protein.ff14SB
    source leaprc.water.tip3p
    pdb = loadpdb {input_pdb}
    solvateOct pdb TIP3PBOX {buffer}
    savepdb pdb {input_pdb}
    quit
    """.format(input_pdb=pdb_out, buffer=buffer)

    with tempfolder():
        parm.save(pdb_out, overwrite=True)
        leap.run(leap_command)
        return parmed.load_file(pdb_out)


def _nab_build(seq, filename, nuc_type='abdna'):
    command = """
    molecule m; 
    m = fd_helix( "{}", "{}", "{}" );
    putpdb("{}", m, "-wwpdb");
    """.format(nuc_type, seq, nuc_type[-3:], filename)
    with tempfolder():
        nab.run(command)
        return parmed.load_file(filename)

def build_adna(seq, filename='nuc.pdb'):
    return _nab_build(seq, filename=filename, nuc_type='adna')

def build_bdna(seq, filename='nuc.pdb'):
    return _nab_build(seq, filename=filename, nuc_type='abdna')

def build_arna(seq, filename='nuc.pdb'):
    return _nab_build(seq, filename=filename, nuc_type='arna')
