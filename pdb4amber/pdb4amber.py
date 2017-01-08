import os
import sys
import subprocess
from itertools import chain
from optparse import OptionParser
import parmed

try:
    from io import StringIO
except ImportError:
    from cStringIO import StringIO

# TODO: include in ParmEd?
from .residue import (RESPROT, RESNA, RESSOLV, 
    RESSUGAR, AMBER_SUPPORTED_RESNAMES
)

__version__ = '1.3'


def assign_his(parm):
    ''' Assign correct name for Histidine based on the atom name

    Parameters
    ----------
    parm : parmed.Structure (or derived)

    Returns
    -------
    parm : updated `parm`
    '''
    amber_his_names = set(['HID', 'HIE' 'HIP'])
    possible_names = set(['HIS', ]) | amber_his_names

    for residue in parm.residues:
        if residue.name in possible_names:
            atom_name_set = set(atom.name for atom in residue.atoms)
            if 'HD1' in atom_name_set and 'HE2' in atom_name_set:
                residue.name = 'HIP'
            elif 'HD1' in atom_name_set and 'HE2' not in atom_name_set:
                residue.name = 'HID'
            elif 'HD1' not in atom_name_set and 'HE2' in atom_name_set:
                residue.name = 'HIE'
            else:
                residue.name = 'HIE'
    return parm


def constph(parm):
    """ Update AS4, GL4, HIP for constph.

    Parameters
    ----------
    parm : parmed.Structure or derived class

    Returns
    -------
    parm : updated `parm`
    """
    for residue in parm.residues:
        if residue.name == 'ASP':
            residue.name = 'AS4'
        elif residue.name == 'GLU':
            residue.name = 'GL4'
        elif residue.name == 'HIS':
            residue.name = 'HIP'
        else:
            pass
    return parm


def find_disulfide(parm):
    """ return set of cys-cys pairs

    Parameters
    ----------
    parm : parmed.Structure (or derived class)

    Returns
    -------
    cys_cys_set : Set[List[int, int]]
    """
    residues = [res for res in parm.residues if 'CYS' in res.name]

    cys_cys_set = set()
    for residue in residues:
        for atom in residue.atoms:
            if 'SG' in atom.name:
                cys_cys_set.add((atom.bonds[0].atom1.residue.idx,
                                 atom.bonds[0].atom2.residue.idx))
    return sorted(cys_cys_set)

def rename_cys_to_cyx(parm, cys_cys_set):
    """ Rename CYS to CYX of having S-S bond.

    Parameters
    ----------
    parm : parmed.Structure (or derived class)
    cys_cys_set : Set[List[int, int]]
    """
    for index in chain.from_iterable(cys_cys_set):
        residue = parm.residues[index]
        residue.name = 'CYX'

def find_non_starndard_resnames(parm):
    return set([res.name for res in parm.residues
                if res.name.strip() not in AMBER_SUPPORTED_RESNAMES])

def find_gaps(parm):
    return []

def find_incomplete(parm):
    return []

def run(arg_pdbout, arg_pdbin,
        arg_nohyd=False,
        arg_dry=False,
        arg_prot=False,
        arg_noter=False,
        arg_constph=False,
        arg_mostpop=False,
        log=None,
        arg_reduce=False,
        arg_model=0,
        arg_elbow=False
        ):
    stderr = sys.stderr
    # if log is not None:
    #     sys.stderr = writer(log)
    filename, extension = os.path.splitext(arg_pdbout)
    if arg_pdbin == 'stdin':
        pdbin = StringIO(sys.stdin.read())
    else:
        pdbin = arg_pdbin

    # optionally run reduce on input file
    if arg_reduce:
        if not hasattr(pdbin, 'read'):
            pdb_fh = open(pdbin, 'r')
        else:
            pdb_fh = pdbin
        try:
            reduce = os.path.join(os.getenv('AMBERHOME', ''),
                                  'bin', 'reduce')
            if not os.path.exists(reduce):
                reduce = 'reduce'
            process = subprocess.Popen([reduce, '-BUILD', '-NUC', '-'], stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate(str.encode(pdb_fh.read()))
            out = out.decode()
            err = err.decode()
            if process.wait():
                print >> sys.stderr, ("REDUCE returned non-zero exit status: "
                                      "See reduce_info.log for more details")
                open('reduce_info.log', 'w').write(err)
            # print out the reduce log even if it worked
            else:
                open('reduce_info.log', 'w').write(err)
            pdbh = StringIO(out)
            # not using load_file since it does not read StringIO
            parm = parmed.read_PDB(pdbh)
        finally:
            pdb_fh.close()
    else:
        if hasattr(pdbin, 'read'):
            # StringIO
            # need to use read_PDB
            parm = parmed.read_PDB(pdbin)
        else:
            parm = parmed.load_file(pdbin)

    # remove hydrogens if option -y is used:==============================
    if arg_nohyd:
        parm.strip('@H=')

    # find non-standard Amber residues:===================================
    #   TODO: why does the following call discard the return array of
    #         non-standard residue names?
    ns_names = find_non_starndard_resnames(parm)

    # write to pdb ifor non-starnd residues
    ns_mask = ':' + ','.join(ns_names)
    if ns_mask != ':':
        parm[ns_mask].save('non_prot.pdb', overwrite=True)
    # ns_names = []
    # if arg_elbow:
    #     ns_names = find_non_starndard_resnames_elbow(parm)

    # keep only protein:==================================================
    if arg_prot:
        parm.strip('!:' + ','.join(RESPROT))

    # remove water if -d option used:=====================================
    if arg_dry:
        water_mask = ':' + ','.join(parmed.residue.WATER_NAMES)
        parm.strip(water_mask)

    # find histidines that might have to be changed:=====================
    if arg_constph:
        constph(parm)
    else:
        assign_his(parm)

    # find possible S-S in the final protein:=============================
    sslist = find_disulfide(parm)
    rename_cys_to_cyx(parm, sslist)

    # find possible gaps:==================================================
    gaplist = find_gaps(parm)

    # count heavy atoms:==================================================
    find_incomplete(parm)

    # =====================================================================
    # make final output to new PDB file
    # =====================================================================
    parm.coordinates = parm.get_coordinates()[arg_model]
    write_kwargs = dict(overwrite=True)
    if arg_pdbout == 'stdout':
        output = StringIO()
        parm.write_pdb(output)
    else:
        output = arg_pdbout
        parm.save(output, overwrite=True)
    if isinstance(output, StringIO):
        output.seek(0)
        print(output.read())
    return ns_names, gaplist, sslist


def main():
    parser = OptionParser(version=__version__)
    parser.add_option("input", nargs='?',
                      help="PDB input file                      (default: stdin)",
    parser.add_option("-i", "--in", metavar="FILE", dest="pdbin",
                      help="PDB input file                      (default: stdin)",
                      default='stdin')
    parser.add_option("-o", "--out", metavar="FILE", dest="pdbout",
                      help="PDB output file                     (default: stdout)",
                      default='stdout')
    parser.add_option("-y", "--nohyd", action="store_true", dest="nohyd",
                      help="remove all hydrogen atoms           (default: no)")
    parser.add_option("-d", "--dry", action="store_true", dest="dry",
                      help="remove all water molecules          (default: no)")
    parser.add_option("-p", "--prot", action="store_true", dest="prot",
                      help="keep only Amber-compatible residues (default: no)")
    parser.add_option("--noter", action="store_true", dest="noter",
                      help="remove TER, MODEL, ENDMDL cards     (default: no)")
    parser.add_option("--constantph", action="store_true", dest="constantph",
                      help="rename GLU,ASP,HIS for constant pH simulation")
    parser.add_option("--most-populous", action="store_true", dest="mostpop",
                      help="keep most populous alt. conf. (default is to keep 'A')")
    parser.add_option("--reduce", action="store_true", dest="reduce",
                      help="Run Reduce first to add hydrogens.  (default: no)")
    parser.add_option("--model", type="int", dest="model", default=0,
                      help="Model to use from a multi-model pdb file (integer).  (default: use all models)")
    (opt, args) = parser.parse_args()

    if opt.pdbin == opt.pdbout:
        raise RuntimeError("The input and output file names cannot be the same!\n")

    run(arg_pdbout=opt.pdbout,
        arg_pdbin=opt.pdbin,
        arg_nohyd=opt.nohyd,
        arg_dry=opt.dry,
        arg_prot=opt.prot,
        arg_noter=opt.noter,
        arg_constph=opt.constantph,
        arg_mostpop=opt.mostpop,
        arg_reduce=opt.reduce,
        arg_model=opt.model)

if __name__ == '__main__':
    main()
