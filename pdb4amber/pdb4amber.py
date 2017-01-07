import os
import sys
from itertools import chain
from optparse import OptionParser
import parmed

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
    if log is not None:
        sys.stderr = writer(log)
    filename, extension = os.path.splitext(arg_pdbout)
    pdbin = arg_pdbin

    # optionally run reduce on input file
    if arg_reduce:
        if arg_pdbin == 'stdin':
            pdbfile = sys.stdin
        else:
            pdbfile = open(arg_pdbin, 'r')
        try:
            reduce = os.path.join(os.getenv('AMBERHOME', ''),
                                  'bin', 'reduce')
            if not os.path.exists(reduce):
                reduce = 'reduce'
            process = subprocess.Popen([reduce, '-BUILD', '-NUC', '-'], stdin=pdbfile,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
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
            parm = parmed.load_file(pdbh)
        finally:
            if pdbfile is not sys.stdin:
                pdbfile.close()
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
    parm[':' + ','.join(ns_names)].save('non_prot.pdb', overwrite=True)
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
    coordinates = parm.get_coordinates()[arg_model]
    parm.write_pdb(arg_pdbout,
                   coordinates=coordinates)
    return ns_names, gaplist, sslist


def main():
    parser = OptionParser(version=__version__)
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
        print("The input and output file names cannot be the same!\n")
        sys.exit(1)

    # Make sure that if we are reading from stdin it's being directed from a pipe
    # or a file. We don't want to wait for user input that will never come.

    if opt.pdbin == 'stdin':
        if os.isatty(sys.stdin.fileno()):
            sys.exit(parser.print_help() or 1)

    run(opt.pdbout, opt.pdbin, opt.nohyd, opt.dry, opt.prot, opt.noter,
        opt.constantph, opt.mostpop, opt.reduce, opt.model)

if __name__ == '__main__':
    main()
