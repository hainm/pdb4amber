import os
import sys
import math
import subprocess
from itertools import chain
import argparse
import parmed

import logging

logger = logging.getLogger('pdb4amber_log')
logger.setLevel(logging.DEBUG)

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

PY3 = sys.version_info[0] == 3
if PY3:
    string_types = str
else:
    string_types = basestring

# TODO: include in ParmEd?
from .residue import (RESPROT, AMBER_SUPPORTED_RESNAMES,
                      HEAVY_ATOM_DICT,
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
            atom_name_set = sorted(set(atom.name for atom in residue.atoms
                                       if atom.atomic_number == 1))
            if set(['HD1', 'HE1', 'HE2']).issubset(atom_name_set):
                residue.name = 'HIP'
            elif 'HD1' in atom_name_set:
                residue.name = 'HID'
            else:
                residue.name = 'HIE'
    return parm


def find_missing_heavy_atoms(parm, heavy_atom_dict=HEAVY_ATOM_DICT):
    residue_collection = []
    for residue in parm.residues:
        if residue.name in heavy_atom_dict:
            n_heavy_atoms = len(set(atom.name for atom in residue.atoms
                                    if atom.atomic_number != 1))
            n_missing = heavy_atom_dict[residue.name] - n_heavy_atoms
            if n_missing > 0:
                logger.warn('{}_{} misses {} heavy atom(s)'.format(
                    residue.name, residue.idx + 1, n_missing))
                residue_collection.append(residue)
    return residue_collection


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

def find_gaps(parm):
    # TODO: doc
    # report original resnum?
    C_atoms = []
    N_atoms = []
    gaplist = []

    #  N.B.: following only finds gaps in protein chains!
    for i, atom in enumerate(parm.atoms):
        if atom.name == 'C' and atom.residue.name in RESPROT:
            C_atoms.append(i)
        if atom.name == 'N' and atom.residue.name in RESPROT:
            N_atoms.append(i)

    nca = len(C_atoms)
    ngaps = 0

    for i in range(nca - 1):
        if parm.atoms[C_atoms[i]].residue.ter:
            continue
        # Changed here to look at the C-N peptide bond distance:
        C_atom = parm.atoms[C_atoms[i]]
        N_atom = parm.atoms[N_atoms[i + 1]]

        dx = float(C_atom.xx) - float(N_atom.xx)
        dy = float(C_atom.xy) - float(N_atom.xy)
        dz = float(C_atom.xz) - float(N_atom.xz)
        gap = math.sqrt(dx * dx + dy * dy + dz * dz)

        if gap > 2.0:
            gaprecord = (gap, C_atom.residue.name, C_atom.residue.idx+1,
                              N_atom.residue.name, N_atom.residue.idx+1)
            gaplist.append(gaprecord)
            ngaps += 1

    if ngaps > 0:
        logger.info("\n---------- Gaps (Renumbered Residues!)")
        cformat = "gap of %lf A between %s %d and %s %d"
        for i, gaprecord in enumerate(gaplist):
            logger.info(cformat % tuple(gaprecord))
    return gaplist


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
                for partner in atom.bond_partners:
                    if partner.residue.name.startswith('CY') and partner.name.startswith('SG'):
                        # use tuple for hashing
                        cys_cys_set.add(tuple(sorted((atom.residue.idx,
                                                      partner.residue.idx))))
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
    ns_names = set()
    for residue in parm.residues:
        if len(residue.name) > 3:
            rname = residue.name[:3]
        else:
            rname = residue.name
        if rname.strip() not in AMBER_SUPPORTED_RESNAMES:
            ns_names.add(rname)
    return ns_names


def find_incomplete(parm):
    return []


def add_hydrogrens(obj):
    ''' Use reduce program to add hydrogen

    Parameters
    ----------
    obj: file object or parmed.Structure or its derived class

    Returns
    -------
    parm : parmed.Structure
    '''
    try:
        if hasattr(obj, 'write_pdb'):
            # assume parmed.Structure
            fileobj = StringIO()
            obj.write_pdb(fileobj)
            fileobj.seek(0)
        else:
            fileobj = obj
            assert hasattr(fileobj, 'read')
        reduce = os.path.join(os.getenv('AMBERHOME', ''), 'bin', 'reduce')
        if not os.path.exists(reduce):
            reduce = 'reduce'
        process = subprocess.Popen([reduce, '-BUILD', '-NUC', '-'], stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate(str.encode(fileobj.read()))
        out = out.decode()
        err = err.decode()
        if process.wait():
            logger.error("REDUCE returned non-zero exit status: "
                         "See reduce_info.log for more details")
            with open('reduce_info.log', 'w') as fh:
                fh.write(err)
        # print out the reduce log even if it worked
        else:
            open('reduce_info.log', 'w').write(err)
        pdbh = StringIO(out)
        # not using load_file since it does not read StringIO
        parm = parmed.read_PDB(pdbh)
    finally:
        fileobj.close()
    return parm


def _write_pdb_to_stringio(parm):
    '''

    Parameters
    ----------
    parm : parmed.Structure or derived class
    '''
    stringio_file = StringIO()
    parm.write_pdb(stringio_file)
    stringio_file.seek(0)
    return stringio_file


def _summary(parm):
    sumdict = dict(has_altlocs=False)

    alt_residues = set()
    chains = set()
    for residue in parm.residues:
        chains.add(residue.chain)
        for atom in residue.atoms:
            if atom.other_locations:
                alt_residues.add(residue)
    # chain
    logger.info('\n----------Chains')
    logger.info('The following (original) chains have been found:')
    for chain_name in sorted(chains):
        logger.info(chain_name)

    # altlocs
    logger.info('\n---------- Alternate Locations (Original Residues!))')
    logger.info('\nThe following residues had alternate locations:')
    if alt_residues:
        sumdict['has_altlocs'] = True
        for residue in sorted(alt_residues):
            logger.info('{}_{}'.format(residue.name, residue.number))
    else:
        logger.info('None')
    return sumdict


def run(arg_pdbout, arg_pdbin,
        arg_nohyd=False,
        arg_dry=False,
        arg_prot=False,
        arg_noter=False,
        arg_constph=False,
        arg_mostpop=False,
        arg_reduce=False,
        arg_model=0,
        arg_elbow=False,
        arg_logfile='pdb4amber.log',
        arg_keep_altlocs=False,
        ):

    # always reset handlers to avoid duplication if run method is called more
    # than once
    logger.handlers = []
    if isinstance(arg_logfile, string_types):
        logfile_handler = logging.FileHandler(arg_logfile)
    elif hasattr(arg_logfile, 'write'):
        logfile_handler = logging.StreamHandler(arg_logfile)
    else:
        raise ValueError(
            "wrong arg_logfile: must be either string or file object")

    logger.addHandler(logfile_handler)
    name = arg_pdbin if not hasattr(
        arg_pdbin, '__name__') else arg_pdbin.__name__
    logger.info("\n==================================================")
    logger.info("Summary of pdb4amber for: %s" % name)
    logger.info("===================================================")

    if arg_pdbin == arg_pdbout:
        raise RuntimeError(
            "The input and output file names cannot be the same!\n")

    filename, extension = os.path.splitext(arg_pdbout)
    if arg_pdbin == 'stdin':
        if PY3:
            pdbin = StringIO(sys.stdin.read())
        else:
            pdbin = sys.stdin
    else:
        pdbin = arg_pdbin

    # optionally run reduce on input file
    if arg_reduce:
        if not hasattr(pdbin, 'read'):
            if isinstance(pdbin, parmed.Structure):
                pdb_fh = _write_pdb_to_stringio(pdbin)
            elif not parmed.formats.PDBFile.id_format(pdbin):
                pdb_fh = _write_pdb_to_stringio(parmed.load_file(pdbin))
            else:
                pdb_fh = open(pdbin, 'r')
        else:
            pdb_fh = pdbin
        parm = add_hydrogrens(pdb_fh)
    else:
        if isinstance(pdbin, parmed.Structure):
            parm = pdbin
        elif hasattr(pdbin, 'read'):
            # StringIO
            # need to use read_PDB
            parm = parmed.read_PDB(pdbin)
        else:
            parm = parmed.load_file(pdbin)

    sumdict = _summary(parm)

    # remove hydrogens if option -y is used:==============================
    if arg_nohyd:
        parm.strip('@H=')

    # find non-standard Amber residues:===================================
    #   TODO: why does the following call discard the return array of
    #         non-standard residue names?
    ns_names = find_non_starndard_resnames(parm)

    ns_mask = ':' + ','.join(ns_names)
    if ns_mask != ':':
        parm[ns_mask].save(filename + '_nonprot.pdb', overwrite=True)
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
    write_kwargs = dict()
    if not arg_keep_altlocs:
        if sumdict['has_altlocs']:
            logger.info('The alternate coordinates have been discarded.')
            if arg_mostpop:
                logger.info('Only the highest occupancy for each atom was kept.')
                write_kwargs = dict(altlocs='occupancy')
            else:
                logger.info('Only the first occurrence for each atom was kept.')
                write_kwargs = dict(altlocs='first')
            # remove altlocs label
        for atom in parm.atoms:
            atom.altloc = ''
    if arg_pdbout == 'stdout':
        output = StringIO()
        parm.write_pdb(output, **write_kwargs)
        output.seek(0)
        print(output.read())
    else:
        output = arg_pdbout
        try:
            parm.save(output, overwrite=True, **write_kwargs)
        except TypeError:
            # mol2 does not accept altloc keyword
            parm.save(output, overwrite=True)
    return ns_names, gaplist, sslist


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs='?',
                        help="PDB input file                      (default: stdin)",)
    parser.add_argument("-i", "--in", metavar="FILE", dest="pdbin",
                        help="PDB input file                      (default: stdin)",
                        default='stdin')
    parser.add_argument("-o", "--out", metavar="FILE", dest="pdbout",
                        help="PDB output file                     (default: stdout)",
                        default='stdout')
    parser.add_argument("-y", "--nohyd", action="store_true", dest="nohyd",
                        help="remove all hydrogen atoms           (default: no)")
    parser.add_argument("-d", "--dry", action="store_true", dest="dry",
                        help="remove all water molecules          (default: no)")
    parser.add_argument("-p", "--prot", action="store_true", dest="prot",
                        help="keep only Amber-compatible residues (default: no)")
    parser.add_argument("--noter", action="store_true", dest="noter",
                        help="remove TER, MODEL, ENDMDL cards     (default: no)")
    parser.add_argument("--constantph", action="store_true", dest="constantph",
                        help="rename GLU,ASP,HIS for constant pH simulation")
    parser.add_argument("--most-populous", action="store_true", dest="mostpop",
                        help="keep most populous alt. conf. (default is to keep 'A')")
    parser.add_argument("--keep-altlocs", action="store_true", dest="keep_altlocs",
                        help="Keep alternative conformations")
    parser.add_argument("--reduce", action="store_true", dest="reduce",
                        help="Run Reduce first to add hydrogens.  (default: no)")
    parser.add_argument("--pdbid", action="store_true", dest="pdbid",
                        help="fetch structure with given pdbid, "
                        "should combined with -i option.\n"
                        "Subjected to change")
    parser.add_argument("--model", type=int, dest="model", default=0,
                        help="Model to use from a multi-model pdb file (integer).  (default: use all models)")
    parser.add_argument("-l", "--logfile", metavar="FILE", dest="logfile",
                        help="log filename", default='pdb4amber.log')
    opt = parser.parse_args()

    # pdbin : {str, file object, parmed.Structure}
    if opt.input is not None:
        pdbin = opt.input
    else:
        pdbin = opt.pdbin

    if opt.pdbid:
        pdbin = parmed.download_PDB(pdbin)

    if opt.pdbin == 'stdin' and opt.input is None:
        if os.isatty(sys.stdin.fileno()):
            parser.print_help()
            sys.exit(0)
    if opt.logfile == 'stderr':
        logfile = sys.stderr
    elif opt.logfile == 'stdout':
        logfile = sys.stdout
    else:
        logfile = opt.logfile

    run(arg_pdbout=opt.pdbout,
        arg_pdbin=pdbin,
        arg_nohyd=opt.nohyd,
        arg_dry=opt.dry,
        arg_prot=opt.prot,
        arg_noter=opt.noter,
        arg_constph=opt.constantph,
        arg_mostpop=opt.mostpop,
        arg_reduce=opt.reduce,
        arg_model=opt.model,
        arg_keep_altlocs=opt.keep_altlocs,
        arg_logfile=logfile)

if __name__ == '__main__':
    main()
