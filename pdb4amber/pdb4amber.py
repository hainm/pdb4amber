import os
import sys
import math
import subprocess
from itertools import chain
import argparse
import parmed
from .compat import StringIO

import logging

logger = logging.getLogger('pdb4amber_log')
logger.setLevel(logging.DEBUG)


PY3 = sys.version_info[0] == 3
if PY3:
    string_types = str
else:
    string_types = basestring

# TODO: include in ParmEd?
from .residue import (RESPROT, AMBER_SUPPORTED_RESNAMES,
                      HEAVY_ATOM_DICT,
                      )

from .utils import tempfolder, amberbin

__version__ = '1.3'


class AmberPDBFixer(object):
    '''

    Parameters
    ----------
    parm : parmed.Structure
    '''
    def __init__(self, parm):
        # TODO: make a copy?
        # Why not now? parm[:] will not correctly assign TER residue
        # self.parm = parm[:]
        self.parm = parm

    def mutate(self, mask_list):
        '''

        Parameters
        ----------
        mask_list: List[Tuple[int, str]]
            [(1, 'ARG'),]
        '''
        idxs = []
        for (idx, resname) in mask_list:
            self.parm.residues[idx].name = resname
            idxs.append(str(idx+1))
        excluded_mask = ':' + ','.join(idxs) + '&!@C,CA,N,O,H'
        self.parm.strip(excluded_mask)
        return self

    def pack(self, mol, n_copies, ig=8888, grid_spacing=0.2):
        ''' add n_copies of mol to AmberPDBFixer

        Parameters
        ----------
        mol : parmed.Structure
        n_copies : number of `mol`
        ig : int
            randome seed
        grid_spacing : float

        Requires
        --------
        AddToBox program
        '''
        add_to_box_exe = amberbin('AddToBox') or 'AddToBox'
        input_pdb = 'input.pdb'
        mol_pdb = 'mol.pdb'
        out_pdb = 'out.pdb'

        with tempfolder():
            mol.save(mol_pdb, overwrite=True)
            self.parm.save(input_pdb, overwrite=True)
            command = [
                add_to_box_exe,
                '-c', input_pdb,
                '-a', mol_pdb,
                '-na', str(n_copies),
                '-IG', str(ig),
                '-G', str(grid_spacing),
                '-o', out_pdb
            ]
            try:
                subprocess.check_output(command, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(e.output.decode())
            self.parm = parmed.load_file(out_pdb)
        return self

    def assign_histidine(self):
        ''' Assign correct name for Histidine based on the atom name
    
        Returns
        -------
        parm : updated `parm`
        '''
        amber_his_names = set(['HID', 'HIE' 'HIP'])
        possible_names = set(['HIS', ]) | amber_his_names
    
        for residue in self.parm.residues:
            if residue.name in possible_names:
                atom_name_set = sorted(set(atom.name for atom in residue.atoms
                                           if atom.atomic_number == 1))
                if set(['HD1', 'HE1', 'HE2']).issubset(atom_name_set):
                    residue.name = 'HIP'
                elif 'HD1' in atom_name_set:
                    residue.name = 'HID'
                elif 'HE2' in atom_name_set:
                    residue.name = 'HIE'
                else:
                    residue.name = 'HIS'
        return self

    def strip(self, mask):
        self.parm.strip(mask)
        return self
    
    def find_missing_heavy_atoms(self, heavy_atom_dict=HEAVY_ATOM_DICT):
        residue_collection = []
        for residue in self.parm.residues:
            if residue.name in heavy_atom_dict:
                n_heavy_atoms = len(set(atom.name for atom in residue.atoms
                                        if atom.atomic_number != 1))
                n_missing = heavy_atom_dict[residue.name] - n_heavy_atoms
                if n_missing > 0:
                    logger.warn('{}_{} misses {} heavy atom(s)'.format(
                        residue.name, residue.idx + 1, n_missing))
                    residue_collection.append(residue)
        return residue_collection

    def add_missing_atoms(self):
        in_pdb = 'in.pdb'
        out_pdb = 'out.pdb'
        with tempfolder():
            self.write_pdb(in_pdb)
            with open('leap.in', 'w') as fh:
                fh.write('source leaprc.protein.ff14SB\n')
                fh.write('source leaprc.dna.bsc1\n')
                fh.write('x = loadpdb {}\n'.format(in_pdb))
                fh.write('savepdb x {}\n'.format(out_pdb))
                fh.write('quit')
            subprocess.check_output('tleap -f leap.in', shell=True)
            self.parm = parmed.load_file(out_pdb)
        return self
    
    
    def constph(self):
        """ Update AS4, GL4, HIP for constph.
    
        Returns
        -------
        parm : updated `parm`
        """
        for residue in self.parm.residues:
            if residue.name == 'ASP':
                residue.name = 'AS4'
            elif residue.name == 'GLU':
                residue.name = 'GL4'
            elif residue.name == 'HIS':
                residue.name = 'HIP'
            else:
                pass
        return self
    
    def find_gaps(self):
        # TODO: doc
        # report original resnum?
        CA_atoms = []
        C_atoms = []
        N_atoms = []
        gaplist = []
        parm = self.parm
    
        #  N.B.: following only finds gaps in protein chains!
        for i, atom in enumerate(parm.atoms):
            if atom.name in ['CA', 'CH3'] and atom.residue.name in RESPROT:
                CA_atoms.append(i)
            if atom.name == 'C' and atom.residue.name in RESPROT:
                C_atoms.append(i)
            if atom.name == 'N' and atom.residue.name in RESPROT:
                N_atoms.append(i)
    
        nca = len(CA_atoms)
        ngaps = 0
    
        for i in range(nca - 1):
            is_ter = parm.atoms[CA_atoms[i]].residue.ter
            if is_ter:
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
    
    
    def find_disulfide(self):
        """ return set of cys-cys pairs
    
        Returns
        -------
        cys_cys_set : Set[List[int, int]]
        """
        residues = [res for res in self.parm.residues if res.name in ['CYS', 'CYX']]
    
        cys_cys_resid_set = set()
        cys_cys_atomidx_set = set()
        for residue in residues:
            for atom in residue.atoms:
                if 'SG' in atom.name:
                    for partner in atom.bond_partners:
                        if (partner.residue.name.startswith('CY') and
                            partner.name.startswith('SG')):
                            # use tuple for hashing
                            cys_cys_resid_set.add(tuple(sorted((atom.residue.idx,
                                                          partner.residue.idx))))
                            cys_cys_atomidx_set.add(tuple(sorted((atom.idx,
                                                          partner.idx))))
        return sorted(cys_cys_resid_set), sorted(cys_cys_atomidx_set)
    
    
    def rename_cys_to_cyx(self, cys_cys_set):
        """ Rename CYS to CYX of having S-S bond.
    
        Parameters
        ----------
        cys_cys_set : Set[List[int, int]]
        """
        for index in chain.from_iterable(cys_cys_set):
            residue = self.parm.residues[index]
            residue.name = 'CYX'
    
    
    def find_non_starndard_resnames(self):
        ns_names = set()
        for residue in self.parm.residues:
            if len(residue.name) > 3:
                rname = residue.name[:3]
            else:
                rname = residue.name
            if rname.strip() not in AMBER_SUPPORTED_RESNAMES:
                ns_names.add(rname)
        return ns_names
    
    
    def add_hydrogen(self):
        ''' Use reduce program to add hydrogen
    
        Parameters
        ----------
        obj: file object or parmed.Structure or its derived class
    
        Returns
        -------
        parm : parmed.Structure
        '''
        try:
            fileobj = StringIO()
            self.write_pdb(fileobj)
            fileobj.seek(0)
            reduce = os.path.join(os.getenv('AMBERHOME', ''), 'bin', 'reduce')
            if not os.path.exists(reduce):
                reduce = 'reduce'
            process = subprocess.Popen([reduce, '-BUILD', '-NUC', '-NOFLIP', '-'], stdin=subprocess.PIPE,
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
            self.parm = parmed.read_PDB(pdbh)
        finally:
            fileobj.close()
        return self

    def visualize(self):
        return self.parm.visualize()

    def write_pdb(self, filename):
        '''

        Parameters
        ----------
        filename : str or file object
        '''
        self.parm.write_pdb(filename)
    
    def _write_renum(self, basename):
        ''' write original and renumbered residue index
        '''

        with open(basename + '_renum.txt', 'w') as fh:
            for residue in self.parm.residues:
                fh.write("%3s %5s    %3s %5s\n" % 
                        (residue.name, residue.number,
                         residue.name, residue.idx))

    def _write_pdb_to_stringio(self, cys_cys_atomidx_set=None,
            disulfide_conect=True, **kwargs):
        stringio_file = StringIO()
        stringio_file_out = StringIO()
        self.parm.write_pdb(stringio_file, **kwargs)
        stringio_file.seek(0)
        lines = stringio_file.readlines()

        # TODO: update ParmEd?
        if disulfide_conect:
            conect_record = [
                    'CONECT%5d%5d\n' % (idx0+1, idx1+1)
                    for (idx0, idx1) in cys_cys_atomidx_set
            ]
            conect_str = ''.join(conect_record)
            lines[-1] = conect_str + 'END\n'
            stringio_file_out.writelines(lines)
            stringio_file_out.seek(0)
            stringio_file = stringio_file_out
        return stringio_file

    def remove_water(self):
        ''' Remove waters and return new `parm` with only waters
        '''
        water_mask = ':' + ','.join(parmed.residue.WATER_NAMES)
        self.parm.strip(water_mask)
        return self

    def _summary(self):
        sumdict = dict(has_altlocs=False)
    
        alt_residues = set()
        chains = set()
        for residue in self.parm.residues:
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
        arg_strip_atom_mask=None,
        arg_mutate_string=None,
        arg_constph=False,
        arg_mostpop=False,
        arg_reduce=False,
        arg_model=0,
        arg_add_missing_atoms=False,
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

    base_filename, extension = os.path.splitext(arg_pdbout)
    if arg_pdbin == 'stdin':
        if PY3:
            pdbin = StringIO(sys.stdin.read())
        else:
            pdbin = sys.stdin
    else:
        pdbin = arg_pdbin

    # optionally run reduce on input file
    if isinstance(pdbin, parmed.Structure):
        parm = pdbin
    elif hasattr(pdbin, 'read'):
        # StringIO
        # need to use read_PDB
        parm = parmed.read_PDB(pdbin)
    else:
        parm = parmed.load_file(pdbin)

    pdbfixer = AmberPDBFixer(parm)

    pdbfixer._write_renum(base_filename)

    if arg_reduce:
        pdbfixer.add_hydrogen()


    sumdict = pdbfixer._summary()

    # remove hydrogens if option -y is used:==============================
    if arg_nohyd:
        pdbfixer.parm.strip('@H=')

    # find non-standard Amber residues:===================================
    #   TODO: why does the following call discard the return array of
    #         non-standard residue names?
    ns_names = pdbfixer.find_non_starndard_resnames()

    ns_mask = ':' + ','.join(ns_names)
    ns_mask_filename = base_filename + '_nonprot.pdb'
    if ns_mask != ':':
        pdbfixer.parm[ns_mask].save(ns_mask_filename, overwrite=True)
    else:
        with open(ns_mask_filename, 'w') as fh:
            fh.write("")

    # if arg_elbow:
    #     ns_names = find_non_starndard_resnames_elbow(parm)

    # keep only protein:==================================================
    if arg_prot:
        pdbfixer.parm.strip('!:' + ','.join(RESPROT))

    # strip atoms with given mask    =====================================
    if arg_strip_atom_mask is not None:
        pdbfixer.parm.strip(arg_strip_atom_mask)

    # remove water if -d option used:=====================================
    if arg_dry:
        water_mask = ':' + ','.join(parmed.residue.WATER_NAMES)
        water_parm = pdbfixer.parm[water_mask]
        pdbfixer.remove_water()
        water_parm.save('{}_water.pdb'.format(base_filename),
                        overwrite=True)
    # find histidines that might have to be changed:=====================
    if arg_constph:
        pdbfixer.constph()
    else:
        pdbfixer.assign_histidine()

    # find possible S-S in the final protein:=============================
    sslist, cys_cys_atomidx_set  = pdbfixer.find_disulfide()
    pdbfixer.rename_cys_to_cyx(sslist)
    with open(base_filename + '_sslink', 'w') as fh:
        for (idx0, idx1) in sslist:
            fh.write('{} {}\n'.format(idx0+1, idx1+1))

    # find possible gaps:==================================================
    gaplist = pdbfixer.find_gaps()

    mask_str_list = []
    if arg_mutate_string is not None:
        # e.g: arg_mutate_str = "3-ALA,4-GLU"
        for mask_str in arg_mutate_string.replace(';', ',').split(','):
            index, resname = mask_str.split('-')
            mask_str_list.append([int(index.strip())-1, resname.strip()])
        pdbfixer.mutate(mask_str_list)

        # mutation will remove all hydrogens
        # add back if using reduce
        if arg_reduce:
            pdbfixer.add_hydrogen()

    # count heavy atoms:==================================================
    pdbfixer.find_missing_heavy_atoms()

    if arg_add_missing_atoms:
        pdbfixer.add_missing_atoms()

    # =====================================================================
    # make final output to new PDB file
    # =====================================================================
    if arg_model >= 0:
        pdbfixer.parm.coordinates = pdbfixer.parm.get_coordinates()[arg_model]

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
        for atom in pdbfixer.parm.atoms:
            atom.altloc = ''
            for oatom in atom.other_locations.values():
                oatom.altloc = ''
    if arg_pdbout in ['stdout', 'stderr'] or arg_pdbout.endswith('.pdb'):
        output = pdbfixer._write_pdb_to_stringio(cys_cys_atomidx_set=cys_cys_atomidx_set,
                disulfide_conect=True,
                **write_kwargs)
        output.seek(0)
        if arg_pdbout in ['stdout', 'stderr']:
             print(output.read())
        else:
            with open(arg_pdbout, 'w') as fh:
                fh.write(output.read())
    else:
        # mol2 does not accept altloc keyword
        pdbfixer.parm.save(arg_pdbout, overwrite=True)
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
    parser.add_argument("-s", "--strip", dest="strip_atom_mask", default=None,
                        help="Strip given atom mask,              (default: no)")
    parser.add_argument("-m", "--mutate", dest="mutation_string", default=None,
                        help="Mutate residue")
    parser.add_argument("-p", "--prot", action="store_true", dest="prot",
                        help="keep only Amber-compatible residues (default: no)")
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
    parser.add_argument("--add-missing-atoms", action="store_true", dest="add_missing_atoms",
                        help="Use tleap to add missing atoms")
    parser.add_argument("--model", type=int, dest="model", default=-1,
                        help="Model to use from a multi-model pdb file (integer).  (default: use all models)")
    parser.add_argument("-l", "--logfile", metavar="FILE", dest="logfile",
                        help="log filename", default='stderr')
    parser.add_argument("-v", "--version", action="store_true", dest="version",
                        help="version")
    opt = parser.parse_args()

    # pdbin : {str, file object, parmed.Structure}
    if opt.version:
        print(__version__)
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
        arg_strip_atom_mask=opt.strip_atom_mask,
        arg_mutate_string=opt.mutation_string,
        arg_prot=opt.prot,
        arg_constph=opt.constantph,
        arg_mostpop=opt.mostpop,
        arg_reduce=opt.reduce,
        arg_model=opt.model-1,
        arg_keep_altlocs=opt.keep_altlocs,
        arg_add_missing_atoms=opt.add_missing_atoms,
        arg_logfile=logfile)

if __name__ == '__main__':
    main()
