import subprocess
import parmed

from .template import default_force_field, leap_template


def run_tleap(parm, ns_names, gaplist, sslist, forcefield_cmd=''):
    # adapted from amber_adaptbx code in phenix
    '''

    Parameters
    ----------
    parm : parmed.Structure
    ns_names : List[str]
    gaplist : List[int]
    sslist : List[Tuple[int, int]]
    forcefield_cmd : str, optional
        If given, use  this for force fields assignment
    '''
    input_pdb = 'x.pdb'
    prmtop = 'x.prmtop'
    rst7 = 'x.rst7'
    # input PDB file
    parm.write_pdb(input_pdb, altlocs='first')

    tleap_input_file = "leap.in"
    f = open(tleap_input_file, "w")

    # box
    box = parm.box
    if box is not None:
        box_info = 'set x box { %s  %s  %s }' % (box[0], box[1], box[2])
    else:
        box_info = ''

    # Now we can assume that we are dealing with AmberTools16:
    more_force_fields = ''
    if not forcefield_cmd:
        for res in ns_names:
            more_force_fields += '%s = loadmol2 %s.mol2\n' % (res, res)
            more_force_fields += 'loadAmberParams %s.frcmod\n' % res

    #  more_leap_cmds 
    more_leap_cmds = ''
    if gaplist:
        for d, res1, resid1, res2, resid2 in gaplist:
            more_leap_cmds += 'deleteBond x.%d.C x.%d.N\n' % (resid1, resid2)

    #  process sslist
    if sslist:
        for resid1, resid2 in sslist:
            more_leap_cmds += 'bond x.%d.SG x.%d.SG\n' % (resid1, resid2)

    leap_string = leap_template.format(
        force_fields=default_force_field,
        more_force_fields=forcefield_cmd,
        box_info=box_info,
        input_pdb=input_pdb,
        prmtop=prmtop,
        rst7=rst7,
        more_leap_cmds=more_leap_cmds)

    f.write(leap_string)
    f.close()

    # strangely tleap appends to the logfile so must delete first
    cmd = ['tleap', '-f', tleap_input_file]
    try:
        output = subprocess.check_output(cmd)
        return parmed.load_file(prmtop, rst7)
    except subprocess.CalledProcessError as e:
        print(e.sdtout)
        raise e
