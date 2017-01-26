import parmed

def run_tleap(parm, ns_names, gaplist, sslist):
    input_pdb = 'x.pdb'
    prmtop = 'x.prmtop'
    rst7 = 'x.rst7'

    tleap_input_file = "leap.in"
    f = open(tleap_input_file, "wb")

    # Now we can assume that we are dealing with AmberTools16:
    f.write('source leaprc.protein.ff14SB\n')
    f.write('source leaprc.DNA.OL15\n')
    f.write('source leaprc.RNA.OL3\n')
    # f.write('source leaprc.GLYCAM_06j-1\n') #un-comment for glycoproteins
    f.write('source leaprc.water.tip3p\n')
    f.write('source leaprc.gaff2\n')
    #  (for the future: have some mechanism for modifying the above list)
    f.write('set default nocenter on\n')

    for res in ns_names:
        f.write('%s = loadmol2 %s.mol2\n' % (res, res))
        f.write('loadAmberParams %s.frcmod\n' % res)

    # input PDB file
    parm.write_pdb(input_pdb)
    f.write('x = loadpdb %s\n' % input_pdb)

    box = parm.box
    f.write('set x box { %s  %s  %s }\n' % (box[0], box[1], box[2]))

    #  process gaplist
    if gaplist:
        for d, res1, resid1, res2, resid2 in gaplist:
            f.write('deleteBond x.%d.C x.%d.N\n' % (resid1, resid2))
    #  process sslist
    if sslist:
        for resid1, resid2 in sslist:
            f.write('bond x.%d.SG x.%d.SG\n' % (resid1, resid2))

    f.write('saveAmberParm x {} {}\n'.format(prmtop, rst7))
    f.write('quit\n')
    f.close()

    # strangely tleap appends to the logfile so must delete first
    cmd = 'tleap -f %s' % tleap_input_file
    try:
        subprocess.check_call(cmd)
        return parmed.load_file(prmtop, rst7)
    except subprocess.CalledProcessError as e:
        return e.sdtout
