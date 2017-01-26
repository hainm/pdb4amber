default_force_field = """
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p
source leaprc.gaff2
"""

leap_template = """
{force_fields}
{more_force_fields}
x = loadpdb {input_pdb}
{box_info}
{more_leap_cmds}
set default PBRaddi mbondi3
set default nocenter on
saveAmberParm x {prmtop} {rst7}
quit
"""
