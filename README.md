[![Build Status](https://travis-ci.org/hainm/pdb4amber.svg?branch=master)](https://travis-ci.org/hainm/pdb4amber)
[![Coverage Status](https://coveralls.io/repos/github/Amber-MD/pdb4amber/badge.svg?branch=master)](https://coveralls.io/github/Amber-MD/pdb4amber?branch=master)

**Please do not use this repo yet**

Install
-------
```bash
python setup.py install
```

Usage
-----
```bash
$ pdb4amber --help

Usage: pdb4amber [options]

Options:
  --version            show program's version number and exit
  -h, --help           show this help message and exit
  -i FILE, --in=FILE   PDB input file                      (default: stdin)
  -o FILE, --out=FILE  PDB output file                     (default: stdout)
  -y, --nohyd          remove all hydrogen atoms           (default: no)
  -d, --dry            remove all water molecules          (default: no)
  -p, --prot           keep only Amber-compatible residues (default: no)
  --noter              remove TER, MODEL, ENDMDL cards     (default: no)
  --constantph         rename GLU,ASP,HIS for constant pH simulation
  --most-populous      keep most populous alt. conf. (default is to keep 'A')
  --reduce             Run Reduce first to add hydrogens.  (default: no)
  --model=MODEL        Model to use from a multi-model pdb file (integer).
                       (default: use all models)

```

Test
----
py.test -vs .
