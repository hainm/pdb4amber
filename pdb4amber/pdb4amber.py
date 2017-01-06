import os
import sys
from optparse import OptionParser

__version__ = '1.3'


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
        print("The input and output file names cannot be the same!\n", file=sys.stderr)
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
