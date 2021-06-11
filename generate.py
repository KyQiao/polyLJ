import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('agg')

# rules for lammps read atomfile
#
# The rules for formatting the file are as follows.
# Each time a set of per-atom values is read, a non-blank line is searched for in the file.
# A comment character “#” can be used anywhere on a line; text starting with the comment character is stripped. Blank lines are skipped.
# The first “word” of a non-blank line, delimited by white-space, is read as the count N of per-atom lines to immediately follow.
# N can be the total number of atoms in the system, or only a subset. The next N lines have the following format

# ID value

# where ID is an atom ID and value is the per-atom numeric value that will be assigned to that atom. IDs can be listed in any order.


def guassian():
    pass


def binaryGaussian(N, small, vs, large, vl, ratio):
    # ratio: small particle: Ntotal
    Ns = int(N*ratio)
    Nb = N-Ns
    xs = np.random.normal(small, vs, Ns)
    xb = np.random.normal(large, vl, Nb)
    return np.arange(1, N+1), np.hstack((xs, xb))


def powerlaw():
    pass


def atomfile(file, ID, charge):
    with open(file, 'w') as f:
        assert len(ID) == len(charge)
        f.write("{} \n".format(len(ID)))
        for i in range(len(ID)):
            f.write("{:<10d}  {:.7f} \n".format(ID[i], charge[i]))
    print("atom file generated")


def testResult(charge, filename):
    plt.hist(charge, bins=100, density=True)
    plt.savefig(filename+'.png')
    print("size figure generated")
    # plt.show()


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("N", help="Number of total particles", type=int)
    parser.add_argument(
        "-s",  help="Size of small particles,default 1", default=1, type=float)
    parser.add_argument(
        "-b",  help="Size of big particles,default 1.4", default=1.4, type=float)
    parser.add_argument(
        "-vs",  help="Variance of small particles, default 0.033", default=0.033, type=float)
    parser.add_argument(
        "-vb",  help="Variance of big particles, default 0.042", default=0.042, type=float)
    parser.add_argument(
        "-r",  help="Ratio of small particles in total number, default 0.5", default=0.5, type=float)
    parser.add_argument("--addfig",  help="Produce figure of distribution",
                        action="store_true", default=False)

    return parser


if __name__ == "__main__":
    import sys
    p = make_parser()
    args = p.parse_args()
    # args = p.parse_args(sys.argv[1:])
    N = args.N
    ID, charge = binaryGaussian(
        N, args.s, args.vs, args.b, args.vb, args.r)
    if args.addfig:
        testResult(charge, str(N)+"atom")
    atomfile(str(N)+"atom.xyz", ID, charge)
