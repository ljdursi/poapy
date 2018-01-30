from __future__ import print_function
try:
    from builtins import zip
    from builtins import str
except ImportError:
    pass


def readfasta(infile):
    labels = []
    sequences = []

    curlabel = None
    cursequence = ""

    def updatelists():
        if len(cursequence) is not 0:
            sequences.append(cursequence)
            if curlabel is not None:
                labels.append(curlabel)
            else:
                labels.append('seq'+str(len(sequences)))

    for line in infile:
        if line[0] == ">":
            updatelists()
            cursequence = ""
            curlabel = line[1:].strip()
        else:
            cursequence += line.strip()

    updatelists()
    return list(zip(labels, sequences))


def writefasta(outfile, entries):
    for label, seq in entries:
        print(">" + label, file=outfile)
        print(seq, file=outfile)


if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    results = readfasta(args.infile)
    for result in results:
        print(result)
