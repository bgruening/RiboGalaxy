"""Write FASTA identifiers and lengths without interpreting the prefix."""

import sys


def chromosome_sizes(fasta, prefix, output):
    name = None
    length = 0
    for line in fasta:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                output.write(f"{prefix}{name}\t{length}\n")
            name = line[1:].split()[0]
            length = 0
        else:
            if name is None:
                raise ValueError("Sequence found before the first FASTA header")
            length += len("".join(line.split()))
    if name is not None:
        output.write(f"{prefix}{name}\t{length}\n")


if __name__ == "__main__":
    with open(sys.argv[1]) as fasta, open(sys.argv[3], "w") as output:
        chromosome_sizes(fasta, sys.argv[2], output)
