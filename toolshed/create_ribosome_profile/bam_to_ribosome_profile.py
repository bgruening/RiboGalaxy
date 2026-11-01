"""Generate zero-based, half-open bedGraph ribosome profiles."""

import os
import sys
import tempfile

from Bio import SeqIO
import pysam


def run_weight_centered(all_reads):
    sequence = {}
    for read in all_reads:
        positions = sorted(read.get_reference_positions())
        if read.query_length < 25 or len(positions) < 25:
            continue
        center = positions[12:-12]
        for asite in center:
            sequence[asite] = sequence.get(asite, 0) + 1.0 / len(center)
    return sequence


def run_offset(all_reads, offset):
    if offset < 0:
        raise ValueError("Offset must be non-negative")
    sequence = {}
    for read in all_reads:
        if read.query_length < 25:
            continue
        positions = sorted(read.get_reference_positions())
        if offset >= len(positions):
            raise ValueError(f"Offset {offset} exceeds aligned length for {read.query_name}")
        asite = positions[-1 - offset] if read.is_reverse else positions[offset]
        sequence[asite] = sequence.get(asite, 0) + 1
    return sequence


def main(bam_path, fasta_path, output, mode="offset", offset=0):
    if mode not in ("offset", "weight"):
        raise ValueError("Unknown profile approach: " + mode)
    with open(fasta_path) as fasta:
        chromosomes = sorted(record.id for record in SeqIO.parse(fasta, "fasta"))
    # Indexes belong to the job, including when the source BAM is read-only.
    with tempfile.TemporaryDirectory(dir=".") as workdir:
        index = os.path.join(workdir, "input.bai")
        pysam.index(bam_path, index)
        with pysam.AlignmentFile(bam_path, "rb", index_filename=index) as alignments:
            if not set(chromosomes).intersection(alignments.references):
                raise ValueError("FASTA identifiers do not match any BAM references")
            with open(output, "w") as bedgraph:
                for chrom in chromosomes:
                    if chrom not in alignments.references:
                        continue
                    reads = alignments.fetch(chrom)
                    if mode == "offset":
                        sequence = run_offset(reads, offset)
                    else:
                        sequence = run_weight_centered(reads)
                    for asite, count in sorted(sequence.items()):
                        bedgraph.write(f"{chrom}\t{asite}\t{asite + 1}\t{count}\n")


if __name__ == "__main__":
    bam_path, offset, source, history, builtin, mode, output = sys.argv[1:]
    main(bam_path, history if source == "history" else builtin,
         output, mode, int(offset))
