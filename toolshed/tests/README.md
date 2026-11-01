The regression suite covers only the local tools listed in
https://github.com/usegalaxy-eu/usegalaxy-eu-tools/pull/1173.

Use Python 3.10 and install `requirements-eu-review.txt`. Put these executables
on PATH for the CLI cases:

- Bowtie 1.2.0 (its Conda package requires a separate Python 3.6 environment).
- samtools 1.13, bedtools 2.30.0, ucsc-bedgraphtobigwig 377 and Cutadapt 3.7
  (these can share a Python 3.9 Conda environment).

Run `python -m pytest -q toolshed/tests/test_eu_review.py` from the repository
root. CLI cases explicitly skip when their executable is unavailable. The
Cheetah command tests flatten newlines as Galaxy does before shell execution.

The small fixtures are synthetic. `reference.fasta` was generated with Python's
`random.Random(1173)` and contains one unique reference and two identical
references to exercise suppression of multiple mappings. Profile fixtures contain
two 28-base reads on opposite strands. The TRIPS counting fixture has two reads
of different lengths and sequences, so losing the final group or using the next
read's sequence changes the asserted result. The annotation fixture contains one
noncoding transcript and checks SQL NULLs and quoted identifiers. Additional generated coding annotations check both strands, exon coordinates,
and the CDS GC slice. Existing large fixtures are retained.
