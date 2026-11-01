# Review of the RiboGalaxy tools in EU PR 1173

Scope: [usegalaxy-eu/usegalaxy-eu-tools#1173](https://github.com/usegalaxy-eu/usegalaxy-eu-tools/pull/1173),
reviewed on 2026-11-01. The installation manifest at head
`18bd60076c9cc5a9ce767018696a1f95a7dcc1c7` lists 19 repositories. This branch
reviews and fixes the 14 implemented in this checkout. The five community
repositories (`sra_tools`, `sam_to_bam`, `fastx_trimmer`, `fastqc`, and
`umi_tools_dedup`) are not vendored here and were not modified or independently
reviewed. `samtools_index`, `samtools_faidx`, and the archived tool bundles are
outside this review.

The requested skill was read from
`../../galaxytools/.agents/skills/galaxy-tool-review/SKILL.md`, together with its
checklist and review evidence. The PR had no submitted reviews, inline threads,
or reported commit statuses when inspected. No GitHub comments were posted.

## Fixed findings

| Tool Shed repository | Local directory | Trigger, consequence, and fix |
| --- | --- | --- |
| `ribogalaxy_bowtie_rrna` | `bowtie_rRNA_removal_wrapper` | Saving unmapped reads changed alignment output away from SAM; disabling mapped output relied on an implicit sentinel. Keep SAM generation consistent and explicitly discard filtered mapped output. Fix paired suppressed-read destinations, preset output flags, and failure propagation. |
| `ribogalaxy_genome_mapping` | `bowtie_genome_wrapper` | Paired suppressed reads were never requested correctly and both destinations referenced mate 1. Fix both mates and preserve quality encoding, orientation, and header settings in the genomic preset. |
| `ribogalaxy_transcriptome_test_new_index` | `bowtie_transcriptome_wrapper` | Fix the same paired suppressed-read and preset-output defects, plus an empty integer default that invalidated the paired mismatch branch. |
| `ribogalaxy_cutadapt` | `cutadapt` | Identical mate identifiers produced colliding symlinks, replacing the first mate with the second; leading dashes could also become options. Use distinct fixed staging names. Correct bzip2 mate extensions and quote the suffix argument. |
| `ribogalaxy_umi_processing` | `UMI_riboseq_tool` | Headers without descriptions crashed; trimming removed an extra biological base; a zero-length 3-prime UMI incorrectly selected the whole sequence. Slice sequence and quality identically, preserve header descriptions, validate lengths, and declare uncompressed output when requested. Align the Biopython pin with the tested 1.79 environment. |
| `ribogalaxy_create_ribosome_profile` | `create_ribosome_profile` | The weight approach returned after one read and returned `None` on empty input. Process every read, preserve total weight, and test both strands. Avoid stale chromosome iterators, reject incompatible references and invalid offsets, index in job-local storage, and emit the actual four-column bedGraph datatype. |
| `ribogalaxy_bedtools_genomcov` | `bedtools_genomecov` | Per-base flags were combined with bedGraph flags and the resulting three-column output was labeled bedGraph. Select one output mode, declare per-base results tabular, and avoid passing the histogram sentinel as `-max 0`. |
| `ribogalaxy_bedgraphtobigwig` | `bedGraphToBigWig` | Execution depended on a bundled platform-specific executable. Declare `ucsc-bedgraphtobigwig=377`, use the resolved executable, sort intervals, quote paths, and document the bedGraph contract. The existing BigWig fixture matches this pinned executable byte for byte. |
| `ribogalaxy_generate_custom_track` | `generate_custom_track` | Unquoted arguments split text and paths. Quote arguments, declare Python and the staged helper, reject multiline attributes, escape quoted track attributes, and repair the duplicated URL scheme. |
| `ribogalaxy_samtools_sort_genome` | `samtools_sort_genome` | The free-text alignment tag was interpolated as an unquoted shell fragment. Quote it and validate the two-character SAM tag syntax. |
| `ribogalaxy_samtools_sort_transcriptome` | `samtools_sort_transcriptome` | Apply the same tag fix. Tests expecting coordinate output now explicitly override this wrapper's name-sort default. |
| `ribogalaxy_get_chrom_sizes_upd` | `get_chrom_sizes` | Prefix text was inserted into awk program text and identifiers were truncated or included descriptions. Replace the pipeline with a staged Python helper that treats the prefix literally, uses the complete FASTA identifier, and includes the last record. |
| `trips_viz_bam_to_sqlite` | `trips_bam_to_sqlite` | The final mapped read group and reads named `read_name` were lost; statistics used the next read's sequence; perfect ambiguous reads were omitted; the description was discarded. Correct these paths and inclusive exon-boundary conversion. Sort in job-local storage, open annotations read-only, and propagate sorting failures. |
| `trips_viz_create_annotation` | `trips_create_annotation` | FASTA/GTF strings were interpolated into SQL; quoted identifiers failed, and invalid annotations could exit successfully. Bind SQL values while retaining genuine SQL NULLs, preserve transcript identifier case, fix reverse-strand exon shifts and CDS GC slicing, use the actual transcriptome-name input, and accept declared GTF/GFF3 datatypes. |

The execution and scientific-result defects above are release-blocking for their
affected paths. Command quoting, dependency/staging declarations, datatype
contracts, and invalid tests are required correctness or reproducibility fixes.
Citations and help improvements accompany the affected tools. Existing tool IDs,
parameter names, reference table names, and output identifiers are preserved;
functional wrapper versions advance. Historical fixtures remain in place, while
broken test definitions use small, auditable replacements. The coverage histogram
fixture now reflects the input: 20 covered bases on chr1 and 500 on chr2, rather
than the old all-zero result produced by `-max 0`. Bowtie tests explicitly enable
the mapping-statistics output before asserting its presence.

All three Bowtie wrappers now explicitly require Python 3.6, compatible with the
pinned Bowtie 1.2.0 Conda package and the helpers' Python 3 syntax. The unused
samtools requirement is removed from these SAM-producing wrappers. No upstream
Bowtie version update is included.

## Validation

- `tests/test_eu_review.py`: 42 passed, none skipped, using Python 3.10,
  Biopython 1.79, pysam 0.19.0, sqlitedict 1.7.0, intervaltree 3.0.2,
  Bowtie 1.2.0, Cutadapt 3.7, bedtools 2.30.0, samtools 1.13, and
  ucsc-bedgraphtobigwig 377. Checks include real CLI execution, Cheetah rendering,
  BAM and SQLite contents, BigWig coordinates/values, and source immutability.
  SQLite integration assertions check the header, schema, and identifiers; unit
  tests query actual row values. This avoids byte-comparison failures caused by
  different SQLite library versions.
- Galaxy wrapper tests: all 37 distinct cases passed across the full run and
  targeted reruns. Helper dependencies were activated with explicit local package
  resolver scripts; the final portable SQLite assertions were also checked with
  Galaxy's assertion verifier.
- Planemo lint: zero errors across all 14 wrappers. The default run reports one
  external bio.tools lookup warning for the existing `bedtools` entry; with only
  `BioToolsValid` skipped, lint passes without warnings. The existing registry
  reference is retained.
- Python compilation and `git diff --check`: passed.

The original local Galaxy checkout could not start because it contains unresolved
merge markers outside this repository. An isolated worktree of its clean HEAD
was used for the wrapper-test attempt. CLI tests use the declared package
versions in separate local environments; this is not BioContainer or Pulsar
execution. No live administrator-provided reference tables were exercised, and
these fixtures do not establish correctness for every annotation convention or
all combinations of advanced mapping parameters.

See [tests/README.md](tests/README.md) for reproduction and fixture details.

## Review basis

- [IUC tool XML and tests](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html)
- [IUC security](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/security.html)
- [IUC remote execution](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/pulsar.html)
- [IUC dependencies](https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/package_xml.html)
- [Bowtie 1.2.0 manual](https://github.com/BenLangmead/bowtie/blob/v1.2.0/MANUAL)
- [Cutadapt 3.7 guide](https://cutadapt.readthedocs.io/en/v3.7/guide.html)
- [samtools 1.13 sort manual](https://www.htslib.org/doc/1.13/samtools-sort.html)
