"""Regression tests for the local tools listed in usegalaxy-eu-tools PR 1173.

Run with pytest; helper dependencies are the versions pinned in the wrappers.
CLI cases also require bowtie 1.2.0, bedtools 2.30.0, samtools 1.13,
Cutadapt 3.7 and ucsc-bedgraphtobigwig 377 on PATH.
"""
import gzip
import importlib.util
import io
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

from Bio import SeqIO
import pysam
import pytest
from sqlitedict import SqliteDict

ROOT = Path(__file__).resolve().parents[1]


def load(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(script, *args, cwd):
    return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)],
                          cwd=cwd, text=True, capture_output=True)


@pytest.mark.parametrize('five,three', [(2, 5), (2, 0), (0, 5), (0, 0)])
@pytest.mark.parametrize('description', ['', ' full header retained'])
def test_umi_bases_and_qualities(five, three, description):
    module = load('UMI_riboseq_tool/UMI.py')
    seq = 'AACCGGTTACGTACGTACGT'
    qual = 'ABCDEFGHIJKLMNOPQRST'
    record = SeqIO.read(io.StringIO(f'@read{description}\n{seq}\n+\n{qual}\n'), 'fastq')
    output = io.StringIO()
    module.process_fastq_record(record, output, five, three)
    end = len(seq) - three
    lines = output.getvalue().splitlines()
    assert lines == [f'@read_{seq[:five] + seq[end:]}{description}',
                     seq[five:end], '+', qual[five:end]]


@pytest.mark.parametrize('compressed', [False, True])
def test_umi_cli_compression(tmp_path, compressed):
    data = ROOT / 'UMI_riboseq_tool/test-data'
    output = tmp_path / 'output'
    result = run('UMI_riboseq_tool/UMI.py', data / ('umi.fastq.gz' if compressed else 'umi.fastq'),
                 output, str(compressed), 2, 5, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    content = gzip.open(output, 'rt').read() if compressed else output.read_text()
    assert content == (data / 'umi.expected.fastq').read_text()


@pytest.mark.parametrize('five,three', [(-1, 2), (30, 5)])
def test_umi_rejects_invalid_lengths(five, three):
    module = load('UMI_riboseq_tool/UMI.py')
    record = SeqIO.read(io.StringIO('@read\nACGT\n+\nIIII\n'), 'fastq')
    with pytest.raises(ValueError):
        module.process_fastq_record(record, io.StringIO(), five, three)


def test_chromosome_sizes_are_literal_and_include_final_record():
    module = load('get_chrom_sizes/calculating_chrom.sizes.py')
    output = io.StringIO()
    name = 'sequence_' + 'x' * 130
    module.chromosome_sizes(io.StringIO(f'>1 description\nAC\nGT\n\n>{name}\nAA\n'),
                            'prefix%$(literal)', output)
    assert output.getvalue() == f'prefix%$(literal)1\t4\nprefix%$(literal){name}\t2\n'


@pytest.mark.parametrize('mode', ['offset', 'weight'])
def test_profile_coordinates_all_reads_and_input_immutability(tmp_path, monkeypatch, mode):
    module = load('create_ribosome_profile/bam_to_ribosome_profile.py')
    data = ROOT / 'create_ribosome_profile/test-data'
    source = tmp_path / 'source with spaces'
    source.mkdir()
    bam = source / 'reads; literal.bam'
    shutil.copyfile(data / 'profile.bam', bam)
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    output = tmp_path / 'profile with spaces.bedgraph'
    monkeypatch.chdir(tmp_path)
    module.main(str(bam), str(data / 'profile.fasta'), str(output), mode, 15)
    assert output.read_text() == (data / f'profile.{mode}.bedgraph').read_text()
    assert {p.name: p.read_bytes() for p in source.iterdir()} == before
    if mode == 'weight':
        assert sum(float(line.split()[3]) for line in output.read_text().splitlines()) == 2


def test_weight_empty_input():
    assert load('create_ribosome_profile/bam_to_ribosome_profile.py').run_weight_centered([]) == {}


@pytest.mark.parametrize('strand,expected', [('+', [9, 11, 19, 21]), ('-', [21, 19, 11, 9])])
def test_trips_exon_boundaries(strand, expected):
    module = load('trips_bam_to_sqlite/bam_to_sqlite.py')
    annotation = {'tx': {'chrom': 'chr1', 'strand': strand, 'exons': [(10, 12), (20, 22)]}}
    assert [module.tran_to_genome('tx', pos, annotation)[1] for pos in [0, 2, 3, 5]] == expected


def test_trips_perfect_ambiguous_read_is_counted():
    module = load('trips_bam_to_sqlite/bam_to_sqlite.py')
    annotation = {name: {'chrom': name, 'strand': '+', 'exons': [(1, 100)]}
                  for name in ['tx1', 'tx2']}
    reads = {'read': [['tx1', 5, [('NM', 0)]], ['tx2', 5, [('NM', 0)]]]}
    result = {}
    assert module.processor(reads, result, annotation, {}, 'A' * 28, {}) == 1
    assert result['tx1']['ambig'] == {28: {5: 1}}
    assert result['tx2']['ambig'] == {28: {5: 1}}


def test_trips_final_read_statistics_description_and_source_immutability(tmp_path, monkeypatch):
    module = load('trips_bam_to_sqlite/bam_to_sqlite.py')
    data = ROOT / 'trips_bam_to_sqlite/test-data'
    source = tmp_path / 'source with spaces'
    source.mkdir()
    bam = source / 'input; literal.bam'
    shutil.copyfile(data / 'small.bam', bam)
    before = bam.read_bytes()
    monkeypatch.chdir(tmp_path)
    output = tmp_path / 'output.sqlite'
    module.process_bam(str(bam), str(data / 'small.annotation.sqlite'), str(output), 'sample with spaces')
    with SqliteDict(str(output)) as result:
        assert result['tx']['unambig'] == {28: {5: 1}, 30: {50: 1}}
        assert result['read_lengths'] == {28: 1, 30: 1}
        assert result['description'] == 'sample with spaces'
        assert result['nuc_counts']['mapped'][28][0]['A'] == 1
        assert result['nuc_counts']['mapped'][30][0]['C'] == 1
        assert result['unambiguous_all_totals']['tx'] == 2
    assert bam.read_bytes() == before
    assert list(source.iterdir()) == [bam]


def test_trips_failed_sort_propagates(tmp_path, monkeypatch):
    module = load('trips_bam_to_sqlite/bam_to_sqlite.py')
    monkeypatch.chdir(tmp_path)
    with pytest.raises(pysam.SamtoolsError):
        module.process_bam(str(tmp_path / 'missing.bam'),
                           str(ROOT / 'trips_bam_to_sqlite/test-data/small.annotation.sqlite'),
                           str(tmp_path / 'output.sqlite'))
    assert not (tmp_path / 'output.sqlite').exists()


def test_annotation_nullable_fields_and_quoted_identifiers(tmp_path):
    data = ROOT / 'trips_create_annotation/test-data'
    annotation = tmp_path / 'input.gtf'
    fasta = tmp_path / 'input.fasta'
    annotation.write_text((data / 'small.gtf').read_text().replace('test_tx', "test'tx"))
    fasta.write_text((data / 'small.fasta').read_text().replace('test_tx', "test'tx"))
    output = tmp_path / 'result.sqlite'
    result = run('trips_create_annotation/create_annotation_sqlite.py', 'test organism', annotation,
                 fasta, 0, "test'tx", 'gene1', output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr + result.stdout
    with sqlite3.connect(output) as db:
        assert db.execute('SELECT transcript, length, cds_start, cds_stop, version FROM transcripts').fetchall() == [
            ("test'tx", 9, None, None, None)]
        assert db.execute('SELECT * FROM exons').fetchall() == [("test'tx", 1, 9)]
        assert db.execute('PRAGMA integrity_check').fetchone() == ('ok',)


def test_annotation_invalid_identifier_fails(tmp_path):
    data = ROOT / 'trips_create_annotation/test-data'
    result = run('trips_create_annotation/create_annotation_sqlite.py', 'test', data / 'small.gtf',
                 data / 'small.fasta', 0, 'missing_identifier', 'gene1', tmp_path / 'output.sqlite', cwd=tmp_path)
    assert result.returncode != 0
    assert not (tmp_path / 'output.sqlite').exists()


def test_custom_track_escaping_and_multiline_rejection(tmp_path):
    module = load('generate_custom_track/construct_custom_track.py')
    output = tmp_path / 'track.txt'
    module.main('bigWig', 'https://example.org/sample.bw', 'sample "one"', 'a description', 'chr1:1-100', output)
    assert 'name="sample \\"one\\""' in output.read_text()
    assert 'https://https://' not in output.read_text()
    with pytest.raises(ValueError):
        module.main('bigWig', 'url', 'sample\ntrack', 'description', 'chr1:1-100', output)


BOWTIE = [('bowtie_genome_wrapper', 'bowtie_genomic_wrapper.py'),
          ('bowtie_transcriptome_wrapper', 'bowtie_transcriptome_wrapper.py'),
          ('bowtie_rRNA_removal_wrapper', 'bowtie_rRNA_tRNA_removal_wrapper.py')]


@pytest.mark.parametrize('folder,script', BOWTIE)
def test_bowtie_paired_suppressed_outputs_with_literal_paths(tmp_path, folder, script):
    if not shutil.which('bowtie'):
        pytest.skip('bowtie CLI unavailable')
    data = ROOT / folder / 'test-data'
    work = tmp_path / 'job with spaces; literal'
    work.mkdir()
    args = ['--threads=1', '--genomeSource=history', '--indexSettings=indexPreSet',
            '--ref=' + str(data / 'reference.fasta'), '--paired=paired', '--params=full',
            '--input1=' + str(data / 'paired_1.fastq'), '--input2=' + str(data / 'paired_2.fastq'),
            '--suppressAlign=1', '--maxInsert=300', '--mateOrient=fr', '--maxMismatches=0',
            '--output=' + str(work / 'output.sam'),
            '--output_suppressed_reads_l=' + str(work / 'left.fastq'),
            '--output_suppressed_reads_r=' + str(work / 'right.fastq')]
    result = run(folder + '/' + script, *args, cwd=work)
    assert result.returncode == 0, result.stderr
    assert (work / 'left.fastq').read_bytes() == (data / 'paired_1.fastq').read_bytes()
    assert (work / 'right.fastq').read_bytes() == (data / 'paired_2.fastq').read_bytes()
    with pysam.AlignmentFile(str(work / 'output.sam'), 'r') as sam:
        assert len(sam.references) == 3


@pytest.mark.parametrize('folder,script', BOWTIE)
def test_bowtie_preset_unmapped_output(tmp_path, folder, script):
    if not shutil.which('bowtie'):
        pytest.skip('bowtie CLI unavailable')
    data = ROOT / folder / 'test-data'
    output = tmp_path / 'unmapped.fastq'
    result = run(folder + '/' + script, '--threads=1', '--genomeSource=history',
                 '--indexSettings=indexPreSet', '--ref=' + str(data / 'reference.fasta'),
                 '--paired=single', '--params=preSet', '--input1=' + str(data / 'reads.fastq'),
                 '--output=' + str(tmp_path / 'output.sam'), '--output_unmapped_reads=' + str(output), cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert [record.id for record in SeqIO.parse(output, 'fastq')] == ['unmapped_read']


def render_command(relative, values):
    from Cheetah.Template import Template
    from lxml import etree
    xml = etree.parse(str(ROOT / relative))
    command = xml.find('command').text
    macro = xml.find('macros/import')
    if macro is not None:
        macros = etree.parse(str((ROOT / relative).parent / macro.text))
        for token in macros.findall('token'):
            command = command.replace(token.get('name'), token.text or '')
    # Galaxy flattens command-template newlines before invoking the shell.
    return " ".join(str(Template(command, searchList=[values])).splitlines())


@pytest.mark.parametrize('zero_based', [False, True])
def test_genomecov_per_base_uses_one_mode_and_tabular_content(tmp_path, zero_based):
    from types import SimpleNamespace as N
    if not shutil.which('bedtools'):
        pytest.skip('bedtools CLI unavailable')
    bed = tmp_path / 'reads.bed'
    bed.write_text('chr1\t1\t3\n')
    sizes = tmp_path / 'genome.sizes'
    sizes.write_text('chr1\t4\n')
    output = tmp_path / 'output'
    values = dict(input_type=N(input_type_select='bed', input=bed,
                              genome_file_opts=N(genome_file_opts_selector='hist', genome=sizes)),
                  report=N(report_select='bg', zero_regions='', scale=1),
                  d='-d', dz='-dz' if zero_based else '', split='', strand='', five='', three='', output=output)
    command = render_command('bedtools_genomecov/genomeCoverageBed.xml', values)
    assert '-bg' not in command and '-bga' not in command
    result = subprocess.run(command, shell=True, cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    expected = 'chr1\t1\t1\nchr1\t2\t1\n' if zero_based else 'chr1\t1\t0\nchr1\t2\t1\nchr1\t3\t1\nchr1\t4\t0\n'
    assert output.read_text() == expected


@pytest.mark.parametrize('folder', ['samtools_sort_genome', 'samtools_sort_transcriptome'])
@pytest.mark.parametrize('key,order', [('', [5, 50]), ('-n', [50, 5])])
def test_samtools_sort_actual_order_and_header(tmp_path, folder, key, order):
    from types import SimpleNamespace as N
    if not shutil.which('samtools'):
        pytest.skip('samtools CLI unavailable')
    # Swap names so coordinate and queryname orders are distinguishable.
    source = tmp_path / 'source.bam'
    data = ROOT / 'trips_bam_to_sqlite/test-data/small.bam'
    with pysam.AlignmentFile(str(data), 'rb') as old, pysam.AlignmentFile(str(source), 'wb', header=old.header) as bam:
        for read in old:
            read.query_name = 'z' if read.reference_start == 5 else 'a'
            bam.write(read)
    output = tmp_path / 'sorted.bam'
    command = render_command(folder + '/samtools_sort.xml', dict(
        prim_key_cond=N(prim_key_select=key), minhash='', input1=source, output1=output))
    result = subprocess.run(command, shell=True, cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        assert [read.reference_start for read in bam] == order
        assert bam.header['HD']['SO'] == ('queryname' if key == '-n' else 'coordinate')


def test_bigwig_conversion_coordinates_and_values(tmp_path):
    if not shutil.which('bedGraphToBigWig'):
        pytest.skip('bedGraphToBigWig CLI unavailable')
    import pyBigWig
    source = tmp_path / 'input.bedgraph'
    source.write_text('chr1\t10\t12\t2\nchr1\t2\t3\t0.5\n')
    sizes = tmp_path / 'genome.sizes'
    sizes.write_text('chr1\t20\n')
    output = tmp_path / 'output.bw'
    command = render_command('bedGraphToBigWig/bedGraphToBigWig.xml', dict(input1=source, input2=sizes, output1=output))
    result = subprocess.run(command, shell=True, cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    with pyBigWig.open(str(output)) as bw:
        assert bw.chroms() == {'chr1': 20}
        assert bw.intervals('chr1') == ((2, 3, 0.5), (10, 12, 2.0))


@pytest.mark.parametrize('bzip', [False, True])
def test_cutadapt_identical_mate_names_keep_both_inputs(tmp_path, bzip):
    import bz2
    from types import SimpleNamespace as N
    if not shutil.which('cutadapt'):
        pytest.skip('cutadapt CLI unavailable')
    class Dataset:
        element_identifier = '-identical mate name'
        def __init__(self, path):
            self.path = path
        def __str__(self):
            return str(self.path)
        def is_of_type(self, *types):
            return ('fastqsanger.bz2' if bzip else 'fastqsanger') in types
    left = '@pair/1\nACGTACGT\n+\nIIIIIIII\n'
    right = '@pair/2\nTGCATGCA\n+\nJJJJJJJJ\n'
    one = tmp_path / 'one.fastq'
    two = tmp_path / 'two.fastq'
    one.write_bytes(bz2.compress(left.encode()) if bzip else left.encode())
    two.write_bytes(bz2.compress(right.encode()) if bzip else right.encode())
    values = dict(library=N(type='paired', input_1=Dataset(one), input_2=Dataset(two),
                           r1=N(adapters=[], anywhere_adapters=[], front_adapters=[], cut=0),
                           r2=N(adapters2=[], anywhere_adapters2=[], front_adapters2=[], cut2=0, quality_cutoff2='')),
                  output_selector=[],
                  adapter_options=N(error_rate=0.1, times=1, overlap=3, no_indels='', match_read_wildcards='', action='trim', revcomp=''),
                  filter_options=N(discard_trimmed='', discard_untrimmed='', minimum_length='', maximum_length='',
                                   length_R2_options=N(length_R2_status='False'), max_n='', pair_filter='', max_expected_errors='', discard_cassava=''),
                  read_mod_options=N(quality_cutoff='0', nextseq_trim='0', trim_n='', strip_suffix='', shorten_options=N(shorten_values='False'),
                                     length_tag='', rename='', zero_cap=''))
    command = render_command('cutadapt/cutadapt.xml', values)
    result = subprocess.run(command, shell=True, cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    def read_output(name):
        return bz2.open(tmp_path / (name + '.fq.bz2'), 'rt').read() if bzip else (tmp_path / (name + '.fq')).read_text()
    assert read_output('out1') == left
    assert read_output('out2') == right


@pytest.mark.parametrize('strand', ['+', '-'])
def test_annotation_coding_coordinates_and_gc_slices(tmp_path, strand):
    data = ROOT / 'trips_create_annotation/test-data'
    annotation = tmp_path / 'coding.gtf'
    attributes = 'gene_id "gene1"; transcript_id "test_tx";'
    annotation.write_text(''.join(
        f'chr1\ttest\t{feature}\t100\t108\t.\t{strand}\t0\t{attributes}\n'
        for feature in ['transcript', 'exon', 'CDS']))
    output = tmp_path / 'result.sqlite'
    result = run('trips_create_annotation/create_annotation_sqlite.py', 'test', annotation,
                 data / 'small.fasta', 0, 'test_tx', 'gene1', output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr + result.stdout
    with sqlite3.connect(output) as db:
        assert db.execute('SELECT exon_start, exon_stop FROM exons').fetchall() == [(100, 108)]
        assert db.execute('SELECT cds_start, cds_stop, strand FROM transcripts').fetchall() == [(1, 9, strand)]
        gc, five_gc, cds_gc, three_gc = db.execute('SELECT gc, five_gc, cds_gc, three_gc FROM transcripts').fetchone()
        assert gc == round(100 / 9, 2)
        assert cds_gc == pytest.approx(gc)
        assert five_gc is None and three_gc is None


def test_trips_read_name_is_a_valid_read_identifier(tmp_path, monkeypatch):
    module = load('trips_bam_to_sqlite/bam_to_sqlite.py')
    data = ROOT / 'trips_bam_to_sqlite/test-data'
    bam = tmp_path / 'read_name.bam'
    with pysam.AlignmentFile(str(data / 'small.bam'), 'rb') as source, pysam.AlignmentFile(str(bam), 'wb', header=source.header) as output:
        read = next(source)
        read.query_name = 'read_name'
        output.write(read)
    monkeypatch.chdir(tmp_path)
    output = tmp_path / 'result.sqlite'
    module.process_bam(str(bam), str(data / 'small.annotation.sqlite'), str(output))
    with SqliteDict(str(output)) as result:
        assert result['tx']['unambig'] == {28: {5: 1}}
        assert result['read_lengths'] == {28: 1}
