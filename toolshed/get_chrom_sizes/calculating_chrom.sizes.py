# input a genome file and return a file genome.chrom.sizes to be associated with the custom build (or just have it as an output to be used later in the history.
# adapted from https://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/
from sys import argv
# python calculating_chrom.sizes.py genome_input.fa output.chrom.sizes
fasta_source = str(argv[1])
prefix = str(argv[2])
genome = str(argv[3])
builtin = str(argv[4])
output = str(argv[5])

# genome = 'test-data/test.fasta'
# output = "test-data/test_chrom.sizes"
if fasta_source == 'builtin':
	genome = builtin

chromSizesoutput = open(output,"w")

records = []
record = False
for line in open(genome, 'r').readlines():
	if line[0] == '>':
		if record:
			records.append(record)
		record = [line.strip("\n").split(' ')[0][1:], 0]

	else:
		sequence = line.strip('\n')
		record[1] += len(sequence)

if record not in records:
	records.append(record)



for seq_record in records:
	if prefix != 'none':
		output_line = f"{prefix}{seq_record[0]}\t{seq_record[1]}\n"
	else:
		output_line = f"{seq_record[0]}\t{seq_record[1]}\n"

	chromSizesoutput.write(output_line)

chromSizesoutput.close()
