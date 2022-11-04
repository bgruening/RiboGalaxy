import gzip
from sys import argv, exit
from itertools import zip_longest
from Bio import SeqIO

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def process_fastq_record(record, output, UMI_5_prime_length=2, UMI_3_prime_length=5):
    '''
    Write UMI to FASTQ header for given biopython record to given open output
     UMI_5_prime_length number of bases to trim from 5' and write to header
     UMI_3_prime_length number of bases to trim from 3' and write to header
     defaults are for McGlincy Ingolia Protocol
    '''
    lines = record.format('fastq').split('\n') # all 4 lines 
    header = lines[0]   
    seq = lines[1]
    sep = lines[2]
    qual = lines[3]

    if (header.startswith('@')):
        trimmed_seq = seq[UMI_5_prime_length:-(UMI_3_prime_length + 1)] # fooprint + barcode
        UMI = seq[0:UMI_5_prime_length]+seq.rstrip()[-UMI_3_prime_length:len(seq)] #7nt in total; 5'NN and last 3'NNNNN
        split_header = header.split(" ")
        new_header = split_header[0]+"_"+UMI+" "+split_header[1]
        new_qual = qual[UMI_5_prime_length:-(UMI_3_prime_length + 1)]
        output.write(new_header+'\n')
        output.write(trimmed_seq+'\n')
        output.write(sep+'\n')
        output.write(new_qual+'\n')   



def UMI_processing(pathToFastaFile, output_path, bool_gzip, UMI_5_prime_length, UMI_3_prime_length):
    
    if bool_gzip:
        output = gzip.open(output_path,"wt")
    else: 
        output = open(output_path, 'w')
    
    if is_gz_file(pathToFastaFile) == True: 
        print ('file is gzipped fastq')
        with gzip.open(pathToFastaFile, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                process_fastq_record(record, output, UMI_5_prime_length, UMI_3_prime_length)

    else: 
        for record in SeqIO.parse(pathToFastaFile, 'fastq'):
            process_fastq_record(record, output, UMI_5_prime_length, UMI_3_prime_length)
            
    output.close()
            

    

def main():
    if len(argv) != 6:
        exit("Usage: 3 arguments required\n1: Path to fasta file \n2: name of output file\n3: string 'True' or 'False' whether to gzip output")

    # Get paths
    pathToFastaFile = argv[1]
    output = argv[2]
    bool_gzip = True if argv[3].lower().capitalize() == 'True' else False
    UMI_5_prime_length = int(argv[4])
    UMI_3_prime_length = int(argv[5])
    UMI_processing(pathToFastaFile, output, bool_gzip, UMI_5_prime_length, UMI_3_prime_length)

if __name__ == "__main__":
    main()
