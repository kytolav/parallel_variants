#!/usr/bin/python
# -*- coding: utf-8 -*-


import argparse
import uuid
import os


VARIANT_OUTPUT_DIRECTORY_BASENAME = 'variants'
TEMPORARY_FILES_DIRECTORY_BASENAME = 'tmp'
NUMBER_OF_CORES = 2
CHUNK_LENGTH = 10000000
PATH_TO_REFERENCE_GENOME_HG19 = '/Users/kytolav/projects/parallel_variants/references/hg19.fa'
PATH_TO_VARSCAN = '/Users/kytolav/tools/VarScan.v2.3.9.jar'
HG19_CHROMOSOME_SIZES = {
    'chr1': 249250621,
    'chr2': 243199373,
    'chr3': 198022430,
    'chr4': 191154276,
    'chr5': 180915260,
    'chr6': 171115067,
    'chr7': 159138663,
    'chrX': 155270560,
    'chr8': 146364022,
    'chr9': 141213431,
    'chr10': 135534747,
    'chr11': 135006516,
    'chr12': 133851895,
    'chr13': 115169878,
    'chr14': 107349540,
    'chr15': 102531392,
    'chr16': 90354753,
    'chr17': 81195210,
    'chr18': 78077248,
    'chr20': 63025520,
    'chrY': 59373566,
    'chr19': 59128983,
    'chr22': 51304566,
    'chr21': 48129895,
    'chrMT': 16571}


def run_analysis(bam_list, chunk_length=CHUNK_LENGTH, n_cores=NUMBER_OF_CORES,
                 out_dir=VARIANT_OUTPUT_DIRECTORY_BASENAME):
    """Launch parallel variant calling for all input files."""
    # Check for required arguments
    if not bam_list:
        raise (Exception("fastq_list parameter cannot be empty."))
    if not out_dir:
        raise (Exception("out_dir parameter cannot be empty."))

    # Looping through all files
    for bam_file in bam_list:

        # Creating tmp directory
        tmp_path = create_tmp_directory(TEMPORARY_FILES_DIRECTORY_BASENAME)

        # Splitting the bam
        bam_chunk_list = split_bam(bam_file, tmp_path, HG19_CHROMOSOME_SIZES, chunk_length)

        # Running variant calling in parallel using simple unix parallel
        call_variants_parallel(bam_chunk_list, tmp_path, n_cores)

        # Gathering results into single VCF
        gather_results(bam_file, tmp_path, VARIANT_OUTPUT_DIRECTORY_BASENAME)


def gather_results(bam_file, tmp_path, output_directory):
    """Gathering separate VCF files into one file the results.
    Note that sorting is currently not implemented."""
    vcf_files = [os.path.join(tmp_path, x) for x in os.listdir(tmp_path) if x.endswith('vcf')]
    output_name = os.path.join(output_directory, os.path.splitext(os.path.basename(bam_file))[0]) + '.vcf'
    cmd = """grep '#' %s > %s; for file in %s; do grep -v '#' $file >> %s; done""" % (vcf_files[0],
                                                                                      output_name,
                                                                                      " ".join(vcf_files),
                                                                                      output_name)
    os.system(cmd)


def call_variants_parallel(bam_chunk_list, tmp_path, n_cores):
    """Calls variants using a parallel VarScan launcher"""
    # Ignoring empty files
    print('Creating mpileup files')
    nonempty_chunks = [x for x in bam_chunk_list if not os.stat(os.path.join(tmp_path, x)).st_size == 0]
    cmd_parallel = "parallel -j %d" % n_cores
    cmd_call = ' "samtools mpileup -B -f %s {} > {}.mpileup"' % PATH_TO_REFERENCE_GENOME_HG19
    cmd_files = " ::: " + " ".join([os.path.join(tmp_path, os.path.basename(x)) for x in nonempty_chunks])
    cmd = cmd_parallel + cmd_call + cmd_files
    os.system(cmd)
    # Listing mpileup files and excluding empty ones
    print('Calling variants')
    mpileup_files = [os.path.join(tmp_path, x) for x in os.listdir(tmp_path) if x.endswith('mpileup')]
    nonempty_mpileups = [x for x in mpileup_files if not os.stat(x).st_size == 0]
    cmd_files = " ::: " + " ".join(nonempty_mpileups)
    cmd_call = ' "java -jar %s mpileup2snp {} --output-vcf 1 > {}.vcf"' % PATH_TO_VARSCAN
    cmd = cmd_parallel + cmd_call + cmd_files
    os.system(cmd)


def split_bam(bam_file, output_dir, chromosome_sizes, chunk_length):
    """Split the BAM file into chunks of specified length. The input BAM file must be sorted."""
    print('Splitting BAM file: %s' % bam_file)
    for key, value in chromosome_sizes.items():
        # Ranges
        ranges_start = range(0, value, chunk_length)
        for start in ranges_start:
            cmd = 'samtools view -h -b %s "%s:%d-%d" > %s/%s_%d_%d.bam' % (bam_file, key, start, start + chunk_length,
                                                                           output_dir, key, start, start + chunk_length)
            os.system(cmd)
    # Listing files
    return [x for x in os.listdir(output_dir) if 'chr' in x]


def bam_line_count(bam_file):
    """Counting lines in a BAM file"""
    cmd = "samtools view %s |wc -l > tmp_line_count.txt" % bam_file
    os.system(cmd)
    with open('tmp_line_count.txt') as f:
        count = "".join(f.readlines()[0].split())
    return int(count)


def create_tmp_directory(path_basename):
    """Creating a tmp directory"""
    random_name = str(uuid.uuid4())

    # Checking that the directory does not exist
    while os.path.exists(os.path.join(path_basename, random_name)):
        # If name exists, generating a new one
        random_name = str(uuid.uuid4())

    # Creating the directory
    new_path = os.path.join(path_basename, random_name)
    print('Created tmp directory: %s' % new_path)
    os.system("mkdir %s " % new_path)
    return new_path


def get_args():
    """Returns command line arguments as dictionary"""
    parser = argparse.ArgumentParser(
        description='Master script for QC report generation.')
    parser.add_argument('--bam_list', help='Space separated paths to FASTQ '
                        'input files. At least one is required.', nargs='+',
                        required=True)
    parser.add_argument('--chunk_length', help=("Length (bp) of the genomic sequence that is used to split"
                                                " input BAM files to chunks. Default: 10,000,000bp."),
                        required=False, default=10000000, type=int)
    parser.add_argument('--n_cores', help='Number of processor cores to use for processing. Default: 2',
                        required=False, default=2, type=int)
    return vars(parser.parse_args())


def main():
    run_analysis(**get_args())


if __name__ == "__main__":
    main()
