#!/usr/bin/python
# -*- coding: utf-8 -*-


import argparse
import uuid
import os
import math


VARIANT_OUTPUT_DIRECTORY_BASENAME = 'variants'
TEMPORARY_FILES_DIRECTORY_BASENAME = 'tmp'
NUMBER_OF_CORES = 2
CHUNK_LENGTH = 10000000
PATH_TO_REFERENCE_GENOME_HG19 = '/data/csb/organisms/homo_sapiens/hg19.fa'
PATH_TO_VARSCAN = '/data/csb/tools/VarScan.v2.3.9.jar'
HG19_CHROMOSOME_SIZES = {
    'chr1':249250621,
    'chr2':243199373,
    'chr3':198022430,
    'chr4':191154276,
    'chr5':180915260,
    'chr6':171115067,
    'chr7':159138663,
    'chrX':155270560,
    'chr8':146364022,
    'chr9':141213431,
    'chr10':135534747,
    'chr11':135006516,
    'chr12':133851895,
    'chr13':115169878,
    'chr14':107349540,
    'chr15':102531392,
    'chr16':90354753,
    'chr17':81195210,
    'chr18':78077248,
    'chr20':63025520,
    'chrY':59373566,
    'chr19':59128983,
    'chr22':51304566,
    'chr21':48129895,
    'chrMT':16571}


def run_analysis(bam_list, out_dir = VARIANT_OUTPUT_DIRECTORY_BASENAME):
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
        bam_chunk_list = split_bam(bam_file, tmp_path, HG19_CHROMOSOME_SIZES, CHUNK_LENGTH)

        # Running variant calling in parallel using simple unix parallel
        call_variants_parallel(bam_chunk_list, tmp_path)

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


def call_variants_parallel(bam_chunk_list, tmp_path):
    """Calls variants using a parallel VarScan launcher"""
    # Ignoring empty files
    nonempty_chunks = [x for x in bam_chunk_list if not os.stat(os.path.join(tmp_path, x)).st_size == 0]
    cmd_parallel = "parallel -j %d" % NUMBER_OF_CORES
    cmd_call = ' "samtools mpileup -B -f %s {} > {}.mpileup"' % (PATH_TO_REFERENCE_GENOME_HG19)
    cmd_files = " ::: " + " ".join([os.path.join(tmp_path, os.path.basename(x)) for x in nonempty_chunks])
    cmd = cmd_parallel + cmd_call + cmd_files
    os.system(cmd)
    # Listing mpileup files and excluding empty ones
    mpileup_files = [os.path.join(tmp_path, x) for x in os.listdir(tmp_path) if x.endswith('mpileup')]
    nonempty_mpileups = [x for x in mpileup_files if not os.stat(x).st_size == 0]
    cmd_files = " ::: " + " ".join(nonempty_mpileups)
    cmd_call = ' "java -jar %s mpileup2snp {} --output-vcf 1 > {}.vcf"' % (PATH_TO_VARSCAN)
    cmd = cmd_parallel + cmd_call + cmd_files
    os.system(cmd)


def split_bam(bam_file, output_dir, chromosome_sizes, chunk_length):
    """Split the BAM file into chunks of specified length. The input BAM file must be sorted."""
    for key, value in chromosome_sizes.items():
        # Ranges
        ranges_start = range(0, value, chunk_length)
        for start in ranges_start:
            cmd = 'samtools view -h -b %s "%s:%d-%d" > %s/%s_%d_%d.bam' % (bam_file, key, start, start + chunk_length,
                                                                           output_dir, key, start, start + chunk_length)
            os.system(cmd)
    # Listing files
    return [x for x in os.listdir(output_dir) if 'chr' in x]


    # Saving the BAM header
    cmd = "samtools view -H %s > %s/header.txt" % (bam_file, output_dir)
    os.system(cmd)
    # Counting lines of the file
    line_count = bam_line_count(bam_file)
    lines_per_chunk = math.ceil(line_count * 1.0 / number_of_chunks)
    # Splitting the file with samtools and unix tools
    cmd = "samtools view %s |split -l %d - %s/chunk_" % (bam_file, lines_per_chunk, output_dir)
    os.system(cmd)
    # Appending header to each chunk
    chunks = [os.path.realpath(x) for x in os.listdir(output_dir) if 'header' not in x]
    for chunk in chunks:
        chunk_base = os.path.basename(chunk)
        cmd = "cat %s/header.txt %s/%s > %s/%s_tmp" % (output_dir, output_dir, chunk_base,
                                                                  output_dir, chunk_base)
        os.system(cmd)
        # Renaming
        cmd = "mv %s/%s_tmp %s/%s" % (output_dir, chunk_base, output_dir, chunk_base)
        os.system(cmd)
    return chunks


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
    os.system("mkdir %s " % os.path.join(path_basename, new_path))
    return new_path


def get_args():
    """Returns command line arguments as dictionary"""
    parser = argparse.ArgumentParser(
        description='Master script for QC report generation.')
    parser.add_argument('--bam_list', help='Space separated paths to FASTQ '
                        'input files. At least one is required.', nargs='+',
                        required=True)
    return vars(parser.parse_args())


def main():
    run_analysis(**get_args())


if __name__ == "__main__":
    main()
