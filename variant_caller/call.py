#!/usr/bin/python
# -*- coding: utf-8 -*-


import argparse
import uuid
import os
import math


VARIANT_OUTPUT_DIRECTORY_BASENAME = 'variants'
TEMPORARY_FILES_DIRECTORY_BASENAME = 'tmp'
NUMBER_OF_CORES = 2
NUMBER_OF_CHUNKS = 10
PATH_TO_REFERENCE_GENOME_HG19 = '/data/csb/organisms/homo_sapiens/hg19.fa'
PATH_TO_VARSCAN = '/data/csb/tools/VarScan.v2.3.9.jar'


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
        bam_chunk_list = split_bam(bam_file, tmp_path, NUMBER_OF_CHUNKS)

        # Running variant calling in parallel using simple unix parallel
        call_variants_parallel(bam_chunk_list, tmp_path)

        # Gathering results into single VCF
        gather_results(bam_file, tmp_path, VARIANT_OUTPUT_DIRECTORY_BASENAME)


def call_variants_parallel(bam_chunk_list, tmp_path):
    """Calls variants using a parallel VarScan launcher"""
    cmd_parallel = "parallel -j %d" % NUMBER_OF_CORES
    cmd_call = ' "samtools mpileup -B -f %s {} |java -jar %s mpileup2snp --output-vcf 1 --output-file {}.vcf"' % (PATH_TO_REFERENCE_GENOME_HG19, PATH_TO_VARSCAN)
    cmd_files = " ::: " + " ".join([os.path.join(tmp_path, os.path.basename(x)) for x in bam_chunk_list])
    cmd = cmd_parallel + cmd_call + cmd_files
    print(cmd)
    #os.system(cmd)


def split_bam(bam_file, output_dir, number_of_chunks):
    """Split the BAM file into a defined number of chunks. The input BAM file must be sorted."""
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
    new_path = os.path.exists(os.path.join(path_basename, random_name))
    os.system("mkdir %s " % new_path)
    return new_path


def get_args():
    """Returns command line arguments as dictionary"""
    parser = argparse.ArgumentParser(
        description='Master script for QC report generation.')
    parser.add_argument('--bam_list', help='Space separated paths to FASTQ '
                        'input files. At least one is required.', nargs='+',
                        required=True)
    parser.add_argument('--out_dir', help='Output root directory under which '
                        'new, results containing folder will be created.',
                        required=True)
    parser.add_argument('--subsample', help='Fraction of reads to be sampled '
                        '[0,1]. Adjust according to input file sizes. Affects '
                        'greatly QC analysis run time.',
                        default=reads.DEFAULT_SUBSAMPLE, type=float)
    parser.add_argument('--rscript_path', help='Path to Rscript.',
                        default=ToolPath.rscript)
    return vars(parser.parse_args())



def main():
    run_analysis(**get_args())


if __name__ == "__main__":
    main()
