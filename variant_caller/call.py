#!/usr/bin/python
# -*- coding: utf-8 -*-


import argparse
import uuid
import os


VARIANT_OUTPUT_DIRECTORY_BASENAME = 'variants'
TEMPORARY_FILES_DIRECTORY_BASENAME = 'tmp'
NUMBER_OF_CORES = 20
NUMBER_OF_CHUNKS = 100


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
        bam_chunk_list = split_bam(bam_file, tmp_path)

        # Running variant calling in parallel using simple unix parallel
        call_variants_parallel(tmp_path)

        # Gathering results into single VCF
        gather_results(bam_file, tmp_path, VARIANT_OUTPUT_DIRECTORY_BASENAME)



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
