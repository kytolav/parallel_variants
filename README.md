# parallel_variants

This repository contains a demo of a parallel variant calling suite using VarScan as the backbone.

## Setup
First, I recommend to use Anaconda to get a stable python environment. Something like
~~~~
conda create -n <name_of_your_environment> python=3.5 pandas numpy scipy flask
~~~~~
should get you started. Activate the environment by typing 
~~~~
source activate <name_of_your_environment>
~~~~~
in terminal.

## Testing variant calling
Clone the repository to a location of your liking using `git clone`. For the testing you are going to need a reference genome and sorted, indexed BAM files. In addition, the following command line utilities must be in your PATH:
* samtools
* parallel
* java
In addition, VarScan2 must be installed. Note that currently only human reference genome hg19 is supported as the chromosome sizes are hard coded into the app. This will be changed in further versions. Before running, edit the constants in the beginning of file "call.py" to fit your environment.

Run the pipeline by typing
~~~~
python <path_to_repo>/variant_caller/call.py --bam_list <list_of_bam_files> --chunk_length 10000000 --n_cores 2
~~~~
By default, this runs the variant calling using VarScan for input BAM files after splitting them into smaller chunks. The results are collected into directory `variants` while all temporary files are collected to `tmp`.

## Tracking execution with web app
The repository contains a simple Flask web app to track the execution of variant calling. Run the web app (in conda environment) by executing script `run_server.sh`, e.g. by
~~~~
run_server.sh 8008
~~~~
This will start a local HTTP server and serves the app at [link](http://localhost:8008 "http://localhost:8008"). Note that the web app by default tracks changes inside the static folder `tmp` and this hard coded folder must be changed if another folder is used for intermediary files.