# -*- coding: utf-8 -*-

import json
import os
from flask import render_template

from . import app

TMP_DIRECTORY_PATH = "/Users/kytolav/projects/parallel_variants/tmp"


@app.route("/", methods=['GET', 'POST'])
def index():
    return render_template('app.html')


@app.route("/table", methods=['POST'])
def post_table_data():
    print('Table query received')
    # Listing the TMP directory
    files = [x for x in os.listdir(TMP_DIRECTORY_PATH) if os.path.isdir(os.path.join(TMP_DIRECTORY_PATH, x))]
    print(files)
    # Creating dict
    response = []
    for i, file in enumerate(files):
        # Counting BAM files
        bam_count = len([x for x in os.listdir(os.path.join(TMP_DIRECTORY_PATH, file)) if x.endswith('bam')])
        mpileup_count = len([x for x in os.listdir(os.path.join(TMP_DIRECTORY_PATH, file)) if x.endswith('mpileup')])
        vcf_count = len([x for x in os.listdir(os.path.join(TMP_DIRECTORY_PATH, file)) if x.endswith('vcf')])
        response.append([i+1, file, bam_count, mpileup_count, vcf_count])

    return json.dumps(response)
