## NOTE: This was created in case analysis of genetic data
##       is done in pandas. Other software like GATK may be
##       more efficient and effective.

"""  Data Reader

This script takes in VCF, BAM, and FASTQ files
and reads them to a pandas data frame.

"""

# Importing libraries
import pandas as pd
import numpy as np
import pysam
import io
from itertools import islice


def read_vcf(fp):
    """
    Reads VCF file
    
    :param fp: String representing file path to VCF file
    """
    
    # Read VCF line by line
    with open(fp, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
        
    # Convert lines to dataframe
    vcf_df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
        
    return vcf_df



def read_bam(fp):
    """
    Reads BAM file
    
    :param fp: String representing file path to BAM file
    """
    
    # Read file, convert to SAM
    imported = pysam.AlignmentFile(fp, mode = 'rb')

    # Store lines from SAM
    sams = []
    for i in imported:
        sams.append(i.to_dict())

    # Convert lines to dataframe
    lines = {}
    for k in sams[0].keys():
          lines[k] = tuple(lines[k] for lines in sams)

    sam_df = pd.DataFrame(lines)
    
    return sam_df



def read_fastq(fp):
    """
    Reads FASTQ file
    
    :param fp: String representing file path to FASTQ file
    """
    
    # Lists to store data frame data
    identifiers, sequences = [], []
    separators, quality_scores = [], []

    # Reads every four lines, stores data
    # NOTE: 1 observation = 4 lines in FASTQ
    with open(fp, "r") as f:
        while True:
            next_4_lines = list(islice(f, 4))

            if not next_4_lines:
                break
            else:
                identifiers.append(next_4_lines[0].replace('\n', ''))
                sequences.append(next_4_lines[1].replace('\n', ''))
                separators.append(next_4_lines[2].replace('\n', ''))
                quality_scores.append(next_4_lines[3].replace('\n', ''))

    fastq_dict = {'identifier': identifiers, 
                  'sequence': sequences, 
                  'separator': separators, 
                  'quality_score':quality_scores}

    fastq_df = pd.DataFrame(fastq_dict)
    
    return fastq_df
