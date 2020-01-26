"""  Data Ingestion

This script retrieves genetic data from the 1000 Genomes
Project through the International Genome Sample Resource's
FTP server. 

"""

# Importing libraries
import pandas as pd
import numpy as np
import shutil
import pysam
import json
import gzip
import io
import os
from itertools import islice
from ftplib import FTP


def get_data(vcf, bam, fastq, **kwargs):
    """
    Main function
    
    Downloads genetic data based on configuration file content.
    
    :param vcf: List containing chromosome values for VCF files
    :param bam: Dictionary containing sample-chromosome pairs
                for BAM files
    :param fastq: List containing sample for FASTQ files
    
    >>> cfg = json.load(open('data-params.json'))
    >>> get_data(**cfg) is None
    True
    """
    
    def run_query(func, data):
        
        if type(data) == dict:
            for key in data.keys():
                func(key, data[key])
        else:
            for arg in data:
                func(arg)
    
    # Download VCF
    run_query(get_vcf, vcf)
    
    # Download BAM
    run_query(get_bam, dict(bam))
        
    # Download FASTQ
    run_query(get_fastq, fastq)
    
    
    
def get_vcf(chr_num):
    """
    Downloads VCF files of a specific chromosome. Chromosomes 
    1-22 are available to download.
    
    :param chr_num: Integer or string representing chromosome to download
    
    >>> get_vcf(21) is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/release/20130502/'
    fname = ('ALL.chr{}.phase3').format(str(chr_num))
        
    # Downloading file in ftp directory
    fnames = download_data(path, 'vcf', fname)
    
    # Unzip file
    unzipper(fnames)
    

    
def get_bam(sample, chr_num):
    """
    Downloads BAM files for a specific individual. Chromosomes
    11 and 20 are available.
    
    :param sample: String identifying individual
    :param chr_num: Integer representing chromosome to download
    
    >>> get_bam("HG00096", 20) is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/phase3/data/{}/alignment/'.format(sample)
    fname = (('{}.chrom{}.ILLUMINA.bwa')
             .format(sample, str(chr_num)))
    
    # Downloading file
    download_data(path, 'bam', fname)
    
    
    
def get_fastq(sample):
    """
    Downloads all FASTQ files for a specific individual.
    
    :param sample: String identifying individual
    
    >>> get_fastq("HG01921") is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/phase3/data/{}/sequence_read/'.format(sample)
    
    # Downloading all sequences in ftp directory to local directory
    fnames = download_data(path, 'fastq')
    
    # Unzip files
    unzipper(fnames)
    


def download_data(path, ftype, fname=None):
    """
    Downloads data at a specifc path within the International 
    Genome Sample Resource's FTP server
    
    :param path: String representing path to directory housing
                 files for download
    :param ftype: String representing file type
    :param fname: List containing file names to download
    """
    
    # Connect to ftp server
    ftp_srv = 'ftp.1000genomes.ebi.ac.uk'
    ftp = FTP(ftp_srv)
    ftp.login()
    
    # Changing directory within ftp
    ftp.cwd(path)
    
    # Searches for files in directory
    all_fnames = []
    ftp.retrlines('NLST', callback=lambda x: all_fnames.append(str(x))) 
    
    if ftype == 'vcf':
        fnames = [f for f in all_fnames if (fname in f) and ('tbi' not in f)]
    
    elif ftype == 'bam':
        fnames = [f for f in all_fnames
                  if (fname in f) and ('bai' not in f) and ('bas' not in f)]
    
    elif ftype == 'fastq':
        fnames = all_fnames
        
    # Creates data directory if it doesnt exist
    if not os.path.exists('data'):
        os.mkdir('data')
    
    # Writes out data to data directory
    for f in fnames:
        file = open('data/'+f, "wb")
        resp = ftp.retrbinary("RETR "+path+f, file.write)
        file.close()
    
    # Disconnects from ftp
    ftp.close()
        
    return fnames

    
    
def unzipper(fnames):
    """
    Unzips files
    
    :param fnames: List containing files to unzip
    """
    
    # Unzipping files in local directory
    for fname in fnames:
        with gzip.open('data/'+fname,'rb') as zipped, open('data/'+fname[:-3], 'wb') as unzipped:
            shutil.copyfileobj(zipped, unzipped)
            
        # Removed old zipped file
        os.remove('data/'+fname) 
              