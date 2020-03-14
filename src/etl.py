"""  Data Ingestion

etl.py retrieves genetic data from the 1000 Genomes
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

    
    
def get_vcf(chr_num, outdir):
    """
    Downloads VCF files of a specific chromosome. Chromosomes 
    1-22 are available to download.
    
    :param chr_num: Integer or string representing chromosome to download
    :param outdir: Directory to write out data
    
    >>> get_vcf(21) is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/release/20130502/'
    fname = ('ALL.chr{}.phase3').format(str(chr_num))
        
    # Downloading file in ftp directory
    fnames = download_data(path, 'vcf', outdir, fname)
    
    # Unzip file
    unzipper(fnames, otudir)
    
    return
    

    
def get_bam(sample, chr_num, outdir):
    """
    Downloads BAM files for a specific individual. Chromosomes
    11 and 20 are available.
    
    :param sample: String identifying individual
    :param chr_num: Integer representing chromosome to download
    :param outdir: Directory to write out data
    
    >>> get_bam("HG00096", 20) is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/phase3/data/{}/alignment/'.format(sample)
    fname = (('{}.chrom{}.ILLUMINA.bwa')
             .format(sample, str(chr_num)))
    
    # Downloading file
    download_data(path, 'bam', outdir, fname)
    
    return
    
    
    
def get_fastq(sample, outdir):
    """
    Downloads all FASTQ files for a specific individual.
    
    :param sample: String identifying individual
    :param outdir: Directory to write out data
    
    >>> get_fastq("HG01921") is None
    True
    """
    
    # Creating ftp path to file
    path = '/vol1/ftp/phase3/data/{}/sequence_read/'.format(sample)
    
    # Downloading all sequences in ftp directory to local directory
    fnames = download_data(path, 'fastq', outdir)
    
    # Unzip files
    unzipper(fnames, outdir)
    
    return
    


def download_data(path, ftype, outdir, fname=None):
    """
    Downloads data at a specifc path within the International 
    Genome Sample Resource's FTP server
    
    :param path: String representing path to directory housing
                 files for download
    :param ftype: String representing file type
    :param fname: List containing file names to download
    :param outdir: Directory to write out data
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
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Writes out data to data directory
    for f in fnames:
        file = open(outdir+f, "wb")
        resp = ftp.retrbinary("RETR "+path+f, file.write)
        file.close()
    
    # Disconnects from ftp
    ftp.close()
        
    return fnames

    
    
def unzipper(fnames, path):
    """
    Unzips files
    
    :param fnames: List containing files to unzip
    :param path: Directory of zipped file
    """
    
    # Unzipping files in local directory
    for fname in fnames:
        with gzip.open(path+fname,'rb') as zipped, open(path+fname[:-3], 'wb') as unzipped:
            shutil.copyfileobj(zipped, unzipped)
            
        # Removed old zipped file
        os.remove(path+fname) 
        
    return
        
        
# ---------------------------------------------------------------------
# Driver Function
# ---------------------------------------------------------------------      
        
        
def get_data(vcf, bam, fastq, outdir, **kwargs):
    """
    Downloads genetic data based on configuration file content.
    
    :param vcf: List containing chromosome values for VCF files
    :param bam: Dictionary containing sample-chromosome pairs
                for BAM files
    :param fastq: List containing sample for FASTQ files
    :param otudir: Directory to write out files
    
    >>> cfg = json.load(open('data-params.json'))
    >>> get_data(**cfg) is None
    True
    """
    
    # Checks if inpath passed in. If so, nothing will
    # be done as only existing test data is needed
    if len(kwargs) > 0:
        return
        
    else:
        def run_query(func, data):

            if type(data) == dict:
                for key in data.keys():
                    func(key, data[key], outdir)
            else:
                for arg in data:
                    func(arg, outdir)

        # Download VCF
        run_query(get_vcf, vcf)

        # Download BAM
        run_query(get_bam, dict(bam))

        # Download FASTQ
        run_query(get_fastq, fastq)
    
    return
              