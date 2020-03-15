
# Clustering the Human Genome

This project will replicate the findings of John Novembre et al. as described in the article "Genes mirror geography within Europe". In this study, genetic variation of 3,000 European individuals, collected by the Population Reference Sample project, are investigated. Interestingly, despite relatively low levels of genetic differences amongst the sample, a strong correlation between genetic variation and geographic distances was found. In fact, clustering and plotting each of these individuals closely resembles the map of Europe as the name of the article suggests. However, given that worldwide genetic variation data is available thanks to the 1000 Genome Project, this project will aim to extend these boundaries to not only Europe but the world.

## Usage Instructions

In order to use the different components of this project, please run `python run.py` along with a target of your choice:

* `clean`: Cleans directory after project run
* `convert`: Converts data files according to convert-params.json
* `data`: Retrieves data files according to data-params.json 'data' key
* `data-test`: Retrieves test file for test run
* `process`: Processes and produces output files
* `test-project`: Tests project – shortcut to running `python run.py data-test process`

## Description of Contents

The project consists of these portions:
```
PROJECT
├── README.md
├── config
│   ├── convert-params.json
│   ├── data-params.json
│   ├── env.json
│   └── test-params.json
├── test
│   └── testdata
├── references
│   └── sample_pop.csv
├── notebooks
│   ├── notebook-resources
│   └── ProjectReport.ipynb
├── requirements.txt
├── run.py
└── src
    ├── conversion.py
    ├── conversion.sh
    ├── etl.py
    ├── process_data.py
    ├── process_data.sh
    └── read_data.py
```

### `src`

* `conversion.py`: Library code to convert between FASTQ, BAM, and FASTQ files.
* `conversion.sh`: Shell script to store BWA, GATK, and SAMTools commands used
                   in 'conversion.py'.
* `etl.py`: Library code that executes tasks useful for getting data from 
            1000 Genomes FTP.
* `process_data.py`: Library code that executes tasks for processing data
                     and generating chromosome cluster plot.
* `process_data.sh`: Shell script to store PLINK2 and VCF merging commands used
                     in 'process_data.py'.
* `read_data.py`: Optional library code to transform BAM, FASTQ,
                  and VCF files into a Pandas dataframe.

### `config`

* `data-params.json`: Common parameters for getting data, serving as
                      inputs to library code.
* `test-params.json`: Parameters for running small process on small
                      test data.
* `convert-parama.json`: Parameters for converting converting between FASTQ,
                         BAM, and VCF files.
* `env.json`: Configuration file containing DockerHub path and ouput filepaths.

### `references`

* `sample_pop.csv`: File containing sample IDs and their corresponding population.
                    Used for population mapping in output plot.

### `notebooks`

* `ProjectReport.ipynb`: Report containing the description, workflow, and findings of this project.

### `test`

* `test`: Files for testing project.