
# Clustering the Human Genome

This project will replicate the findings of John Novembre et al. as described in the article "Genes mirror geography within Europe". In this study, genetic variation of 3,000 European individuals, collected by the Population Reference Sample project, are investigated. Interestingly, despite relatively low levels of genetic differences amongst the sample, a strong correlation between genetic variation and geographic distances was found. In fact, clustering and plotting each of these individuals closely resembles the map of Europe as the name of the article suggests. However, given that worldwide genetic variation data is available thanks to the 1000 Genome Project, this project will aim to extend these boundaries to not only Europe but the world.

## Usage Instructions

* Description of targets and using `run.py`

## Description of Contents

The project consists of these portions:
```
PROJECT
├── README.md
├── config
│   ├── data-params.json
│   └── test-params.json
├── test
│   └── testdata
├── references
│   └── sample_pop.csv
├── notebooks
│   └── .gitkeep
├── requirements.txt
├── run.py
└── src
    ├── etl.py
    ├── process_data.py
    ├── read_data.py
```

### `src`

* `etl.py`: Library code that executes tasks useful for getting data.
* `process_data.py`: Library code that executes tasks for processing data
                     and generating chromosome cluster plot
* `read_data.py`: Optional library code useful for importing BAM, FASTQ,
                  and VCF files into a Pandas dataframe

### `config`

* `data-params.json`: Common parameters for getting data, serving as
  inputs to library code.
  
* `test-params.json`: parameters for running small process on small
  test data.

### `references`

* Data Dictionaries, references to external sources

### `notebooks`

* Jupyter notebooks for *analyses*
  - notebooks are not for data processing; they should import code
    from `src`.