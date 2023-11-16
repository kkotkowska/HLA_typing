# thesis-HLA_typing

## Introduction
Pipeline for HLA typing based on Nanopore sequencing data.

## Pipeline summary
1. Demultiplex based on primers
2. Index database files ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
3. Map against the database ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
4. Map the best hits again ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
5. Gather results into a csv

## Quick start
1. Install [Docker](https://docs.docker.com/engine/installation/), or [Anaconda](https://conda.io/miniconda.html)
    For docker - activate the deamon, and build an image
    ```bash
    docker build -t hla-typing .
    ```
    For conda - make a virtual environment 
    ```bash
    conda env create -f hla-env.yml
    ```
2. Download the pipeline
3. Start running your analysis!
