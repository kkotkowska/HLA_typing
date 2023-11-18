# HLA_typing

## Introduction
Pipeline for HLA typing based on Nanopore sequencing data.

## Pipeline summary
1. Demultiplex based on primers
2. Index database files ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
3. Map against the database ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
4. Map the best hits again ([KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/))
5. Gather results into a csv file

## Quick start
1. Install [Docker](https://docs.docker.com/engine/installation/), or [Anaconda](https://conda.io/miniconda.html).
2. Download the pipeline
    ```bash
    git clone git@github.com:kkotkowska/HLA_typing.git
    ```
    Download the database from [IMGT-HLA](https://www.ebi.ac.uk/ipd/imgt/hla/download/) and place in the data folder
    ```bash

    wget -r -np -nH --cut-dirs 3 ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/ && mv imgt data
    ```

    For docker - activate the deamon, and build an image
    ```bash
    docker build -t hla-typing .

    ```
    run with "-profile docker"

    For conda - make a virtual environment 
    ```bash
    conda env create -f hla-env.yml
    conda activate hla-env
    ```
    download KMA in the project directory. (Make sure you have a c compiler - for linux - gcc)
    ```bash
    sudo apt-get install gcc
    ```
    ```bash
    git clone https://bitbucket.org/genomicepidemiology/kma.git
    cd kma && make
    mv kma ../bin/
    ```
3. Start running your analysis!
    Run without parameters for a test.
    ```bash
    nextflow run workflow.nf 
    ```
