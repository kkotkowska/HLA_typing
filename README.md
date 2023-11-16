# thesis-HLA_typing

## Introduction
Pipeline for HLA typing based on Nanopore sequencing data.

## Pipeline summary
1. Demultiplex based on primers
2. Index database files (['KMA'](https://bitbucket.org/genomicepidemiology/kma/src/master/))
3. Map against the database (['KMA'](https://bitbucket.org/genomicepidemiology/kma/src/master/))
4. Map the best hits again (['KMA'](https://bitbucket.org/genomicepidemiology/kma/src/master/))
5. Gather results into a csv

## Quick start
1. Install ['Docker']((https://docs.docker.com/engine/installation/), or ['anaconda']((https://conda.io/miniconda.html)
2. Download the pipeline
3. Start running your analysis!
