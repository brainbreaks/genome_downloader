=================
genome_downloader
=================

Introduction
============

This package can be used to download commonly used reference genomes and annotations associated with
them. Currently **hg19**, **hg19** and **hg19** genomes are supported

Downloading a genome
====================

    python -m genome_downloader <genome> <destination>

Genome directory structure
==========================
Each genome is downloaded into <destination>/<genome> folder (e.g. ``python -m genome_downloader mm10 genomes``
will download data into **genomes/mm10** folder)

.. code-block::

    <destination>
    └── <genome>
        │   annotation
        │   ├── ChromInfo.txt               - Chromosome sizes
        │   ├── cytoBand.txt                - Cytoband information
        │   ├── <genome>.chrom.sizes        - Chromosome sizes
        │   ├── <genome>.ncbiRefSeq.gtf.gz  - Genome annotations with putative transcripts
        │   ├── <genome>.refGene.gtf.gz     - Curated genome annotations
        │   └── ucsc_repeatmasker.tsv       - Repeatmasker's repeats coordinates
        ├── <genome>.fa                     - Reference genome
        ├── <genome>.fai                    - Reference genome index
        ├── <genome>.*.bt2                  - Bowtie2 index files
        └── <genome>.rev.*.bt2              - Bowtie2 index files
