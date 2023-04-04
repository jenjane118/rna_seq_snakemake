__author__ = "Jennifer J Stiens"
__copyright__ = "Copyright 2022, Jennifer J Stiens"
__email__ = "j.j.stiens@gmail.com"
__license__ = "MIT"


#snakefile for downloading files
#use sra-tools to download fastq files from single-end data and sanity check report to make sure reads look ok

configfile: "config.yaml"

rule all:
    input:
        expand("files/{accession}.fastq.gz", accession=config["accession"]),
        expand(["files/{accession}_sc.txt"], accession=config["accession"])

rule fasterq_single:
    output:
        "files/{accession}.fastq"
    params:
        # optional extra arguments
        extra=""
    threads: 6  # defaults to 6
    shell:
        "fasterq-dump {wildcards.accession} -O files"

rule compress_files:
    input:
        "files/{accession}.fastq"
    output:
        "files/{accession}.fastq.gz"
    shell:
        "gzip {input}"

rule sanity_check:
    input:
        "files/{accession}_1.fastq.gz",
        "files/{accession}_2.fastq.gz"
    output:
        "files/{accession}_sc.txt"
    shell:
         "for file in {input}; do zcat $file | wc -l; done >> {output}"
