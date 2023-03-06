__author__ = "Jennifer J Stiens"
__copyright__ = "Copyright 2022, Jennifer J Stiens"
__email__ = "j.j.stiens@gmail.com"
__license__ = "MIT"

#run fastq for filtering and trimming reads
configfile: "config.yaml"

rule all:
    input:
        expand(["trimmed/trimmed_{sample}_1.fastq.gz", "trimmed/trimmed_{sample}_2.fastq.gz"], sample=config["samples"])

rule fastp_pe:
    input:
        #sample=[config["PATH"] + "ncbi/files/{sample}_1.fastq.gz", config["PATH"] + "ncbi/files/{sample}_2.fastq.gz"]
        sample=["files/{sample}_1.fastq.gz", "files/{sample}_2.fastq.gz"]
    output:
        trimmed=["trimmed/trimmed_{sample}_1.fastq.gz", "trimmed/trimmed_{sample}_2.fastq.gz"],
        #unpaired in single file
        unpaired="trimmed/{sample}.singletons.fastq.gz",
        failed="trimmed/{sample}.failed.fastq.gz",
        html="trimmed/{sample}_fastp.html",
        json="trimmed/{sample}_fastp.json"
    log:
        "logs/fastp_trim/{sample}.log"
    threads: 2
    wrapper:
        "v1.14.1/bio/fastp"

