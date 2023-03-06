__author__ = "Jennifer J Stiens"
__copyright__ = "Copyright 2022, Jennifer J Stiens"
__email__ = "j.j.stiens@gmail.com"
__license__ = "MIT"

#run fastq for filtering and trimming reads for single end
# run from project directory

configfile: "config.yaml"

rule all:
    input:
        expand("trimmed/se/{sample}.fastq.gz", sample=config["single_samples"])

rule fastp_se:
    input:
        sample=config["PATH"] + "ncbi/files/se/{sample}.fastq.gz"
    output:
        trimmed="trimmed/se/{sample}.fastq.gz",
        failed="trimmed/se/{sample}.failed.fastq.gz",
        html="report/se/{sample}.html",
        json="report/se/{sample}.json"
    log:
        "logs/fastp/se/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA",
        extra=""
    threads: 1
    shell:
        "fastp -i {input.sample} -o {output.trimmed} "
        "--failed_out {output.failed} {params.adapters} "
	"-j {output.json} -h {output.html}"
    #wrapper:
    #    "v1.17.2/bio/fastp"

