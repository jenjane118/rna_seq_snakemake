configfile: "config.yaml"

rule all:
    input:
        expand(["files/{accession}_1.fastq.gz", "files/{accession}_2.fastq.gz"], accession=config["accession"]),
	    expand(["files/{accession}_sc.txt"], accession=config["accession"])

# rule get_fastq_pe:
#     output:
#         "data/{accession}_1.fastq",
#         "data/{accession}_2.fastq"
#     params:
#         # optional extra arguments
#         extra=""
#     threads: 6  # defaults to 6
#     wrapper:
#         "0.59.1/bio/sra-tools/fasterq-dump"


rule fasterq_pe:
    output:
        "files/{accession}_1.fastq",
        "files/{accession}_2.fastq"
    params:
        # optional extra arguments
        extra=""
    threads: 6  # defaults to 6
    shell:
        "fasterq-dump {wildcards.accession} -O files"

rule compress_files:
    input:
        "files/{accession}_{read}.fastq"
    output:
        "files/{accession}_{read}.fastq.gz"
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
