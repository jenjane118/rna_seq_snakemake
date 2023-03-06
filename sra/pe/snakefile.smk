configfile: "config.yaml"

rule all:
    input:
        expand(["files/{accession}_1.fastq.gz", "files/{accession}_2.fastq.gz"], accession=config["accession"]),
	"files/sc.txt"

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
    output:
        "files/sc.txt"
    shell:
        "for file in `ls files/*.fastq.gz`; do zcat $file | wc -l; done > {output}"

    
