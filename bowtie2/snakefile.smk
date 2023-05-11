#map eukaryotic genomes with splice aware aligner, bowtie2

configfile: "config.yaml"

rule all:
    input:
        expand("sorted_reads/flagstat_{sample_name}.txt", sample=config["samples"]),
	    expand("sorted_reads/{sample_name}.bam.bai", sample=config["samples"])

rule bowtie2_build:
    input:
        ref=config["genomeFile"],
        #name=config["genomeName"],
    output:
        multiext(
            "config["genomeName"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra="config["genomeName"],  # optional parameters
    threads: 8,
    shell:
        "bowtie2-build --threads {threads} {input.ref} {params.extra} 2> {log}"

rule bowtie2:
    input:
        #read1="trimmed/trimmed_{sample}_1.fastq.gz", 
        #read2="trimmed/trimmed_{sample}_2.fastq.gz",
        # requires reads to be in .1.fastq format
        sample=["trimmed/trimmed_{sample_name}.1.fastq", "trimmed/trimmed_{sample_name}.2.fastq"],
        idx=multiext(
            "index/" + config["genomeName"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"
        ),
        #name=config["genomeName"],
    output:
        "mapped_reads/{sample_name}.bam",
    log:
        "logs/bowtie2/{sample_name}.log",
    params:
        extra="--fast-local --no-unal",  # optional parameters
    threads: 8  # Use at least two threads
    #shell: 
    #    "(bowtie2 --threads {threads} -x {input.name} -1 {input.read1} -2 {input.read2} {params.extra} | "
    #    "samtools view -Sb - > {output}) 2> {log}"
    wrapper:
        "v1.25.0/bio/bowtie2/align"


rule samtools_sort:
    input:
        "mapped_reads/{sample_name}.bam"
    output:
        protected("sorted_reads/{sample_name}.bam")
    log:
        "logs/samtools_sort/{sample_name}.log"
    shell:
        "(samtools sort -T sorted_reads/{wildcards.sample_name} "
        "-O bam {input} > {output}) 2> {log}"

rule samtools_index:
    input:
        "sorted_reads/{sample_name}.bam"
    output:
        "sorted_reads/{sample_name}.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_flagstat:
    input:
        "sorted_reads/{sample_name}.bam"
    output:
        "sorted_reads/flagstat_{sample_name}.txt"
    shell:
        "samtools flagstat {input} > {output}"