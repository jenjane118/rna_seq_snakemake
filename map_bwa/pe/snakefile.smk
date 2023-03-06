#snakemake snakefile map with bwa-mem, sort with samtools, mapping stats with samtools flagstat

configfile: "config.yaml"

rule all:
    input:
        expand("sorted_reads/flagstat_{sample}.txt", sample=config["samples"]),
	expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])    

#def get_bwa_map_input_fastqs(wildcards):
#    return config["samples"][wildcards.sample]

rule bwa_map_pe:
    input: 
        fa=config["genomeFile"],
        reads=["trimmed/trimmed_{sample}_1.fastq.gz", "trimmed/trimmed_{sample}_2.fastq.gz"]
    output: temp("mapped_reads/{sample}.bam")
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'", #this is read-group header if needed
    threads: 8
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    log:
        "logs/bwa_mem/{sample}.log"
    shell: 
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    log:
        "logs/samtools/{sample}.log"
    shell:
        "(samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}) 2> {log}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_flagstat:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/flagstat_{sample}.txt"
    shell:
        "samtools flagstat {input} > {output}"
