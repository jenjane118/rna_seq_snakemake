#snakemake rules for calculating coverage .bigwig files from .bam files
#requires deeptools in conda env
#usage: snakemake --cores {num cores} -s bam_coverage.smk

configfile: "config.yaml"

rule all:
    input:
        expand("covg_bigwigs/{sample}_{dir}.bw", sample=config["samples"], dir=["fwd", "rev"])    

rule bam_coverage: 
    input:
        bamfile="sorted_reads/{sample}.bam"
    output:
        fwfile="covg_bigwigs/{sample}_fwd.bw",
        rvfile="covg_bigwigs/{sample}_rev.bw"
    log:
        "logs/bamCov/{sample}.log"
    threads: 4
    run: 
        shell("bamCoverage -b {input} -o {output.fwfile} -of bigwig --filterRNAstrand forward --normalizeUsing RPKM -p 8 --binSize 10 --extendReads > {log}")
        shell("bamCoverage -b {input} -o {output.rvfile} -of bigwig --filterRNAstrand reverse --normalizeUsing RPKM -p 8 --binSize 10 --extendReads > {log}")
