#run fastqc and multiqc on paired end fastq files

configfile: "config.yaml"

rule all:
    input:
        expand("fastqc_outputs/{sample}_{read}_fastqc.html", sample=config["samples"], read=[1,2])
 
rule run_fastqc: 
    input:
        "trimmed/trimmed_{sample}_{read}.fastq.gz"
    output:
        html="fastqc_outputs/{sample}_{read}_fastqc.html",
        zip_files="fastqc_outputs/{sample}_{read}_fastqc.zip"
        #done=touch("fastqc_outputs/{sample}_{read}_done.txt")  #this made when all done
    threads: 8
    shell: 
        "fastqc {input} -o fastqc_outputs"
    #wrapper: "0.31.1/bio/fastqc"

#rule run_multiqc:
#    input:
#        expand("fastqc_outputs/{sample}_{read}_fastqc.html", sample=config["samples"], read=[1,2])
        #"run_fastq_done.txt" #requires run_fastqc is all done first?
#    output:
#        "fastqc_outputs/multiqc_report.html"
#    shell:
#        "multiqc -f {input}"
    #wrapper: "0.31.1/bio/multiqc"
