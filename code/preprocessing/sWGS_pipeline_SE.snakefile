# Sample list
IDS, = glob_wildcards("{sample}_R1.fastq.gz")

# Import config and make report
configfile: "config.yaml"
report: "report/workflow.rst"

rule all:
    input:
        "multiqc_report.html"

rule mapping_bwa:
    input:
        fq = "{sample}_R1.fastq.gz"
    output:
        "mapped_reads/{sample}_R1.bam"
    threads: 8
    params:
        walltime="2:00:00",
        ppn = 10,
        library = config["library"],
        genome = config["humangenome"],
        outdir = "mapped_reads"
    shell:
        "ml purge; ml FastQC/0.11.5-Java-1.8.0_74 ; "
        "fastqc {input.fq} --outdir {params.outdir} ; "
        "ml purge; ml BWA/0.7.17-intel-2018a; ml SAMtools/1.8-intel-2018a; "
        "bwa mem -t {threads} {params.genome} {input.fq} | samtools sort -O bam -o {output}"

rule deduplicate:
    input:
        "mapped_reads/{sample}_R1.bam"
    output:
        bam = "nodups/{sample}_nodups.bam",
        metrics = "nodups/{sample}_dup_metrics.txt"
    params:
        walltime="0:45:00",
        ppn = 6,
    shell:
        "ml purge && ml load picard/2.1.1-Java-1.8.0_74; "
        "java -jar $EBROOTPICARD/picard.jar MarkDuplicates "
        "I={input} "
        "O={output.bam} "
        "M={output.metrics}  "
        "VALIDATION_STRINGENCY=SILENT "
        "REMOVE_DUPLICATES=true; "
        "ml purge && ml SAMtools/1.6-intel-2017b; "
        "samtools index {output.bam}"

rule convertnpz:
    input:
        "nodups/{sample}_nodups.bam"
    output:
        "convertnpz/{sample}.npz"
    params:
        walltime="0:20:00",
        ppn = 1,
        wiseX = config["wisex"],
        binsize = config["wisexkb"] * 1000,
        wiseXref = config["wisexref"]
    shell:
        "ml purge && ml R/3.4.4-intel-2018a-X11-20180131 && ml Python/3.6.4-intel-2018a; "
        "{params.wiseX} convert --binsize {params.binsize} {input} {output}"

rule plotnpz:
    input:
        "convertnpz/{sample}.npz"
    output:
        bed = "plotnpz/{sample}_bins.bed"
    params:
        walltime="0:20:00",
        ppn = 1,
        wiseX = config["wisex"],
        binsize = config["wisexkb"] * 1000,
        wiseXref = config["wisexref"],
        blacklist = config["blacklist"],
        outputid = "plotnpz/{sample}",
    shell:
        "ml purge && ml R/3.4.4-intel-2018a-X11-20180131 && ml Python/3.6.4-intel-2018a; "
        "{params.wiseX} predict {input} {params.wiseXref} {params.outputid} --blacklist {params.blacklist} --plot --bed"

rule hmmcopy:
    input:
        "nodups/{sample}_nodups.bam"
    output:
        "wigfiles/{sample}.wig"
    params:
        walltime="0:20:00",
        ppn = 1,
        hmmcopy = config["hmmcopy"],
        binsize = config["ichorkb"] * 1000
    shell:
        "{params.hmmcopy} --window {params.binsize} --quality 20 "
        "--chromosome 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' "
        "{input} > {output}"

rule ichorCNA:
    input:
        tum = "wigfiles/{sample}.wig"
    output:
        directory("ichorCNA/{sample}/")
    params:
        walltime="0:40:00",
        ppn = 1,
        runichorCNA = config["runichorCNA"],
        ichorCNAfolder = config["ichorCNAfolder"],
        PoNichorCNA = config["PoNichorCNA"],
        kb = config["ichorkb"],
        id = "{sample}"
    shell:
        "ml purge && ml  R/3.4.4-intel-2018a-X11-20180131; "
        "mkdir -p {output} ; "
        "Rscript {params.runichorCNA} --id {params.id} "
        "--WIG {input.tum} "
        "--gcWig {params.ichorCNAfolder}/inst/extdata/gc_hg38_{params.kb}kb.wig "
        "--mapWig {params.ichorCNAfolder}/inst/extdata/map_hg38_{params.kb}kb.wig "
        "--centromere {params.ichorCNAfolder}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt "
        "--normalPanel {params.PoNichorCNA} "
        "--chrs 'c(1:22)' "
        " --chrTrain 'c(1:22)' "
        "--outDir {output} --scStates 'c(1,3)' "
        "--txnE 0.9999 --txnStrength 10000 "
        "--normal 'c(0.2,0.35,0.5,0.65,0.8)' --maxCN 5 "

rule multiqc:
    input:
        expand("nodups/{sample}_nodups.bam", sample = IDS),
        expand(directory("ichorCNA/{sample}/"), sample = IDS),
        expand("plotnpz/{sample}_bins.bed", sample = IDS)
    output:
        report("multiqc_report.html")
    params:
        walltime="0:20:00",
        ppn = 1
    shell:
        "ml purge && ml MultiQC/1.2-intel-2017b-Python-2.7.14; "
        "multiqc -f ."
