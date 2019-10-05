# Parallelisation options
import multiprocessing
import sys
cpuCount = (multiprocessing.cpu_count() - 2)

# Sample list
IDS, = glob_wildcards("{sample}.fastq.gz")

# Import config and make report
configfile: "config.yaml"
report: "report/workflow.rst"

rule all:
    input:
        "multiqc_report.html"

rule adaptortrimming:
    input:
        "{sample}.fastq.gz"
    params:
        fastqc="fastqc_before_trimming/{sample}_fastqc.html"
    output:
        temp("trimmed_reads/{sample}_trimmed.fq.gz")
    shell:
        "ml purge && ml Trim_Galore/0.6.0-foss-2018b-Python-3.6.6; ml FastQC/0.11.8-Java-1.8;"
        "trim_galore {input} --fastqc --gzip {config[rrbs]} --three_prime_clip_R1 1 -o ./trimmed_reads"

rule mapping_bismark:
    input:
        genome = config["humangenome"],
        fq = "trimmed_reads/{sample}_trimmed.fq.gz"
    output:
        temp("mapped_reads/{sample}_trimmed_bismark_bt2.bam"),
        temp(directory("bismarkTemp_{sample}"))
    threads: 4
    params:
        tempdir = "bismarkTemp_{sample}"
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b; "
        "bismark --nucleotide_coverage --temp_dir {params.tempdir} --multicore {threads} {input.genome} {input.fq} -o ./mapped_reads/ ;"

rule sort_samtools:
    input:
        bam = "mapped_reads/{sample}_trimmed_bismark_bt2.bam"
    output:
        bam = "sorted_reads/{sample}_trimmed_bismark_bt2.sorted.bam",
        bai = "sorted_reads/{sample}_trimmed_bismark_bt2.sorted.bam.bai"
    threads: 4
    shell:
        "ml purge && ml SAMtools/1.9-intel-2018b; "
        "samtools sort -@ {threads} {input.bam} -o {output.bam}; "
        "samtools index {output.bam} {output.bai} "

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38":
    rule hsmetrics:
        input:
            bam = "sorted_reads/{sample}_trimmed_bismark_bt2.sorted.bam",
            bai = "sorted_reads/{sample}_trimmed_bismark_bt2.sorted.bam.bai"
        output:
            txt = "{sample}_hs_metrics.txt"
        threads: 1
        params:
            genomepath = config["humangenomefa"],
            targetpath = config["target"],
        shell:
            "ml purge && ml picard/2.18.27-Java-1.8; "
            "java -jar $EBROOTPICARD/picard.jar CollectHsMetrics "
            "I={input.bam} "
            "O={output.txt} "
            "R={params.genomepath} "
            "TARGET_INTERVALS={params.targetpath} "
            "BAIT_INTERVALS={params.targetpath} "
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37 or GRCh38)")

rule mapping_to_lambda_bismark:
    input:
        genome = config["lambdagenome"],
        fq = "trimmed_reads/{sample}_trimmed.fq.gz"
    output:
        protected("mapped_reads_lambda/{sample}_lambda_SE_report.txt"),
        temp("mapped_reads_lambda/{sample}_trimmed_bismark_bt2.bam"),
        temp(directory("bismarkTempLambda_{sample}"))
    threads: 2
    params:
        tempdir = "bismarkTempLambda_{sample}"
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "bismark --multicore {threads} --temp_dir {params.tempdir} {input.genome} {input.fq} -o ./mapped_reads_lambda/ ;"
        "cd mapped_reads_lambda/; "
        "rename _trimmed_bismark_bt2_SE_report.txt _lambda_SE_report.txt *txt"

rule extract_methylation_calls:
    input:
        "mapped_reads/{sample}_trimmed_bismark_bt2.bam"
    output:
        "methylation_extractor_output/{sample}_trimmed_bismark_bt2.bismark.cov.gz"
    threads: 6
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b; "
        "bismark_methylation_extractor --bedGraph --gzip --multicore {threads} {input} -o ./methylation_extractor_output/"

rule deduplicate:
    input:
        "mapped_reads/{sample}_trimmed_bismark_bt2.bam"
    output:
        "dedupl_reads/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "deduplicate_bismark -s --bam {input} --output_dir ./dedupl_reads"

rule extract_methylation_calls_deduplicated:
    input:
        "dedupl_reads/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        "methylation_extractor_output/{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
    threads: 5
    params:
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "bismark_methylation_extractor --bedGraph --gzip --multicore {threads} {input} -o ./methylation_extractor_output/"

rule rename:
    input:
        expand("methylation_extractor_output/{sample}_trimmed_bismark_bt2.{dupl}bismark.cov.gz", sample=IDS, dupl = ['deduplicated.', ''])
    output:
        expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=IDS, dupl = ['dedupl', ''])
    shell:
        "cd methylation_extractor_output/ ;"
        "rename _trimmed_bismark_bt2.bismark.cov.gz .cov.gz *_trimmed_bismark_bt2.bismark.cov.gz; "
        "rename _trimmed_bismark_bt2.deduplicated.bismark.cov.gz dedupl.cov.gz *_trimmed_bismark_bt2.deduplicated.bismark.cov.gz;"
        "rm *txt.gz"

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38":
    rule multiqc:
        input:
            expand("{sample}_hs_metrics.txt", sample = IDS),
            expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=IDS, dupl = ['dedupl', '']),
            expand("mapped_reads_lambda/{sample}_lambda_SE_report.txt", sample=IDS)
        output:
            report("multiqc_report.html")
        shell:
            "ml purge && ml Python/3.6.6-intel-2018b; "
            "multiqc -f ."
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37 or GRCh38)")
