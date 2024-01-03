import os

SUFFIX = '.reads.bam'
suffix_length = len(SUFFIX)
SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
# SAMPLES = ["bam"]
print(SAMPLES)

REF = "/storage1/fs1/hprc/Active/xzhuo/ref/hg38.fa"
configfile: "config.yaml"
I = config["sample"]
rule all:
    input:
        I + ".model.combined.bed",
        I + ".count.combined.bed"

rule extracthifi:
    input:
        bam = "{sample}.reads.bam"
    output:
        bam = "{sample}.hifi.bam"
    threads:
        4
    shell:
        "extracthifi -j 4 {input.bam} {output.bam}"

rule jasmine:
    input:
        bam = "{sample}.hifi.bam"
    output:
        bam = "{sample}.5mc.bam"
    threads:
        16
    shell:
        "jasmine -j {threads} {input.bam} {output.bam}"

rule fofn:
    input:
        expand("{sample}.5mc.bam", sample=SAMPLES)
    output:
        I + ".5mc.fofn"
    shell:
        "ls {input} > {output}"

rule pbmm2:
    input:
        fofn = I + ".5mc.fofn"
        ref = REF
    output:
        bam = I + ".5mc.hg38.bam",
        bai = I + ".5mc.hg38.bam.bai"
    threads:
        16
    shell:
        "pbmm2 align {input.ref} {input.fofn} {output.bam} --preset CCS --sort -j {threads} -m 2G"

rule pbCpGtools:
    input:
        bam = I + ".5mc.hg38.bam"
    output:
        model = I + ".model.combined.bed",
        counts = I + ".count.combined.bed"
    threads:
        8
    params:
        model = "/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        prefix_model = I +".model",
        prefix_count = I + ".count"
    run:
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_model} --model {params.model} --threads {threads})")
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_count} --pileup-mode count --threads {threads}")
