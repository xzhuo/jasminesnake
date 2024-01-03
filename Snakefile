import os

SUFFIX = '.hifi.bam'
suffix_length = len(SUFFIX)
SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
# SAMPLES = ["bam"]
print(SAMPLES)

REF = "/storage1/fs1/hprc/Active/xzhuo/ref/hg38.fa"

rule all:
    input:
        expand("{sample}.model.combined.bed", sample=SAMPLES),
        expand("{sample}.count.combined.bed", sample=SAMPLES)

rule jasmine:
    input: 
        bam = "{sample}.hifi.bam"
    output:
        bam = "{sample}.5mc.bam"
    threads:
        16
    shell:
        "jasmine -j {threads} {input.bam} {output.bam}"

rule pbmm2:
    input:
        bam = "{sample}.5mc.bam",
        ref = REF
    output:
        bam = "{sample}.5mc.hg38.bam",
        bai = "{sample}.5mc.hg38.bam.bai"
    threads:
        16
    shell:
        "pbmm2 align {input.ref} {input.bam} {output.bam} --preset CCS --sort -j {threads} -m 2G"

rule pbCpGtools:
    input:
        bam = "{sample}.5mc.hg38.bam"
    output:
        model = "{sample}.model.combined.bed",
        counts = "{sample}.count.combined.bed"
    threads:
        8
    params:
        model = "/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        prefix_model = "{sample}.model",
        prefix_count = "{sample}.count"
    run:
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_model} --model {params.model} --threads {threads})")
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_count} --pileup-mode count --threads {threads}")
