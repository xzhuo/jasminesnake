import os

configfile: "config.yaml"
REF = config["ref"]
FASTA = config["fa"]
tmpdir=config["tmpdir"]
I = config["sample"]
SUFFIX = config["suffix"]
INIT = config["init"]  # one of ["reads", "hifi", "5mc"]
suffix_length = len(SUFFIX)
FILES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))

print(FILES)
if SUFFIX != INIT + ".bam":
    for f in FILEs:
        os.symlink("{file}" + SUFFIX, "{file}" + INIT + ".bam")

rule all:
    input:
        I + ".5mc." + REF + ".bam",
        I + ".5mc." + REF + ".model.combined.bed",
        I + ".5mc." + REF + ".count.combined.bed"

rule extracthifi:
    input:
        bam = "{file}.reads.bam"
    output:
        bam = "{file}.hifi.bam"
    threads:
        4
    shell:
        "extracthifi -j 4 {input.bam} {output.bam}"

rule jasmine:
    input:
        bam = "{file}.hifi.bam"
    output:
        bam = "{file}.5mc.bam",
        log = "{file}.5mc.jasmine.log"
    threads:
        16
    params:
        np = 2,
    shell:
        "jasmine -j {threads} --min-passes {params.np} --log-file {output.log} {input.bam} {output.bam}"

rule fofn:
    input:
        expand("{file}.5mc.bam", file=FILES)
    output:
        I + ".5mc.fofn"
    shell:
        "ls {input} > {output}"

rule pbmm2:
    input:
        fofn = I + ".5mc.fofn",
        fasta = FASTA
    output:
        bam = I + ".5mc." + REF + ".bam",
        bai = I + ".5mc." + REF + ".bam.bai"
    threads:
        32
    shell:
        "pbmm2 align {input.fasta} {input.fofn} {output.bam} --preset CCS --sort -j {threads} -m 2G"

rule pbCpGtools:
    input:
        bam = I + ".5mc." + REF + ".bam"
    output:
        model = I + ".5mc." + REF + ".model.combined.bed",
        counts = I + ".5mc." + REF + ".count.combined.bed"
    threads:
        32
    params:
        model = "/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        prefix_model = I + ".5mc." + REF + ".model",
        prefix_count = I + ".5mc." + REF + ".count"
    run:
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_model} --model {params.model} --threads {threads}")
        shell("aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_count} --pileup-mode count --threads {threads}")
