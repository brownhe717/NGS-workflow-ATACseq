# =========================
# Bigwig creation and normalization
# =========================

wildcard_constraints:
    frag_size = "small|large|total"

# -------------------------
# Individual sample bigwigs
# -------------------------
rule make_bigwigs_ind:
    input:
        bam="results/aligned_reads/split_fragments/{sample}_{frag_size}.bam",
        bai="results/aligned_reads/split_fragments/{sample}_{frag_size}.bam.bai"
    output:
        temp("results/bigwigs/coverage/individual/{sample}_{frag_size}.bw")
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_ind"]
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

# -------------------------
# Merged bigwigs
# -------------------------
rule merge_bam:
    input:
        get_bam_merge
    output:
        temp("results/aligned_reads/merged/{sample_group}_{frag_size}.bam")
    threads: 8
    wrapper:
        "v1.1.0/bio/samtools/merge"

rule samtools_index_merged:
    input:
        "results/aligned_reads/merged/{sample}_{frag_size}.bam"
    output:
        temp("results/aligned_reads/merged/{sample}_{frag_size}.bam.bai")
    log:
        "logs/samtools_index/{sample}_{frag_size}.log"
    threads: 4
    wrapper:
        "v1.1.0/bio/samtools/index"

rule make_bigwigs_merged:
    input:
        bam="results/aligned_reads/merged/{sample}_{frag_size}.bam",
        bai="results/aligned_reads/merged/{sample}_{frag_size}.bam.bai"
    output:
        temp("results/bigwigs/coverage/merged/{sample}_{frag_size}.bw")
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_merged"]
    threads: 8
    resources:
        mem_mb=32000
    shell:
        "bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

# -------------------------
# Z-score normalization
# -------------------------
rule zscore_normalize_ind_bigwigs:
    input:
        "results/bigwigs/coverage/individual/{sample}_{frag_size}.bw"
    output:
        "results/bigwigs/zscore_normalized/individual/{sample}_{frag_size}.bw"
    conda:
        "../envs/zscore_normalize_bw.yaml"
    threads: 1
    resources:
        mem_mb=32000
    script:
        "../scripts/zscore_normalize_bw.R"

rule zscore_normalize_merged_bigwigs:
    input:
        "results/bigwigs/coverage/merged/{sample}_{frag_size}.bw"
    output:
        "results/bigwigs/zscore_normalized/merged/{sample}_{frag_size}.bw"
    conda:
        "../envs/zscore_normalize_bw.yaml"
    threads: 1
    resources:
        mem_mb=32000
    script:
        "../scripts/zscore_normalize_bw.R"
