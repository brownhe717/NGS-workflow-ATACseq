import os

HYBRID_ATAC_ENV = "../envs/hybrid_atac.yaml"


rule bowtie2_index_concat:
    input:
        fasta=config["ref_concat"]["genome"]["fasta"]
    output:
        expand(
            config["ref_concat"]["bowtie2_index"]["prefix"] + ".{i}.bt2",
            i=range(1, 5)
        ),
        expand(
            config["ref_concat"]["bowtie2_index"]["prefix"] + ".rev.{i}.bt2",
            i=range(1, 3)
        )
    conda:
        HYBRID_ATAC_ENV
    params:
        prefix=config["ref_concat"]["bowtie2_index"]["prefix"],
        index_dir=lambda wc: os.path.dirname(
            config["ref_concat"]["bowtie2_index"]["prefix"]
        )
    log:
        "logs/hybrid_atac/bowtie2_index_concat.log"
    threads: 8
    shell:
        r"""
        mkdir -p {params.index_dir} logs/hybrid_atac

        bowtie2-build \
            --threads {threads} \
            {input.fasta} \
            {params.prefix} \
            > {log} 2>&1
        """


rule bowtie2_align_concat:
    input:
        idx=rules.bowtie2_index_concat.output,
        fq1=lambda wc: units.loc[wc.sample_name, "fq1"],
        fq2=lambda wc: units.loc[wc.sample_name, "fq2"]
    output:
        bam="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam",
        bai="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    params:
        idx_prefix=config["ref_concat"]["bowtie2_index"]["prefix"],
        extra=lambda wc: (
            config["params"]["bowtie2_align_concat"]
            + f" -k {config['hybrid_atac']['bowtie2_k']}"
        )
    log:
        "logs/hybrid_atac/bowtie2_align/{sample_name}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned logs/hybrid_atac/bowtie2_align

        bowtie2 \
            -x {params.idx_prefix} \
            -1 {input.fq1} \
            -2 {input.fq2} \
            -p {threads} \
            {params.extra} \
            2> {log} | \
        samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index {output.bam} 2>> {log}
        """


rule filter_unique_concat_atac_bam:
    input:
        bam="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam",
        bai="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam.bai"
    output:
        bam="results/hybrid_atac/aligned_unique/{sample_name}.concat.unique.sorted.bam",
        bai="results/hybrid_atac/aligned_unique/{sample_name}.concat.unique.sorted.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    params:
        min_mapq=config["hybrid_atac"]["unique_bam_min_mapq"]
    log:
        "logs/hybrid_atac/filter_unique/{sample_name}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned_unique logs/hybrid_atac/filter_unique

        samtools view \
            -h \
            -q {params.min_mapq} \
            {input.bam} \
            2> {log} \
            | awk '$0 ~ /^@/ || $0 !~ /XS:i:/' \
            | samtools sort \
                -@ {threads} \
                -o {output.bam} \
                2>> {log}

        samtools index {output.bam} 2>> {log}

        echo "Created unique-only ATAC BAM from {input.bam}" >> {log}
        """


rule filter_noMT_multi_bam:
    input:
        bam="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam",
        bai="results/hybrid_atac/aligned/{sample_name}.concat.multimapper.sorted.bam.bai"
    output:
        bam="results/hybrid_atac/aligned_noMT/multimapper_inclusive/{sample_name}.concat.multimapper.noMT.sorted.bam",
        bai="results/hybrid_atac/aligned_noMT/multimapper_inclusive/{sample_name}.concat.multimapper.noMT.sorted.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    params:
        keep=config["hybrid_atac"]["keep_major_chroms_concat"]
    log:
        "logs/hybrid_atac/filter_noMT_multi/{sample_name}.log"
    threads: 4
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned_noMT/multimapper_inclusive logs/hybrid_atac/filter_noMT_multi

        samtools view \
            -h \
            {input.bam} \
            $(cat {params.keep}) \
            2> {log} | \
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            2>> {log}

        samtools index {output.bam} 2>> {log}

        echo "Created noMT multimapper-inclusive BAM from {input.bam}" >> {log}
        """


rule filter_noMT_unique_bam:
    input:
        bam="results/hybrid_atac/aligned_unique/{sample_name}.concat.unique.sorted.bam",
        bai="results/hybrid_atac/aligned_unique/{sample_name}.concat.unique.sorted.bam.bai"
    output:
        bam="results/hybrid_atac/aligned_noMT/unique_only/{sample_name}.concat.unique.noMT.sorted.bam",
        bai="results/hybrid_atac/aligned_noMT/unique_only/{sample_name}.concat.unique.noMT.sorted.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    params:
        keep=config["hybrid_atac"]["keep_major_chroms_concat"]
    log:
        "logs/hybrid_atac/filter_noMT_unique/{sample_name}.log"
    threads: 4
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned_noMT/unique_only logs/hybrid_atac/filter_noMT_unique

        samtools view \
            -h \
            {input.bam} \
            $(cat {params.keep}) \
            2> {log} | \
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            2>> {log}

        samtools index {output.bam} 2>> {log}

        echo "Created noMT unique-only BAM from {input.bam}" >> {log}
        """


HYBRID_ATAC_SAMPLE_GROUPS = sorted(samples["condition"].unique())


def get_hybrid_atac_noMT_multi_bams_by_group(wildcards):
    group_samples = samples.query("condition == @wildcards.sample_group").index
    return expand(
        "results/hybrid_atac/aligned_noMT/multimapper_inclusive/{sample_name}.concat.multimapper.noMT.sorted.bam",
        sample_name=group_samples
    )


def get_hybrid_atac_noMT_unique_bams_by_group(wildcards):
    group_samples = samples.query("condition == @wildcards.sample_group").index
    return expand(
        "results/hybrid_atac/aligned_noMT/unique_only/{sample_name}.concat.unique.noMT.sorted.bam",
        sample_name=group_samples
    )


rule merge_hybrid_atac_noMT_multi_bam:
    input:
        get_hybrid_atac_noMT_multi_bams_by_group
    output:
        bam="results/hybrid_atac/aligned_merged_noMT/multimapper_inclusive/{sample_group}.noMT.bam",
        bai="results/hybrid_atac/aligned_merged_noMT/multimapper_inclusive/{sample_group}.noMT.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    log:
        "logs/hybrid_atac/merge_noMT_multi/{sample_group}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned_merged_noMT/multimapper_inclusive logs/hybrid_atac/merge_noMT_multi

        samtools merge \
            -@ {threads} \
            -f \
            {output.bam} \
            {input} \
            2> {log}

        samtools index {output.bam} 2>> {log}
        """


rule merge_hybrid_atac_noMT_unique_bam:
    input:
        get_hybrid_atac_noMT_unique_bams_by_group
    output:
        bam="results/hybrid_atac/aligned_merged_noMT/unique_only/{sample_group}.noMT.bam",
        bai="results/hybrid_atac/aligned_merged_noMT/unique_only/{sample_group}.noMT.bam.bai"
    conda:
        HYBRID_ATAC_ENV
    log:
        "logs/hybrid_atac/merge_noMT_unique/{sample_group}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/aligned_merged_noMT/unique_only logs/hybrid_atac/merge_noMT_unique

        samtools merge \
            -@ {threads} \
            -f \
            {output.bam} \
            {input} \
            2> {log}

        samtools index {output.bam} 2>> {log}
        """


rule make_bigwigs_hybrid_atac_noMT_multi:
    input:
        bam="results/hybrid_atac/aligned_noMT/multimapper_inclusive/{sample_name}.concat.multimapper.noMT.sorted.bam",
        bai="results/hybrid_atac/aligned_noMT/multimapper_inclusive/{sample_name}.concat.multimapper.noMT.sorted.bam.bai"
    output:
        "results/hybrid_atac/bigwigs_noMT/multimapper_inclusive/{sample_name}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_ind"]
    log:
        "logs/hybrid_atac/bigwigs_noMT_multi/{sample_name}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/bigwigs_noMT/multimapper_inclusive logs/hybrid_atac/bigwigs_noMT_multi

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule make_bigwigs_hybrid_atac_noMT_unique:
    input:
        bam="results/hybrid_atac/aligned_noMT/unique_only/{sample_name}.concat.unique.noMT.sorted.bam",
        bai="results/hybrid_atac/aligned_noMT/unique_only/{sample_name}.concat.unique.noMT.sorted.bam.bai"
    output:
        "results/hybrid_atac/bigwigs_noMT/unique_only/{sample_name}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_ind"]
    log:
        "logs/hybrid_atac/bigwigs_noMT_unique/{sample_name}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/bigwigs_noMT/unique_only logs/hybrid_atac/bigwigs_noMT_unique

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule make_bigwigs_hybrid_atac_noMT_multi_merged:
    input:
        bam="results/hybrid_atac/aligned_merged_noMT/multimapper_inclusive/{sample_group}.noMT.bam",
        bai="results/hybrid_atac/aligned_merged_noMT/multimapper_inclusive/{sample_group}.noMT.bam.bai"
    output:
        "results/hybrid_atac/bigwigs_noMT/multimapper_inclusive_merged/{sample_group}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_merged"]
    log:
        "logs/hybrid_atac/bigwigs_noMT_multi_merged/{sample_group}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/bigwigs_noMT/multimapper_inclusive_merged logs/hybrid_atac/bigwigs_noMT_multi_merged

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule make_bigwigs_hybrid_atac_noMT_unique_merged:
    input:
        bam="results/hybrid_atac/aligned_merged_noMT/unique_only/{sample_group}.noMT.bam",
        bai="results/hybrid_atac/aligned_merged_noMT/unique_only/{sample_group}.noMT.bam.bai"
    output:
        "results/hybrid_atac/bigwigs_noMT/unique_only_merged/{sample_group}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_merged"]
    log:
        "logs/hybrid_atac/bigwigs_noMT_unique_merged/{sample_group}.log"
    threads: 2
    shell:
        r"""
        mkdir -p results/hybrid_atac/bigwigs_noMT/unique_only_merged logs/hybrid_atac/bigwigs_noMT_unique_merged

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule zscore_normalize_hybrid_atac_noMT_multi_ind_bigwigs:
    input:
        "results/hybrid_atac/bigwigs_noMT/multimapper_inclusive/{sample_name}.bw"
    output:
        "results/hybrid_atac/bigwigs_zscore_noMT/multimapper_inclusive/individual/{sample_name}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_atac_noMT_unique_ind_bigwigs:
    input:
        "results/hybrid_atac/bigwigs_noMT/unique_only/{sample_name}.bw"
    output:
        "results/hybrid_atac/bigwigs_zscore_noMT/unique_only/individual/{sample_name}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_atac_noMT_multi_merged_bigwigs:
    input:
        "results/hybrid_atac/bigwigs_noMT/multimapper_inclusive_merged/{sample_group}.bw"
    output:
        "results/hybrid_atac/bigwigs_zscore_noMT/multimapper_inclusive/merged/{sample_group}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_atac_noMT_unique_merged_bigwigs:
    input:
        "results/hybrid_atac/bigwigs_noMT/unique_only_merged/{sample_group}.bw"
    output:
        "results/hybrid_atac/bigwigs_zscore_noMT/unique_only/merged/{sample_group}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"
