rule samtools_idxstats_unfiltered:
    input:
        bam="results/aligned_reads/dedup/{sample}.dedup.bam",  # ← changed from sorted to dedup
        idx="results/aligned_reads/dedup/{sample}.dedup.bam.bai"
    output:
        "results/aligned_reads/stats/{sample}_unfiltered.idxstats"
    log:
        "logs/samtools/idxstats/{sample}.log"
    wrapper:
        "v1.1.0/bio/samtools/idxstats"


if config["filter_chroms"]:
    rule filter_multireads:
        input:
            "results/aligned_reads/dedup/{sample}.dedup.bam"  # ← changed from sorted to dedup
        output:
            temp("results/aligned_reads/unireads/{sample}.bam")
        log:
            "logs/filter_multireads/{sample}.log"
        params:
            extra=config["params"]["filter_multireads"]
        threads: 8
        wrapper:
            "v1.1.0/bio/sambamba/view"

    rule samtools_index_unireads:
        input:
            "results/aligned_reads/unireads/{sample}.bam"
        output:
            temp("results/aligned_reads/unireads/{sample}.bam.bai")
        log:
            "logs/samtools_index/{sample}.log"
        params:
            ""
        threads: 4
        wrapper:
            "v1.1.0/bio/samtools/index"

    rule samtools_idxstats_unireads:
        input:
            bam="results/aligned_reads/unireads/{sample}.bam",
            idx="results/aligned_reads/unireads/{sample}.bam.bai"
        output:
            "results/aligned_reads/stats/{sample}_unireads.idxstats"
        log:
            "logs/samtools/idxstats/{sample}.log"
        wrapper:
            "v1.1.0/bio/samtools/idxstats"

    rule filter_chroms:
        input:
            bam="results/aligned_reads/unireads/{sample}.bam",
            keep_chroms="resources/keep_chroms.bed"
        output:
            "results/aligned_reads/filtered/{sample}.bam"
        log:
            "logs/filter_chroms/{sample}.log"
        conda:
            "../envs/samtools.yaml"
        shell:
            "samtools view -bh -L {input.keep_chroms} --output-fmt BAM -o {output} {input.bam} 2>> {log}"

else:
    rule filter_multireads:
        input:
            "results/aligned_reads/dedup/{sample}.dedup.bam"  # ← changed from sorted to dedup
        output:
            "results/aligned_reads/filtered/{sample}.bam"
        log:
            "logs/filter_multireads/{sample}.log"
        params:
            extra=config["params"]["filter_multireads"]
        threads: 8
        wrapper:
            "v1.1.0/bio/sambamba/view"


rule samtools_index_filtered:
    input:
        "results/aligned_reads/filtered/{sample}.bam"
    output:
        "results/aligned_reads/filtered/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        ""
    threads: 4
    wrapper:
        "v1.1.0/bio/samtools/index"


rule samtools_idxstats_filtered:
    input:
        bam="results/aligned_reads/filtered/{sample}.bam",
        idx="results/aligned_reads/filtered/{sample}.bam.bai"
    output:
        "results/aligned_reads/stats/{sample}_filtered.idxstats"
    log:
        "logs/samtools/idxstats/{sample}.log"
    wrapper:
        "v1.1.0/bio/samtools/idxstats"
