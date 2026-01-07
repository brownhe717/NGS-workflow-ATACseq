rule get_sra_se:
	output:
		temp("data/sra/se/{accession}.fastq.gz"),
	conda:
		"../envs/sratools.yaml"
	log:
		"logs/get_sra/{accession}.log",
	threads: 2
	script:
			"../scripts/fasterq-dump.py"
                
rule get_sra_pe:
	output:
		temp("data/sra/pe/{accession}_1.fastq.gz"),
		temp("data/sra/pe/{accession}_2.fastq.gz"),
	conda:
		"../envs/sratools.yaml"
	log:
		"logs/get_sra/{accession}.log",
	threads: 2
	script:
			"../scripts/fasterq-dump.py"
        
rule merge_fastqs:
    input:
        get_fq_merge,
    output:
        temp("data/merged/{sample}_{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|1|2",
    shell:
        "cat {input} > {output} 2> {log}"
        
rule bowtie2_align:
	input:
		sample=["data/trimmed/{sample}_1.fastq.gz","data/trimmed/{sample}_2.fastq.gz"],
		index=rules.bowtie2_index.output
	output:
		temp("results/aligned_reads/mapped/{sample}.bam")
	log:
		"logs/bowtie2/{sample}.log"
	params:
		index=lambda w, input: input.index[0].split('.')[0],  # prefix of reference genome index (built with bowtie2-build)
		extra=config["params"]["bowtie2_align"]   # optional parameters
	threads: 8  # Use at least two threads
	wrapper:
		"v1.1.0/bio/bowtie2/align"

rule summarize_hybrid_mapping:
    input:
        bam="results/aligned_reads/mapped/{sample}.bam"
    output:
        "results/qc/hybrid_mapping/{sample}_mapping_summary.tsv"
    log:
        "logs/qc/hybrid_mapping_{sample}.log"
    conda:
        "../envs/pysam.yaml"
    shell:
        """
        python workflow/scripts/summarize_hybrid_mapping.py \
            {input.bam} {output} > {log} 2>&1
        """

rule samtools_sort:
    input:
       "results/aligned_reads/mapped/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted/{sample}.bam")
    log:
        "logs/samtools_sort/{sample}.log"
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v1.1.0/bio/samtools/sort"
        
rule samtools_index_aligned:
    input:
        "results/aligned_reads/sorted/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted/{sample}.bam.bai")
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"

rule split_hybrid_bam:
    """
    Split sorted BAM into parental BAMs (Dmel vs Dsim) and write a list of reads mapping to both.
    """
    input:
        bam="results/aligned_reads/sorted/{sample}.bam"
    output:
        Dmel_bam="results/aligned_reads/sorted/{sample}_Dmel.bam",
        Dsim_bam="results/aligned_reads/sorted/{sample}_Dsim.bam",
        multi_log="logs/bowtie2/{sample}_multi_parent_reads.txt"
    log:
        "logs/bowtie2/{sample}_split.log"
    threads: 4
    params:
        dmel_prefix=lambda w: config["hybrid_chrom_prefixes"]["Dmel"],
        dsim_prefix=lambda w: config["hybrid_chrom_prefixes"]["Dsim"],
        enabled=lambda w: config["params"]["split_hybrid_bam"]["enable"]
    shell:
        r"""
        set -euo pipefail

        if [ "{params.enabled}" != "True" ] && [ "{params.enabled}" != "true" ] && [ "{params.enabled}" != "1" ]; then
          echo "Hybrid BAM splitting disabled in config; skipping {wildcards.sample}" > {log}
          # Still create empty outputs so downstream doesn't crash if it expects them:
          : > {output.multi_log}
          samtools view -H {input.bam} | samtools view -b -o {output.Dmel_bam} -
          samtools view -H {input.bam} | samtools view -b -o {output.Dsim_bam} -
          exit 0
        fi

        # Extract Dmel and Dsim reads by reference name prefix
        samtools view -b {input.bam} "{params.dmel_prefix}*" > {output.Dmel_bam} 2> {log}
        samtools view -b {input.bam} "{params.dsim_prefix}*" > {output.Dsim_bam} 2>> {log}

        # Identify read names that appear in both parental BAMs
        samtools view {output.Dmel_bam} | cut -f1 | sort > Dmel_reads.txt
        samtools view {output.Dsim_bam} | cut -f1 | sort > Dsim_reads.txt
        comm -12 Dmel_reads.txt Dsim_reads.txt > {output.multi_log}

        rm -f Dmel_reads.txt Dsim_reads.txt
        """
