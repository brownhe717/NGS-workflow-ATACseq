# Get SRA data
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

# Alignment        
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

#Sorting, Fixmate, and Deduplication
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

rule samtools_sort_by_name:
    """Sort BAM by read name (-n) before fixmate."""
    input:
        "results/aligned_reads/mapped/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted_name/{sample}.bam")
    log:
        "logs/samtools_sort_by_name/{sample}.log"
    threads: 8
    shell:
        "samtools sort -n -@ {threads} -o {output} {input} 2> {log}"

rule samtools_fixmate:
    """Add mate information to name-sorted BAM."""
    input:
        "results/aligned_reads/sorted_name/{sample}.bam"
    output:
        temp("results/aligned_reads/fixmate/{sample}.bam")
    log:
        "logs/samtools_fixmate/{sample}.log"
    threads: 4
    shell:
        "samtools fixmate -m {input} {output} 2> {log}"

rule samtools_sort_position:
    """Sort BAM by coordinate after fixmate."""
    input:
        "results/aligned_reads/fixmate/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted_position/{sample}.bam")
    log:
        "logs/samtools_sort_position/{sample}.log"
    threads: 8
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2> {log}"

rule samtools_markdup:
    """Remove duplicates from coordinate-sorted BAM."""
    input:
        "results/aligned_reads/sorted_position/{sample}.bam"
    output:
        "results/aligned_reads/dedup/{sample}.dedup.bam"
    log:
        "logs/samtools_markdup/{sample}.log"
    threads: 8
    shell:
        "samtools markdup -r -@ {threads} {input} {output} 2> {log}"

rule samtools_index_dedup:
    """Index the final deduplicated BAM."""
    input:
        "results/aligned_reads/dedup/{sample}.dedup.bam"
    output:
        "results/aligned_reads/dedup/{sample}.dedup.bam.bai"
    log:
        "logs/samtools_index_dedup/{sample}.log"
    wrapper:
        "v1.1.0/bio/samtools/index"

#Index deduped alignment        
rule samtools_index_aligned:
    input:
        "results/aligned_reads/dedup/{sample}.dedup.bam"
    output:
        temp("results/aligned_reads/dedup/{sample}.dedup.bam.bai")
    log:
        "logs/samtools_index_dedup/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"
