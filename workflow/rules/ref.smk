# ----------------------------
# Reference genome (LOCAL FASTA)
# ----------------------------
# Expect config:
# ref_genome:
#   fasta: hybrid_genome.fasta   # local path on execute node
#
# This will stage it into:
#   resources/ref_genome.fasta.gz
# then "rename" to:
#   resources/genome.fasta.gz

rule get_ref_genome:
    """
    Stage a local FASTA reference into resources/ as a gzipped file.
    If the input is already .gz, just copy it. Otherwise gzip it.
    """
    input:
        ref=lambda w: config["ref_genome"]["fasta"]
    output:
        temp("resources/ref_genome.fasta.gz")
    log:
        "logs/get_ref_genome.log"
    conda:
        "../envs/seqkit.yaml"  # provides bgzip/gzip on most installs; if not, use coreutils env
    cache: True
    shell:
        r"""
        set -euo pipefail

        REF="{input.ref}"

        if [ ! -f "$REF" ]; then
            echo "ERROR: Reference FASTA not found: $REF" >&2
            exit 1
        fi

        mkdir -p resources logs

        # If already gzipped, copy; else gzip
        if [[ "$REF" == *.gz ]]; then
            cp -f "$REF" {output} 2> {log}
        else
            gzip -c "$REF" > {output} 2> {log}
        fi
        """

rule rename_genome:
    """
    Keep the original behavior: move staged reference to resources/genome.fasta.gz
    """
    input:
        "resources/ref_genome.fasta.gz"
    output:
        temp("resources/genome.fasta.gz")
    log:
        "logs/rename_genome.log"
    cache: True
    shell:
        r"""
        set -euo pipefail
        mv {input} {output} 2> {log}
        """

# ----------------------------
# Optional: define chromosome keep BED
# ----------------------------
if config.get("filter_chroms", False):
    rule define_keep_chroms:
        input:
            genome="resources/genome.fasta.gz",
            keep_chroms=config["keep_chroms"]
        output:
            "resources/keep_chroms.bed"
        log:
            "logs/define_keep_chroms.log"
        conda:
            "../envs/seqkit.yaml"
        cache: True
        shell:
            r"""
            set -euo pipefail
            seqkit grep -f {input.keep_chroms} {input.genome} \
              | seqkit fx2tab -nil \
              | awk -v OFS='\t' '{{print $1, 1, $2}}' > {output} 2> {log}
            """

# ----------------------------
# Bowtie2 index
# ----------------------------
rule bowtie2_index:
    input:
        reference="resources/genome.fasta.gz"
    output:
        multiext(
            "resources/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        )
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "v1.1.0/bio/bowtie2/build"
