"""
Germline Variant Calling Pipeline (GATK4 best practices)
Snakemake workflow that reproduces the pipeline used to generate the
results committed to this repository.

Run with:
    snakemake --cores 1 -p

Requires GATK4, BWA, samtools, and bcftools on PATH (see environment.yml).
"""

REF = "data/reference/chr21.fasta"
R1 = "data/raw/sample_R1.fastq.gz"
R2 = "data/raw/sample_R2.fastq.gz"
DBSNP = "data/known_sites/dbsnp_138.hg38.vcf.gz"
MILLS = "data/known_sites/mills_and_1000G.indels.hg38.vcf.gz"
SAMPLE = "sample1"

rule all:
    input:
        "results/vcf/final_pass.vcf.gz",
        "results/tables/bcftools_stats.txt",
        "results/figures/variant_counts.png",
        "results/figures/quality_distribution.png",
        "results/figures/filter_outcome.png",
        "results/figures/depth_distribution.png"

rule index_reference:
    input: REF
    output:
        REF + ".bwt",
        REF + ".fai",
        REF.replace(".fasta", ".dict")
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """

rule index_known_sites:
    input: vcf="{prefix}.vcf.gz"
    output: "{prefix}.vcf.gz.tbi"
    shell: "gatk IndexFeatureFile -I {input.vcf}"

rule align_reads:
    input:
        ref=REF,
        idx=REF + ".bwt",
        r1=R1,
        r2=R2
    output:
        bam="results/bam/sorted.bam",
        bai="results/bam/sorted.bam.bai"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:lib1".format(sample=SAMPLE)
    shell:
        """
        bwa mem -t 1 -R '{params.rg}' {input.ref} {input.r1} {input.r2} \
            | samtools sort -o {output.bam} -
        samtools index {output.bam}
        """

rule mark_duplicates:
    input: "results/bam/sorted.bam"
    output:
        bam="results/bam/dedup.bam",
        metrics="results/tables/dup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}
        samtools index {output.bam}
        """

rule base_recalibrator:
    input:
        bam="results/bam/dedup.bam",
        ref=REF,
        dbsnp=DBSNP,
        dbsnp_idx=DBSNP + ".tbi",
        mills=MILLS,
        mills_idx=MILLS + ".tbi"
    output:
        "results/tables/recal_data.table"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} \
            --known-sites {input.dbsnp} --known-sites {input.mills} \
            -O {output}
        """

rule apply_bqsr:
    input:
        bam="results/bam/dedup.bam",
        ref=REF,
        table="results/tables/recal_data.table"
    output:
        bam="results/bam/recalibrated.bam"
    shell:
        """
        gatk ApplyBQSR -I {input.bam} -R {input.ref} \
            --bqsr-recal-file {input.table} -O {output.bam}
        samtools index {output.bam}
        """

rule haplotype_caller:
    input:
        bam="results/bam/recalibrated.bam",
        ref=REF
    output:
        vcf="results/vcf/raw_variants.vcf.gz"
    shell:
        "gatk HaplotypeCaller -I {input.bam} -R {input.ref} -O {output.vcf}"

rule select_snps:
    input: vcf="results/vcf/raw_variants.vcf.gz", ref=REF
    output: "results/vcf/raw_snps.vcf.gz"
    shell: "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP -O {output}"

rule select_indels:
    input: vcf="results/vcf/raw_variants.vcf.gz", ref=REF
    output: "results/vcf/raw_indels.vcf.gz"
    shell: "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include INDEL -O {output}"

rule filter_snps:
    input: vcf="results/vcf/raw_snps.vcf.gz", ref=REF
    output: "results/vcf/filtered_snps.vcf.gz"
    shell:
        """
        gatk VariantFiltration -R {input.ref} -V {input.vcf} \
            --filter-expression "QD < 2.0" --filter-name QD2 \
            --filter-expression "FS > 60.0" --filter-name FS60 \
            --filter-expression "MQ < 40.0" --filter-name MQ40 \
            --filter-expression "MQRankSum < -12.5" --filter-name MQRankSum-12.5 \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRankSum-8 \
            --filter-expression "SOR > 3.0" --filter-name SOR3 \
            -O {output}
        """

rule filter_indels:
    input: vcf="results/vcf/raw_indels.vcf.gz", ref=REF
    output: "results/vcf/filtered_indels.vcf.gz"
    shell:
        """
        gatk VariantFiltration -R {input.ref} -V {input.vcf} \
            --filter-expression "QD < 2.0" --filter-name QD2 \
            --filter-expression "FS > 200.0" --filter-name FS200 \
            --filter-expression "SOR > 10.0" --filter-name SOR10 \
            -O {output}
        """

rule merge_filtered:
    input:
        snps="results/vcf/filtered_snps.vcf.gz",
        indels="results/vcf/filtered_indels.vcf.gz"
    output:
        "results/vcf/final_filtered.vcf.gz"
    shell:
        "gatk MergeVcfs -I {input.snps} -I {input.indels} -O {output}"

rule pass_only:
    input: "results/vcf/final_filtered.vcf.gz"
    output: "results/vcf/final_pass.vcf.gz"
    shell:
        """
        bcftools view -f PASS {input} -O z -o {output}
        bcftools index -t {output}
        """

rule bcftools_stats:
    input: "results/vcf/final_pass.vcf.gz"
    output: "results/tables/bcftools_stats.txt"
    shell: "bcftools stats {input} > {output}"

rule figures:
    input:
        raw="results/vcf/raw_variants.vcf.gz",
        final="results/vcf/final_filtered.vcf.gz"
    output:
        "results/figures/variant_counts.png",
        "results/figures/quality_distribution.png",
        "results/figures/filter_outcome.png",
        "results/figures/depth_distribution.png"
    shell:
        "python3 scripts/generate_figures.py"
