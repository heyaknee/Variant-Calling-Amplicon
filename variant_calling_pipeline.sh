#!/usr/bin/env bash
set -euo pipefail

 ── CONFIGURATION ────────────────────────────────────────────────────────── 
THREADS=40

# Reference files
REF_DIR="$HOME/reference_data"
REF="$REF_DIR/hg38.fa"
BED="$REF_DIR/targets.bed"
DICT="${REF%.fa}.dict"
FAI="${REF}.fai"
DBSNP="$REF_DIR/dbsnp.vcf.gz"
SNPEFF_DB="hg38"

# Input directories (Paired-end FASTQs)
CONTROL_DIR="$HOME/data/control_samples"
CASE_DIR="$HOME/data/case_samples"

# Output directories
PROJECT="$HOME/variant_analysis"
LOGS="$PROJECT/logs"
TRIMMED="$PROJECT/trimmed"
ALIGN_BAM="$PROJECT/alignment/bam"
ALIGN_BAI="$PROJECT/alignment/bai"
COV_DIR="$PROJECT/coverage"
GVCF_DIR="$PROJECT/gvcf"
VCF_DIR="$PROJECT/vcf"
ANNOT_DIR="$PROJECT/annotation"

mkdir -p "$LOGS" "$TRIMMED" "$ALIGN_BAM" "$ALIGN_BAI" "$COV_DIR" "$GVCF_DIR" "$VCF_DIR" "$ANNOT_DIR"

echo "[`date`] Starting Variant Calling Pipeline" | tee "$LOGS/pipeline.log"

 ── STEP 0: REFERENCE INDEXING ───────────────────────────────────────────── 
echo "[`date`]  Checking reference indexes..." | tee -a "$LOGS/pipeline.log"
[[ -f "$FAI" ]] || samtools faidx "$REF"
[[ -f "$DICT" ]] || picard CreateSequenceDictionary R="$REF" O="$DICT"
[[ -f "${REF}.bwt" ]] || bwa index -a bwtsw "$REF"

 ── FUNCTION: PROCESS SAMPLE DIRECTORY ───────────────────────────────────── 
process_samples() {
  local SAMPLE_DIR="$1"
  local CONDITION="$2"

  echo "[`date`] Processing $CONDITION samples..." | tee -a "$LOGS/pipeline.log"

  for R1 in "$SAMPLE_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
    R2="$SAMPLE_DIR/${SAMPLE}_R2_001.fastq.gz"
    OUT1="$TRIMMED/${SAMPLE}_R1.trimmed.fastq.gz"
    OUT2="$TRIMMED/${SAMPLE}_R2.trimmed.fastq.gz"
    LOG="$LOGS/${SAMPLE}.fastp.log"

    # STEP 1: Trimming
    fastp -i "$R1" -I "$R2" -o "$OUT1" -O "$OUT2" \
      --detect_adapter_for_pe --cut_front --cut_tail \
      --qualified_quality_phred 20 --length_required 50 \
      --thread "$THREADS" > "$LOG" 2>&1

    echo "[`date`]  Trimmed: $SAMPLE" >> "$LOGS/pipeline.log"

    # STEP 2: Alignment
    BAM_UNSORTED="$ALIGN_BAM/${SAMPLE}.unsorted.bam"
    BAM_SORTED="$ALIGN_BAM/${SAMPLE}.bam"
    BAM_TMP="$ALIGN_BAM/${SAMPLE}.tmp.bam"
    LOG_ALIGN="$LOGS/${SAMPLE}.align.log"

    echo "[`date`] 🔧 Aligning $SAMPLE" | tee -a "$LOGS/pipeline.log"

    bwa mem -t "$THREADS" "$REF" "$OUT1" "$OUT2" 2>>"$LOG_ALIGN" | \
      samtools view -@ "$THREADS" -Sb - > "$BAM_UNSORTED"

    samtools sort -@ "$THREADS" -o "$BAM_SORTED" "$BAM_UNSORTED"
    rm -f "$BAM_UNSORTED"

    gatk AddOrReplaceReadGroups -I "$BAM_SORTED" -O "$BAM_TMP" \
      -LB LIBRARY -PL ILLUMINA -PU "$SAMPLE" -SM "$SAMPLE" \
      --CREATE_INDEX true -SO coordinate 2>>"$LOG_ALIGN"

    mv "$BAM_TMP" "$BAM_SORTED"

    if [[ -f "${BAM_TMP}.bai" ]]; then
      mv "${BAM_TMP}.bai" "$ALIGN_BAI/${SAMPLE}.bam.bai"
    elif [[ -f "${BAM_SORTED}.bai" ]]; then
      mv "${BAM_SORTED}.bai" "$ALIGN_BAI/${SAMPLE}.bam.bai"
    else
      samtools index "$BAM_SORTED"
      cp "${BAM_SORTED}.bai" "$ALIGN_BAI/${SAMPLE}.bam.bai"
    fi

    echo "[`date`]  Aligned and indexed: $SAMPLE" >> "$LOGS/pipeline.log"

    # STEP 3: Coverage
    bedtools coverage -f 0.35 -a "$BED" -b "$BAM_SORTED" > "$COV_DIR/${SAMPLE}_coverage.txt"
    echo "[`date`] Coverage calculated: $SAMPLE" >> "$LOGS/pipeline.log"
  done
}

 ── PROCESS CONTROL AND CASE SAMPLES ────────────────────────────────────────
process_samples "$CONTROL_DIR" "Control"
process_samples "$CASE_DIR" "Case"

 ── STEP 4: GATK HAPLOTYPECALLER ──────────────────────────────────────────── 
echo "[`date`]  Running GATK HaplotypeCaller..." | tee -a "$LOGS/pipeline.log"
for BAM in "$ALIGN_BAM"/*.bam; do
  SAMPLE=$(basename "$BAM" .bam)
  OUT_GVCF="$GVCF_DIR/${SAMPLE}.g.vcf.gz"
  LOG_GVCF="$LOGS/${SAMPLE}.gvcf.log"

  gatk --java-options "-Xmx16G" HaplotypeCaller \
    -R "$REF" -I "$BAM" -O "$OUT_GVCF" \
    -ERC GVCF -L "$BED" \
    --native-pair-hmm-threads "$THREADS" 2>>"$LOG_GVCF"

  echo "[`date`]  GVCF created: $SAMPLE" >> "$LOGS/pipeline.log"
done

 ── STEP 5: GATK JOINT GENOTYPING ─────────────────────────────────────────── 
echo "[`date`]  Performing joint genotyping..." | tee -a "$LOGS/pipeline.log"
COHORT_GVCF="$VCF_DIR/cohort.g.vcf.gz"
COHORT_VCF="$VCF_DIR/cohort.vcf.gz"

GVCF_LIST=( "$GVCF_DIR"/*.g.vcf.gz )
gatk CombineGVCFs -R "$REF" "${GVCF_LIST[@]/#/-V }" -O "$COHORT_GVCF" 2>>"$LOGS/combine.log"

gatk GenotypeGVCFs -R "$REF" \
  -V "$COHORT_GVCF" \
  --dbsnp "$DBSNP" \
  -O "$COHORT_VCF" 2>>"$LOGS/genotype.log"

 ── DONE ──────────────────────────────────────────────────────────────────── 
echo "[`date`]  COmpleted!" | tee -a "$LOGS/pipeline.log"
