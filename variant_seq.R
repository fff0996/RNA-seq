#!/usr/bin/env Rscript
# HGVSc/HGVSp 기반 WT/MUT 서열 윈도우 추출 파이프라인
# Author: 혜인's bioinformatics team
# Date: 2025-01-19

library(data.table)
library(biomaRt)
library(stringr)
library(Biostrings)

# ============================================================================
# 1. Ensembl에서 CDS/Peptide 시퀀스 가져오기
# ============================================================================
cat("Step 1: Fetching sequences from Ensembl...\n")
mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", GRCh=38)
tx_ids <- unique(x2$Transcript_ID)

cds_df <- getSequence(id=tx_ids, type="ensembl_transcript_id", seqType="coding", mart=mart)
pep_df <- getSequence(id=tx_ids, type="ensembl_transcript_id", seqType="peptide", mart=mart)

seq_dt <- merge(
  as.data.table(cds_df)[, .(Transcript_ID=ensembl_transcript_id, cds_seq=coding)],
  as.data.table(pep_df)[, .(Transcript_ID=ensembl_transcript_id, pep_seq=peptide)],
  by="Transcript_ID", all=TRUE
)

cat(sprintf("  - %d transcripts fetched\n", nrow(seq_dt)))

# ============================================================================
# 2. CDS 변이 적용 함수 (SNP/Del/Ins/Frameshift 대응)
# ============================================================================
apply_hgvsc <- function(cds, hgvsc) {
  if (is.na(cds) | is.na(hgvsc) | cds == "Sequence unavailable") {
    return(NA_character_)
  }
  
  tryCatch({
    # Frameshift, nonsense, splice는 일단 원본 반환 (추후 개선 가능)
    if (str_detect(hgvsc, "fs|\\*|splice|dup")) {
      return(cds)
    }
    
    # SNP: c.100A>G
    if (str_detect(hgvsc, ">")) {
      m <- str_match(hgvsc, "c\\.(\\d+)([ACGT])>([ACGT])")
      if (!any(is.na(m))) {
        pos <- as.integer(m[2])
        alt <- m[4]
        if (pos <= nchar(cds)) {
          substr(cds, pos, pos) <- alt
          return(cds)
        }
      }
    }
    
    # Deletion: c.100_102del or c.100del
    if (str_detect(hgvsc, "del")) {
      m <- str_match(hgvsc, "c\\.(\\d+)_?(\\d+)?del")
      if (!any(is.na(m[1:2]))) {
        s <- as.integer(m[2])
        e <- if(is.na(m[3])) s else as.integer(m[3])
        if (s <= nchar(cds) && e <= nchar(cds)) {
          return(paste0(substr(cds, 1, s-1), substr(cds, e+1, nchar(cds))))
        }
      }
    }
    
    # Insertion: c.100_101insATT
    if (str_detect(hgvsc, "ins")) {
      m <- str_match(hgvsc, "c\\.(\\d+)_(\\d+)ins([ACGT]+)")
      if (!any(is.na(m[1:3]))) {
        s <- as.integer(m[2])
        ins_seq <- m[4]
        if (s <= nchar(cds)) {
          return(paste0(substr(cds, 1, s), ins_seq, substr(cds, s+1, nchar(cds))))
        }
      }
    }
    
    return(NA_character_)
  }, error = function(e) {
    warning(sprintf("Failed to apply HGVSc '%s': %s", hgvsc, e$message))
    return(NA_character_)
  })
}

# ============================================================================
# 3. 번역 함수 (stop codon 이후 제거)
# ============================================================================
translate_cds <- function(dna) {
  if (is.na(dna) || nchar(dna) < 3) return(NA_character_)
  
  tryCatch({
    pep <- as.character(translate(DNAString(dna), if.fuzzy.codon="solve"))
    # Stop codon(*) 이후 모두 제거
    pep_clean <- str_extract(pep, "^[^*]+\\*?")
    return(pep_clean)
  }, error = function(e) {
    return(NA_character_)
  })
}

# ============================================================================
# 4. 윈도우 슬라이싱 함수
# ============================================================================
slice_window <- function(seq, center, flank) {
  if (is.na(seq) || is.na(center) || center < 1) return(NA_character_)
  
  s <- max(1, center - flank)
  e <- min(nchar(seq), center + flank)
  
  if (s > nchar(seq)) return(NA_character_)
  
  substr(seq, s, e)
}

# ============================================================================
# 5. 메인 파이프라인
# ============================================================================
cat("Step 2: Merging with input data...\n")

# **중요**: x2에 pep_seq 컬럼이 있으면 충돌 방지를 위해 삭제
if ("pep_seq" %in% colnames(x2)) {
  x2[, pep_seq := NULL]
}

x3 <- merge(x2, seq_dt, by="Transcript_ID", all.x=TRUE)

cat("Step 3: Applying mutations to CDS...\n")
x3[, cds_mut := mapply(apply_hgvsc, cds_seq, HGVSc, SIMPLIFY=TRUE)]

cat("Step 4: Translating mutant CDS to peptide...\n")
x3[, pep_mut := sapply(cds_mut, translate_cds)]

cat("Step 5: Extracting position info and creating windows...\n")
# HGVSc/HGVSp에서 위치 추출
x3[, cds_pos := as.integer(str_extract(HGVSc, "\\d+"))]
x3[, aa_pos  := as.integer(str_extract(HGVSp_Short, "\\d+"))]

# 61nt CDS 윈도우 (±30)
x3[, WT_CDS_61  := mapply(slice_window, cds_seq, cds_pos, MoreArgs=list(flank=30))]
x3[, MUT_CDS_61 := mapply(slice_window, cds_mut, cds_pos, MoreArgs=list(flank=30))]

# 21aa Peptide 윈도우 (±10)
x3[, WT_PEP_21  := mapply(slice_window, pep_seq, aa_pos, MoreArgs=list(flank=10))]
x3[, MUT_PEP_21 := mapply(slice_window, pep_mut, aa_pos, MoreArgs=list(flank=10))]

# ============================================================================
# 6. QC 체크 (선택사항)
# ============================================================================
cat("Step 6: Quality control checks...\n")

n_total <- nrow(x3)
n_cds_ok <- sum(!is.na(x3$cds_mut))
n_pep_ok <- sum(!is.na(x3$pep_mut))

cat(sprintf("  - Total variants: %d\n", n_total))
cat(sprintf("  - CDS mutations applied: %d (%.1f%%)\n", n_cds_ok, 100*n_cds_ok/n_total))
cat(sprintf("  - Peptide translations: %d (%.1f%%)\n", n_pep_ok, 100*n_pep_ok/n_total))

# 번역 길이 체크
x3[, cds_len := nchar(cds_seq)]
x3[, expected_pep_len := cds_len %/% 3]
x3[, actual_pep_len := nchar(pep_seq)]

weird_translations <- x3[abs(expected_pep_len - actual_pep_len) > 2]
if (nrow(weird_translations) > 0) {
  cat(sprintf("  ⚠️  Warning: %d transcripts have unexpected peptide lengths\n", 
              nrow(weird_translations)))
}

# ============================================================================
# 7. 최종 결과 출력
# ============================================================================
final_output <- x3[, .(
  Hugo_Symbol, 
  Transcript_ID, 
  HGVSc, 
  HGVSp_Short, 
  WT_CDS_61,   # Wild-type CDS 61nt window
  MUT_CDS_61,  # Mutant CDS 61nt window
  WT_PEP_21,   # Wild-type Peptide 21aa window
  MUT_PEP_21   # Mutant Peptide 21aa window
)]

cat("\nStep 7: Done! Preview of results:\n")
print(head(final_output, 10))

# 파일로 저장
output_file <- "hgvs_windows_output.tsv"
fwrite(final_output, output_file, sep="\t", quote=FALSE)
cat(sprintf("\n✅ Results saved to: %s\n", output_file))

# 결과 통계
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Complete windows (all 4 columns filled): %d (%.1f%%)\n",
            sum(complete.cases(final_output[, .(WT_CDS_61, MUT_CDS_61, WT_PEP_21, MUT_PEP_21)])),
            100*mean(complete.cases(final_output[, .(WT_CDS_61, MUT_CDS_61, WT_PEP_21, MUT_PEP_21)]))))

