# 필수 라이브러리 (Biostrings 없으면 BiocManager::install("Biostrings") 필요)
library(data.table)
library(biomaRt)
library(stringr)
library(Biostrings)

# 1. Mart 설정 및 데이터 가져오기 (기존 코드 유지)
mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", GRCh=38)
tx_ids <- unique(x2$Transcript_ID)

cds_df <- getSequence(id=tx_ids, type="ensembl_transcript_id", seqType="coding", mart=mart)
pep_df <- getSequence(id=tx_ids, type="ensembl_transcript_id", seqType="peptide", mart=mart)

seq_dt <- merge(
  as.data.table(cds_df)[, .(Transcript_ID=ensembl_transcript_id, cds_seq=coding)],
  as.data.table(pep_df)[, .(Transcript_ID=ensembl_transcript_id, pep_seq=peptide)],
  by="Transcript_ID", all=TRUE
)

# 2. CDS 변이 적용 함수 (SNP, Del, Ins 대응)
apply_hgvsc <- function(cds, hgvsc) {
  if (is.na(cds) | cds == "Sequence unavailable") return(NA_character_)
  
  tryCatch({
    # SNP: c.100A>G
    if (str_detect(hgvsc, ">")) {
      m <- str_match(hgvsc, "c\\.(\\d+)([ACGT])>([ACGT])")
      pos <- as.integer(m[2]); alt <- m[4]
      substr(cds, pos, pos) <- alt
      return(cds)
    }
    # Del: c.100_102del
    if (str_detect(hgvsc, "del")) {
      m <- str_match(hgvsc, "c\\.(\\d+)_?(\\d+)?del")
      s <- as.integer(m[2]); e <- if(is.na(m[3])) s else as.integer(m[3])
      return(paste0(substr(cds, 1, s-1), substr(cds, e+1, nchar(cds))))
    }
    # Ins: c.100_101insATT
    if (str_detect(hgvsc, "ins")) {
      m <- str_match(hgvsc, "c\\.(\\d+)_(\\d+)ins([ACGT]+)")
      s <- as.integer(m[2]); ins_seq <- m[4]
      return(paste0(substr(cds, 1, s), ins_seq, substr(cds, s+1, nchar(cds))))
    }
    return(NA_character_)
  }, error = function(e) return(NA_character_))
}

# 3. 데이터 병합 및 메인 처리
x3 <- merge(x2, seq_dt, by="Transcript_ID", all.x=TRUE)

# (1) Mutant CDS 생성
x3[, cds_mut := mapply(apply_hgvsc, cds_seq, HGVSc)]

# (2) Mutant Peptide 생성 (실제 번역!)
# Biostrings의 translate 함수를 사용하여 DNA -> AA 변환
x3[, pep_mut := sapply(cds_mut, function(dna) {
  if (is.na(dna)) return(NA_character_)
  # DNAString 객체로 변환 후 번역 (Stop codon은 *로 표기)
  as.character(translate(DNAString(dna), if.fuzzy.codon="solve"))
})]

# 4. 윈도우 추출 (중심점 기준 자르기)
slice_window <- function(seq, center, flank) {
  if (is.na(seq) | is.na(center)) return(NA_character_)
  s <- max(1, center - flank)
  e <- min(nchar(seq), center + flank)
  substr(seq, s, e)
}

# 위치 정보 파싱
x3[, cds_pos := as.integer(str_extract(HGVSc, "\\d+"))]
x3[, aa_pos  := as.integer(str_extract(HGVSp_Short, "\\d+"))]

# 최종 윈도우 생성
x3[, WT_CDS_61  := mapply(slice_window, cds_seq, cds_pos, MoreArgs=list(flank=30))]
x3[, MUT_CDS_61 := mapply(slice_window, cds_mut, cds_pos, MoreArgs=list(flank=30))]
#x3[, WT_PEP_21  := mapply(slice_window, pep_seq, aa_pos,  MoreArgs=list(flank=10))]
#x3[, MUT_PEP_21 := mapply(slice_window, pep_mut, aa_pos,  MoreArgs=list(flank=10))]
# 1. 중복된 컬럼 정리 (pep_seq.y를 pep_seq로 변경)
setnames(x3, "pep_seq.y", "pep_seq", skip_absent=TRUE)

# 2. (선택사항) 불필요해진 pep_seq.x는 삭제
if("pep_seq.x" %in% colnames(x3)) x3[, pep_seq.x := NULL]

# 3. 이제 다시 윈도우 계산 실행
x3[, WT_PEP_21 := mapply(slice_window, pep_seq, aa_pos, MoreArgs=list(flank=10))]

# 4. MUT_PEP_21도 혹시 안 되어 있다면 다시 실행
x3[, MUT_PEP_21 := mapply(slice_window, pep_mut, aa_pos, MoreArgs=list(flank=10))]

# 확인
head(x3[, .(Transcript_ID, HGVSp_Short, WT_PEP_21, MUT_PEP_21)])
           
# 결과 확인 및 저장
out <- x3[, .(Hugo_Symbol, Transcript_ID, HGVSc, HGVSp_Short, 
              WT_CDS_61, MUT_CDS_61, WT_PEP_21, MUT_PEP_21)]

fwrite(out, "final_variant_windows.tsv", sep="\t")
