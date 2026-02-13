

tpm <- function(raw_counts, gene_lengths) {
  x <- raw_counts * 1e3 / gene_lengths        # RPK
  return(t(t(x) * 1e6 / colSums(x)))          # normalize by library size
}


rownames(d) <- d$Geneid
lengths <- d$Length
counts <- d[, 7:ncol(d), drop = FALSE]

# TPM 계산
tpm_mat <- tpm(counts, lengths)

# log2(TPM + 1)
log_tpm <- log2(tpm_mat + 1)

for ( i in 1:length(m2.file.list) ){
d <- read.table(paste("/home/fff0996/MUTATION2/RNA/",m2.file.list[i],sep=""),sep="\t",header=T)
rownames(d) <- d$Geneid
lengths <- d$Length
counts <- d[, 7:ncol(d), drop = FALSE]
tpm_mat <- tpm(counts, lengths)
log_tpm <- log2(tpm_mat + 1)
if(i == 1){
df <- log_tpm
}
else{
df <- cbind(df,log_tpm)
}
}



#library(DESeq2) #raw count
library(edgeR) 
library(limma) #FPKM, TPM

# 데이터 불러오기
data <- read.csv(“gene_expressiong.csv”, row.names = 1, header = TRUE, check.names=F)

# 표현형 데이터 정보 불러오기 
# 유전자 발현데이터 샘플과 일치해야함
# 일부 customizing 필요
metadata <- read.csv(“metadata.csv”)

# FPKM이라면 TPM변화 추천 (샘플 간 비교면 TPM변환 일반적)#fpkmToTpm <- function(fpkm)
#{
  #  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
#}

# FPKM 데이터를 TPM으로 변환
#Tpm_data <- fpkmToTpm(data)

# log2 적용
data_log2 <- log2(data + 1) 

# 디자인 매트릭스 생성 정규화
data_log_norm <- normalizeBetweenArrays(data_log, method="quantile")

# 차이를 볼 그룹
grp <- metadata$group
design <- model.matrix(~0 + grp)
colnames(design) <- c("APOBEC","non_APOBEC")
# 그룹이 3개 이상일 때,
#colnames(design) <- c("group1","group2","group3")

# Contrast 설정 (APOBEC vs non-APOBEC)
contrast.matrix <- makeContrasts(APOBEC - non_APOBEC, levels=design)

# 그룹이 3개 이상일 때,
contrast_matrix <- makeContrasts(
    Group2_vs_Group1 = group2 - group1,
    Group3_vs_Group1 = group3 - group1,
    Group3_vs_Group2 = group3 - group2,
    levels = design
)

# voom 변환
v <- voom(data_log2, design, plot=TRUE)

# Linear Model 피팅
fit <- lmFit(data_log_norm, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#DEG 추출
DEG_results <- topTable(fit2, number=Inf, adjust="fdr")

# 그룹이 3개 이상일 때,
# Group2 vs Group1
#res_2_vs_1 <- topTable(fit2, coef="Group2_vs_Group1", adjust="BH", number=Inf)

# Group3 vs Group1
#res_3_vs_1 <- topTable(fit2, coef="Group3_vs_Group1", adjust="BH", number=Inf)

# Group3 vs Group2
#res_3_vs_2 <- topTable(fit2, coef="Group3_vs_Group2", adjust="BH", number=Inf)

# 발현도 차이가 통계적으로 유의했던 유전자 구별 (p-value < 0.05 & logFC > 1)
DEG_results$sig <- "Not significance"
DEG_results[DEG_results$adj.P.Val < 0.05 & DEG_results$logFC > 1,]$sig <- "UP"
DEG_results[DEG_results$adj.P.Val < 0.05 & DEG_results$logFC < -1,]$sig <- "DOWN"

# 결과 볼테이노 플랏으로 그리기
library(ggplot2)

ggplot(DEG_results, aes(x=logFC, y=-log10(adj.P.Val), color=sig)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "Not significance" = "gray")) +
  labs(title = "Volcano Plot of DEG Analysis",
       x = "Log2 Fold Change",
       y = "-Log10(P-value)") +
  theme_minimal()
       
# PCA그리기        
pca_res <- prcomp(t(data_log_norm), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- group_info$APOBEC_group  # 예: "APOBEC" or "non_APOBEC"

 ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA: APOBEC vs non-APOBEC")



###GSVA

library(GSVA)
library(msigdbr)

# gene set 가져오기 (Hallmark)
msig <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# expression matrix 만들기
expr_mat <- as.matrix(m[, -1])
rownames(expr_mat) <- m$gene_symbol

# GSVA 실행
gsva_res <- gsva(expr_mat, gene_sets, method = "gsva")
# ssGSEA 하고 싶으면 method = "ssgsea"



###GSEA
library(DESeq2)
library(clusterProfiler)
library(msigdbr)

# 1. count matrix 준비 (raw count, gene_symbol x samples)
count_mat <- as.matrix(m[, -1])
rownames(count_mat) <- m$gene_symbol

# 중복 gene symbol 합치기
count_mat <- round(count_mat)  # DESeq2는 integer 필요

# 2. sample info (여기 수정)
condition <- factor(c("shNC","shNC","shTP53","shTP53",
                      "xeno","xeno","xeno","xeno","xeno","xeno","xeno","xeno","xeno","xeno"))
coldata <- data.frame(condition = condition, row.names = colnames(count_mat))

# 3. DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

# 비교 그룹 지정 (treat vs control)
res <- results(dds, contrast = c("condition", "shTP53", "shNC"))
res <- as.data.frame(res)
res$gene_symbol <- rownames(res)

# NA 제거
res <- res[!is.na(res$log2FoldChange) & !is.na(res$pvalue), ]

# 4. ranked gene list
gene_list <- res$log2FoldChange
names(gene_list) <- res$gene_symbol
gene_list <- sort(gene_list, decreasing = TRUE)

# 5. gene set
msig <- msigdbr(species = "Homo sapiens", category = "H")
msig_t2g <- msig[, c("gs_name", "gene_symbol")]

# 6. GSEA
gsea_res <- GSEA(geneList = gene_list,
                 TERM2GENE = msig_t2g,
                 pvalueCutoff = 0.05)

# 7. 시각화
dotplot(gsea_res, showCategory = 20)
gseaplot2(gsea_res, geneSetID = 1:3)
ridgeplot(gsea_res, showCategory = 20)

# 8. 결과 저장
write.csv(as.data.frame(gsea_res), "gsea_results.csv", row.names = FALSE)
