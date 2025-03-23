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


