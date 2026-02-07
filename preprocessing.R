# Load data
UROMOL_TaLG.teachingcohort <- readRDS("./data/UROMOL_TaLG.teachingcohort.rds")
knowles_matched_TaLG_final <- readRDS("./data/knowles_matched_TaLG_final.rds")

# Remove samples with no Recurrence or RFS_time
UROMOL_TaLG.teachingcohort = dplyr::filter(UROMOL_TaLG.teachingcohort, !is.na(Recurrence))
UROMOL_TaLG.teachingcohort = dplyr::filter(UROMOL_TaLG.teachingcohort, !is.na(RFS_time))

# Get differentially expressed genes
library(limma)

expr = t(UROMOL_TaLG.teachingcohort$exprs)
group <- factor(UROMOL_TaLG.teachingcohort$Recurrence)
design <- model.matrix(~ group)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

res <- topTable(fit, coef=2, number=Inf, sort.by="P")
sig <- subset(res, P.Value < 0.05 & abs(logFC) > 1)

# Save data as csv
write.csv(UROMOL_TaLG.teachingcohort, "./data/UROMOL.csv", quote = FALSE, row.names = FALSE)
write.csv(knowles_matched_TaLG_final, "./data/knowles.csv", quote = FALSE, row.names = FALSE)
write.csv(UROMOL_TaLG.teachingcohort$exprs[,rownames(sig)], "./data/UROMOL_sigExpression.csv", quote = FALSE, row.names = FALSE)
