library(limma)
load(file = 'datasets.Rdata')

res <- lapply(DAT, function(x) {
  m <- Biobase::exprs(x)
  md <- Biobase::pData(x)
  platform <- x@annotation
  design <- cbind(Control = 1, Tumor = md$Group == 'd')
  fit <- lmFit(as.data.frame(m), design)
  fit <- eBayes(fit)
  y <- topTable(fit, number = 'Inf', adjust.method = 'BH')
  # summary(decideTests(fit))
  entrez <- mapIds(x = eval(parse(text = paste0(platform, '.db'))),
                   keys = rownames(y),
                   keytype = 'PROBEID', 
                   column = 'ENTREZID')
  geneList <- y$adj.P.Val
  names(geneList) <- entrez
  geneList <- -log10(geneList[geneList < 0.05])
  geneList <- geneList[!duplicated(names(geneList))]
  # sig_entrez <- entrez[y$adj.P.Val < 0.05]
  kk <- gseKEGG(gene = geneList,
                organism = 'hsa',
                nPerm = 10000,
                minGSSize = 10,
                pvalueCutoff = 0.05,
                verbose = FALSE)
  return(head(kk, 10)$Description)
})

geneList <- y$adj.P.Val
names(geneList) <- entrez
geneList <- -log10(geneList[geneList < 0.05])
geneList <- geneList[!duplicated(names(geneList))]
kk <- gseKEGG(gene = geneList,
              organism = 'hsa',
              nPerm = 10000,
              minGSSize = 10,
              pvalueCutoff = 0.05,
              verbose = FALSE)
head(kk, 10)$Description

x <- DAT[[5]]
m <- Biobase::exprs(x)
colSums(m)
md <- Biobase::pData(x)
platform <- x@annotation
design <- model.matrix(~ md$Group)
fit <- lmFit(as.data.frame(m), design)
fit <- eBayes(fit)
y <- topTable(fit, number = 'Inf', adjust.method = 'BH')
sum(y$adj.P.Val < 0.05)
# summary(decideTests(fit))
entrez <- mapIds(x = eval(parse(text = paste0(platform, '.db'))),
                 keys = rownames(y),
                 keytype = 'PROBEID', 
                 column = 'ENTREZID')
sig_entrez <- entrez[y$adj.P.Val < 0.05]
kk <- enrichKEGG(gene = sig_entrez,
                 organism = 'hsa', 
                 pvalueCutoff = 0.05)
head(kk, 10)$Description

library(ggplot2)
library(cowplot)
library(dplyr)
ggplots <- lapply(DAT, function(x) {
  m <- Biobase::exprs(x)
  md <- Biobase::pData(x)
  pca <- prcomp(t(m), scale. = FALSE)
  pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  data_pca <- tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Tissue = md$Group)
  gg <- ggplot(data = data_pca, mapping = aes(x = PC1, y = PC2, color = Tissue)) +
    geom_point(size = 3) + ggtitle('PCA plot') +
    xlab(paste0('PC1, (', pv[1], ')%')) + 
    ylab(paste0('PC2, (', pv[2], ')%'))
})
plot_grid(plotlist = ggplots, ncol = 1)

# edgering
library(edgeR)
er <- DGEList(counts = m, group = md$Group)
er <- calcNormFactors(er, method = 'RLE')
er <- estimateDisp(er, model.matrix(~ md$Group))
gfit <- glmFit(er, model.matrix(~ md$Group), dispersion = pmax(er$tagwise.dispersion,
                                                               er$trended.dispersion,
                                                               er$common.dispersion))
lrt_group <- glmLRT(gfit, 2)
z <- topTags(lrt_group, n = Inf, adjust.method = 'BH', sort.by = 'none')
sum(z$table$FDR < 0.05)
entrez <- mapIds(x = eval(parse(text = paste0(platform, '.db'))),
                 keys = rownames(z$table),
                 keytype = 'PROBEID', 
                 column = 'ENTREZID')
sig_entrez <- entrez[z$table$FDR < 0.05]
kk <- enrichKEGG(gene = sig_entrez,
                 organism = 'hsa', 
                 pvalueCutoff = 0.05)
head(kk, 10)$Description

# new
n <- 13
m <- as.data.frame(exprs(DAT[[n]]))
md <- pData(DAT[[n]])
platform <- DAT[[n]]@annotation
plotMDS(m, col = ifelse(md$Group == 'c', 'blue', 'red'), pch = 19)
plotRLDF(y = m, design = model.matrix(~ md$Group), pch = 19,
         col = ifelse(md$Group == 'c', 'blue', 'red'))
anno <- AnnotationDbi::select(x = hgu133plus2.db, keys = rownames(m),
                      columns = c('ENTREZID'),
                      keytype = 'PROBEID')
anno <- anno[!is.na(anno$ENTREZID), ]
anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- summarise(anno_grouped, no_of_matches = n_distinct(ENTREZID))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
anno <- anno[!anno$PROBEID %in% anno_filtered$PROBEID, ]
length(unique(anno$ENTREZID))
m$PROBEID <- rownames(m)
m <- left_join(anno, m, by = "PROBEID")
rownames(m) <- m$PROBEID
m$PROBEID <- NULL
m$ENTREZID <- NULL
design <- model.matrix(~ md$Group)
fit <- lmFit(as.data.frame(m), design)
fit <- eBayes(fit)
y <- topTable(fit, number = 'Inf', adjust.method = 'BH')
y$PROBEID <- rownames(y)
y <- left_join(anno, y, "PROBEID")
sigenes <- unique(y$ENTREZID[y$adj.P.Val < 0.05 & y$logFC > 2])
kk <- enrichKEGG(gene = sigenes,
                 organism = 'hsa', 
                 pvalueCutoff = 0.05)
head(kk, 10)$Description
library(DOSE)
edo <- enrichDO(gene = sigenes, ont = 'DO', pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                universe = unique(y$ENTREZID), minGSSize = 5, maxGSSize = 500,
                qvalueCutoff = 0.2, readable = FALSE)
head(edo)$Description

# permutations
ranked <- y %>% 
  group_by(ENTREZID) %>% 
  summarise(q_value = min(adj.P.Val), logFC = max(logFC)) %>% 
  filter(q_value < 0.05, logFC > 2) %>% 
  arrange(desc(logFC))
geneList <- ranked$logFC
names(geneList) <- ranked$ENTREZID
kk2 <- gseKEGG(geneList = geneList,
               organism = 'hsa',
               nPerm = 10000,
               minGSSize = 5,
               pvalueCutoff = 0.05,
               verbose = FALSE)
head(kk2)$Description
edo2 <- gseDO(geneList = geneList,
              minGSSize = 5,
              pvalueCutoff = 0.2, 
              pAdjustMethod = 'BH',
              verbose = FALSE)
head(edo2)$Description

kegres <- kegga(de = sigenes, universe = unique(y$ENTREZID), species = 'Hs', FDR = 0.05)
head(kegres %>% arrange(P.DE), 10)
write.table(x = sigenes, file = 'lol.txt', row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# new try
library(Biobase)
library(hgu133plus2.db)
library(hgu133a.db)
library(dplyr)
library(clusterProfiler)
library(DOSE)
n <- 3
expr_mat <- as.data.frame(exprs(ma_data[[n]]))
dim(expr_mat)
f <- factor(pData(ma_data[[n]])$Group)
# plotMDS(expr_mat, col = ifelse(f == 'c', 'blue', 'red'), pch = 19)

anno <- AnnotationDbi::select(x = hgu133plus2.db, keys = rownames(expr_mat),
                              columns = c('ENTREZID'),
                              keytype = 'PROBEID')
anno <- anno[!is.na(anno$ENTREZID), ]
anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- summarise(anno_grouped, no_of_matches = n_distinct(ENTREZID))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
anno <- anno[!anno$PROBEID %in% anno_filtered$PROBEID, ]
length(unique(anno$ENTREZID))
expr_mat$PROBEID <- rownames(expr_mat)
expr_mat <- left_join(anno, expr_mat, by = "PROBEID")
rownames(expr_mat) <- expr_mat$PROBEID
expr_mat$PROBEID <- NULL
expr_mat$ENTREZID <- NULL

design <- model.matrix(~ 0 + f)
colnames(design) <- c('control', 'disease')
design

expr_fit <- lmFit(object = expr_mat, design = design)
expr_fit$coefficients[1:10, ]

ctst_mat <- makeContrasts(disease-control, levels = design)
expr_ctst_fit <- contrasts.fit(fit = expr_fit, contrasts = ctst_mat)
expr_fit_eb <- eBayes(fit = expr_ctst_fit)
names(expr_fit_eb)
expr_fit_eb$coefficients[1:10, ]
expr_fit_eb$p.value[1:10, ]
volcanoplot(expr_fit_eb)
de <- decideTests(object = expr_fit_eb, method = 'global', adjust.method = 'BH',
            p.value = 0.05, lfc = 1)
ups <- unique(anno$ENTREZID[anno$PROBEID %in% rownames(de[de == 1, ])])
downs <- unique(anno$ENTREZID[anno$PROBEID %in% rownames(de[de == -1, ])])
# x <- kegga(de = c(ups, downs), universe = unique(anno$ENTREZID), species = 'Hs')
# x$q <- p.adjust(x$P.DE, 'BH')
# x <- arrange(x, q)
# x[grep('cancer', x$Pathway), ]
# edo <- enrichDO(gene = ups, ont = 'DO', pvalueCutoff = 0.05, pAdjustMethod = 'BH',
#                 universe = unique(anno$ENTREZID), minGSSize = 5, maxGSSize = 500,
#                 qvalueCutoff = 0.2, readable = FALSE)
# head(edo)$Description
kk <- enrichKEGG(gene = ups,
           organism = 'hsa', 
           pvalueCutoff = 0.05)
kk@result[grepl('cancer', kk@result$Description) | grepl('carcinoma', kk@result$Description),
          c(2, 7)]
pData(ma_data[[9]])  

glm.control(maxit = 100)

glmFitting <- function(ExObj) {
  # Применяет линейную модель (lmFit из limma) к данным
  # Добавляет в ExObj саму модель
  f <- ExObj$ph == 'd'
  temat <- as_tibble(t(ExObj$emat))
  pvals <- sapply(temat, function(x) return(anova(glm.nb(f ~ x))$`Pr(>Chi)`[2]))
  names(pvals) <- ExObj$genes$ENTREZID
  return(p.adjust(pvals, 'BH'))
}

esetGFit <- function(eset, threshold = 0) {
  return(
    glmFitting(
      filteredExpressionObj(
        filterByIntensity(eset = eset, threshold = threshold)
      )
    )
  )
}

y <- esetGFit(norm_data$`Dataset #1`)
