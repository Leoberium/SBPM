# Загрузим файл с ридами, полученный на предыдущем этапе
raw_counts <- read.delim(file = 'counts.out', header = FALSE, comment.char = '#',
                         stringsAsFactors = FALSE)
samples <- c('B14.5', 'B15.5', 'B17.5', 'B20', 'B34',
             'C14.5', 'C15.5', 'C17.5', 'C20', 'C34')
age <- as.numeric(substring(samples, 2))
tissue <- as.factor(substring(samples, 1, 1))
colnames(raw_counts) <- c('gene', samples)

# Отфильтруем неэкспрессирующиеся гены
keep <- rowSums(raw_counts[, 2:11]) >= 1
counts <- raw_counts[keep, ]

# Загрузим библиотеку edgeR и создадим объект DGEList, нормализуем данные
library(edgeR)
rownames(counts) <- counts$gene
counts$gene <- NULL
er <- DGEList(counts = as.matrix(counts),
              group = tissue)
er <- calcNormFactors(er, method = 'RLE')

# Оценим дисперсионные параметры
formula_er <- ~ tissue + age
dm <- model.matrix(formula_er)
er <- estimateGLMCommonDisp(er, dm)
er <- estimateGLMTrendedDisp(er, dm)
er <- estimateGLMTagwiseDisp(er, dm)
strict.disp <- pmax(er$tagwise.dispersion,
                    er$trended.dispersion,
                    er$common.dispersion)
plotBCV(er)

# Поcтроим GLM-модель
gfit <- glmFit(er, dm, dispersion = strict.disp)
lrt_tissue <- glmLRT(gfit, 2)
lrt_age <- glmLRT(gfit, 3)

# Скорректированные p-value
library(ggplot2)
res_tissue <- topTags(lrt_tissue, n = Inf, adjust.method = 'BH', sort.by = 'none')
qplot(res_tissue$table$PValue, bins = 100,
      xlab = 'Adjusted P-value', ylab = 'Frequency',
      main = 'By tissue')
res_age <- topTags(lrt_age, n = Inf, adjust.method = 'BH', sort.by = 'none')
qplot(res_age$table$PValue, bins = 100,
      xlab = 'Adjusted P-value', ylab = 'Frequency',
      main = 'By age')

# Межтканевые различия не менее, чем в 2 раза
t2f <- res_tissue$table$logFC >= 1 | res_tissue$table$logFC <= -1
sum(t2f)  # 237
# Гены со значимым изменением между тканями
st <- t2f & (res_tissue$table$PValue < 0.05)
sum(st)  # 73 гена
rownames(res_tissue$table)[st]
# Гены со значимым изменением по возрасту
sa <- res_age$table$PValue < 0.05
sum(sa) # 311 генов
rownames(res_age$table)[sa]

# Скластеризуем гены значимые хотя бы по одному фактору
significant <- (st | sa)
sum(significant)  # 352 гена
norm_counts <- cpm(er) # нормализованные каунты
sig_norm_counts <- norm_counts[significant, ] # векторы экспрессии по значимым генам
corr_data <- cor(t(sig_norm_counts), method = 'spearman')
h_data <- hclust(as.dist(1 - corr_data))
cl <- cutree(h_data, k = 6)
plot(h_data, hang = -1, cex = 0.75, xlab = 'Gene')
rect.hclust(tree = h_data, k = 6, border = 1:6, cluster = cl)

# Z-score от возраста по каждому кластеру
z <- t(scale(t(sig_norm_counts)))
all(rownames(z) == names(cl))
sel <- tissue == 'B'
par(mfrow = c(3, 2))
legy <- c(0.75, 0, 1.0, 0.25, 1.0, 0)
for (i in 1:6) {
  v <- colMeans(z[cl == i, ])
  title_ <- paste('Cluster', i)
  plot(x = age[sel], y = v[sel], type='b', col = 'red',
       xlab = 'Age', ylab = 'Average Z-Score',
       main = title_, pch = 19,
       ylim = c(min(v), max(v)))
  lines(x = age[!sel], y = v[!sel], type = 'b', col = 'blue',
        pch = 19)
  legend(x = 30, y = legy[i], legend = c('B', 'C'), fill = c('red', 'blue'))
}
# Ткань B - красный, видимо это кора
# Ткань C - синий, видимо это мозжечок (cerebellum)

