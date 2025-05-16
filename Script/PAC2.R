## ----error=FALSE,message=FALSE,warning=FALSE,include=FALSE,echo=FALSE------------------------------------

library(dplyr)
library(tidyverse)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
## Lectura de de la matriu de counts.
counts <- read.csv("Dades/GSE161731_counts.csv", row.names = 1)
## Lectura de la matriu de metadades. 
metadades <- read.csv("Dades/GSE161731_counts_Key.csv")

## S'han de canviar una mica els noms de les columnes de "counts" per a poder "emparellar les dades
colnames(counts) <-  gsub("\\.", "_", colnames(counts))
colnames(counts) <-  gsub("X", "", colnames(counts))

#De les metadades canviem els signes extranys per a "_"
metadades <- metadades %>%
  mutate(across(where(is.character), ~ str_replace_all(., "[\\s\\-/]", "_")))


## Si fem la intersecció entre les columnes de "counts" i rna_id de les metadades
## veiem que ens quadren els 198 individus que hi tocaria haver
## Ens quedem només amb les metadades i counts de les mostres en comú
## entre la matriu de expressió i les metadades

counts_metadades <- intersect(colnames(counts), metadades$rna_id)
metadades <- metadades[metadades$rna_id %in%counts_metadades,]
counts <- counts[,counts_metadades]

str(metadades)
#Abanas de crear el SummarizedExperiment el que farem serà definir bé les covariables
## la variable edat té un >89, que s'haurà de considerar NA, ja que no tenim l'edat i no la podem 
## tractar com a numèrica si posem >89. El Gender, Race, cohort, time_since_onset i hospitalized seran factors
metadades <- metadades %>%
  mutate(
    age = as.numeric(age),                  # edat com a numèrica
    gender = as.factor(gender),                   # sexe com a factor
    cohort = as.factor(cohort),
    race = as.factor(race),
    time_since_onset = as.factor(time_since_onset),
    hospitalized = as.factor(hospitalized),
    batch = as.factor(batch)
  )
# Converteix metadata en DataFrame per a SummarizedExperiment
colData <- as.data.frame(metadades)
rownames(colData) <- colData$rna_id


## Obtenim la base de dades Anotacions dels gens
gens <- genes(EnsDb.Hsapiens.v86, return.type = "GRanges")

# Alinea gens (files) entre matriu i anotació
## Mirem quins gens de la matriu "counts" trobem anotats en la base de dades
## de Ensembl i ens quedem només amb els que tinguin anotació (coordenades)

matriu_gens <- intersect(rownames(counts), names(gens))
counts <- counts[matriu_gens, ]
gens <- gens[matriu_gens]


## Creem l'objecte SummarizedExperiment
library(SummarizedExperiment)
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = gens,
  colData = colData
)

## Ara volem treballar només amb les cohorts COVID19, BActerial, i Healthy
coldata <- as.data.frame(colData(se))
coldata <- coldata %>%
  dplyr::filter(cohort %in% c("COVID_19", "Bacterial", "healthy"))

## Volem eliminar els individus repetits
coldata <- coldata %>%
  distinct(subject_id, .keep_all = TRUE)

## Treballem només amb el se que no té individus repetits i amb els filtres de cohort utilitzats
## De fet podríem ficar com a colnames de se el id del subjecte.
se <- se[, colnames(se) %in% coldata$rna_id]
nom_sub <- coldata$subject_id[match(colnames(se), coldata$rna_id)]
colnames(se) <- nom_sub


## Segons l'enunciat, seleccionem una mostra de 75 individus de manera aleatoria
myseed <- sum(utf8ToInt("joanmartinezsancho")) 
set.seed(myseed)
mostres_random <- sample(colnames(se), 75)
se_75 <- se[,mostres_random]




## ----echo=FALSE, error=FALSE,message=FALSE,warning=FALSE,include=FALSE-----------------------------------
library(edgeR)
## Expressem els contatges com CPMs es a dir, counts per milió. 
## Podem guardar els contatjes "bruts" per si mai es volen recuperar
counts.CPM <- cpm(assay(se_75))

assays(se_75)$countCPM <- counts.CPM


## Un cop tenim les dades com CPMs podem filtrar. Com a mètode de filtratge utilitzarel el filtre que 
## realitz l'article de referència que hem escollit per obtenir les dades
## Genes with counts per million greater than 1 in fewer than 20% of samples were dropped along with three samples with a high proportion of lowly expressed reads.
## És a dir, s'eliminem els gens amb recomptes per milló superior a 1 en menys del 20% de les mostres.
## Per tant, ens quedem amb aquells gens que tenen recomptes per milió sperior a 1 en almenys el 20% de mostres. Per tant, en el nostre cas ens quedaríem amb aquells gens que tenen un recompte per milió superior a 1 en almenys 15 (0.2*75mostres) mostres

gens_keep <- rowSums(counts.CPM>1) >= (0.2*ncol(se_75))

# Eliminar 3 mostres amb més baixa expressió, els quals són les mostres amb id 896282, 5835FC i 9C7138
low_expr_samples <- order(colSums(counts.CPM>1))[1:3]
#colSums(counts.CPM>1)[c(20,43,54)]

## Tenim el nou SummarizedExperiment amb les dades filtrades
se_75_filtrat <- se_75[gens_keep,-low_expr_samples]


## Normalització
## Per fer-ho primer guardarem les dades de recomptes a un objecte dgeList
dgeObj <- DGEList(counts = assays(se_75_filtrat)$countCPM, samples = colData(se_75_filtrat), genes = rowData(se_75_filtrat),group = se_75_filtrat$cohort,lib.size = colSums(assays(se_75_filtrat)$countCPM), remove.zeros = FALSE, norm.factors = rep(1,ncol(assays(se_75_filtrat)$countCPM)))

#show(dgeObj)

##Además de estandarizar los contajes, es importante eliminar otros sesgos de composición entre librerías. Esto puede hacerse aplicando la normalización por el método TMM que genera un conjunto de factores de normalización, tal que producto de estos factores y los tamaños de librería (el número de secuencias de cada muestra) definen el tamaño efectivo de dichas muestras, es decir el peso real que se les asignará en las comparaciones posteriores.
## Així com indica l'apartat de mètodes de l'article de referència realitzarem una normalització pel mètode TMM i llavors aplicarem el logaritme en base 2

dgeObj_norm <- calcNormFactors(dgeObj, method = "TMM")
head(dgeObj_norm$samples, 10)

## Aquestes transformacions cerquen compensar la mida diferent de les llibreríes o la diferent composició
## d'aquestes, pero les distribucions en cada mostra són asimètriques

boxplot(dgeObj_norm$counts, col = dgeObj_norm$samples$group, las = 2, cex.axis = 0.7,
    main = "Contajes normalizados",ylim = c(0, 10000))

log2count_norm <- cpm(dgeObj_norm, log = TRUE)

orden <- order(dgeObj_norm$samples$group)
log2count_ordenat <- log2count_norm[, orden]
cols_ordenats <- dgeObj_norm$samples$group[orden]
grups_ordenats <- dgeObj_norm$samples$group[orden]
##grups_ordenats <- dgeObj_norm$samples$subject_id[orden]
boxplot(log2count_ordenat,
        col = cols_ordenats,
        las = 2,
        cex.axis = 0.7,
        main = "Contajes normalizados (log2)",
        names = grups_ordenats)

assays(se_75_filtrat)$log2count_nom <- log2count_norm

ordered_sample_names <- rep(colnames(log2count_ordenat), each = nrow(log2count_ordenat))
group_order <- rep(c(rep("Bacterial", 19), rep("COVID_19", 40), rep("Healthy", 13)),
                   each = nrow(log2count_norm))

d <- tibble(
  y = as.numeric(log2count_ordenat),
  x = factor(rep(colnames(log2count_ordenat), each = nrow(log2count_ordenat))),
  co = factor(c(rep("Bacterial", each = nrow(log2count_norm) * 19),
                rep("COVID_19", each = nrow(log2count_norm) * 40),
                rep("Healthy", each = nrow(log2count_norm) * 13)))
)

# Reordena x segons co
d$x <- factor(d$x, levels = unique(d$x[order(d$co)]))

gg1 <- ggplot(d,aes(x,y,fill=co)) + geom_boxplot() +
  xlab('') + ylab(expression(bold('Normalized'~'counts'~(log[2])))) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.5, size=10),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=15),
        axis.title = element_text(face='bold'),
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        legend.position = 'bottom',
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"))+
  ggtitle("Comptatges Normalitzats i transformats")



## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----

boxplot(dgeObj_norm$counts, col = dgeObj_norm$samples$group, las = 2, cex.axis = 0.7,
    main = "Contatges Normalitzats sense transformació",ylim = c(0, 10000))

gg1


## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----
## Per poder visualitzar millor, el que podem fer es definir el nom de les mostres com
## COVID19_1, COVID19_2, etc, Bacterial_1, Bacterial_2, etc...
cohorts <- as.character(colData(se_75_filtrat)$cohort)
# Crear un vector nuu amb números únics enumerats per cohort
# Agrupar per cohort i enumerar cada element
cohort_counts <- ave(cohorts, cohorts, FUN = seq_along)
new_names <- paste0(cohorts, "_", cohort_counts)

colnames(se_75_filtrat) <- new_names

sampleDists <- dist(t(assays(se_75_filtrat)$log2count_nom))
#round(sampleDists, 1)

#Visualització heatmap
library(factoextra)
fviz_dist(sampleDists,
          show_labels = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6))


## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----

library(ggdendro)
library(dendextend)

## Calculem l'agrupament jeràrquic
hc <- hclust(sampleDists,method = 'ward.D2')

dend <- as.dendrogram(hc)

cohort_colors <- c(
  "Bacterial" = "#E78181",  # vermell suau
  "COVID_19"  = "#4DAF4A",  # verd mitjà
  "healthy"   = "#80B1D3"   # blau clar
)
# Assignar colors per cohort
ordered_cohorts <- cohorts[match(labels(dend), colnames(se_75_filtrat))]

label_colors <- scales::col_factor(palette = "Dark2", domain = unique(ordered_cohorts))(ordered_cohorts)

# 5. Assignem els colors als textos del dendrograma
dend <- dend %>%
  set("labels_col", cohort_colors[cohorts[match(labels(dend), new_names)]]) %>%
  set("labels_cex", 0.6)

# 6. Dibuixem el dendrograma
plot(dend, main = "Dendrograma amb etiquetes per cohort")
legend("topright", legend = unique(ordered_cohorts),
       fill = unique(cohort_colors), border = NA, bty = "n", cex = 0.8)


## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10,include=FALSE----

library(limma)
library(RColorBrewer)

mds <- limma::plotMDS(assays(se_75_filtrat)$log2count_nom, main = "Status", cex = 0.7)



## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----


# 4. Fer el plot amb colors per grup
plot(mds$x, mds$y,
     col = cohort_colors[cohorts],
     pch = 16,
     cex = 1.3,
     xlab = paste("Dimensió 1 (", round(mds$var.explained[1]*100, 1), "%)", sep = ""),
     ylab = paste("Dimensió 2 (", round(mds$var.explained[2]*100, 1), "%)", sep = ""),
     main = "MDS plot per cohort")
legend("topright", legend = names(cohort_colors),
       col = cohort_colors, pch = 16, cex = 0.8)
text(mds$x, mds$y,
     labels = colnames(se_75_filtrat),
     pos = 4,  # posició: 3 = sobre del punt
     cex = 0.6,
     col = cohort_colors[cohorts])


## ----error=FALSE,message=FALSE,warning=FALSE, include=FALSE----------------------------------------------

mostres_eliminar <- c("Bacterial_2", "Bacterial_5", "Bacterial_7", "COVID_19_2", "COVID_19_3", "COVID_19_29")
noms_mostres_eliminar <- colData(se_75_filtrat)$subject_id[which(colnames(se_75_filtrat) %in% mostres_eliminar)]
se_75_filtrat <- se_75_filtrat[, !colnames(se_75_filtrat) %in% mostres_eliminar]



## ----echo=FALSE,error=FALSE,message=FALSE,warning=FALSE--------------------------------------------------
library(compareGroups)
t <- compareGroups(cohort~age+gender+race+time_since_onset+hospitalized,data=as.data.frame(colData(se_75_filtrat)),method=NA)
t <- createTable(t)
export2md(t,caption = "Comparació Bivariada")




## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----


## Treballarem amb l'objecte dgeObje_norm que havíem creat al principi, per fer-ho haurem
## d'eliminar les 6 mostres que hem eliminat del objecte summarizedExperiment
dgeObj_norm <- dgeObj_norm[, !(colnames(dgeObj_norm) %in% noms_mostres_eliminar)]


group <- relevel(factor(colData(se_75_filtrat)$cohort), ref = "healthy")
age <- colData(se_75_filtrat)$age
gender <- factor(colData(se_75_filtrat)$gender)

## Creem la matriu de disseny
design <- model.matrix(~ 0+group + age + gender)
#colnames(design)
head(design)
## Estem interesants en les diferències entre els grups, necessitem especificar quines comparacions volem dur a terme.
cont.matrix <- makeContrasts(BacterialVsHealthy=groupBacterial - grouphealthy, 
                             COVID19VsHealthy=groupCOVID_19 - grouphealthy, levels=design)
cont.matrix

## voom+lima
library(edgeR)
library(limma)
set.seed(myseed)
library(kableExtra)
sample(c("edgeR", "voom+limma", "DESeq2"), size = 1)
## Ens ha tocat el mètode voom+lima
## Transformem les dades amb voom
voomObj <- voom(dgeObj_norm, design)

## Realitzem el model amb la matriu de disseny, ens quedarà ajustat per edat i sexe
## A més acabem el procés amb la regularització del estimador del error utilitzant la funció
## eBayes
fit <- lmFit(voomObj)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## Per Bacterial vs Healty
res_bact <- topTable(fit2,coef=1,sort.by="p", number=nrow(fit2))
res_bact$significant <- with(res_bact, abs(logFC) > 1.5 & adj.P.Val < 0.05)
## Com a threshold posem el límit en el pvalo <0.05 i abs(logFC) > 1.5
##head(res_bact)

# Total de gens analitzats
total_genes <- nrow(res_bact)

# Gens diferencialment expressats
num_DEGs <- sum(res_bact$significant)

# Sobreexpressats en Bacterial (logFC > 1.5)
num_upregulated <- sum(res_bact$logFC > 1.5 & res_bact$adj.P.Val < 0.05)

# Infraexpressats en Bacterial (logFC < -1.5)
num_downregulated <- sum(res_bact$logFC < -1.5 & res_bact$adj.P.Val < 0.05)

# Gens no diferencialment expressats
num_non_DEGs <- total_genes - num_DEGs

# Mostrar taula resum
summary_table <- data.frame(
  `Estat dels gens` = c(
    "Gens diferencialment expressats (|logFC| > 1.5 & FDR < 0.05)",
    "  - Sobreexpressats en Bacterial (logFC > 1.5)",
    "  - Infraexpressats en Bacterial (logFC < -1.5)",
    "Gens no diferencialment expressats",
    "Total de gens analitzats"
  ),
  `Nombre de gens` = c(
    num_DEGs,
    num_upregulated,
    num_downregulated,
    num_non_DEGs,
    total_genes
  )
)

kable(summary_table, caption = "Resum dels gens diferencialment expressats en la comparació Bacterial vs Healthy") %>%
  kable_styling(latex_options = c("striped", "hold_position"))


library(ggrepel)
ggplot(res_bact, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(res_bact, significant),
                  aes(label = symbol), size = 3, max.overlaps = 15) +
  labs(title = "Volcano Plot: Bacterial vs Healthy",
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~FDR)) +
  theme_minimal() +
  theme(legend.position = "none")




## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----


## COVID19 vs Healty
res_covid <- topTable(fit2,coef=2,sort.by="p", number=nrow(fit2))
res_covid$significant <- with(res_covid, abs(logFC) > 1.5 & adj.P.Val < 0.05)

# Total de gens analitzats
total_genes <- nrow(res_covid)

# Gens diferencialment expressats
num_DEGs <- sum(res_covid$significant)

# Sobreexpressats en Bacterial (logFC > 1.5)
num_upregulated <- sum(res_covid$logFC > 1.5 & res_covid$adj.P.Val < 0.05)

# Infraexpressats en Bacterial (logFC < -1.5)
num_downregulated <- sum(res_covid$logFC < -1.5 & res_covid$adj.P.Val < 0.05)

# Gens no diferencialment expressats
num_non_DEGs <- total_genes - num_DEGs

# Mostrar taula resum
summary_table <- data.frame(
  `Estat dels gens` = c(
    "Gens diferencialment expressats (|logFC| > 1.5 & FDR < 0.05)",
    "  - Sobreexpressats en COVID19 (logFC > 1.5)",
    "  - Infraexpressats en COVID19 (logFC < -1.5)",
    "Gens no diferencialment expressats",
    "Total de gens analitzats"
  ),
  `Nombre de gens` = c(
    num_DEGs,
    num_upregulated,
    num_downregulated,
    num_non_DEGs,
    total_genes
  )
)

kable(summary_table, caption = "Resum dels gens diferencialment expressats en la comparació COVID19 vs Healthy") %>%
  kable_styling(latex_options = c("striped", "hold_position"))

ggplot(res_covid, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(res_covid, significant),
                  aes(label = symbol), size = 3, max.overlaps = 15) +
  labs(title = "Volcano Plot: COVID19 vs Healthy",
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~FDR)) +
  theme_minimal()+
  theme(legend.position = "none")



## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----

##Cream una taula similar a la que crea decideTests però amb els pvalors ajustats

fdr_cutoff <- 0.05
logFC_cutoff <- 1.5
result_matrix <- sapply(colnames(fit2$coefficients), function(contrast) {
  tt <- topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
  sign <- ifelse(tt$adj.P.Val < fdr_cutoff & tt$logFC > logFC_cutoff, 1,
                 ifelse(tt$adj.P.Val < fdr_cutoff & tt$logFC < -logFC_cutoff, -1, 0))
  return(sign)
})
## Aquesta taula ens dona 1 si el gene esta sobreespressat, 0 si no es significatiu i -1 si està 
## Infraexpressat
rownames(result_matrix) <- rownames(fit2$coefficients)
result_matrix <- as.data.frame(result_matrix)

summa.fit <- data.frame(BacterialvsHealthy=c(table(result_matrix$BacterialVsHealthy)[names(table(result_matrix$BacterialVsHealthy))==-1],table(result_matrix$BacterialVsHealthy)[names(table(result_matrix$BacterialVsHealthy))==0],table(result_matrix$BacterialVsHealthy)[names(table(result_matrix$BacterialVsHealthy))==1]),COVID19VsHealthy=c(table(result_matrix$COVID19VsHealthy)[names(table(result_matrix$COVID19VsHealthy))==-1],table(result_matrix$COVID19VsHealthy)[names(table(result_matrix$COVID19VsHealthy))==0],table(result_matrix$COVID19VsHealthy)[names(table(result_matrix$COVID19VsHealthy))==1]))


vc<- vennCounts(result_matrix)
vennDiagram(vc, include=c("Up", "Down"),
    counts.col=c("red", "blue"),
    circle.col = c("red", "blue", "green3"), cex=c(1,1,1))





## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")

res_covid <- topTable(fit2,coef=2,sort.by="p", number=nrow(fit2))


library(clusterProfiler)
library(org.Hs.eg.db)
topTab<- res_covid

allEntrezs <- as.vector(as.character(topTab$entrezid))
selectedEntrezsUP <- as.vector(as.character(subset(topTab, (logFC> 1.5) & (adj.P.Val < 0.05))$entrezid))
#length(allEntrezs)
#length(selectedEntrezsUP)

ego_up <- enrichGO(gene         = selectedEntrezsUP,
                   universe     = allEntrezs,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05,
                   readable     = TRUE)

#head(ego_up)
ego_results <- data.frame(ego_up)
write.csv(ego_results, "clusterProfiler_ORAresults_UpGO.csv")


## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----
head(ego_up)
dotplot(ego_up, showCategory=9, font.size=6)


## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----

cnetplot(ego_up, font.size=1) 



## ----echo=FALSE, error=FALSE,message=FALSE,fig.align='center',warning=FALSE,fig.height=5,fig.width=10----
library(edgeR)

# Eliminar les mostres atípiques
dgeObj_edge <- dgeObj_norm[, !(colnames(dgeObj_norm) %in% noms_mostres_eliminar)]

# Factors d’interès
group <- relevel(factor(colData(se_75_filtrat)$cohort), ref = "healthy")
age <- colData(se_75_filtrat)$age
gender <- factor(colData(se_75_filtrat)$gender)

# Matriu de disseny (ajustat per edat i sexe)
design <- model.matrix(~ 0 + group + age + gender)
colnames(design) <- make.names(colnames(design))

# Crear objecte DGEList amb la informació de disseny
dge <- estimateDisp(dgeObj_edge, design)

# Ajust del model amb GLM
fit <- glmQLFit(dge, design)

# Definició dels contrastos
cont.matrix <- makeContrasts(
  BacterialVsHealthy = groupBacterial - grouphealthy,
  COVID19VsHealthy = groupCOVID_19 - grouphealthy,
  levels = design
)

# Test de quasi-likelihood per Bacterial vs Healthy
qlf_bact <- glmQLFTest(fit, contrast = cont.matrix[, "BacterialVsHealthy"])

# Resultats ordenats i anotació de significatius
res_bact_edge <- topTags(qlf_bact, n = nrow(dge$counts))$table
res_bact_edge <- subset(res_bact_edge, abs(logFC) > 1.5 & FDR < 0.05)

res_bact <- topTable(fit2,coef=1,sort.by="p", number=nrow(fit2))
res_bact <- subset(res_bact, abs(logFC) > 1.5 & adj.P.Val < 0.05)


topGenes_voom <- rownames(res_bact) 
topGenes_edge <- rownames(res_bact_edge) 
library(ggvenn)
x = list(LimmaVoom = topGenes_voom, edgeR = topGenes_edge)
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"), stroke_size = 0.5, set_name_size = 3)+
  ggtitle("Comparacio Gens Diferencialment Expressats Voom i Edge Bacterial")


## ----error=FALSE,message=FALSE,warning=FALSE, include=FALSE----------------------------------------------
library(edgeR)

# Eliminar les mostres atípiques
dgeObj_edge <- dgeObj_norm[, !(colnames(dgeObj_norm) %in% noms_mostres_eliminar)]

# Factors d’interès
group <- relevel(factor(colData(se_75_filtrat)$cohort), ref = "healthy")
age <- colData(se_75_filtrat)$age
gender <- factor(colData(se_75_filtrat)$gender)

# Matriu de disseny (ajustat per edat i sexe)
design <- model.matrix(~ 0 + group + age + gender)
colnames(design) <- make.names(colnames(design))

# Crear objecte DGEList amb la informació de disseny
dge <- estimateDisp(dgeObj_edge, design)

# Ajust del model amb GLM
fit <- glmQLFit(dge, design)

# Definició dels contrastos
cont.matrix <- makeContrasts(
  BacterialVsHealthy = groupBacterial - grouphealthy,
  COVID19VsHealthy = groupCOVID_19 - grouphealthy,
  levels = design
)

# Test de quasi-likelihood per Bacterial vs Healthy
qlf_bact <- glmQLFTest(fit, contrast = cont.matrix[, "COVID19VsHealthy"])

# Resultats ordenats i anotació de significatius
res_cov_edge <- topTags(qlf_bact, n = nrow(dge$counts))$table
res_cov_edge <- subset(res_cov_edge, (logFC> 1.5) & (FDR < 0.05))

res_covid <- topTable(fit2,coef=2,sort.by="p", number=nrow(fit2))
res_covid <- subset(res_covid, abs(logFC) > 1.5 & adj.P.Val < 0.05)


topGenes_voom <- rownames(res_covid) 
topGenes_edge <- rownames(res_cov_edge) 
library(ggvenn)
x = list(LimmaVoom = topGenes_voom, edgeR = topGenes_edge)
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"), stroke_size = 0.5, set_name_size = 3)+
  ggtitle("Comparacio Gens Diferencialment Expressats Voom i Edge COVID")


## ----echo=FALSE,eval=T,message=FALSE,error=FALSE,warning=FALSE-------------------------------------------
knitr::purl("PAC2.Rmd")

