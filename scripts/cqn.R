setwd("~/OneDrive - Bayer/2025/publish/DE")
#install.packages("matrixStats")
#install.packages("BiocManager")
#BiocManager::install("DESeq2",update = F)
#install.packages("devtools")
#devtools::install_github("mhahsler/rBLAST",type = "source")
#BiocManager::install("cqn",update = F)
#BiocManager::install('EnhancedVolcano',update=F,dependencies=T)
library(cqn);library(EnhancedVolcano)
load("./Data_prep.RData")

# Subset WG, NC, SC batch 2 fams from experimental data and raw counts #####
expt_dat_b2 <- expt_dat[which(expt_dat$batch ==2 & expt_dat$State %in% c("WG","NC","SC")),]
expt_dat_b2$State <- as.factor(expt_dat_b2$State)
raw_counts_b2 <- raw_counts[,which(expt_dat$batch ==2 & expt_dat$State %in% c("WG","NC","SC"))]
head(rownames(raw_counts_b2))

raw_counts_b2 <- apply(raw_counts_b2,2,as.integer)
rownames(raw_counts_b2) <- rownames(raw_counts)

# Create design grouping for DESeq #####  
TIP_FAMS <- which(!expt_dat_b2$State %in% c("WG"))
des <- rep("WG",length(expt_dat_b2$State))
des[TIP_FAMS] <- "TIP"
expt_dat_b2$group <- as.factor(des)
rm(des)

expt_dat_b2 <- droplevels.data.frame(expt_dat_b2)

# Initial DESeq PCA plotting #####
dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts_b2,
                                      colData = expt_dat_b2,
                                      design = ~ group)
vs_log = DESeq2::vst(dds,nsub = 20000,blind = T,fitType = "parametric")

DESeq2::plotPCA(vs_log,intgroup = "group") # Note outliers bottom left
DESeq2::plotPCA(vs_log,intgroup = "fam_id") # Note outliers bottom left

# Idenitfy outliers and remove ####
cc <- DESeq2::plotPCA(vs_log,intgroup = "fam_id")
table(as.character(expt_dat_b2$fam_id))
cc$data$fam_id[which(cc$data[,2] < -20)]
table(cc$data$fam_id[which(cc$data[,1] > 0 & cc$data[,1] < 11)])
cc$data[,4][which(cc$data[,2] < -15)]
rownames(cc)[which(cc$data[,1] < -10)]
rm(cc,vs_log,TIP_FAMS)

# Fam N08061 is an outlier, remove this family #
raw_counts_b2 <- raw_counts_b2[,-which(expt_dat_b2$fam_id == "N08061")]
expt_dat_b2 <- expt_dat_b2[-which(expt_dat_b2$fam_id == "N08061"),]
expt_dat_b2 <- droplevels.data.frame(expt_dat_b2)
table(as.character(expt_dat_b2$fam_id))
dim(raw_counts_b2)

# Recreate DESeq PCA plots after dropping N08061 #####
dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts_b2,
                                      colData = expt_dat_b2,
                                      design = ~ group)
vs_log = DESeq2::vst(dds,nsub = 20000,blind = T,fitType = "parametric")
DESeq2::plotPCA(vs_log,intgroup = "group")
DESeq2::plotPCA(vs_log,intgroup = "fam_id")

# Load Evi 78k transcript lengths and gc content #####
pita_lengths <- read.table("./78.3txpts.len.gc.final")
head(pita_lengths)
identical(rownames(raw_counts_b2),pita_lengths[,1])
#TRUE
sum_cts <- apply(raw_counts_b2,2,sum) # get sum of counts for cqn #

# Run cqn and pass offset to deseq data #####
cqnObject <- cqn(counts = raw_counts_b2,lengths = pita_lengths$V2,x = pita_lengths$V3,sizeFactors = sum_cts,verbose=T)
#uncomment below for plotting of cqn object #
#par(mfrow=c(1,2))
#cqnplot(cqnObject, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
#cqnplot(cqnObject, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
cqnOffset <- cqnObject$glm.offset # retrieve offest
cqnNormFactors <- exp(cqnOffset) #exponentiate and apply as cqnNormFactors

dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts_b2,
                                      colData = expt_dat_b2,
                                      design = ~  group)

DESeq2::normalizationFactors(dds) <- cqnNormFactors


# Run DESeq2 and results ####
dds <- DESeq2::DESeq(dds,parallel = T,fitType = "parametric",minReplicatesForReplace = Inf)

res_cqn <- DESeq2::results(dds, contrast=c("group", "WG", "TIP"),alpha = .05,parallel = T,lfcThreshold = 0)
summary(res_cqn)

# Generate lfc shrinkage for volcano plot ####
res_cqn_lfc_shrink <- DESeq2::lfcShrink(dds,coef = "group_WG_vs_TIP", res=res_cqn,parallel = T,type = "normal")

EnhancedVolcano(res_cqn_lfc_shrink,labCol = NA,FCcutoff = 0,
                lab = rownames(res_cqn_lfc_shrink),
                x = 'log2FoldChange',
                y = 'pvalue',xlim = c(-3,3))
#plot(metadata(res_cqn)$filterNumRej, 
#     type="b", ylab="number of rejections",
#     xlab="quantiles of filter") + lines(metadata(res_cqn)$lo.fit, col="red") +abline(v=metadata(res_cqn)$filterTheta)

# Subset the upregeulated and downregulated transcripts ####
up_regulated_cqn <- res_cqn[which(res_cqn$padj <= .05 & (res_cqn$log2FoldChange) > 0),]
down_regulated_cqn <- res_cqn[which(res_cqn$padj <= .05 & (res_cqn$log2FoldChange) < 0),]
summary(up_regulated_cqn$log2FoldChange)
summary(down_regulated_cqn$log2FoldChange)

#sig_down_reg <- res_cqn[down_regulated_cqn,]
#sig_up_reg <- res_cqn[up_regulated_cqn,]
#summary(sig_up_reg$log2FoldChange)
#summary(sig_down_reg$log2FoldChange)

up_reg_all <- rownames(up_regulated_cqn)
down_reg_all <- rownames(down_regulated_cqn)


## Load all pita blast hits to arabidopsis ####
all_pita_blast <- read.table("./all_pita_blast.txt",header = T)
head(all_pita_blast)

all_pita_blast_arab <- (sub("\\..*", "",all_pita_blast$SubjectID)) # strip off the last two characters from arabidopsis transcripts to convert to genes #
length(c(table(all_pita_blast_arab)))
#17603 - the 78k transcripts mapped to 17603 unique arabidopsis genes #
arab_occurances <- (c(table(all_pita_blast_arab))) # <- Contains how many times a particular gene was found
#length(which(arab_occurances ==6))
hist(as.numeric(all_pita_blast$Perc.Ident)) # histogram of the identity found to each arabidopsis transcript
summary(as.numeric(all_pita_blast$Perc.Ident)) 
summary(as.numeric(all_pita_blast$E))
hist(as.numeric(all_pita_blast$E))

# Subset all of the pita blast hits to include only those with "low" E values ####
all_pita_blast_hq <- all_pita_blast[which(as.numeric(all_pita_blast$E) < 10e-5),]
dim(all_pita_blast_hq)
all_pita_blast_arab <- (sub("\\..*", "",all_pita_blast_hq$SubjectID))
length(c(table(all_pita_blast_arab)))
#12736 - the 57,996k transcripts mapped to 12736 unique arabidopsis genes given "low" E value#
arab_occurances <- (c(table(all_pita_blast_arab)))
summary(arab_occurances)
length(which(arab_occurances ==5))
summary(as.numeric(all_pita_blast_hq$Perc.Ident))
summary(as.numeric(all_pita_blast_hq$E))

### Load network connectivity #####
#matching pine to spruce is much easier than pine to arabidopis transcripts
#how many of the 54 pine transcripts showed up in the araibdopsis genes
NetworkConnectivity <- read_excel("NetworkConnectivity.xlsx", 
                                  sheet = "PineContigsSimilarToSpruceDEepi")
nrow(NetworkConnectivity)
in_up_nc <- up_reg_all[which(up_reg_all %in% NetworkConnectivity$Evi78k_PineContigID)]
in_down_nc <- down_reg_all[which(down_reg_all %in% NetworkConnectivity$Evi78k_PineContigID)]
length(c(in_up_nc,in_down_nc))/nrow(NetworkConnectivity)
in_up_down_nc_blast_hq <- all_pita_blast_hq[which(all_pita_blast_hq$QueryID %in% c(in_up_nc,in_down_nc)),]

# Identify which of the upregulated pine transcripts overlap with sig blast hits #####
up_reg_blast <- all_pita_blast[which(all_pita_blast$QueryID %in% up_reg_all),]
up_reg_blast <- up_reg_blast[which(as.numeric(up_reg_blast$E) < 10e-5),]
up_gene_hits <- sub("\\..*", "",up_reg_blast$SubjectID)
#length(which(c(table(up_gene_hits))== 3))
#c(table(up_gene_hits))[which.max(c(table(up_gene_hits)))]
up_gene_unique_hits <- unique(up_gene_hits)
length(up_gene_unique_hits)
# Identify which of the downregulated pine transcripts overlap with sig blast hits #####
down_reg_blast <- all_pita_blast[which(all_pita_blast$QueryID %in% down_reg_all),]
down_reg_blast <- down_reg_blast[which(as.numeric(down_reg_blast$E) < 10e-5),]
#dim(down_reg_blast)
down_gene_hits <- sub("\\..*", "",down_reg_blast$SubjectID)
#length(which(c(table(down_gene_hits))== 3))
#c(table(down_gene_hits))[which.max(c(table(down_gene_hits)))]
down_gene_unique_hits <- unique(down_gene_hits)
length(down_gene_unique_hits)

# Keep only unique arabidopsis genes that were in the unique up & down regulated set ####
u_genes <- names(which(c(table(c(up_gene_unique_hits,down_gene_unique_hits))) ==1))
down_gene_unique_single_hits <- down_gene_unique_hits[which(down_gene_unique_hits %in% u_genes)]
length(down_gene_unique_single_hits)

up_gene_unique_single_hits <- up_gene_unique_hits[which(up_gene_unique_hits %in% u_genes)]
length(up_gene_unique_single_hits)
#save.image("./Batch2_Results_analysis.RData",compress=T)
#write.table(as.matrix(down_gene_unique_single_hits),file = "./down_gene_hits.txt",col.names = F,row.names = F,quote = F)
#write.table(as.matrix(up_gene_unique_single_hits),file = "./up_gene_hits.txt",col.names = F,row.names = F,quote = F)
all_pita_blast_nc_arab <- (sub("\\..*", "",in_up_down_nc_blast_hq$SubjectID))
unique(all_pita_blast_nc_arab)
length(which(down_gene_unique_single_hits %in% unique(all_pita_blast_nc_arab)))
length(which(up_gene_unique_single_hits %in% unique(all_pita_blast_nc_arab)))
# Use ensemble for goseq ####
ensembl_arabidopsis <- biomaRt::useEnsemblGenomes(biomart="plants_mart",dataset="athaliana_eg_gene")
attributes = biomaRt::listAttributes(ensembl_arabidopsis)

gene_ids <- unlist(read.table("/mnt/aribdopsis/ath.geneIDs.txt")[,1])
test <- biomaRt::getBM(c("ensembl_gene_id","go_id"), filters="ensembl_gene_id",values = gene_ids,mart=ensembl_arabidopsis)
length(unique(test$ensembl_gene_id))
length(unique(gene_ids))
length(unique(test$go_id))
head(test)
dim(test)
cdna_lengths <- read.table("/mnt/aribdopsis/ath.gene.len.gc.txt")
library(tidyverse)
gene_median_length <- cdna_lengths %>%
  group_by(V2) %>%
  summarise(x = median(V3))
gene_names <- unlist(gene_median_length[,1])
gene_lengths <- unlist(gene_median_length[,2])
names(gene_lengths) <- gene_names


# create empty vector of 0's for all genes
ex_dat <- rep(0,length(gene_lengths))
# set the unique down regulated genes to 1 showing DE
ex_dat[which(gene_names %in% down_gene_regulated_hits)] <- 1
names(ex_dat) <- gene_names
# Run go_seq null p and label rownames as gene names
test_ex <- goseq::nullp(DEgenes = ex_dat,bias.data = gene_lengths)

# For each gene create a list of which GO categories they have
list_go <- test %>%
  group_by(ensembl_gene_id) %>%
  summarise(x = list(go_id))

# Convert this to a list
elist <- lapply(1:27655,function(x){
  as.vector(unlist(list_go[x,2]))
})
length(elist)
names(elist) <- unlist(list_go[,1])
rm(list_go)
GO.wall.down <- goseq::goseq(pwf = test_ex,gene2cat = elist,test.cats = "GO:BP",use_genes_without_cat = F,method = "Wallenius")
# conduct p-value adjustments
GO.wall.down$over_represented_pvalue_fdr <- p.adjust(GO.wall.down$over_represented_pvalue,method = "BH")
GO.wall.down$under_represented_pvalue_fdr <- p.adjust(GO.wall.down$under_represented_pvalue,
                                                      method="BH")
enriched.GO.down=GO.wall.down[which(GO.wall.down$over_represented_pvalue_fdr <.05 | GO.wall.down$under_represented_pvalue_fdr < .05),]
go.down.p =GO.wall.down[which(GO.wall.down$over_represented_pvalue <.05 | GO.wall.down$under_represented_pvalue < .05),]


# create empty vector of 0's for all genes
ex_dat <- rep(0,length(gene_lengths))

# set the unique down regulated genes to 1 showing DE
ex_dat[which(gene_names %in% up_gene_regulated_hits)] <- 1
names(ex_dat) <- gene_names
# Run go_seq null p and label rownames as gene names
test_ex <- goseq::nullp(DEgenes = ex_dat,bias.data = gene_lengths)

# For each gene create a list of which GO categories they have
list_go <- test %>%
  group_by(ensembl_gene_id) %>%
  summarise(x = list(go_id))

# Convert this to a list
elist <- lapply(1:27655,function(x){
  as.vector(unlist(list_go[x,2]))
})
length(elist)
names(elist) <- unlist(list_go[,1])
rm(list_go)

GO.wall.up <- goseq::goseq(pwf = test_ex,gene2cat = elist,test.cats = "GO:BP",use_genes_without_cat = F,method = "Wallenius")
# conduct p-value adjustments
GO.wall.up$over_represented_pvalue_fdr <- p.adjust(GO.wall.up$over_represented_pvalue,
                                                   method="BH")
GO.wall.up$under_represented_pvalue_fdr <- p.adjust(GO.wall.up$under_represented_pvalue,
                                                    method="BH")
enriched.GO.up=GO.wall.up[which(GO.wall.up$over_represented_pvalue_fdr <.05 | GO.wall.up$under_represented_pvalue_fdr < .05),]
go.up.p =GO.wall.up[which(GO.wall.up$over_represented_pvalue <.05 | GO.wall.up$under_represented_pvalue < .05),]

#####
enriched.GO.down=GO.wall.down[which(GO.wall.down$over_represented_pvalue_fdr <.05 ),]
dim(enriched.GO.down)
table(enriched.GO.down$ontology)
enriched.GO.up=GO.wall.up[which(GO.wall.up$over_represented_pvalue_fdr <.05),]
dim(enriched.GO.up)
table(enriched.GO.up$ontology)

table(c(table(c(enriched.GO.down$category,enriched.GO.up$category))))
enriched.GO.down[which(enriched.GO.down$category %in% enriched.GO.up$category),]

####
enriched.GO.down=GO.wall.down[which(GO.wall.down$over_represented_pvalue_fdr <.1 ),]
dim(enriched.GO.down)
table(enriched.GO.down$ontology)
enriched.GO.up=GO.wall.up[which(GO.wall.up$over_represented_pvalue_fdr <.1),]
dim(enriched.GO.up)
table(enriched.GO.up$ontology)

table(c(table(c(enriched.GO.down$category,enriched.GO.up$category))))
enriched.GO.down[which(enriched.GO.down$category %in% enriched.GO.up$category),]
