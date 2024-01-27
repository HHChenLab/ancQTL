library(DESeq2)
library(data.table)
library(gap)
library(peer)

##
ind_IDs = read.table($ID_list, header=F, stringsAsFactors=F)
rnaseq = lapply(ind_IDs, function(x) fread( $ind.RNASEQC.gene_reads.gct, header=T, skip=2, stringsAsFactors=F, data.table=>
#geneinfo = read.table( paste0('/vgipiper04/CCHC/local_ancestry/replicate/rnaseqc/',mus[1,1],'.gene_reads.gct'), header=T, skip=2, stringsAsFactors=F)[,c(1,2)]

rnacount = do.call(ind_IDs, lapply(mus$RRID, function(x) fread(paste0(x,'RNASeQC.gene_reads.gct'), header=T, skip=2, stringsAsFactors=F, data.table=F)[,3]))
rownames(rnacount) = geneinfo$Name
colnames(rnacount) = ind_IDs

rnatpm = do.call( ind_IDs, lapply(mus$RRID, function(x) fread(paste0(x,'RNASeQC.gene_tpm.gct'), header=T, skip=2, stringsAsFactors=F, data.table=F)[,3]))
rownames(rnatpm) = geneinfo$Name
colnames(rnatpm) = ind_IDs

count_perc = apply(rnacount, MARGIN=1, function(x) sum(sapply(x, function(y) ifelse(y>=6, 1,0)))/473)
tpm_perc = apply(rnatpm, MARGIN=1, function(x) sum(sapply(x, function(y) ifelse(y>=0.1, 1,0)))/473)
length( intersect(names(count_perc)[which(count_perc>=0.2)], names(tpm_perc)[which(tpm_perc>=0.2)]))
##[1] 27071
passed_gene = intersect(names(count_perc)[which(count_perc>=0.2)], names(tpm_perc)[which(tpm_perc>=0.2)])
rnacount_passed = rnacount[passed_gene, ]
dim(rnacount_passed)
#[1] 27071   473
write.table(passed_gene, 'passed_gene.txt', col.names=F, row.names=F,  sep='\t', quote=F)
write.table(rnacount_passed, 'raw_readcounts.txt', col.names=T, row.names=T, sep='\t', quote=F)

deseq = estimateSizeFactorsForMatrix(rnacount_passed)
write.table(cchc_deseq, 'DEseq2_scale.txt', col.names=F, row.names=T, sep='\t', quote=F)

col = as.data.frame(cbind(rrid = colnames(rnacount), cond = 1))
RNA_deseq2 = DESeqDataSetFromMatrix(countData = rnacount_passed, colData = col, design = ~ 1)
RNA_deseq2 = estimateSizeFactors(RNA_deseq2)
RNA_deseq2 <- estimateDispersions(RNA_deseq2)
RNA_deseq_norm = counts(RNA_deseq2, normalized=T)

RNA_deseq_int = t(apply(RNA_deseq_norm, MARGIN=1, invnormal))

model = PEER()
PEER_setPhenoMean(model,as.matrix(t(RNA_deseq_int)))
PEER_setNk(model,60)
PEER_update(model)
##Converged (var(residuals)) after 81 iterations
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
rownames(factors) = rownames(RNA_deseq_int)
colnames(factors) = sapply(1:60, function(x) paste0('PEER',x))
colnames(weights) = sapply(1:60, function(x) paste0('PEER',x))
rownames(weights) = colnames(RNA_deseq_int)
rownames(residuals) = rownames(RNA_deseq_int)
colnames(residuals) = colnames(RNA_deseq_int)

write.table( factors, 'DEseq2_peer_factors.txt', col.names=T, row.names=T, sep='\t', quote=T)
write.table( precision, 'DEeq2_peer_precision.txt', col.names=T, row.names=T, sep='\t', quote=T)
write.table( weights, 'DEseq2_peer_weights.txt', col.names=T, row.names=T, sep='\t', quote=T)
write.table( residuals, 'DEseq2_peer_residuals.txt', col.names=T, row.names=T, sep='\t', quote=T)
