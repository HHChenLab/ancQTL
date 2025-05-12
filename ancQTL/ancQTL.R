library(mixqtl)
library(data.table)
library(metap)
library(multcomp)

ancqtl = function(geno1, geno2, anc1, anc2, y1, y2, ytotal, lib_size, cov_offset = NULL,passed_n = 30, trc_cutoff = 20, asc_cutoff = 5, weight_cap = 100, asc_cap = 5000) {
	if(is.null(cov_offset)) {
		cov_offset = rep(0, length(lib_size))
	}
	# prepare X
	com_snp = (apply(rbind(geno1, geno2), MARGIN=2, mean)>0.05 & apply(rbind(geno1, geno2), MARGIN=2, mean)<0.95 & apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)>0.05 & apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)<0.95 & apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)>0.05 & apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)<0.95)
	h1 = geno1[,com_snp]
	h1[is.na(h1)] = 0.5
	h2 = geno2[,com_snp]
	h2[is.na(h2)] = 0.5
	Xasc = h1 - h2
	Xtrc = (h1 + h2) / 2
	Aasc = anc1 - anc2
	Atrc = (anc1+anc2)/20
	XAasc = anc1*h1 - anc2*h2
	XAtrc = (anc1*h1 + anc2*h2)/2
	passed_ind = y1 >= asc_cutoff & y2 >= asc_cutoff & y1 <= asc_cap & y2 <= asc_cap & !is.na(Aasc)
#	print(table(com_snp))
	ancqtl_out = sapply( colnames(h1), function(i){
	asc = ancqtl_asc(y1, y2, Xasc[,i], Aasc, XAasc[,i], asc_cutoff = asc_cutoff, weight_cap = weight_cap, asc_cap = asc_cap)
	trc = ancqtl_trc(ytotal, lib_size, Xtrc[,i], Atrc, XAtrc[,i], cov_offset, trc_cutoff = trc_cutoff)
	c(asc, trc)
	})
	ancqtl_out = as.data.frame(t(ancqtl_out))
	ancqtl_out$anc1.n = length(which(anc1==0)) + length(which(anc2==0))
	ancqtl_out$anc1.af = apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)[com_snp]
	ancqtl_out$anc2.n = length(which(anc1==1)) + length(which(anc2==1))
	ancqtl_out$anc2.af = apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)[com_snp]
	na_ind = is.na(ytotal) | ytotal < trc_cutoff | is.na(as.numeric(Atrc))
	ancqtl_out$trc.anc1.n = length(which(anc1==0 & !na_ind)) + length(which(anc2==0 & !na_ind))
	ancqtl_out$trc.anc1.af = apply(rbind(geno1[which(anc1==0 & !na_ind),], geno2[which(anc2==0 & !na_ind),]), MARGIN=2, mean)[com_snp]
	ancqtl_out$trc.anc2.n = length(which(anc1==1 & !na_ind)) + length(which(anc2==1 & !na_ind))
	ancqtl_out$trc.anc2.af = apply(rbind(geno1[which(anc1==1 & !na_ind),], geno2[which(anc2==1 & !na_ind),]), MARGIN=2, mean)[com_snp]
	passed_ind = y1 >= asc_cutoff & y2 >= asc_cutoff & y1 <= asc_cap & y2 <= asc_cap & !is.na(Aasc)
	if(table(passed_ind)['TRUE']>passed_n){
	geno1 = geno1[passed_ind,com_snp]
	geno2 = geno2[passed_ind,com_snp]
	anc1 = anc1[passed_ind]
	anc2 = anc2[passed_ind]
	ancqtl_out$asc.anc1.n = sum(anc1==0) + sum(anc2==0)
	ancqtl_out$asc.anc1.af = apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)
	ancqtl_out$asc.anc2.n = sum(anc1==1) + sum(anc2==1)
	ancqtl_out$asc.anc2.af = apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)
	} else{
	ancqtl_out$asc.anc1.n = NA
	ancqtl_out$asc.anc1.af = NA
	ancqtl_out$asc.anc2.n = NA
	ancqtl_out$asc.anc2.af = NA
	}
	rownames(ancqtl_out) = colnames(h1)
#	recom_snp = (apply(rbind(geno1, geno2), MARGIN=2, mean)>0.05 & apply(rbind(geno1, geno2), MARGIN=2, mean)<0.95 & apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)>0.05 & apply(rbind(geno1[which(anc1==0),], geno2[which(anc2==0),]), MARGIN=2, mean)<0.95 & apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)>0.05 & apply(rbind(geno1[which(anc1==1),], geno2[which(anc2==1),]), MARGIN=2, mean)<0.95)
#	ancqtl_out = ancqtl_out[recom_snp,]
	anc1_meta = t(apply(ancqtl_out[, c('trc.snp.beta','asc.snp.beta','trc.snp.se','asc.snp.se')], MARGIN =1 , function(o) mix_meta(o[1], o[2], o[3], o[4])))
	colnames(anc1_meta) = c('meta.anc1.beta','meta.anc1.se','meta.anc1.z','meta.anc1.p')
	anc2_meta = t(apply(ancqtl_out[, c('trc.anc2.beta','asc.anc2.beta','trc.anc2.se','asc.anc2.se')], MARGIN =1 , function(o) mix_meta(o[1], o[2], o[3],o[4])))
	colnames(anc2_meta) = c('meta.anc2.beta','meta.anc2.se','meta.anc2.z','meta.anc2.p')
	ancqtl_out = cbind(ancqtl_out, anc1_meta)
	ancqtl_out = cbind(ancqtl_out, anc2_meta)
	anctest = t(apply(ancqtl_out[,c('asc.anc1.n','trc.anc1.n','meta.anc1.beta','meta.anc1.se','asc.anc2.n','trc.anc2.n','meta.anc2.beta','meta.anc2.se')], MARGIN = 1, function(t) anc_t(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8])))
	colnames(anctest) = c('meta.ANC.t','meta.ANC.p')
	ancqtl_out = cbind(ancqtl_out, anctest)
#	print(ancqtl_out[1:2,])
	ancqtl_out$meta.model.p = apply( ancqtl_out[,c('trc.model.p','asc.model.p')], MARGIN = 1, function(p) ifelse( any(is.na(p)), NA, invchisq(as.numeric(p), c(2,2))$p))
	ancqtl_out$asc.anc1.fdr = p.adjust(ancqtl_out$asc.anc1.p, method = 'BH')
	ancqtl_out$asc.anc2.fdr = p.adjust(ancqtl_out$asc.anc2.p, method = 'BH')
	ancqtl_out$asc.model.fdr = p.adjust(ancqtl_out$asc.model.p, method = 'BH')
	ancqtl_out$asc.ANC.fdr = p.adjust(ancqtl_out$asc.ANC.p, method = 'BH')
	ancqtl_out$trc.anc1.fdr = p.adjust(ancqtl_out$trc.anc1.p, method = 'BH')
	ancqtl_out$trc.anc2.fdr = p.adjust(ancqtl_out$trc.anc2.p, method = 'BH')
	ancqtl_out$trc.model.fdr = p.adjust(ancqtl_out$trc.model.p, method = 'BH')
	ancqtl_out$trc.ANC.fdr = p.adjust(ancqtl_out$trc.ANC.p, method = 'BH')
	ancqtl_out$meta.anc1.fdr = p.adjust(ancqtl_out$meta.anc1.p, method='BH')
	ancqtl_out$meta.anc2.fdr = p.adjust(ancqtl_out$meta.anc2.p, method='BH')
	ancqtl_out$meta.ANC.fdr = p.adjust(ancqtl_out$meta.ANC.p, method='BH')
	ancqtl_out$meta.model.fdr = p.adjust(ancqtl_out$meta.model.p, method='BH')
#	return(ancqtl_out)
#	print(ancqtl_out[1:2,])
#	print(dim(ancqtl_out))
	return(ancqtl_out[,c('anc1.n','anc1.af','anc2.n','anc2.af','asc.n','asc.anc1.n','asc.anc1.af','asc.anc2.n','asc.anc2.af','asc.snp.beta','asc.snp.se','asc.snp.t','asc.snp.p','asc.anc.beta','asc.anc.se','asc.anc.t','asc.anc.p','asc.ancsnp.beta','asc.ancsnp.se','asc.ancsnp.t','asc.ancsnp.p','asc.anc1.beta','asc.anc1.se','asc.anc1.t','asc.anc1.p','asc.anc1.fdr','asc.anc2.beta','asc.anc2.se','asc.anc2.t','asc.anc2.p','asc.anc2.fdr','asc.ANC.t','asc.ANC.p','asc.ANC.fdr','asc.model.p','asc.model.fdr','trc.n','trc.anc1.n','trc.anc1.af','trc.anc2.n','trc.anc2.af','trc.snp.beta','trc.snp.se','trc.snp.t','trc.snp.p','trc.anc.beta','trc.anc.se','trc.anc.t','trc.anc.p','trc.ancsnp.beta','trc.ancsnp.se','trc.ancsnp.t','trc.ancsnp.p','trc.anc1.beta','trc.anc1.se','trc.anc1.t','trc.anc1.p','trc.anc1.fdr','trc.anc2.beta','trc.anc2.se','trc.anc2.t','trc.anc2.p','trc.anc2.fdr','trc.ANC.t','trc.ANC.p','trc.ANC.fdr','trc.model.p','trc.model.fdr','meta.anc1.beta','meta.anc1.se','meta.anc1.z','meta.anc1.p','meta.anc1.fdr','meta.anc2.beta','meta.anc2.se','meta.anc2.z','meta.anc2.p','meta.anc2.fdr','meta.ANC.t','meta.ANC.p','meta.ANC.fdr','meta.model.p','meta.model.fdr')])
}

mix_meta = function(b1, b2, s1, s2) {
	if(any(is.na(c(b1,b2)))){
	return(c(NA,NA,NA,NA))
	} else {
	w1 = 1 / s1^2
	w2 = 1 / s2^2
	b = (w1 * b1 + w2 * b2) / (w1 + w2)
	s = sqrt(1 / (w1 + w2))
	z = b / s
	p = exp(pnorm(abs(z), log.p=T, lower.tail=F)) * 2
	return(c(b, s, z, p))
	}
}

anc_t = function(n11, n12, m1, s1, n21, n22, m2, s2){
	if(any(is.na(c(m1, m2)))){
		return(c(NA, NA))
	} else {
	tryCatch( {
	n1 = min(n11, n12, na.rm=T)
	n2 = min(n21, n22, na.rm=T)
	anct = (((m1 - m2))/sqrt((s1^2) + (s2^2)))
	ancp = 2*pt( abs(anct), (n1 + n2 - 2), lower.tail =F)
	return(c(t.stat = anct, pval = ancp))
	},
	error = function(e) {
	return(c(NA, NA))
	})
}}

ancqtl_asc = function(asc1, asc2, x, anc, xanc, passed_n = 20, asc_cutoff = 5, weight_cap = 100, asc_cap = 5000) {
	passed_ind = asc1 >= asc_cutoff & asc2 >= asc_cutoff & asc1 <= asc_cap & asc2 <= asc_cap & !is.na(anc)
	x = x[passed_ind]
	anc = anc[passed_ind]
	xanc = xanc[passed_ind]
	asc1 = asc1[passed_ind]
	asc2 = asc2[passed_ind]
	asc = log(asc1 / asc2)
	sample_size = table(passed_ind)['TRUE']
	if(table(passed_ind)['TRUE'] < passed_n){
#		print(paste0('no enough N for asc model: n=', sample_size))
		ascout = c(rep(NA, times = 21), sample_size, NA, NA)
		names(ascout) = c('asc.snp.beta','asc.snp.se','asc.snp.t','asc.snp.p', 'asc.anc.beta','asc.anc.se','asc.anc.t','asc.anc.p','asc.ancsnp.beta','asc.ancsnp.se','asc.ancsnp.t','asc.ancsnp.p','asc.anc1.beta','asc.anc1.se','asc.anc1.t','asc.anc1.p','asc.anc2.beta','asc.anc2.se','asc.anc2.t','asc.anc2.p','asc.model.p','asc.n','asc.ANC.t', 'asc.ANC.p')
		return(ascout)
	} else {
	weights = 1 / (1 / asc1 + 1 / asc2) 
	weight_cap = min(weight_cap, floor(sample_size / 10))
	weight_cutoff = min(weights) * weight_cap
	weights[weights > weight_cutoff] = weight_cutoff
	asc = diag(sqrt(weights)) %*% asc
	x = diag(sqrt(weights)) %*% x
	anc = diag(sqrt(weights)) %*% anc
	xanc = diag(sqrt(weights)) %*% xanc
	lm_asc = lm(asc ~ x + anc + xanc)
#	print(summary(lm_asc))
#	print(is.na(c(summary(lm_asc)$coefficients[c('x', 'anc', 'xanc'),'Estimate'])))
	if(length(summary(lm_asc)$coefficients[,'Estimate'])<4){
#		print(summary(lm_asc))
		anc2_sum = c(NA, NA, NA, NA)
		ANC_sum = c(NA, NA)
		model.pval = NA
	} else {
		anc2_glht = glht(lm_asc, linfct = 'x + xanc = 0')
		anc2_sum = c(summary(anc2_glht)$test$coefficients, summary(anc2_glht)$test$sigma, summary(anc2_glht)$test$tstat, summary(anc2_glht)$test$pvalues[1])
#		print(summary(anc2_glht))
		ANC_glht = glht(lm_asc, linfct = 'xanc = 0')
#	print(summary(ANC_glht))
		ANC_sum = c(summary(ANC_glht)$test$tstat, summary(ANC_glht)$test$pvalues[1])
		model.pval = anova(lm_asc, lm(asc ~ x))[[6]][2]
	}
	ascout = tryCatch({
		c(summary(lm_asc)$coefficients['x',], summary(lm_asc)$coefficients['anc',], summary(lm_asc)$coefficients['xanc',], summary(lm_asc)$coefficients['x',], anc2_sum, model.pval, sample_size, ANC_sum)
	 }, error = function(e){
		return(rep(NA, times=24))
	 })
#	print(ascout)
	names(ascout) = c('asc.snp.beta','asc.snp.se','asc.snp.t','asc.snp.p', 'asc.anc.beta','asc.anc.se','asc.anc.t','asc.anc.p','asc.ancsnp.beta','asc.ancsnp.se','asc.ancsnp.t','asc.ancsnp.p','asc.anc1.beta','asc.anc1.se','asc.anc1.t','asc.anc1.p','asc.anc2.beta','asc.anc2.se','asc.anc2.t','asc.anc2.p','asc.model.p','asc.n','asc.ANC.t', 'asc.ANC.p')
	return(ascout)
}}

ancqtl_trc = function(trc, lib_size, x, anc, xanc, cov, trc_cutoff = 20) {
	trc_in = trc
	trc = log(trc / 2 / lib_size) - cov
	na_ind = is.na(trc) | trc_in < trc_cutoff | is.na(as.numeric(anc))
	x = x[!na_ind]
	anc = anc[!na_ind]
	trc = trc[!na_ind]
	xanc = xanc[!na_ind]
	sample_size = length(trc)
	lm_trc = lm(trc ~ x + anc +xanc)
	if(length(summary(lm_trc)$coefficients[,'Estimate'])<4){
		anc2_sum = c(NA, NA, NA, NA)
		ANC_sum = c(NA, NA)
		model.pval = NA
	} else {
		anc2_glht = glht(lm_trc, linfct = 'x + xanc = 0')
		anc2_sum = c(summary(anc2_glht)$test$coefficients, summary(anc2_glht)$test$sigma, summary(anc2_glht)$test$tstat, summary(anc2_glht)$test$pvalues[1])
		ANC_glht = glht(lm_trc, linfct = 'xanc = 0')
		ANC_sum = c(summary(ANC_glht)$test$tstat, summary(ANC_glht)$test$pvalues[1])
		model.pval = anova(lm_trc, lm(trc ~ x))[[6]][2]
	}
	trcout = tryCatch({
		c(summary(lm_trc)$coefficients['x',], summary(lm_trc)$coefficients['anc',], summary(lm_trc)$coefficients['xanc',],summary(lm_trc)$coefficients['x',] ,anc2_sum, model.pval, sample_size, ANC_sum)
	}, error = function(e){
		return(rep(NA, times=24))
	})
	names(trcout) = c('trc.snp.beta','trc.snp.se','trc.snp.t','trc.snp.p', 'trc.anc.beta','trc.anc.se','trc.anc.t','trc.anc.p','trc.ancsnp.beta','trc.ancsnp.se','trc.ancsnp.t','trc.ancsnp.p','trc.anc1.beta','trc.anc1.se','trc.anc1.t','trc.anc1.p','trc.anc2.beta','trc.anc2.se','trc.anc2.t','trc.anc2.p','trc.model.p','trc.n','trc.ANC.t','trc.ANC.p')
	return(trcout)
}

