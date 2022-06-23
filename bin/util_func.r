
library(gridExtra)
library(ggplot2)

##########
require(scales)
cubic_trans = function() trans_new("cubic", function(x) x^(1/3), function(x) x^(3))
## cubic transformation for ggplot2

getMin <- function(x) min(x, na.rm = TRUE)

getMax <- function(x) max(x, na.rm = TRUE)

getSum <- function(x) sum(x, na.rm = TRUE)

getMedian <- function(x) median(x, na.rm = TRUE)



rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

calcLR <- function(x, y) {
    # Calculate log likelihood ratio statistic value
    y = scale(y)
    x = scale(x, center = TRUE, scale = FALSE)
    n = length(y);
    yty = t(y) %*% y;
    xtx = t(x) %*% x;
    xty = t(x) %*% y;
    return(n * (log(yty) - log(yty - xty * xty / xtx)));
} # Function to calculate log likelihood ratio of y ~ x

calc_pval_chisq <- function(LRT, lambda, df){
    return(1 - pchisq(LRT, df=df, ncp = lambda))
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
}

panel.outlier <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    points(x[col ==2], y[col ==2], pch = pch, col = col[col ==2], bg = bg, cex = cex)
}

na_mean <- function(x){
    return(mean(x, na.rm = TRUE))
}

Inverse_Normalize <- function(x){
    x.rank <- rank(x)
    x.prob <- (x.rank - 0.5) / length(x)
    x.norm <- qnorm(x.prob)
    return(x.norm)
}

rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

count_na <- function(x){return(sum(is.na(x)))}

Impute <- function(x){
	# replace NA by sample mean
    meanx = mean(x, na.rm = TRUE)
    x[is.na(x)] = meanx
    return(x)
}

standardize <- function(x){
    meanx = mean(x, na.rm = TRUE)
    sdx = sd(x, na.rm = TRUE)
    x[is.na(x)] = meanx
    if(sdx > 0){
        x = (x - meanx) / sdx
        }else{
            x = (x - meanx)
        }
     return(x)
}


EffectSize <- function(N_case, N_control){
	return(4/(1/N_case + 1/N_control))
}

calc_p_var <- function(p_case, p_control, N_case, N_control){
	p = (p_case * N_case + p_control * N_control) / (N_case + N_control)

	return( p * (1-p) * (1 / (2 * N_case) + 1/(2 * N_control) ) )
}

calc_zscore_dich <- function(p_case, p_control, N_case=1000, N_control=1000){
	return((p_case - p_control) / sqrt(calc_p_var(p_case, p_control, N_case, N_control)))
}

get_p_value <- function(zscore){
	return(2 * (1 - pnorm(abs(zscore))))
}

calc_power <- function(C, mu){
	return( 1 - pnorm(C - mu) + pnorm(-C - mu) )
}

calc_power_ncp <- function(alpha, lambda, df = 1){
	C = qchisq(1 - alpha, df = df, ncp = 0)
	return(1 - pchisq(C, df = df, ncp = lambda))
}

sqrtMat <- function(a){
	a.eig <- eigen(a)
	a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
	return(a.sqrt)
}

calcMAF <- function(x) {
    n = sum(!is.na(x))
    if(n >0){
        maf  = sum(x, na.rm = TRUE) / (2 * sum(!is.na(x))) # unfold maf
    }else{
        maf = NA # all samples has NA genotypes
    }
	return(maf)
}

calcXmean_all <- function(xmean_vec, n_vec){
        n_weight = n_vec / sum(n_vec)
        xmean = apply(t(xmean_vec) * n_weight, 2, sum)
        return(xmean)
}

getMAFw_adj <- function(maf_vec, maf_vec_fitted, n_vec){
    maf_adj = maf_vec - maf_vec_fitted
    maf_adj_avg = calcXmean_all(maf_adj, n_vec)

    n_weight = n_vec / sum(n_vec)
    w_mat = t(maf_vec * (1 - maf_vec) + 2 * (maf_adj - maf_adj_avg)^2 ) * n_weight
    w = apply(w_mat, 2, sum)

    w[w>0] = 1 / sqrt(w[w>0])
    w = w / sum(w, na.rm = TRUE)

    return(w)
}

getAVGw_adj <- function(xvar_mat, maf_vec_adj, n_vec){

    n_weight = n_vec / sum(n_vec)
    w_mat = t( xvar_mat + 2 * maf_vec_adj^2 ) * n_weight
    w = apply(w_mat, 2, sum)
    
    w[w>0] = 1 / sqrt(w[w>0])
    w = w / sum(w, na.rm = TRUE)

    return(w)
}



getAVGw <- function(maf_vec, n_vec){

    n_weight = n_vec / sum(n_vec)
    w_mat = t(maf_vec * (1 - maf_vec) ) * n_weight
    w = apply(w_mat, 2, sum)

    w[w>0] = 1 / sqrt(w[w>0])
    w = w / sum(w, na.rm = TRUE)

    return(w)
}

getVw <- function(V, order = 1/2){
    v_vec = diag(V)
    w = rep(0, length(v_vec))
    w[v_vec >0] = 1 / (v_vec[v_vec>0] ^ order)
    w = w / sum(w)
    return(w)
}

getMAFw <- function(maf){
    w =  1 / sqrt( maf * (1 - maf) )
    w[maf == 0 | maf == 1] <- 0
    w <- w / sum(w, na.rm = TRUE)
    return(w)
}

getEqw <- function(maf){
    pos_maf = (maf > 0) & (maf < 1 )
    w = rep(0, length(maf))
    w[pos_maf] = 1 / sum(pos_maf)
    return(w)
}

getBetaw <- function(maf, a = 1, b = 25, f0 = 0.0008333333){
    # f0 equivalent to 4 MAC in the sample
    w = rep(0, length(maf))
    maf[(maf >0) & (maf < f0)] = f0
    maf[(maf >0) & (maf > (1-f0))] = 1-f0

    w[maf>0] = dbeta(maf[maf>0], a, b)
    w = w / sum(w)
    return(w)
}


calcMAF_all <- function(maf_vec, n_vec){
	K = length(n_vec)
	m = dim(maf_vec)[1]
	maf <- mac <- rep(0, m)
    for(i in 1:m){
        n_effect = 0
        for(k in 1:K){
            if(! is.na(maf_vec[i, k]) ){
                mac[i] = mac[i] + n_vec[k] * maf_vec[i, k]
                n_effect = n_effect + n_vec[k]
            }
            
        }
        maf[i] = mac[i] / n_effect 
    } 
	return(maf)
}

Logit <- function(p){log(p / (1-p))}


getSKAT_pval <- function(Q, Psi, b = 10, method = "lc2"){
    if(method == "skat"){
        pval = Get_PValue_SKAT(Psi, Q)$p.value 
    }else if(method == "lc2"){
        lambda = eigen(Psi, only.values = TRUE)$values
        pval_lc_out = system(paste("/net/fantasia/home/yjingj/lib/lc2-1.0/lc2_s_linux ", paste(lambda, collapse=" "), "-b", b, "-x", Q), intern=TRUE)
        pval = as.numeric(unlist( strsplit(pval_lc_out[[1]], split="=") )[2])
    }else if(method == "davies"){
        lambda = eigen(Psi, only.values = TRUE)$values
        pval = davies(Q, lambda, acc = 1e-20)$Qq
    }else{
        stop("method need to be 'lc2', 'davies', or 'skat' ... ") 
    }
    return(pval)
}

count_reject <- function(pval, alpha){
        return( sum(pval < alpha, na.rm = TRUE) )
}

calc_typeI <- function(pval, alpha){
    m = sum(!is.na(pval))
    if(m>0){
        typeI_err = sum(pval < alpha, na.rm = TRUE) / m
    }else{
        typeI_err = NA
        print("input pvalues are all NAs!")
    }
    return(typeI_err)
}


GetR2 <- function(model){
    R2<- 1-(model$deviance/model$null.deviance)
    return(R2)
}

GetGC <- function(pval_vec){
    x = qchisq((1-pval_vec), df = 1)
    return(median(x, na.rm = TRUE) / qchisq(0.5, df =1)) 
}


qq <- function(pvector, title="Quantile-quantile plot of p-values", ypos=0.8, xpos=1.5, size = 28, spartan=F) {
    pvector = pvector[!is.na(pvector)]
    n = length(pvector)
    lambda_gc = GetGC(pvector)
    pvector_uniq = unique(pvector)

    o = -log10(sort(pvector_uniq))
    e = -log10( seq(1/(n), n/(n+1), length.out = length(o)) )

    qqplot=qplot(e, o, xlim=c(0,max(e)), ylim=c(0,max(o[o<Inf]))) + 
        geom_abline(intercept=0,slope=1, col="red") +
        labs(title=title) + xlab(expression(Expected~~-log[10](italic(p)))) +
        ylab(expression(Observed~~-log[10](italic(p))))+ 
        annotate("text", label = paste("lambda[GC] == ", round(lambda_gc, 4)), parse = TRUE, x = xpos, y = max(o)*ypos, size = 8) +
        theme_grey(base_size = size)
    # ggsave("/net/fantasia/home/yjingj/METAL/results/qq_adj.pdf")

    return(qqplot)
}

myqq <- function(pvector, title="Quantile-quantile plot of p-values", add=FALSE, colour = "blue", lty = 1) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    if(!add) {
        plot(e, o, type = "l", main = title, xlim=c(0,max(e)), ylim=c(0,max(o[o<Inf])), col = colour, lty = lty, lwd = 2, xlab = "expectation", ylab = "observation")
        abline(0, 1)
    }else{
        points(e, o, type = "l", col = colour, lty = lty, lwd = 2)
    }
}



