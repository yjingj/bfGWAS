####### Functions for logistic model
options(stringsAsFactors=F)
library(data.table)
library(ggplot2)

LoadEMdata <- function(filename, header = FALSE){
    paramdata = fread(filename, sep = "\t", header = header)
    setnames(paramdata, c("ID", "chr", "bp", "ref", "alt", "maf", "func", "beta", "pi", "Zscore", "SE_beta", "LRT", "pval_LRT", "rank"))
    setkey(paramdata, "ID")
    return(paramdata)
}

LoadEMhyp <- function(filename, header=FALSE){
    EMhyp_data <- read.table(filename, sep = "\t", header = header)
    n_type <- (dim(EMhyp_data)[2] - 3)/4
    temp_names <- c("Iter", "h", "loglike")
    for(i in 1:n_type){
        temp_names <- c(temp_names, paste(c("pi", "pi_se", "sigma2", "sigma2_se"), (i-1), sep="_"))
    }
    colnames(EMhyp_data)<- temp_names
    return(EMhyp_data)
}

getCI <- function(est, est_se, alpha){
    z_alpha <- -qnorm((1-alpha)/2)
    
    pi_low <- est - z_alpha * est_se
    if(pi_low < 0) pi_low = 0
    
    pi_up <- est + z_alpha * est_se
    return(c(pi_low, pi_up))
}

CItable <- function(v, n_type, alpha, funcgroup=NA){

    temp_table <- data.frame(pi = NA, pi_lcl = NA, pi_ucl = NA, v = NA, v_lcl = NA, v_ucl = NA, funcgroup= funcgroup)
    
    for(i in 1:n_type){
        pi_hat <- v[paste("pi", (i-1), sep="_")]
        pi_se <- v[paste("pi_se", (i-1), sep="_")]

        sigma2_hat <- v[paste("sigma2", (i-1), sep="_")]
        sigma2_se <- v[paste("sigma2_se", (i-1), sep="_")]

        temp_table[i, 1:6] <- c(pi_hat, getCI(pi_hat, pi_se, alpha), sigma2_hat, getCI(sigma2_hat, sigma2_se, alpha))
    }
    return(temp_table)
}

PlotCI_groupPP <- function(hyp_table, pdfname="", size = 28, tit="", wid=10, scale = TRUE){

    p = ggplot(hyp_table, aes(y=pi, x = funcgroup, colour=funcgroup)) + 
        geom_point(size = 4)+ guides (colour=FALSE) + 
        geom_errorbar(aes(ymax=pi_ucl, ymin=pi_lcl), width=0.5, size = 1.5) + 
        labs(title = tit, x = NULL, y = "Causal Probability") + 
        theme_grey(base_size = size) + 
        theme(text = element_text(size=size), axis.text.x = element_text(angle = 30, hjust = 0.8), plot.margin=unit(c(4,0,0,4),"pt"))

    if(scale){trans = "cubic"} else {trans = "identity"}

    p = p + scale_y_continuous(trans="cubic", label = scientific_format())
 
    ggsave(pdfname, plot = p, width = wid)
}


PlotCI_groupESvar <- function(hyp_table = AMDgwas_CI_table, pdfname="", size = 14, tit = "", wid=10, scale = TRUE){

    p = ggplot(hyp_table, aes(y=v, x = funcgroup, colour=funcgroup)) + 
        geom_point(size = 4)+ guides (colour=FALSE) + 
        geom_errorbar(aes(ymax=v_ucl, ymin=v_lcl), width=0.5, size = 1.5) + 
        labs(title = tit, x = NULL, y = "Effect-size Variance") + 
        theme_grey(base_size = size) + 
        theme(text = element_text(size=size), axis.text.x = element_text(angle = 30, hjust = 0.8), plot.margin=unit(c(4,0,0,4),"pt"))

    if(scale){trans = "cubic"} else {trans = "identity"}

    p = p + scale_y_continuous(trans="cubic", label = scientific_format())
 
    ggsave(pdfname, plot = p, width = wid)

}

####### Compare group-wise estimates #######
comp_groupEst <- function(est = group_pip, est_se = pip_se, n_vec, i = 1, conf = 0.95){

    est_1 = est[i]
    w = n_vec / sum(n_vec)
    est_all = sum(est * w )

    var_1 = est_se[i]^2
    var_all = sum(est_se^2 * w^2)

    log_ratio = unlist ( log(est_1 / est_all) )
    log_ratio_se = sqrt(var_1 / est_1^2 + var_all / est_all^2)

    pvalue = 2 * (1 - pnorm( abs(log_ratio), 0, log_ratio_se))
    print(c(exp(log_ratio), pvalue, est_1) )

    z = abs(qnorm((1 - conf)/2))
    log_ratio_lcl = log_ratio - z * log_ratio_se
    log_ratio_ucl = log_ratio + z * log_ratio_se

    output = unname( unlist(c(log_ratio_lcl, log_ratio, log_ratio_ucl)) )
    return(output)

}

###### plot comparison ratioes
PlotRatio <- function(comp_dat, tit ="", pdfname = "", size = 28, ymode = 1, wid = 8, scale = TRUE){

    ymin = min(comp_dat[, 1:3], na.rm = TRUE)
    ymax = max(comp_dat[, 1:3], na.rm = TRUE)

    dodge=position_dodge(width=0.9)
    if(ymode == 1){
        p = ggplot(comp_dat, aes(y=log_ratio, x = funcgroup, colour=funcgroup)) + 
        geom_point(size = 4) + 
        geom_errorbar(aes(ymax=log_ratio_ucl, ymin=log_ratio_lcl), width=0.5, 
                        size = 1.5) + 
        guides(colour=FALSE) + 
        labs(title = tit, x = NULL, y = expression(paste(pi[q], "/", pi[avg])) ) + 
        theme_grey(base_size = size) +
        theme(text = element_text(size=size), axis.text.x = element_text(angle=30, hjust = 0.8), plot.margin=unit(c(4,0,0,4),"pt"))
    }else{
        p = ggplot(comp_dat, aes(y=log_ratio, x = funcgroup, colour=funcgroup)) + 
        geom_point(size=4) + guides(colour=FALSE) + 
        geom_errorbar(aes(ymax=log_ratio_ucl, ymin=log_ratio_lcl), width=0.5, 
                        size = 1.5) + 
        labs(title = tit, x = NULL, y = expression(paste(sigma[q] ^2, "/", sigma[avg] ^2))) +
        theme_grey(base_size = size) +
        theme(text = element_text(size=size), axis.text.x = element_text(angle=30, hjust = 0.8), 
            plot.margin=unit(c(4,0,0,4),"pt"))
    }

    if(scale){trans = "cubic"} else {trans = "identity"}
    p = p + scale_y_continuous(breaks = c(1, floor(seq(ymin, ymax, length.out = 4))), trans=trans) 

    ggsave(pdfname,plot = p, width = wid)
}

### Load Annotation file
LoadAnnodata <- function(filename, header = FALSE){
    Annodata = fread(filename, sep = "\t", header = header)
    n_anno = dim(Annodata)[2] - 3
    setnames(Annodata, c("ID", "chr", "bp", paste("anno", 1:n_anno, sep="_")))
    return(Annodata)
}

getPi <- function(lambda, A){
	lambda <- matrix(lambda, ncol = 1)
	return(1/(1 + exp(-A %*% lambda)))
}


loglike_lambda <- function(lambda, lambda_bar, r, A, S){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	lambda <- lambda - lambda_bar
	l = sum(r * log(pi_pre) + (1-r) * log(1 - pi_pre)) - 0.5 * t(lambda) %*% S %*% (lambda)
	return(l)
}

loglike_lambda_g <- function(lambda, lambda_bar, r, A, S){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	lambda <- lambda - lambda_bar
	g = t(A) %*% (r - pi_pre) - S %*% lambda
	return(g)
}

### Hessian matrix, i.e., Fisher's information
FisherInfo_lambda <- function(lambda, A){
	pi_pre <- 1/(1 + exp(-A %*% lambda))
	h=diag(0, dim(A)[2])
	for(i in 1:dim(A)[1]){
		h = h + pi_pre[i] * (1 - pi_pre[i]) * (A[i, ] %*% t(A[i, ]))
	}
	return(h)
}

###########
weight <- function(x){
	return(x/sum(x))
}

loglike_v <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	l = sum(r * (0.5 * log(tau / s) - 0.5 * tau * beta^2/s ) ) + sum(-(a+1) * log(v) - b / v)
	return(l)
}

loglike_v_g <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	g = t(W) %*% (r * (-0.5/s + 0.5*tau*(beta/s)^2)) + (-(a+1) / v + b/v^2)
	return(g)
}


FisherInfo_v <- function(v, r, beta, tau, W, a, b){
	s = (W) %*% v
	h = diag(0, length(v))

	for(i in which(r>0)){
		s2 = s[i]^2
		s3 = s2 * s[i]
		h = h + r[i] * (0.5 / s2 - tau * beta[i]^2/s3) * ( W[i, ] %*% t(W[i, ]) )
	}
	# h = t(W) %*% diag(as.vector(r * (0.5 / s^2 - tau * beta^2/s^3))) %*% W 
		#+ diag(as.vector((a+1)/v^2 - 2*b/v^3))
	return(-h)
}

######################## log prior functions
logprior_sigma <- function(a, b, x){ return(-(1+a) * log(x) - b/x) }

logprior_pi <- function(a, b, x){ return((a-1) * log(x) + (b-1) * log(1 - x)) }

logprior_lambda <- function(s, x){ return(- 0.5 * x^2 / s) }


############## Functions for non-overlap annotation moel ########
####### Functions to calculate MLE estimate and CIs

Est_pi <- function(m, p, a, b){
	pi_hat = (sum(m) + a - 1.0) / (p + a + b - 2.0)
	if(pi_hat <= 0 || pi_hat > 1){
		pi_hat = a/(a+b)
	}
	return(pi_hat)
}

CI_pi <- function(m, m_low, m_upper, p, a, b){
	pi_hat = Est_pi(m, p, a, b)
	if(pi_hat * (1-pi_hat) / (p + a + b - 2.0) > 0){
		  se_pi = sqrt(pi_hat * (1-pi_hat) / (p + a + b - 2.0))
		}else{se_pi=0}
	pi_low_fisher = pi_hat - 1.645 * se_pi
	pi_upper_fisher = pi_hat + 1.645 * se_pi
	pi_low_mcmc = Est_pi(m_low, p, a, b)
	pi_upper_mcmc = Est_pi(m_upper, p, a, b)
	return(c(pi_hat, se_pi, pi_low_fisher, pi_upper_fisher, pi_low_mcmc, pi_upper_mcmc))
}

CI_fish_pi <- function(m, p, a, b){
	pi_hat = (sum(m) + a - 1.0) / (p + a + b - 2.0)
	if(pi_hat <= 0 || pi_hat > 1){
		pi_hat = a/(a+b)
		se_pi = 0
	}else{
		if(pi_hat * (1-pi_hat) / (p + a + b - 2.0) > 0){
		  se_pi = sqrt(pi_hat * (1-pi_hat) / (p + a + b - 2.0))
		}else{se_pi=0}
	}
	return(c(pi_hat, se_pi))
}

# Without MCMC CI
Est_sigma <- function(beta_est, pi_est, tau, a, b){
	sigma2_hat = (sum(beta_est^2 * pi_est) * tau + 2 * b) / (sum(pi_est) + 2 * (a + 1))
	return(sigma2_hat)
}

CI_sigma <- function(beta_est, pi_est, tau, a, b){
	sigma2_hat = Est_sigma(beta_est, pi_est, tau, a, b)
	se_sigma2 = sigma2_hat * sqrt(1/(sum(pi_est) * tau - sum(pi_est)/2 - (a+1) + 2*b/sigma2_hat))
	sigma2_low_fisher = sigma2_hat - 1.645 * se_sigma2
	sigma2_upper_fisher = sigma2_hat + 1.645 * se_sigma2
	return(c(sigma2_hat, se_sigma2, sigma2_low_fisher, sigma2_upper_fisher, NA, NA))
}

Est_sigma2 <- function(sigma2, m, tau, a, b){
	sigma2_hat = (sum(sigma2 * m) * tau + 2 * b) / (sum(m) + 2 * (a + 1))
	return(sigma2_hat)
}

CI_fish_sigma2 <- function(sigma2, m, tau, a, b){
	sigma2_hat = Est_sigma2(sigma2, m, tau, a, b)
	if( (sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat) < 0){
		se_sigma2=0
	}else{
		se_sigma2 = sigma2_hat * sqrt(1/(sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat))
	}
	return(c(sigma2_hat, se_sigma2))
}

# With MCMC CI

CI_sigma2 <- function(sigma2, m, sigma2_low, sigma2_upper, tau, a, b){
	sigma2_hat = Est_sigma2(sigma2, m, tau, a, b)
	se_sigma2 = sigma2_hat * sqrt(1/(sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat))
	sigma2_low_fisher = sigma2_hat - 1.645 * se_sigma2
	sigma2_upper_fisher = sigma2_hat + 1.645 * se_sigma2
	sigma2_low_mcmc = Est_sigma2(sigma2_low, m, tau, a, b)
	sigma2_upper_mcmc = Est_sigma2(sigma2_upper, m, tau, a, b)
	return(c(sigma2_hat, se_sigma2, sigma2_low_fisher, sigma2_upper_fisher, sigma2_low_mcmc, sigma2_upper_mcmc))
}


CI_LnRatio <- function(a, b, se_a, se_b){
	Ln_r_hat = log(a/b)
	se_Ln_r = sqrt(se_a^2 / a^2 + se_b^2 / b^2)
	Ln_r_low = Ln_r_hat - 1.645 * se_Ln_r
	Ln_r_upper = Ln_r_hat + 1.645 * se_Ln_r
	return(c(Ln_r_hat, se_Ln_r, Ln_r_low, Ln_r_upper, NA, NA))
}