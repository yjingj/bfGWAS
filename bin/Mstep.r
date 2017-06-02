Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)
source("/net/fantasia/home/yjingj/GIT/bfGWAS/bin/R_funcs.r")

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
# print(args)

hypfile=args[[1]]
k=as.numeric(args[[2]])
pp = as.numeric(args[[3]])
abgamma = as.numeric(args[[4]])
EM_result_file = args[[5]]
hypcurrent_file=args[[6]]

# abgamma=0.1
a_gamma <- b_gamma <- abgamma;
print(paste("a_gamma=b_gamma = ", abgamma))

ptm <- proc.time()

########### Load data ....

# hypfile="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/hypval.current 

hypdata = read.table(hypfile, sep="\t", header=FALSE)
n_type = (dim(hypdata)[2] - 4)/4
print(paste(" Total Annotation categories : ", n_type))

temp_col_names <- c("block", "loglike", "GV", "rv")
for(i in 1:n_type){
	temp_col_names <- c(temp_col_names, 
			paste(c("n", "G", "m", "sigma2"), (i-1), sep = "_"))
}

colnames(hypdata) <-  temp_col_names

########### Update hyper parameter values
rv = mean(hypdata[, "rv"])
tau = 1.0 / rv
pve = sum(hypdata[, "GV"])

prehyp <- read.table(hypcurrent_file, header=TRUE)
print("hyper parameter values before MCMC: ")
print(prehyp)

######### Set hierarchical parameter values
n_vec = rep(0, n_type)
for(i in 1:n_type){
	 n_vec[i] <- sum(hypdata[, paste("n", (i-1), sep="_")])
}

#### updating hyper pi and sigma2 values for each group
hypcurrent <- NULL
hypmat <- NULL

for(i in 1:n_type){
	# print(i)
	if(n_vec[i] > 0){
		a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
		}else{a_beta=1; b_beta = 1e6 - 1;}
	
	m_temp = hypdata[, paste("m", (i-1), sep="_")]

	pi_temp = CI_fish_pi(m_temp, n_vec[i], a_beta, b_beta)

	sigma2_temp = CI_fish_sigma2(hypdata[, paste("sigma2", (i-1), sep="_")], m_temp, tau, a_gamma, b_gamma)

	hypcurrent <- c(hypcurrent, pi_temp, sigma2_temp)
	# print(cbind(pi_temp, sigma2_temp))
	hypmat <- rbind(hypmat, c(pi_temp[1], sigma2_temp[1]))
}

########## Write out updated hyper parameter values
colnames(hypmat) <- c("pi", "sigma2")
print("hyper parameter values updates after MCMC: ")
print(hypmat)
write.table(format(hypmat, scientific=TRUE), file=hypcurrent_file, 
	quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

#### Summarize log-likelihood
loglike_total = sum(hypdata$loglike)

for(i in 1:n_type){
	if(sum(prehyp[i, ]>0)==2){

		if(n_vec[i] > 0){
			a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
		}else{a_beta=1; b_beta = 1e6 - 1;}

		loglike_total = loglike_total + 
				logprior_pi(a_beta, b_beta, prehyp[i, 1]) +
    			logprior_sigma(a_gamma, b_gamma, prehyp[i, 2]) 
    }else{
    	print("pre-hyper-parameter <= 0... ")
    }
}

########## Write out updated hyper parameter values and se to EM_result_file
# EM_result_file="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/Eoutput/EM_result.txt"
hypcurrent = c(pve, loglike_total, hypcurrent)
hypcurrent <- format(hypcurrent, scientific = TRUE)
print("write to hypcurrent file with hyper parameter values after MCMC: ")
print(c(k, hypcurrent))
write.table(matrix(c(k, hypcurrent), nrow=1), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)

print("EM step time cost (in minutes) : ")
print((proc.time() - ptm)/60)








