options(stringsAsFactors=FALSE)
library(data.table)


##### load all hits genotype data from 1KG

VCF_Data = fread("/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/vcfs/causalSNP.vcf", sep="\t", header=TRUE)

dim(VCF_Data)

GenoData = t(data.frame(VCF_Data[, -(1:9), with=FALSE]))
n = 2504
p = 15

GenoMat = matrix(0, nrow = 2504, ncol = 15)
snpID = VCF_Data$ID
sampleID = names(VCF_Data)[-(1:9)]

for(i in 1:n){
	for(j in 1:p){
		GenoMat[i, j] = sum( as.numeric ( unlist(strsplit(GenoData[i, j], "[|]")) ) )
	}
}
x_mean = apply(GenoMat, 2, mean)

###### load Odds Ratios #########

OR = read.table(file = "/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/vcfs/causalSNP_OR.txt", header = TRUE)
beta = log(OR$OR) * 2

###### Simulate AMD phenotypes by logistic regression

prev = 0.5 # for those > 75 years old

a0 = log(prev/(1-prev)) - sum(x_mean * beta) 

logit_y = a0 + GenoMat %*% beta

pr = 1 / (1 + exp(-logit_y)); 


Ysimu = rbinom(n, 1, pr) # generate cases

n_case = sum(Ysimu)
n_control = n - n_case

print(c(n_case, n_control))

########### Test model ###########

##### fit linear model
fit1 = glm(Ysimu ~ GenoMat)
summary(fit1)
1 - (fit1$deviance/fit1$null.deviance) # r2 = 22.2%

##### fit logistic model
fit2 = glm(Ysimu ~ GenoMat, family = binomial)
summary(fit2)

########## Write phenotypes ############
write.table(data.frame(sampleID, Ysimu), file = "/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/phenoAMD_1KG.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)













