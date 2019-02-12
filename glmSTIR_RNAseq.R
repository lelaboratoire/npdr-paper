
library(devtools)
install_github("insilico/glmSTIR")
library(glmSTIR)

### RNA-Seq data
load("unfilteredMostafavi.Rdata")
dim(unfiltered.mostafavi)  # 915 subjects x (15230 genes + 1 class)
class.lab <- "class"
sum(colnames(unfiltered.mostafavi)==class.lab) # check for class
rnaseq.mdd.phenotype <- unfiltered.mostafavi[, class.lab]
# 463 MDD, 452 HC

## Match the RNAseq subject ids with the GWAS data to get the sex covariate 
rnaseq.subject.ids.df <- data.frame(ids = row.names(unfiltered.mostafavi))  
# need to match these in MDD_LD.fam file to get sex covariate

### GWAS data, which contains sex covariate
#library(snpStats)
#mdd.gwas.data <- read.plink(bed="MDD_LD")

splitFN <- function(elt, string_vector){
  unlist(strsplit(string_vector[elt],"_"))[3]
}

fam.data <- read.table("MDD_LD.fam", header=F)
colnames(fam.data) <- c("pedigree","member", "father", "mother", "sex", "affected")
gwas.subject.LONG.ids <- as.character(fam.data$pedigree)
sex.column <- fam.data$sex
  
#gwas.subject.LONG.ids <- mdd.gwas.data$fam$pedigree
#sex.column <- mdd.gwas.data$fam$sex
gwas.subject.ids <- sapply(1:length(gwas.subject.LONG.ids), function(x){splitFN(x,gwas.subject.LONG.ids)})
gwas.sex.df <- data.frame(ids = gwas.subject.ids, sex=sex.column)

#create sex covariate column vector that has same subject id order as RNA-Seq data
sex.merged <- merge(rnaseq.subject.ids.df,gwas.sex.df,by="ids")
# sex.vec or sex.covar.mat will work with glmSTIR
sex.vec <- sex.merged$sex
sex.covar.mat <- as.matrix(sex.merged$sex, ncol=1)
colnames(sex.covar.mat) <- "sex"  # not required by glmSTIR to add column name
# 274 male, 641 female

# assocation between MDD and sex
chisq.test(rnaseq.mdd.phenotype,sex.vec)
table(rnaseq.mdd.phenotype,sex.vec)

### Filter RNA-Seq and run glmSTIR

unfiltered.predictors.mat <- unfiltered.mostafavi[, - which(colnames(unfiltered.mostafavi) == class.lab)]


geneLowVarianceFilter <- function(dataMatrix, percentile=0.5) {
  variances <- apply(as.matrix(dataMatrix), 2, var)
  threshold <- quantile(variances, c(percentile))
  #hist(variances)
  #abline(v = threshold)
  # remove variable columns with lowest percentile variance
  mask <- apply(dataMatrix, 2, function(x) var(x) > threshold)
  fdata <- dataMatrix[, mask]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

# strict filter so it finishes in a couple minutes
# use .5 in real analysis, but it will take a while (a day?)
pct <- 0.7 # .9, 1,524 genes, .75 3,809 vars, .98 306vars 4.7mins, .7 4570vars
filter <- geneLowVarianceFilter(unfiltered.predictors.mat, pct)
#filtered.mostafavi <- filter$fdata
filtered.mostafavi.df <- data.frame(filter$fdata, class = rnaseq.mdd.phenotype)
#unfiltered.mostafavi <- data.frame(predictors.mat, class = pheno.class)
dim(filtered.mostafavi.df)


########################################## Univariate
class.idx <- length(colnames(filtered.mostafavi.df))
colnames(filtered.mostafavi.df)[class.idx]

rnaseq.univariate <- univariateRegression(outcome="class", dataset=filtered.mostafavi.df, regression.type="glm")
rnaseq.univariate[rnaseq.univariate[,3]<.05,]
write.csv(rnaseq.univariate[1:20,],"ResultsRNA/univar.mdd.csv")
rnaseq.univariate[1:40,]

rnaseq.univariate.sexadjust <- univariateRegression(outcome="class", dataset=filtered.mostafavi.df, 
                                                    regression.type="glm", covars=as.factor(sex.vec))
rnaseq.univariate.sexadjust[rnaseq.univariate.sexadjust[,3]<.05,]
rnaseq.univariate.sexadjust[1:10,]
write.csv(rnaseq.univariate.sexadjust[1:20,],"ResultsRNA/univar.mdd.sexadjust.csv")

# univariate with sex as outcome
rnaseq.sextest <- univariateRegression(outcome=as.factor(sex.vec), dataset=filtered.mostafavi.df[,-class.idx], regression.type="glm")
rnaseq.sextest[rnaseq.sextest[,3]<.05,]
write.csv(rnaseq.sextest[rnaseq.sextest[,3]<.05,],"ResultsRNA/univar.sexassoc.csv")
write.csv(rnaseq.sextest[1:50,],"ResultsRNA/univar.sexassoc.csv")

# compare list of genes for no adjust, sex-adjusted, and sex-associated
topx <- 20
univ.compare <- cbind(rownames(rnaseq.univariate)[1:topx], rownames(rnaseq.univariate.sexadjust)[1:topx],
      rownames(rnaseq.sextest)[1:topx]); colnames(univ.compare) <-c("no-adjust", "sex-adjust","sex-assoc")
univ.compare

#cbind(rnaseq.univariate[1:topx,],
#      rownames(rnaseq.univariate.sexadjust)[1:topx], rnaseq.univariate.sexadjust[1:topx,],
#      rownames(rnaseq.sextest)[1:topx], rnaseq.sextest[1:topx,])


################################### run glmstir, no covariate adjustment
start_time <- Sys.time()
glmstir.mdd.rnaseq.results <- glmSTIR("class", filtered.mostafavi.df, regression.type="glm", attr.diff.type="numeric-abs",
                                      nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                      fdr.method="bonferroni")
end_time <- Sys.time()
end_time - start_time  # about 5 min for pct=.98, 306vars, 18min for pct=.9 and 1524 vars

glmstir.mdd.rnaseq.results[1:10,]  # top 10, p-value sorted
#pval.adj    pval.attr beta.attr   beta.0       pval.0
#UTY      1.123919e-12 3.684982e-15 -7.865203 11.06278 1.901189e-28
#RPS4Y1   3.090827e-11 1.013386e-13 -7.439146 10.71743 8.431560e-27
#PRKY     3.098578e-11 1.015927e-13 -7.438815 10.72223 8.005256e-27
#KDM5D    7.089015e-11 2.324267e-13 -7.328682 10.65013 1.741177e-26
#EIF1AY   7.169368e-10 2.350613e-12 -7.011925 10.42237 1.960046e-25
#TTTY15   2.649630e-09 8.687310e-12 -6.826728 10.29496 7.424777e-25
#USP9Y    8.428545e-09 2.763457e-11 -6.658654 10.15141 3.266155e-24
#CYorf15A 1.483237e-08 4.863073e-11 -6.575069 10.10700 5.143328e-24
#DDX3Y    1.849974e-08 6.065488e-11 -6.542115 10.06183 8.147345e-24
#CYorf15B 2.475517e-08 8.116450e-11 -6.498426 10.03785 1.039183e-23

glmstir.mdd.fdr.idx <- glmstir.mdd.rnaseq.results[,1]<.05
glmstir.mdd.rnaseq.positives <- row.names(glmstir.mdd.rnaseq.results[glmstir.mdd.fdr.idx,]) # reSTIR p.adj<.05
glmstir.mdd.rnaseq.positives
glmstir.mdd.rnaseq.results[glmstir.mdd.fdr.idx,]
write.csv(glmstir.mdd.rnaseq.results[glmstir.mdd.fdr.idx,],"ResultsRNA/NPDR.mdd.csv")


library(glmnet)
######################### VECTOR MEMORY EXHAUSTED (LIMIT REACHED)
##### Run glmnetSTIR, penalized glmSTIR, no covar
glmnetSTIR.mdd.results <- glmSTIR("class", filtered.mostafavi.df, regression.type="glmnet", attr.diff.type="numeric-abs",
                                 nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5,
                                 glmnet.alpha=1, glmnet.family="binomial",
                                 fdr.method="bonferroni", verbose=T)
# attributes with glmSTIR adjusted p-value less than .05 
glmnetSTIR.mdd.results.mat <- as.matrix(glmnetSTIR.mdd.results)
# .05 regression coefficient threshold is arbitrary
# not sure why glment did not force zeros
# Finds more interactions than regular glmnet, but not nearly as good as regular glmSTIR
nonzero.glmnetSTIR.mdd.mask <- abs(glmnetSTIR.mdd.results.mat[,1])>.05  
as.matrix(glmnetSTIR.mdd.results.mat[nonzero.glmnetSTIR.mdd.mask,],ncol=1)

## run glmstir with covariates
start_time <- Sys.time()
# sex.vec or sex.covar.mat will work as the covariate variable.
# works with multiple covariates. 
glmstir.mdd.sexadj1.rnaseq.results <- glmSTIR("class", filtered.mostafavi.df, regression.type="glm", 
                                      attr.diff.type="numeric-abs", 
                                      nbd.method="multisurf", 
                                      nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                      covars=sex.vec,                   # works with sex.covar.mat as well
                                      covar.diff.type="match-mismatch", # for categorical covar like sex 
                                      fdr.method="bonferroni", verbose=TRUE)
end_time <- Sys.time()
end_time - start_time

# When you adjust for sex, these obviously male-related confounders
# are no longer significant: USP9Y, UTY, PRKY.
glmstir.mdd.sexadj1.rnaseq.results[1:20,]  # top 10, p-value sorted
#              pval.adj    pval.attr beta.attr     beta.0       pval.0
#TTTY10    4.271299e-10 1.400426e-12  7.084009 -10.369698 3.406002e-25
#LOC220594 3.977758e-05 1.304183e-07 -5.278257  -4.643482 3.425865e-06
#NEBL      2.879843e-04 9.442108e-07  4.902923  -9.064953 1.246584e-19
#FAM22A    4.053387e-04 1.328979e-06 -4.835370  -4.725455 2.296010e-06
#HTR6      2.981664e-03 9.775946e-06 -4.422071  -5.414826 6.134819e-08
#CHGB      3.596614e-03 1.179218e-05 -4.381395  -4.703197 2.561187e-06
#HBG2      6.686266e-03 2.192218e-05  4.244356  -8.827761 1.067923e-18
#HLA.L     6.742079e-03 2.210518e-05 -4.242492  -4.898121 9.675759e-07
#HLA.DQA2  9.171833e-03 3.007158e-05 -4.172924  -5.204169 1.948665e-07
#RPS4Y2    1.962924e-02 6.435816e-05  3.996234  -7.980825 1.453586e-15

glmstir.mdd.sex1.fdr.idx <- glmstir.mdd.sexadj1.rnaseq.results$pval.adj<.05
glmstir.mdd.sexadj1.rnaseq.positives <- row.names(glmstir.mdd.sexadj1.rnaseq.results[glmstir.mdd.sex1.fdr.idx,]) # reSTIR p.adj<.05
glmstir.mdd.sexadj1.rnaseq.positives
glmstir.mdd.sexadj1.rnaseq.results[glmstir.mdd.sex1.fdr.idx,]
write.csv(glmstir.mdd.sexadj1.rnaseq.results[glmstir.mdd.sex1.fdr.idx,],"ResultsRNA/NPDR.mdd.sexadjust.csv")

##### residualize?
#r=glm(trait~1, offset=c*t, family=binomial(link=probit))$residuals
#fit=lm(r~genotypes) 
#plink --bfile data --linear --covar covars.txt --adjust --out data

# sex-associated by glmSTIR
start_time <- Sys.time()
class.idx <- length(colnames(filtered.mostafavi.df))
glmstir.mdd.sexassoc.rnaseq.results <- glmSTIR(as.factor(sex.vec), filtered.mostafavi.df[,-class.idx], 
                                      regression.type="glm", attr.diff.type="numeric-abs",
                                      nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                      fdr.method="bonferroni")
end_time <- Sys.time()
end_time - start_time  # about 5 min for pct=.98, 306vars

sexassoc.idx <- glmstir.mdd.sexassoc.rnaseq.results$pval.adj<.05
length(sexassoc.idx)
glmstir.mdd.sexassoc.rnaseq.results[1:20,] 

topx <- 20
stir.sex.compare <- cbind(rownames(glmstir.mdd.rnaseq.results)[1:topx], 
                      rownames(glmstir.mdd.sexadj.rnaseq.results)[1:topx],
                      rownames(glmstir.mdd.sexassoc.rnaseq.results)[1:topx])
colnames(stir.sex.compare) <-c("no-adjust", "sex-adjust","sex-assoc")
stir.sex.compare

### remove sex-associated/confounded genes from distance matrix and neighborhood calculation 
glmstir.mdd.sexadj2.rnaseq <- glmSTIR("class", filtered.mostafavi.df, regression.type="glm", 
                                       attr.diff.type="numeric-abs", 
                                       nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                       covars=sex.vec,                   # works with sex.covar.mat as well
                                       covar.diff.type="match-mismatch", # for categorical covar like sex
                                       #rm.attr.from.dist=c(), # length=0
                                       #rm.attr.from.dist = rownames(glmstir.mdd.sexassoc.rnaseq.results)[1:topx],
                                       rm.attr.from.dist = rownames(rnaseq.sextest[1:40,]),
                                       fdr.method="fdr", verbose=TRUE)
rownames(glmstir.mdd.sexadj2.rnaseq[glmstir.mdd.sexadj2.rnaseq$pval.adj<.05,]) 
stir.sex.compare2 <- cbind(rownames(glmstir.mdd.rnaseq.results)[1:topx], 
                           rownames(glmstir.mdd.sexassoc.rnaseq.results)[1:topx],
                           rownames(glmstir.mdd.sexadj.rnaseq.results)[1:topx],
                           rownames(glmstir.mdd.sexadj2.rnaseq)[1:topx]
                          )
colnames(stir.sex.compare2) <-c("no-adjust", "sex-assoc","sex-adjust1","sex-adjust2")
stir.sex.compare2

library(randomForest)
ranfor.fit <- randomForest(as.factor(class) ~ ., data = filtered.mostafavi.df) 
rf.importance <- importance(ranfor.fit)  # variable 3 best
x<-sort(rf.importance, decreasing=T, index.return=T)
x$ix
rownames(rf.importance)[x$ix]
rf.sorted <- cbind(rownames(rf.importance)[x$ix],rf.importance[x$ix])
rf.sorted[1:50,]
write.csv(rf.sorted[1:50,],"ResultsRNA/randomForest.mdd.csv")


library(privateEC)
cncv.rnaseq <- consensus_nestedCV(train.ds = filtered.mostafavi.df, 
                                        validation.ds =  filtered.mostafavi.df, 
                                        label = "class",
                                        method.model = "classification",
                                        is.simulated = TRUE,
                                        ncv_folds = c(10, 10),
                                        param.tune = FALSE,
                                        learning_method = "rf", 
                                        importance.algorithm = "ReliefFequalK",
                                        relief.k.method = "k_half_sigma",     # surf k
                                        num_tree = 500,
                                        verbose = F)

cat("\n Train Accuracy [",cncv.rnaseq$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv.rnaseq$Validation,"]\n")
cat("\n Selected Features \n [",cncv.rnaseq$Features,"]\n")
cat("\n Elapsed Time [",cncv.rnaseq$Elapsed,"]\n")

rownames(glmstir.mdd.rnaseq.results[glmstir.mdd.rnaseq.results$pval.adj<.05,])
