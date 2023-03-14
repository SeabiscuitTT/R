# THIS SCRIPT EVALUATES MCMC CONVERGENCE OF HMSC MODELS OF THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# The preliminaries are as in script S1
# 开头部分如脚本S1所示
# Change the working directory if necessary
# wd = "C:/HMSC_book/R-scripts/Section_7_9_fungi"
# setwd(wd)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
library(Hmsc)
set.seed(1)

# We first read in the object models that includes all six fitted models
# We recommed running the script first with setting thin = 1
# Then rerun it with thin = 10, thin = 100, ... to examine how the convergence statistics improve
# 我们首先读入包括所有六个拟合模型的对象模型。
# 我们建议先运行脚本，设置Thin=1。
# 然后使用Thin=10，Thin=100重新运行它，...。要检查收敛统计数据如何改进，请执行以下操作

nChains = 2
samples = 100
thin = 1 # try with thin = 1, thin = 10, thin = 100, etc.
filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We evaluate here MCMC convergence only for the $\beta$ and $\Omega$ parameters, and only in terms of the potential scale reduction factor.
# While we could evaluate MCMC convergence also in terms of effective sample size, we note that usually examining the potential scale reduction factor is more critical: if convergence is satisfactory in terms of it,
# then it is usually satisfactory also in terms of effective sample size
# 我们在此仅针对$\beta$和$\omga$参数评估MCMC收敛，并且仅根据潜在的规模缩减系数进行评估。
# 虽然我们也可以根据有效样本量来评估MCMC收敛，但我们注意到，通常检查潜在的规模缩减因子更为关键：如果从它的角度看收敛是令人满意的，
# 那么就有效样本量而言，通常也是令人满意的

par(mfrow=c(3,2))
for (i in 1:3){
  mpost = convertToCodaObject(models[[i]][[2]])
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  psrf.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
  mymain = switch (i, "Lognormal Poisson", "Probit", "Normal")
  myxlab1 = switch (i, "", "", "psrf (beta)")
  myxlab2 = switch (i, "", "", "psrf (Omega)")
  hist(psrf.beta, main=mymain, xlab = myxlab1)
  hist(psrf.omega, main=mymain, xlab = myxlab2)
}

# With these data, good MCMC convergence is the most difficult to obtain for the lognormal Poisson model.
# That is, esepcially for the lognormal Poisson model a large value of thin is needed to make the potential scale reduction factors to go close to one.
# While usually the normal model shows the best MCMC convergence, here it behaves essentially identically with the probit model.
# This is because now the normal model has a large fraction of missing data, which makes MCMC convergence more challenging.
# 有了这些数据，对于对数正态泊松模型来说，良好的MCMC收敛性是最难获得的。
# 也就是说，特别是对于对数正态泊松模型，需要一个很大的Thin值才能使潜在的缩尺因子接近于1。
# 虽然通常正态模型表现出最佳的MCMC收敛，但在这里它的行为本质上与概率模型相同。
# 这是因为现在的正常模型有很大一部分数据丢失，这使得MCMC收敛更具挑战性。