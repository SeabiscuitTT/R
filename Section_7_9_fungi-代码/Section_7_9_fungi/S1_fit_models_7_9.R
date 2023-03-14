# THIS SCRIPT CONSTRUCTS AND FITS HMSC MODELS FOR THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# We first set the working directory to the directory that has all the files for this case study
# We recommend that you create always the subdirectories "data" and "models" to the main directory
# If running the script in another computer, all that needs to be changed in the script below is the main directory
# 我们首先将工作目录设置为包含此案例研究的所有文件的目录。
# 我们建议您在主目录中始终创建子目录\“data\”和\“model\”
# 如果在另一台计算机上运行该脚本，则需要在下面的脚本中更改的只是主目录

# wd = "C:/HMSC_book/R-scripts/Section_7_9_fungi"
# setwd(wd)
wd = "C:/Users/Google/Documents/R/Joint Species Distribution Modelling in R/section_7_9_fungi_2020_05_31/Section_7_9_fungi"
setwd(wd)
localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# We load the Hmsc package, and set the random seed, so that all results from this script are reproducable
# 我们加载hMSC包，并设置随机种子，这样该脚本的所有结果都是可重现的
library(Hmsc)
set.seed(1)

# We next read in the data and have a look at it
# 接下来我们读入数据并进行查看
data = read.csv(file=file.path(data.directory, "data.csv"),stringsAsFactors=TRUE)
n = dim(data)[1]
head(data[,1:6])

# For each of n=99 logs, the decay class (DC) of the log is classified as 1-4, with 1 represening freshly fallen logs and 4 much decayed logs.
# 对于n=99个原木中的每一个，原木的衰变级别(DC)被分类为1-4级，其中1个代表刚倒下的原木，4个较多腐烂的原木。
# The column readcount is the total number of sequences obtained for each sample, which number can be viewed to represent observation effort and is commonly called the sequencing depth.
# The Column Readcount是每个样本获得的序列总数，可以查看该数字以表示观察工作，通常称为测序深度。
# The remaing columns include the sequence counts for each species.
# 剩余栏包括每个物种的顺序计数。
# The sequences were identified (in a data pre-processing step) using the probabilistic taxonomic placement algorithm (PROTAX-Fungi, for reference see the book section).
# 使用概率分类放置算法(PROTAX-Fungi，参考本书部分)识别序列(在数据预处理步骤中)。
# The sequences that could not be assigned with at least 50% probability to any species are pooled in the category `unk`.
# 不能以至少50%的概率将其分配给任何物种的序列被汇集在“未知”类别中。
# This category contains typically the majority of the samples. This is because species-level identification is challenging e.g. due to incompleteness of the reference databases.
# 此类别通常包含大部分样本。这是因为物种一级的鉴定具有挑战性，例如由于参考数据库的不完整。
# The remaining columns consist of those 413 fungal species that were detected in the data at least once with  at least 50% identification probability. 
# 其余列包括在数据中至少检测到一次且识别概率至少为50%的413种真菌。
# Using 50% as a threshold means that many of the species-level names may be wrong, in which the correct species name likely to belong to some closely related species. 
# 使用50%作为阈值意味着许多物种级别的名称可能是错误的，其中正确的物种名称可能属于一些密切相关的物种。

# We first organise the data. We construct the dataframe XData that includes the decay class and the total sequence count, as these will be used as fixed effects.
# 我们首先整理数据。我们构造包含衰减类和总序列计数的数据帧扩展数据，因为它们将用作固定效果。
# We then construct the dataframe YData of the species data, including the unknown class.
# 然后我们构造物种数据的数据帧YData，包括未知类。

XData = data.frame(DC = as.factor(data$DC), readcount = data$readcount)
YData = data[,4:dim(data)[2]]

# We are primarily interested in co-occurrences among the species. Detecting co-occurrences requires a substantial amount of information, for which reason the above script selects only those species that were found in at least 10 logs.
# This leaves us with 31 species, as you can verify by checking dim(YData) before and after evaluating the script below
# 我们主要对物种间的共生现象感兴趣。检测重现需要大量信息，因此上面的脚本只选择那些在至少10个日志中发现的物种。
# 这给我们留下了31个物种，您可以在评估下面的脚本之前和之后选中dim(YData)进行验证

sel.sp = colSums(YData>0)>=10
YData = YData[,sel.sp]

# We next explore the raw data by plotting histograms of species prevalence (P) and abundance (A).
# While in the model we will use the raw sequence counts, here we illustrate the data as relative abundance, and thus below we normalize abundance by "sum(YData)".
# Further, the relative abundances involve many orders of magnitude, some being e.g. 0.0001 and other as 0.1. 
# For this reason, below we take the log-10 transform of the relative abundances when making the histogram.
# This "spreads" the variation so that it more visible. 
# The value of -4 corresponds to relative abundance 0.0001, 
# the value of -3 to relative abundance 0.001, 
# the value of -2 to relative abundance 0.01,
# the value of -1 to relative abundance 0.1, 
# and the value of 0 to relative abudance of 1 (which would mean that the species dominates the entire sample).
# 接下来，我们通过绘制物种流行率(P)和丰度(A)的直方图来探索原始数据。
# 而在模型中，我们将使用原始顺序计数，这里我们将数据说明为相对丰度，因此下面我们通过\“SUM(YData)\”对丰度进行归一化。
# 此外，相对丰度涉及许多数量级，有些为0.0001，另一些为0.1。
# 为此，下面我们在制作直方图时对相对丰度进行log-10变换。
# 这\“传播\”变化，使其更明显。
# -4的值对应于相对丰度0.0001，
# 值-3对应于相对丰度0.001，
# 值-2对应于相对丰度0.01，
# -1的值对应于相对丰度0.1，
# 和值0对应于相对丰度1(这意味着该物种在整个样本中占主导地位)。

P = colMeans(YData>0)
A = colSums(YData)/sum(YData)

par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "Abundance")

# As usually, there is great variation among the species, some being common and others rare.
# By exploring the variables P and A, you can see e.g. that excluding the unknown class that is present in all sampling units,
# the most prevalent species are  Fomitopsis pinicola (prevalence=0.54), Capronia kleinmondensis (0.42) and Phellophilus nigrolimitatus (0.37).
# In terms of abundance, the unknown class represents the vast majority (86%) of all the sequences included in our analyses.
# The most abundant species are F. pinicola (5%), P. nigrolimitatus (1.8%) and A. serialis (1.5%). 
# 与往常一样，物种之间有很大的差异，一些是常见的，另一些是稀有的。
# 通过研究变量P和A，您可以看到，例如，排除所有采样单位中存在的未知类别，
# 最常见的种类是松果壳霉(患病率=0.54)、克氏毛壳虫(0.42)和黑线姬鼠(0.37)。
# 就丰度而言，未知类别代表了我们分析中包括的所有序列的绝大多数(86%)。
# 数量最多的种类是松材线虫(5%)、黑线线虫(1.8%)和串珠线虫(1.5%)。

# We will fit six different HMSC models to these data, consisting of three model types and two choices of the explanatory variables. 
# Below, and in the book, all of these models are fitted in a single loop.
# However, it is not recommended to construct multiple models straight away in a loop, as that comes with the risk of not checking each of the models carefully. 
# For example, with these data the simplest starting point would be to fit a presence-absence model with decay class and log-transformed read-count as explanator variables,and the random effect of the log as a random effect. 
# This model can be defined as follows:
# 我们将为这些数据拟合六种不同的hMSC模型，包括三种模型类型和两种解释变量选择。
# 下面，在这本书中，所有这些型号都安装在一个循环中。
# 但是，不建议直接在循环中构建多个模型，因为这会带来不仔细检查每个模型的风险。
# 例如，利用这些数据，最简单的起点将是以衰减类和对数变换的读取计数作为解释变量，并将日志的随机效应作为随机效应来拟合存在-缺席模型。
# 该模型定义如下：

Y = as.matrix(YData)

# For this model We truncate the data to presence-absence
# 对于此模型，我们将数据截断为存在-缺席

Y = Y>0

# Even if the model definition HMSC accepts logical values (TRUE/FALSE) as the response matrix Y for presence-absence data, some of the post-processing functions assume that Y has numerical 0/1 values. 
# So it is safer to convert Y to numerical either as Y = 1*Y or as
# 即使模型定义hMSC接受逻辑值(真/假)作为存在-缺席数据的响应矩阵Y，一些后处理函数假定Y具有数值0/1值。
# 因此，当Y=1*Y或AS时，将Y转换为数值型会更安全

Y = apply(Y,MARGIN = 2,FUN = as.numeric)

XFormula = ~DC + log(readcount)

studyDesign = data.frame(sample = data$LogID)
rL = HmscRandomLevel(units = studyDesign$sample)
m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
         studyDesign = studyDesign, ranLevels = list(sample = rL),
         distr="probit")

# To understand what was done with the above script, let us look at the objects:
# The study design simply lists the id of each log.
# 要了解上面的脚本做了什么，我们来看一下对象：
# 研究设计简单地列出了每个日志的ID。

head(studyDesign)

# The motivation for including the studyDesign when defining the model object m was that we wished to include a random effect at the log level. 
# The random effects are included by the ranLevels option of Hmsc(...)
# The list includes all levels at which a random effect is to be included. Here the random effect is included as "sample = rL".
# The left-hand side (sample) corresponds to the relevant column of the studyDesign, for which the log-id column (the only column) was named as sample.
# The right-hand side (rL) describes the structure of random effect. To see the structure, evaluate
# 定义模型对象m时包含StudyDesign的动机是我们希望在日志级别包含随机效果。
# 随机效果包含在hMSC(...)的ranLeveles选项中。
# 该列表包括要包含随机效果的所有级别。这里，随机效果包含为\“Sample=RL\”。
# 左侧(Sample)对应StudyDesign的相关列，log-id列(唯一一列)命名为Sample。
# 右手边(RL)描述随机效果的结构。要查看结构，请评估

rL

# This gives the information that "Hmsc random level object with 99 units. Spatial dimensionality is 0 and number of covariates is 0."
# Thus, rL is an "unstructured" random effect, with no reference to space or covariates. This is the simplest type of a random effect.
# There are 99 units, as each log is a separate unit. 
# We note in passing that often the random effects are not defined at the sampling unit level, but at a higher hierarchical level that includes many sampling units. 
# Let us assume e.g. that the studyDesign would include another column describing the plot to which the log would belong, and that we would like to set up another random effect at the level of the plot.
# In such a case, this should NOT be defined as
# rL.plot = HmscRandomLevel(units = studyDesign$plot), 
# but rather as rL.plot = HmscRandomLevel(units = levels(studyDesign$plot)). 
# The levels() picks each level (here each plot) only once, and thus they should not be repeated when defining the random effect strurcture.
# This was a side remark that is not relevant to this case study, see examples on hieararchical study designs for more information on this topic. 
# 这给出了这样的信息：\“hMSC随机级别对象，99个单位。空间维度为0，协变量个数为0。\”
# 因此，RL是一种\“非结构化\”随机效果，与空间或协变量无关。这是最简单的随机效果类型。
# 有99个单位，因为每个日志都是一个独立的单位。
# 顺便说一句，我们注意到随机效果通常不是在采样单元级别定义的，而是在包括许多采样单元的更高层次级别定义的。
# 让我们假设，例如，StudyDesign将包括另一栏，描述日志将属于的曲线图，并且我们想要在曲线图级别设置另一随机效果。
# 在这种情况下，这不应定义为。
# rL.lot=HmscRandomLevel(单位=StudyDesign$PLOT)，
# 而不是rL.lot=HmscRandomLevel(单位=级别(StudyDesign$Plot))。
# Levels()只选择每个级别(这里是每个绘图)一次，因此在定义随机效果结构时不应重复。
# 这是一个与此案例研究无关的附带评论，有关此主题的更多信息，请参阅关于分层研究设计的示例。

# When defining the object m, we imported the data as Y = 1*(Y>0) to truncate it to presence-absence
# Corresponding to this, we assumed the probit link function
#定义对象m时，我们将数据导入为Y=1*(Y&gt;0)，将其截断为存在-不存在。
#与此相对应，我们假设概率位链接函数

# It is always a good idea to look at how the XData and XFormula are translate to the design matrix X. 
# These can be seen as e.g. by evaluating the following
#查看XDATA和XFormula如何转换为设计矩阵X总是一个好主意。
#这些可视为例如通过评估以下各项

m$XFormula
head(m$XData)
head(m$X)

# We observed that DC is treated as factor (as it should), with DC=1 being the baseline level and thus not visible in the matrix X
# Further, log(readcount) is treated as continuous covariate (as it should)
#我们观察到DC被视为因子(理应如此)，DC=1是基线水平，因此在矩阵X中不可见。
#此外，log(Readcount)被视为连续协变量(理应如此)

# It is always a good idea to look at the model object as well:
#同时查看模型对象始终是个好主意：

m

# This should give "Hmsc object with 99 sampling units, 31 species, 5 covariates, 1 traits and 1 random levels"
# Note that "5 covariates" corresponds to the columns of the X matrix. We have two environmental variables,
# but with adding the intercept and expanind the factor of DC into three columns, there are in total 5 "covariates"
# Note further that there is "1 trait" even if we did not define any traits. To see why this is the case, evaluate
#这应提供\“hMSC对象，具有99个采样单位、31个物种、5个协变量、1个特性和1个随机水平\”
#请注意，\“5个协变量\”对应于X矩阵的列。我们有两个环境变量，
#但加上截距并将DC因子展开为三列，总共有5个\“协变量\”
#请进一步注意，即使我们没有定义任何特性，也存在\“1个特性\”。要了解为什么会出现这种情况，请评估

head(m$Tr)

# This shows that the trait matrix contains the intercept, which is the "1 trait". 
# The intercept models the mean response of the species to the environmental covariates.

# Above we exemplified how to construct a single model, and we recommend the users of Hmsc to typically define their models one by one.
# But when the collection of the models is "ready", it is often convenient to treat all of them within a single loop.
# In this way, e.g. the code for model fitting or evaluating MCMC convergence needs to be written only once.
# Thus, i the scripts below, the object models will include all six models that were discussed in the book chapter, organised as a list of lists, 
# so the models[[i]][[j]] is the model type i=1,2,3 for the choice of explanatory variables j=1,2.
#上面我们举例说明了如何构建单个模型，我们建议hMSC的用户通常逐个定义自己的模型。
#但当模型集合\“就绪\”时，在单个循环中处理所有模型通常是很方便的。
#这样，例如，只需编写一次模型拟合或评估MCMC收敛性的代码。
#因此，在下面的脚本中，对象模型将包括本书章节中讨论的所有六个模型，这些模型被组织为列表列表。
#所以模型[[i]][[j]]对于解释变量j=1，2的选择是模型类型i=1，2，3。

# The model type i=1 is a lognormal Poisson model that is fitted to the sequence count data as they are.
# So this model will simultaneously account both for whether the species is present or not as well as how abundant it is when it is present.
# The model types i=2 and i=3 form together a hurdle model that separates presence-absence variation from abundance variation. 
# Thus, model i=2 is a probit model that is fitted to sequences counts truncated to presence-absence, whereas model i=3 is a normal model that is fitted to log-transformed sequence counts conditional on presence.
# This means that for i=3, sampling units where the species is not present are shown in the Y matrix as missing data (NA) rather than as zeros.
# 模型类型i=1是对数正态泊松模型，该模型按原样适用于顺序计数数据。
# 所以这个模型将同时考虑物种是否存在，以及当物种存在时物种的丰富程度。
# 模型类型i=2和i=3一起形成了将存在-缺席变化与丰度变化分开的栅栏模型。
# 因此，模型i=2是适合于截断到存在-不存在的序列计数的概率模型，而模型i=3是适合于以存在为条件的对数变换的序列计数的正常模型。
# 这意味着对于i=3，物种不存在的采样单元在Y矩阵中显示为缺失数据(NA)而不是零。

# Concerning the choices on the explanatory variables, we always include log-transformed sequencing depth, as this variable measures the total observation effort.
# Sequencing depth is the only variable included with j=1, and thus in these models we will estimate raw co-occurrences.
# With j=2, we additionally include the categorical variable of the decay class, and thus in these models we will estimate residual co-occurrences.
# We note that there are also many other properties of the log than the decay class that are likely to influence fungal occurrences, such as the diameter of the log.
# However, we ignore the other variables, as the data was specifically sampled to include much variation in decay class but less variation in other properties such as diameter.
#关于解释变量的选择，我们总是包括对数变换测序深度，因为该变量衡量的是总观测工作量。
#测序深度是j=1包含的唯一变量，因此在这些模型中，我们将估计原始的共现。
#当j=2时，我们另外包括衰变类的分类变量，因此在这些模型中我们将估计剩余的共现。
#我们注意到，除了腐烂级别之外，原木还有许多其他属性可能会影响真菌的发生，例如原木的直径。
#然而，我们忽略了其他变量，因为数据是专门抽样的，以包括衰变级别的较大变化，但其他属性(如直径)的变化较小。





# In all models, we also a random effect at the sampling unit level. 
# The random effect models associations among the species, which is what we are primarily interested about.
#在所有模型中，我们也在抽样单位水平上进行随机效应。
#随机效应模拟了物种之间的联系，这是我们主要感兴趣的。

models = list()
for (i in 1:3){
  Y = as.matrix(YData)
  if (i==2) {Y = 1*(Y>0)}
  if (i==3) {
    Y[Y==0] = NA
    Y = log(Y)
  }
  tmp = list()
  for (j in 1:2){
    XFormula = switch(j, ~1 + log(readcount), ~DC + log(readcount))
    m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
             studyDesign = studyDesign, ranLevels = list(sample = rL),
             distr=switch(i,"lognormal poisson","probit","normal"),
             YScale = TRUE)
    tmp[[j]] = m
  }
  models[[i]] = tmp
}

# In the above script, we have used the option Yscale = TRUE to scale the response data to zero mean and unit variance.
# It is discussed in more detail in Section 8.3 of the book, this scaling influences only the normal model,
# and it is done to make the default priors of Hmsc compatible with the data.
#在上面的脚本中，我们使用选项Yscale=true将响应数据缩放到零均值和单位方差。
#在本书的第8.3节中进行了更详细的讨论，此缩放仅影响法线模型，
#并且这样做是为了使hMSC的默认优先级与数据兼容。

# We will fit each of the six models so that we store 100 posterior samples for each of two chains
# We note that for more "final" results, one might wish to have e.g. 1000 samples for each of four chains
#我们将拟合这六个模型中的每一个，以便为两个链中的每一个存储100个后验样本。
#我们注意到，对于更多的\“最终\”结果，人们可能希望四条链中的每条链都有1000个样本

nChains = 2
samples = 100

# We next loop over both the model types as well as the selections of explanatory variables to fit all the six models.
# After fitting all models, we save the models object (including the six fitted model objects) to a file
# Loading the fitted models then serves as the starting point for exploring the results
# The script runs over a loop where thin is first 1, then 10, then 100, and so on
# Thin is the thinning parameter of the MCMC chain.
# The transient (also called burn-in) is set to 50*thin
# When thin = 1, there will be 50 burn-in and 100 actual iterations. All actual iterations are stored.
# When thin = 10, there will be 500 burn-in and 1000 actual iterations. The actual iterations are thinned by 10, so 100 are stored.
# When thin = 100, there will be 5000 burn-in and 10000 actual iterations. The actual iterations are thinned by 100, so 100 are stored.
# A long MCMC chain is needed to achieve convergence
# Thinning is applied to avoid storing model objects of very large size
# Even if in the end thin = 1000 is required to achieve converge, We recommend to run the loop thin = 1, 10, 100, 1000
# This is for several reasons.
# First of all, it will not be known beforehand how much thinning is needed to achieve satisfactory convergence
# Second, thin = 1 will run very fast, whereas thin = 1000 will take very long (1000 times longer)
# After thin = 1 is completed, it is already possible to develop all the remaining scripts that explore the fitted model
# When exploring the fitted model, often one realizes changes that need to be made, even if the fitting has not converged
# Third, running the model fitting for thin = 1, 10, 100, 1000 does not take much longer than running it just for thin = 1000 (it takes ca. 12% longer)
# Thus, in summary, running the model fitting for thin = 1, 10, 100, 1000 typically saves a lot of time,
# as it allows one to proceed fast in writing (and revising) all the scripts that are needed from defining the model to producing the result tables and figures
# The idea is not to run the entire loop in one go, as that would take a lot of time. Just run thin = 1, and then move to develop the next scripts. 
# You may then leave the remaining part of the loop (e.g. thin = 10, 100, 1000) to run e.g. overnight
# 接下来，我们将遍历模型类型以及解释变量的选择，以适应所有六个模型。
# 拟合完所有模型后，我们将Models对象(包括6个拟合的模型对象)保存到一个文件中。
# 然后加载拟合的模型作为探索结果的起点。
# 脚本运行一个循环，其中Thin首先是1，然后是10，然后是100，依此类推。
# Thin是MCMC链的细化参数。
# 瞬变(也称为老化)设置为50*Thin。
# 当Thin=1时，将有50次老化和100次实际迭代。所有实际小版本都会存储。
# 当Thin=10时，将有500次老化和1000次实际迭代。实际迭代减少了10次，因此存储了100次。
# 当Thin=100时，将有5,000次老化和10000次实际迭代。实际迭代减少了100次，因此存储了100次。
# 需要一条较长的MCMC链条才能实现收敛。
# 应用细化可避免存储非常大的模型对象。
# 即使最终需要Thin=1000才能实现聚合，我们也建议运行循环Thin=1、10、100、1000。
# 这有几个原因。
# 首先，要达到令人满意的收敛，需要进行多少细化，这一点事先并不为人所知。
# 其次，Thin=1会运行得非常快，而Thin=1000会花费很长时间(长1000倍)。
# Thin=1完成后，已经可以开发探索拟合模型的所有剩余脚本。
# 在探索拟合的模型时，人们通常会意识到需要进行的更改，即使配件尚未收敛。
# 第三，运行适合Thin=1、10、100、1000的模型并不比仅运行Thin=1000花费的时间长(大约多12%)。
# 因此，总而言之，运行适合Thin=1、10、100、1000的模型通常节省了大量时间，
# 因为它允许您快速编写(和修改)从定义模型到生成结果表格和图形所需的所有脚本。
# 我们的想法是不要一次运行整个循环，因为这会花费很多时间。只需运行Thin=1，然后移动以开发下一个脚本。
# 然后您可以让循环的其余部分(例如Thin=10、100、1000)运行，例如通宵运行

for (thin in c(1,10)){
  transient = 50*thin
  for (i in 1:3){
    for (j in 1:2){
      cat("model = ",i, ", modeltype = ",j,"\n",sep="")
      models[[i]][[j]] = sampleMcmc(models[[i]][[j]], thin = thin, samples = samples, transient = transient,
                                    nChains = nChains, nParallel = nChains, initPar = if(i==3) {NULL} else {"fixed effects"})
    }
  }
  filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(models,file=filename)
}

# 只运行一次Thin=1
for (thin in c(1)){
  transient = 50*thin
  for (i in 1:3){
    for (j in 1:2){
      cat("model = ",i, ", modeltype = ",j,"\n",sep="")
      models[[i]][[j]] = sampleMcmc(models[[i]][[j]], thin = thin, samples = samples, transient = transient,
                                    nChains = nChains, nParallel = nChains, initPar = if(i==3) {NULL} else {"fixed effects"})
    }
  }
  filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(models,file=filename)
}



# MCMC convergence can be difficult to achieve especially in those models that are not base on normal distribution
# For this reason, in the script above we initialize the "lognormal poisson" (i=1) and "probit" (i=2) models with initPar="fixed effects", 
# with which option the MCMC chains are not started from locations randomized from the prior but from a maximum likelihood solution to the fixed-effects part of the model
# MCMC收敛可能很难实现，特别是在那些不是基于正态分布的模型中。
# 因此，在上面的脚本中，我们使用initPar=\“Fixed Effects\”来初始化\“logNormal Poisson\”(i=1)和\“Probit\”(i=2)模型，
# 利用该选项，MCMC链不是从从先前随机的位置开始的，而是从模型的固定效果部分的最大似然解开始的





