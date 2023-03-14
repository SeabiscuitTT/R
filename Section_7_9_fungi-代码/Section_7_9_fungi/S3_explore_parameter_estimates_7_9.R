# THIS SCRIPT EXPLORES THE PARAMETER ESTIMATES OF HMSC MODELS OF THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# The preliminaries are as in script S1

# Change the working directory if necessary
wd = "C:/Users/Google/Documents/R/Joint Species Distribution Modelling in R/section_7_9_fungi_2020_05_31/Section_7_9_fungi"
setwd(wd)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
library(Hmsc)
set.seed(1)

# We first read in the object models that includes all six fitted models
# You may change the thin parameter to the highest thin for which you have fitted the model, to get as reliable results as possible
#我们首先读入包括所有六个拟合模型的对象模型。
#您可以将Thin参数更改为适合模型的最高Thin，以获得尽可能可靠的结果

nChains = 2
samples = 100
thin = 1
filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We first plot the responses of the species to the environmental covariates in the models that include both sequencing depth and decay class (the second index is 2 in m = models[[i]][[2]]).

# By setting the value of i the lognormal Poisson model on sequence counts (i=1), the probit model on presence-absence data (i=2), and normal model on log-transformed sequence count conditional on presence (i=3). 
# Note that only the last plot may be shown in your screen, so you need to scroll back with the images to see the other ones
#我们首先在同时包含排序深度和衰减级的模型中绘制物种对环境协变量的响应(m=模型[[i]][[2]]中的第二个指数是2)。

#通过设置i的值，即关于顺序计数的对数正态泊松模型(i=1)、关于存在-缺席数据的概率模型(i=2)、
#以及关于有条件存在的对数变换后的顺序计数的正态模型(i=3)。
#请注意，您的屏幕上可能只会显示最后一张图，所以您需要回滚图片才能看到其他的图

i=2
m = models[[i]][[2]]
postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, param="Support", supportLevel = 0.95, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))


#如果绘图不显示或有问题请用这段代码关闭画图，重新运行上段代码
dev.off()


# In these heatmaps of estimated species niches. Red (respectively, blue) colour shows parameters that are estimated to be positive (respectively, negative) with at least 0.95 posterior probability.
# The intensity of the colour refers to the posterior mean estimate of the parameter. 
#在这些估计物种生态位的热图中。红色(分别为蓝色)表示估计为正(分别为负)且后验概率至少为0.95的参数。
#颜色的强度是指参数的后验平均估计值。

# For all model types, some species are estimated to respond positively to sequencing depth,
# whereas no species are estimated negatively to it. 
# This can be expected, as having more sequences in general means having more sequences for the focal species.
# As decay class is a categorical variable, so the responses of the species to that are somewhat hard to interpret from the beta-plot and best visualised through gradient plots.
# We next construct a gradient plot over decay classes to the presence-absence model (i=2) 
# For simplicity, we call this model(models[[2]][[2]]) as m
#对于所有模式类型，一些物种估计对测序深度有积极反应，
#而没有任何物种被估计为负面的。
#这是可以预料到的，因为拥有更多的序列通常意味着有更多的局部物种的序列。
#由于衰变类是一个分类变量，所以物种对此的反应很难从β图中解释出来，而通过梯度图表现得最好。
#接下来，我们构造一个关于衰变类的梯度图，以得到存在-不存在模型(i=2)。
#为简单起见，我们将此模型(Models[[2]][[2]])称为m

m = models[[2]][[2]]

# We explore the gradient over decay class, and thus we set focalVariable = "DC"
# We wish to normalize readcount to its mean over all the data, and thus we set non.focalVariables = list("readcount"=list(1))
# See help for the function constructGradient to see the options for non-focal variables
#我们探索衰减类上的渐变，因此我们设置FocalVariable=\“DC\”
#我们希望将所有数据的重新计数标准化为其平均值，因此我们设置了Non.foalVariables=list(\“readcount\”=list(1))。
#查看函数structGradient的帮助以查看非焦点变量的选项

Gradient = constructGradient(m,focalVariable = "DC",
                             non.focalVariables = list("readcount"=list(1)))

# It is usually a good idea to look at the Gradient object that has been constructed, to ensure that it is as it should be
# Thus, you may evaluate
#查看已构造的渐变对象通常是一个好主意，以确保它是应该的。
#因此，您可以评估

Gradient$XDataNew

# The most obvious part of this object is Gradient$XDataNew. This shows that predictions are to be made for four sampling units, that have decay classes 1,2,3 and 4. 
# All have exactly the same readcount to make them otherwise identical except the focal variable, i.e. the readcount
# The Gradient object has also information about the study design and the random effects. 
# The predictions are to be done for a "new_unit", meaning that the random effect has not been estimated for the unit but is randomized from its distribution.
#此对象最明显的部分是渐变$XDataNew。这表明要对衰变级别为1、2、3和4的四个采样单元进行预测。
#所有的读取计数都完全相同，以使它们在其他方面完全相同，除了焦点变量，即读取计数。
#渐变对象还包含有关研究设计和随机效果的信息。
#要对\“new_unit\”进行预测，这意味着尚未估计该单元的随机效应，但根据其分布随机化。
# We next make predicted species communities over this gradient
#接下来我们在这个梯度上预测物种群落
predY = predict(m,Gradient = Gradient, expected = TRUE)

# Let's explore this object:
#让我们来探索一下这个对象：
class(predY)

# It is a list...
#这是一份名单...
length(predY)

# ..of length 200. This is because the predictions are done for each posterior sample, and we stored 100 samples for each of the two chains, thus 200 in total.
#..长度为200。这是因为预测是针对每个后验样本进行的，我们为两个链中的每一个存储了100个样本，因此总共存储了200个样本。
dim(predY[[1]])

# each prediction has a dimension of 4 x 31, as the predictions are made for 4 sampling units (as explained above) for each of the 31 species.
#每个预测的维度为4x31，因为预测是针对31个物种中的每个物种的4个采样单位(如上所述)进行的。
head(predY[[1]])

# Examining the predictions shows that they are numbers between 0 and 1. 
# This is because we are exploring a presence-absence model, and hence the model predicts occurrence probabilities. 
# Note that when making the predictions, we set expected = TRUE, meaning that we asked for probabilities, not occurrences.
# If setting expected = FALSE, the predictions would be 0 or 1, and thus they would involve also the Bernoulli randomness around the occurrence probabilities.
# We next visualize the predictions with the plotGradient function, either for individual species
# (measure="Y", index selecting the species) or species richness (measure="S"). 
# With the option showData, one can decide whether the raw data is plotted as well. 
# With the option jigger, one can avoid overlapping points. 
# With the option prob, one can choose the credible intervals to be plotted 
#检查预测显示它们是介于0和1之间的数字。
#这是因为我们正在探索一个在场-缺席模型，因此该模型预测了发生的概率。
#请注意，在进行预测时，我们设置了Expect=True，这意味着我们要求的是概率，而不是事件。
#如果设置Expect=False，则预测将为0或1，因此它们还会涉及发生概率周围的伯努利随机性。
#接下来，我们使用plotGradient函数将预测可视化，无论是针对单个物种。
#(MEASURE=\“Y\”，物种选择指数)或物种丰富度(MEASURE=\“S\”)。
#使用选项showData，可以决定是否也绘制原始数据。
#使用选项jigger，可以避免重叠点。
#使用选项prob，可以选择要绘制的可信区间
prob = c(0.25,0.5,0.75)
plotGradient(m, Gradient, pred=predY, measure="Y", index=23, showData = TRUE, jigger = 0.15, prob = prob)
plotGradient(m, Gradient, pred=predY, measure="Y", index=16, showData = TRUE, jigger = 0.15, prob = prob)
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, jigger = 0.15, prob = prob)

# We observe that the overall species richness increases over the decay classes, 
# but e.g. the occurrences probabilities of the two species Phellopilus nigrolimitatus and Heterobasidion parviporum are lowest in the last decay class.
#我们观察到，总体物种丰富度随着腐烂级别的增加而增加，
#但例如，在最后一个腐烂类中，黑腹黄姑鱼和细小异尖线虫的出现概率最低。(2)在最后一类腐烂类中，黑腹扁平线虫和细小异尖线虫的出现概率最低。
# We next plot the species association networks. To do so, we need the corrplot package
#接下来我们绘制物种协会网络图。要做到这一点，我们需要corrlot包
library(corrplot)

# We plot the associations here only for the presence-absence model (i=2) for the version that does not account for decay class in the fixed effects (j=1)
#我们在此仅为不考虑固定效果(j=1)中衰变级别的版本的存在-缺席模型(i=2)绘制关联。
m = models[[2]][[1]]

# To plot the association network for any of the other models, simply pick another part of the models object
# We first use computeAssociations to derive the correlations from the object m
#要绘制任何其他模型的关联网络，只需选择Models对象的另一部分。
#我们首先使用ComputeAssociations从对象m导出相关性
OmegaCor = computeAssociations(m)

# Let's explore this object
#让我们探索一下这个物体
class(OmegaCor)

# It is a list..
#这是一份名单..
length(OmegaCor)

# ...of length 1. This is because we included only one random effect to the model, i.e. that at the sampling unit level
#...长度为1。这是因为我们只在模型中包含了一个随机效果，即在采样单元级别
class(OmegaCor[[1]])

# Also the first element of this list is a list. 
# It is a named list with posterior mean (mean) values and the posterior supports for the correlations being positive (support). 
# These can be explore numerically as follows:
#此列表的第一个元素也是列表。
#这是一个命名列表，其后验均值(Mean)值和相关性的后验支持度为正(支持度)。
#这些可以按以下数字进行探索：
head(OmegaCor[[1]]$mean)
head(OmegaCor[[1]]$support)

# We next plot the network of assocications, colouring only those that are either positive or negative with at least 95% posterior probability
#接下来，我们绘制关联网络，只给那些具有至少95%的后验概率的正负组合上色
supportLevel = 0.95
supportLevel
toPlot = ((OmegaCor[[1]]$support>supportLevel) + (OmegaCor[[1]]$support<(1-supportLevel))>0)* OmegaCor[[1]]$mean
corrplot(toPlot, method = "color", col=colorRampPalette(c("blue","white","red"))(200))

# With m = models[[2]][[1]] you will visualize the "raw associations" of the presence-absence model
# Repeating the above with m = models[[2]][[2]] you will visualize the "residual associations" of the presence-absence model
# There are many more raw association than residual associations. 
# This is typically to be expected, as the raw associations are also generated by
# differential responses of the species to the environmental conditions, here the decay class.
#with m=Models[[2]][[1]]您将可视化存在-缺席模型的\“原始关联\”
#使用m=Models[[2]][[2]]重复上述步骤，您将看到在场-缺席模型的\“残差关联\”
#原始关联比残留关联多得多。
#这通常是意料之中的，因为原始关联也是由。
#物种对环境条件的不同反应，这里指的是腐烂类。

# We may also visualize the association network as a model-based ordination using the biPlot function
# Let us first do that e.g. for the raw associations of the presence-absence model 
#我们还可以使用biPlot函数将关联网络可视化为基于模型的排序。
#让我们首先这样做，例如，对于在场-缺席模型的原始关联
m = models[[2]][[1]]
biPlot(m,etaPost = getPostEstimate(m,"Eta"),lambdaPost = getPostEstimate(m,"Lambda"), colVar = 1)

# The circles are the sampling units and the triangles the dots. 
# The sampling units are coloured accoring the decay class 
# because we selected colVar = 1 and decay class is the first column of XData
# 圆圈是取样单位，三角形是点。
# 采样单位的颜色是根据腐烂等级来确定的，因为我们选择了colVar = 1，而腐烂等级是XData的第一列
head(m$XData)

# Note that even if decay class is part of XData, in this model it is not part of XFormula,
# and thus not part of the X-matrix:
# 注意，即使腐烂等级是XData的一部分，在这个模型中它也不是XFormula的一部分。
# 因此也不是X矩阵的一部分。
head(m$X)

# Thus decay class is not part of the model, even if we have included information about it in the XData
# 因此腐烂等级不是模型的一部分，即使我们在XData中包含了它的信息。


# Let us do a similar plot for residual associations:
# 让我们为残差关联做一个类似的图。
m = models[[2]][[2]]
biPlot(m,etaPost = getPostEstimate(m,"Eta"),lambdaPost = getPostEstimate(m,"Lambda"), colVar = 1)
