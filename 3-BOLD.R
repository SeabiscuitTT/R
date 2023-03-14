# 20221104 WY 
# 使用R下载常用数据库的分布数据及其他
# 下面主要是BOLD http://www.boldsystems.org/ 
# 下面代码除了从BOLD下载分布数据，还可以下载序列

rm(list = ls()) #清除环境变量
gc()


library(bold)
library(openxlsx)
library(tidyr)
library(beepr)

# 分布数据----------------】-----

# 科-----
# a <- bold_specimens(taxon='Nolidae') #瘤蛾科测试 不限制地理范围时间会比较长
a <- bold_specimens(taxon='Nolidae',geo='China') #限制样本采集地为中国，不限制默认为全世界
beepr::beep(8)#因为BOLD比较慢（相比于GBIF），这边设置一个程序运行完的提示音
Location <- unique(a[c(22,47,48)]) # 47,48是纬度和经度，去重可能会有一行是空值
good <- complete.cases(Location) #找到非空值,这里都是逻辑值
Location_family <- Location[good, ] #提取非空值所有行，就是所有非重复坐标
# 保存数据
#如果物种名一栏是空的，表明只定到了属，没定到种，可根据分析目的合理筛选
write.xlsx(Location_family,'Location_family.xlsx')
# 作图
# 请使用【物种批量出图.R】脚本作图

# 属-----------------
# b <- bold_specimens(taxon='Meganola') #瘤蛾科 Meganola属 测试
b <- bold_specimens(taxon='Meganola',geo='China')
beepr::beep(sound = "mario")
Location <- unique(b[c(22,47,48)]) # 47,48是纬度和经度，去重可能会有一行是空值
good <- complete.cases(Location) #找到非空值,这里都是逻辑值
Location_genus <- Location[good, ] #提取非空值所有行，就是所有非重复坐标
# 保存数据
#如果属名一栏是空的，表明只定到了科，没定到属，可根据分析目的合理筛选
write.xlsx(Location_genus,'Location_genus.xlsx')
# 作图
# 请使用【物种批量出图.R】脚本作图

# 种-----------------
# 如果是一个物种可以向上面科属类似，把“taxon=”后面的内容换成某个物种就行
# 但是一般我们需要批量查询多个物种，可以采用下面的循环
# species_list.xlsx 存入需要问询的物种，格式参考实例文件
# GBIF代码是使用的文本，这边是表格，随意。
species <- as.vector(unlist(read.xlsx('species_list.xlsx',colNames = F)))

Location_sp <- data.frame()#新建一个空数据框，用于分别存放循环中产生的坐标数据

# i="Theretra alecto" 
for(i in species){
  # 如果问询的物种在库内没有记录，循环结束之后会有报错，不用在意
  a <- bold_specimens(taxon = i) # bold_specimens返回的是一个data.frame
  if(nrow(a) < 1){next}else{
    # 坐标信息
    Location <- unique(a[c(22,47,48)]) # 47,48是纬度和经度，去重可能会有一行是空值
    good <- complete.cases(Location) #找到非空值
    Location <- Location[good, ] #提取非空值所有行，就是所有非重复坐标
    Location_sp <- rbind(Location_sp,Location)
    
    print(i)  # 看看运行到第几个了（对于异常终止比较有用）
  }
}
beepr::beep(8)

# 保存数据
write.xlsx(da_Location,'Location_sp.xlsx')

# 序列获取---------------------】------
rm(list = ls())
gc()
# 某个物种的所有序列--------------------
species = 'Theretra alecto' #你需要问询的物种
a <- bold_seq(taxon = species) # bold_seq返回的是一个包含4个list(id name gene sequence)的list
beep(8)
length(a) # a包含的list的个数，下面一行命令会用到
b <- as.data.frame(matrix(unlist(a), 
                          nrow=length(a), #nrow是根据a中包含的list的个数定的
                          byrow=T)) 
# 1和3列是重复的，把第一列删除
b <- b[,-1]
# 把每行合并起来，新生成一列放在b后面，方便后续的使用
names(b) <- c('name','gene','sequence') # 先把数据框重命名
c <- tidyr::unite(b,'name_gene',name,gene) #remove = FALSE 参数则不会删除原来的数据列
names(c) <- c('ID','sequence') # 重命名一下列名
d <- tidyr::unite(c,'ID,sequence',ID,sequence,sep = ",")#这边用,隔开,方便后续按照,替换为换行符做成.fasta格式
e <- data.frame(matrix(rep('>',length(a)), nrow=length(a), byrow=T))
f <- cbind(e,d)
names(f) <- c('start','ID_sequence')
g <- tidyr::unite(f,'start,ID_sequence',start,ID_sequence,sep = "")#起到的效果就是在d的每一行开头添加一个>
names(g) <- c('sequences')
filename <- paste0('2_','某个物种的所有序列_',species,".fas")
out <- write.table(g,filename,row.names = F,col.names = F)

# 接下来需要在文本编辑器中把所有引号都替换为空，把英文的逗号替换为换行符
# 注意：物种拉丁名中间是空格，要注意。

# 多个物种的所有序列-----------------
rm(list = ls())
species <- as.vector(unlist(read.xlsx('species_list.xlsx',colNames = F)))
da <- data.frame()   #新建一个空数据框，用于存放循环中产生的da
for(i in species) {
  species = i #你需要问询的物种
  a <- bold_seq(taxon = species) # bold_seq返回的是一个包含4个list(id name gene sequence)的list
  length(a) # a包含的list的个数，下面一行命令会用到
  if(length(a)==0){next}else{ #如果需要问询的物种没有序列，则length(a)=0，则跳过该物种
    b <- as.data.frame(matrix(unlist(a), nrow=length(a), byrow=T)) #nrow是根据a中包含的list的个数定的
    # 1和3列是重复的，把第一列删除
    b <- b[,-1]
    # 把每行合并起来，新生成一列放在b后面，方便后续的使用
    names(b) <- c('name','gene','sequence') # 先把数据框重命名
    c <- tidyr::unite(b,'name_gene',name,gene) #remove = FALSE 参数则不会删除原来的数据列
    names(c) <- c('ID','sequence') # 重命名一下列名
    d <- tidyr::unite(c,'ID,sequence',ID,sequence,sep = ",")#这边用,隔开,方便后续按照,替换为换行符做成.fasta格式
    e <- data.frame(matrix(rep('>',length(a)), nrow=length(a), byrow=T))
    f <- cbind(e,d)
    names(f) <- c('start','ID_sequence')
    g <- tidyr::unite(f,'start,ID_sequence',start,ID_sequence,sep = "")#起到的效果就是在d的每一行开头添加一个>
    names(g) <- c('sequences')
    da <- rbind(da,g)
  }
}
out <- write.table(da,'2_多个物种所有条码.fas',row.names = F,col.names = F)
# 接下来需要在文本编辑器中把所有引号都替换为空(就删删除)，把英文的逗号替换为换行符
# 注意：物种拉丁名中间是空格。
# 如果需要特定属的，只需要把物种名换成属名就行,那每一个小的list就是一个物种。

# 多个物种各选1条条码------------------------------------------------------------------

# 某个物种的第一条序列（也就是选一条序列来代表这个物种）,当然输入肯定不只是一个物种
# 那就在上面代码的基础上
rm(list = ls())
species <- as.vector(unlist(read.xlsx('species_list.xlsx',colNames = F)))
da <- data.frame()   #新建一个空数据框，用于存放循环中产生的da
for(i in species) {
  species = i #你需要问询的物种
  s <- bold_seq(taxon = species)
  if(length(s)==0){next}else{ #如果需要问询的物种没有序列，则length(s)=0，则跳过该物种
    a <- bold_seq(taxon = species)[[1]] # bold_seq返回的是一个包含4个list(id name gene sequence)的list
    b <- as.data.frame(matrix(unlist(a), nrow=1, byrow=T)) 
    b <- b[,-1]
    # 把每行合并起来，新生成一列放在b后面，方便后续的使用
    names(b) <- c('name','gene','sequence') # 先把数据框重命名
    c <- tidyr::unite(b,'name_gene',name,gene) #remove = FALSE 参数则不会删除原来的数据列
    names(c) <- c('ID','sequence') # 重命名一下列名
    d <- tidyr::unite(c,'ID,sequence',ID,sequence,sep = ",")#这边用,隔开,方便后续按照,替换为换行符做成.fasta格式
    e <- data.frame(matrix(rep('>',length(a)), nrow=length(a), byrow=T))
    f <- cbind(e,d)
    names(f) <- c('start','ID_sequence')
    g <- tidyr::unite(f,'start,ID_sequence',start,ID_sequence,sep = "")#起到的效果就是在d的每一行开头添加一个>
    names(g) <- c('sequences')
    da <- rbind(da,g)
  }
}
beep(8)
da <- unique(da)
out <- write.table(da,'2_多个物种每种1条条码.fas',row.names = F,col.names = F)
# 接下来需要在文本编辑器中把所有引号都替换为空，把英文的逗号替换为换行符
# 注意：物种拉丁名中间是空格，要注意。

# 多个物种挑选固定数量的条码------------------------------------------------------------------
rm(list = ls())
species <- as.vector(unlist(read.xlsx('species_list.xlsx',colNames = F)))
da <- data.frame()   #新建一个空数据框，用于存放循环中产生的da
n = 10 #这边是你需要的每个物种重复条码的数量，这边以1为例
for(i in species) {
  s <- bold_seq(taxon = i) # bold_seq返回的是一个包含4个list(id name gene sequence)的list
  length(s) # a包含的list的个数，下面一行命令会用到
  if(length(s)==0){next      #如果需要问询的物种无条码,则跳过
  }else if(length(s) < n){ #如果需要问询的物种一共的条码数都不到n,则全取
    a <- bold_seq(taxon = i) # bold_seq返回的是一个包含4个list(id name gene sequence)的list
    length(a) # a包含的list的个数，下面一行命令会用到
    b <- as.data.frame(matrix(unlist(a), nrow=length(a), byrow=T)) #nrow是根据a中包含的list的个数定的
    # 1和3列是重复的，把第一列删除
    b <- b[,-1]
    # 把每行合并起来，新生成一列放在b后面，方便后续的使用
    names(b) <- c('name','gene','sequence') # 先把数据框重命名
    c <- tidyr::unite(b,'name_gene',name,gene) #remove = FALSE 参数则不会删除原来的数据列
    names(c) <- c('ID','sequence') # 重命名一下列名
    d <- tidyr::unite(c,'ID,sequence',ID,sequence,sep = ",")#这边用,隔开,方便后续按照,替换为换行符做成.fasta格式
    e <- data.frame(matrix(rep('>',length(a)), nrow=length(a), byrow=T))
    f <- cbind(e,d)
    names(f) <- c('start','ID_sequence')
    g <- tidyr::unite(f,'start,ID_sequence',start,ID_sequence,sep = "")#起到的效果就是在d的每一行开头添加一个>
    names(g) <- c('sequences')
    da <- rbind(da,g)
  }else{ 
    a <- bold_seq(taxon = i)[c(1:n)]# list取子集
    b <- as.data.frame(matrix(unlist(a), nrow=n, byrow=T)) 
    # 1和3列是重复的，把第一列删除
    b <- b[,-1]
    # 把每行合并起来，新生成一列放在b后面，方便后续的使用
    names(b) <- c('name','gene','sequence') # 先把数据框重命名
    c <- tidyr::unite(b,'name_gene',name,gene) #remove = FALSE 参数则不会删除原来的数据列
    names(c) <- c('ID','sequence') # 重命名一下列名
    d <- tidyr::unite(c,'ID,sequence',ID,sequence,sep = ",")#这边用,隔开,方便后续按照,替换为换行符做成.fasta格式
    e <- data.frame(matrix(rep('>',length(a)), nrow=length(a), byrow=T))
    f <- cbind(e,d)
    names(f) <- c('start','ID_sequence')
    g <- tidyr::unite(f,'start,ID_sequence',start,ID_sequence,sep = "")#起到的效果就是在d的每一行开头添加一个>
    names(g) <- c('sequences')
    da <- rbind(da,g)
  }
}
beep(8)
write.table(da,'2_多个物种每种固定数量条码.fas',row.names = F,col.names = F)


