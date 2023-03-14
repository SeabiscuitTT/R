# 20221104 WY 
# 使用R下载常用数据库的分布数据及其他
# 下面主要是GBIF https://www.gbif.org/ 

rm(list = ls()) #清除环境变量
gc()

#导包，若没有请自行安装
library(rgbif) 
library(readr)
library(openxlsx)
library(beepr)  #代码完成后发出提示音 

# step1 结合分类需求查询物种的key--------
# 如果需要总科级的数据，直接像科、属、种，
# 这样下载是无法直接操作的，
# 建议列出总科下面的科，用科水平代码下载。

## Search for a family----------
a <- name_lookup(query="Nolidae", rank="family",limit=99999) #瘤蛾科
# 会返回所有这个科物种数据的list
#从list里面挑选出物种的key
key <- a[["data"]][["key"]] 

## Search for a genus---------
b <- name_lookup(query="Theretra", rank="genus",limit=99999) #斜纹天蛾属
key <- b[["data"]][["key"]] 

## Search for a species-------
c <- name_lookup('Parasa consocia',rank="species",limit=99999) #褐边绿刺蛾
key <- c[["data"]][["key"]] 

## Search for many species---------
# 物种过多不方便如上写成字符向量
# 将需要查询的物种名逐行写入到TXT文件，
# 不需要表头
f <- read.table("species_list.txt",
                sep = "\t", header = F) 
#将读取的拉丁名结果转换成向量
d <- as.vector(unlist(f[1])) 
#新建一个空列表用于存放循环中产生的key
ls <- list()   
# 下面这个循环主要是用来产生key,并保存到list中
for ( i in d ) {
  print(i)
  s <- print(name_backbone(name = i)$speciesKey)
  ls <- list(ls,s)
}
#将所有产生的key，转化为一个向量
key <- as.vector(unlist(ls)) 

# step2 用上面得到的key下载需要的分布数据--------
res <- occ_download(
  pred_in("taxonKey", key),   #key就是上面得到的所有key
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',#描述保存标本的发生记录
                             'HUMAN_OBSERVATION', #多人观测记录
                             'OBSERVATION', #描述观察结果的发生记录
                             'LIVING_SPECIMEN',#活体样本
                             'MACHINE_OBSERVATION', #由机器进行的观察结果
                             'LITERATURE')), #文献记录
  # pred_in("country",c("CN","TW")), #注意台湾数据需要单独写地名,目前是中国,如果需要全世界的分布数据,就把这行注释
  pred_within(gbif_bbox2wkt(minx=118, miny=27, maxx=123, maxy=31.25)),#目标范围，可以缩小检索范围
  pred("hasCoordinate", TRUE), #只挑选有坐标的记录
  pred("hasGeospatialIssue", FALSE),
  # pred_gte("year", 1990), #可以限定数据开始的年份，这边不需要
  format = "SIMPLE_CSV",  #结果文件的形式
  user="wangying", pwd="wy211994", email="1187468386@qq.com") #你的账户信息
occ_download_meta(res)  #这步是向API提交问询，需要重复多次，直到Status: SUCCEEDED
# z <- occ_download_get(res,overwrite=T)   # 把数据下载到本地(覆盖写入)
z <- occ_download_get(res)   # 把数据下载到本地(不进行覆盖写入)
beepr::beep(8) # 上一步有时会比较慢，这选中上面和这两步一起运行，那么上面的执行完了就会有提示
df <- occ_download_import(z) # 载入数据

colnames(df) #查看所有列名，找到表示经纬度的列名
df$decimalLongitude  #粗略查看维度
df$decimalLatitude   #粗略查看经度
#挑出物种及经纬度信息
df1 <- df[,c('taxonKey','speciesKey','species','scientificName','verbatimScientificName',
             'decimalLongitude','decimalLatitude')]
head(df1) #查看结果前6行，看看结构
dim(df1) #产看记录的总数
df2 <- unique(df1) #去重
dim(df2) #产看去重后的总数,如果问询的是科或者属水平的数据，speciesKey可能是NA,也就是没办法定到种
# df2 =na.omit(df2) #如果不想要没定到种的数据，可以运行这一行
# dim(df2)
head(df2) #查看结果前6行，看看结构

# step3 数据清洗--------
# 上面简单去重
# 下面使用CoordinateCleaner包删除有问题的分布记录
library(CoordinateCleaner)
df3 <- subset(df2,decimalLongitude !=0 & decimalLatitude != 0)# 挑选出所有经纬度都是不是0的行
dim(df3)

# 检测记录点坐标是否围绕首都、国家的中心，是否落入海洋，为零，
# 或在饲养动物的博物馆(机构)周围。
problem_records<- clean_coordinates(df3, lon = "decimalLongitude", lat = "decimalLatitude",
                                    species = "Species",
                                    tests = c("capitals", "centroids","equal","institutions","zeros", "seas"))

# 边检查边安装一些库文件
summary(problem_records) # 会显示有问题的分布点

# 从数据集中排除所有有问题的记录
data_clean <- df3[which(problem_records$.summary== "TRUE"),]
str(data_clean)

# 保存数据---------------
write.xlsx(data_clean,"data_clean.xlsx")