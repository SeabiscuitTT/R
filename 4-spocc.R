# 20221104 WY 
# 使用R下载常用数据库的分布数据及其他
# 下面主要是使用spocc集中下载某些数据库的分布数据

rm(list = ls()) #清除环境变量
gc()

# 从GBIF下载某种植物在中国发布的位置信息---------
library(spocc)
library(dplyr)

# 从gbif获取指定物种的信息,from后面可以接9种指定的数据库，包括动植物。
df <- occ(query = 'Theretra alecto',
          from = c('gbif','inat'),
          has_coords = TRUE,#只要有坐标的数据
          limit = 9999 )  # 50 此处只是先试试，看有没有记录（提高速度）
df1 <- occ2df(df) #提取坐标信息
# 虽然上面参数设置里面我们参数设置了只需要有坐标的
# 但是，结果里面还是有存在坐标数据列是NA的情况
# 剔除那些坐标位置为NA的行
out <- na.omit(df1[,1:4]) #只选1：4列，不要包含后面的，否则只要时间里面也有NA，那么就算前面数据都齐全也会被删除

# 保存数据
write.xlsx(out,'4_spocc.xlsx')
