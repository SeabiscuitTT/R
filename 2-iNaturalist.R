# 20221104 WY 
# 使用R下载常用数据库的分布数据及其他
# 下面主要是iNaturalist https://www.inaturalist.org/ 
# 下面代码主要比较了iNaturalist 和 GBIF

rm(list = ls()) #清除环境变量
gc()

# 安装
# install.packages("rinat")

# 导入
library(rinat)

#默认只返回100个观测记录
# cai <- get_inat_obs(query = "Pieris rapae") 

# maxresults 参数设置记录的数量最高不超过1万,数量多了会比较慢
# cai <- get_inat_obs(query = "Pieris rapae",maxresults=9999) 

# iNaturalist数据下载------------------------------
# 使用bounds参数限制地理范围（框，不是行政区域）
# minx=118, miny=27, maxx=123, maxy=31.25     
# 南纬、西经、北纬和东经
bounds <- c(27, 118, 31.25, 123) # 浙江全省
cai_inat <- get_inat_obs(query = "Pieris rapae", #菜粉蝶为例
                    maxresults = 9999,
                    bounds = bounds)
# 查看数据列名
colnames(cai_inat)
dim(cai_inat)    # 去重前 41 37
unique(cai_inat) #去重一下
dim(cai_inat)    # 去重后 41 37 没有重复记录

# 建议质控
# 【num_identification_agreements】列会显示有多少人对这条记录进行了鉴定
# 【quality_grade】列里面会有几个等级：research\needs_id\casual,
# research是两个以上人对一个观察鉴定相同且到种级，是可信度最高的，
# needs id是没有到种或只有一个人鉴定的，
# casual  是缺少信息或不符合正规观察的要求，可信度最低

library(dplyr)
# 挑选出研究级别的数据
cai_research <- filter(cai_inat,quality_grade == 'research') 

# 对比一下gbif------------------------------
library(rgbif)
library(beepr)
cai_gbif <- name_lookup(query="Pieris rapae", rank="species",limit=99999)
key <- cai_gbif[["data"]][["key"]]
res <- occ_download(
  pred_in("taxonKey", key),   #key就是上面得到的所有key
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',#描述保存标本的发生记录
                             'HUMAN_OBSERVATION', #多人观测记录
                             'OBSERVATION', #描述观察结果的发生记录
                             'LIVING_SPECIMEN',#活体样本
                             'MACHINE_OBSERVATION', #由机器进行的观察结果
                             'LITERATURE')), #文献记录
  # pred_in("country",  c("CN","TW")), #注意台湾数据需要单独写地名,目前是中国,如果需要全世界的分布数据,就把这行注释
  pred_within(gbif_bbox2wkt(minx=118, miny=27, maxx=123, maxy=31.25)),#目标范围
  pred("hasCoordinate", TRUE), #只挑选有坐标的记录
  pred("hasGeospatialIssue", FALSE),
  # pred_gte("year", 1970), #限定数据开始的年份，为了匹配Worldclim上1970-2000的环境图层
  format = "SIMPLE_CSV",  #结果文件的形式
  user="wangying", pwd="wy211994", email="1187468386@qq.com") #账户信息
occ_download_meta(res)  #这步是向API提交问询，需要重复多次，直到Status: SUCCEEDED
z <- occ_download_get(res)   # 把数据下载到本地
beepr::beep(8) # 上一步有时会比较慢，这选中上面和这两步一起运行，那么上面的执行完了就会有提示
df <- occ_download_import(z) # 载入数据


# 挑选出gbif中来源是iNaturalist的数据
cai_gbif <- unique(df[,c('species','decimalLongitude','decimalLatitude',
                         'institutionCode')])
cai_gbif <- cai_gbif[-16,-4]
# institutionCode == 'iNaturalist'

# 取并集-----------------------------
names(cai_gbif) <- c('species','lon','lat')
cai_gbif$from <- 'gbif_iNaturalist'
cai_inat_u <- unique(cai_inat[,c(1,6,5,33)])
names(cai_inat_u) <- c('species','lon','lat','from')
# df3 <- unique(rbind(df1,df2))
bind <- rbind(cai_gbif,cai_inat_u)

out <- filter(bind,from %in% c('gbif_iNaturalist','research'))
out1 <- unique(out[,1:3]) #这里需要注意，有些看着一样，但是由于小数点保留的问题，其实是不一样的。
# 但是这个不一样是数据库之间传输的时候造成的，
# 而不是原始搜集数据时就是非常近的。
# 通过设置统一的小数点位数进行过滤
# 统一保留3位
out2 <- data.frame(out$species,round(out$lon,3),round(out$lat,3))
out_final <- unique(out2)
write.csv(out_final,'组合数据.csv')
