# 20221104 WY 
# 使用R下载常用数据库的分布数据及其他
# 下面主要是对于搜集到的分布数据就行物种名标准化、分布数据清洗以及可视化


rm(list = ls()) #清除环境变量
gc()

# 物种名标准化---------------------
library(taxize) # 本流程使用的主要R包
library(openxlsx)

# 读取需要问询的物种名
data <- read.xlsx('species_list.xlsx')
name <- unique(data[,1]) # 去重一下避免重复

# 密钥前后一个空格都不能有
# 这是NCBI分配的密钥，可以去NCBI网页申请
# 没有密钥也没关系，就是慢一点
options(ENTREZ_KEY='你的密钥') 
ids_NCBI <- get_ids(name, db = 'ncbi')

# 得到物种在NCBI的taxid
out_id <- data.frame(unlist(ids_NCBI)) 

# 建立输入的物种名和taxid的对应关系，后面要用
out_id <- data.frame(rownames(out_id),out_id[,1]) 
colnames(out_id) <- c("input.name","taxid")

# 基于NCBI taxid 查询分类
calssf <- classification(unlist(ids_NCBI),db = 'ncbi')

# 写个循环处理一下上面得到的结果
d0 <- data.frame(rank=c("superfamily","family","genus","species"))
# i=1 #测试用的
for(i in 1:length(calssf)){
  if(is.na(calssf[i])){next}else{
    #倒数5行，确保想要的数据都在里面
    d1 <- calssf[[i]][-(1:(length(calssf[[i]]$rank)-5)),1:2] 
    #修改每个物种列名
    names(d1) <- c(names(calssf)[i],'rank') 
    d0 <- merge(d0,d1,by='rank',all = TRUE) 
  }
}

# 进一步处理循环的结果，如果没有到想要的，就增加倒数的行数
# 根据不同数据酌情调整
d <- d0[c(1:3),] #1到3行是科、属、种
d1 <- t(d) #转置
d2 <- data.frame(rownames(d1),d1) #让原始的taxid显示出来，下面好据此合并
colnames(d2)=c('taxid','NCBI.family','NCBI.genus','NCBI.species')
d2 <- d2[-1,] #去掉第一行

# 合并前面提交的物种名和后面的结果
out <- merge(out_id,d2,by='taxid',all = T)
write.xlsx(out,'物种名录校对.xlsx')
# 对输入的物种名和NCBI返回的结果比较，查看同物异名
# 后续操作在excel里面完成更方便

# NCBI返回结果与输入不符的情况有3种：
# 1.同物异名
# 2.分类变动（变了属之类）
# 3.NCBI数据库没有匹配到（建议直接保留之前的命名）


# 分布数据清洗--------------------
rm(list = ls()) #清除环境变量
gc()

# 使用 CoordinateCleaner 包删除有问题的分布记录
library(CoordinateCleaner)
df <- read.xlsx('4_spocc.xlsx')
df1 <- subset(df,longitude !=0 & latitude != 0)# 挑选出所有经纬度都是不是0的行
dim(df1)
# 为了防止坐标数据是文本，这边强制转化一下
df2 <- data.frame(df1$name,
                  as.numeric(df1$longitude),
                  as.numeric(df1$latitude))
names(df2) <- c('name','longitude','latitude')
# 检测记录点坐标是否围绕首都、国家的中心，是否落入海洋，为零，
# 或在饲养动物的博物馆(机构)周围。
problem_records<- clean_coordinates(df2, lon = "longitude", lat = "latitude",
                                    species = "name",
                                    tests = c("capitals", "centroids","equal","institutions","zeros", "seas"))

# 边检查边安装一些库文件
summary(problem_records) # 会显示有问题的分布点

# 从数据集中排除所有有问题的记录
data_clean <- df2[which(problem_records$.summary== "TRUE"),]
str(data_clean)

# 保存清洗后的数据
write.xlsx(data_clean,'4_spocc_clean.xlsx')

# 批量可视化-------------------------
# 定义作图函数(不单独保存文件，用于拼图,PPT内展示)-----------------------------------------
library(ggplot2)
world <- map_data("world")
plot_world <- function(x) {
  df <- subset(data_clean,name == x)
  ggplot(world)+
    geom_polygon(aes(x=long,y=lat,group=group),fill='white',colour='black')+
    coord_quickmap()+ #为地图设置合适的纵横比
    geom_point(data = df,aes(longitude,latitude),
               shape=16,colour='red',size=0.8)+
    ggtitle(x)+
    theme(plot.title = element_text(hjust = .5),
          panel.grid = element_blank(),
          axis.title = element_blank())
}
# 利用purrr的map函数循环出图，利用cowplot::plot_grid()函数排列图
sp <- as.vector(unlist(unique(data_clean$name)))
world_plots <- purrr::map(sp, plot_world) #批量绘图,默认的是7*7的尺寸
# 下面这行是用来组合图片的
p <- cowplot::plot_grid(plotlist = world_plots,labels = LETTERS[1:length(sp)])

# 保存组合的图片
pdf("批量分布绘图.pdf",width=15,height=6)
p
dev.off()

