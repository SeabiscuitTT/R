# 20220321
# 此代码结合最新的SDM流程进行修改
# 调整了稀疏化和拆分测试集和训练集的顺序
# 结合参考文献和视频优化流程 修改了交叉验证的内容
# 【这个版本适合阈值一般的物种建模】
# 值得注意的是，如果选用bolckCV来挑选测试集和训练集的时候可能会出现
# 测试集为空的情况
# 故，批量操作的时候还是选择使用随机的方法产生训练集和测试集。



# 分布数据和气候数据已完成
# 算法默认使用Maxent

rm(list = ls()) #清除环境变量
gc()

# 0 所有样本分布点的情况-------------------------------------------

#导入数据(全部)
input <- read.xlsx('data/!Sphingidae_HMD_clean_3597_置换过物种名.xlsx')[,2:3]
China <- st_read("mapdata/China.json") #导入中国地图json文件

df1 <- st_as_sf(input,coords = c("lon", "lat"),crs = 4326)#转换格式

name_p1 <- paste0('EMN/',species_name,'/',species_name,'_all_points.pdf')
pdf(name_p1,width =12 ,height = 10)

# 先画全国地图主体部分把坐标点打上
ChinaMap <- ggplot() +
  geom_sf(data = China, colour = "#525252",fill="white")+ #画地图
  geom_sf(data = df1,shape=20,colour='red',size=1)+ #注意这边打点不是geom_point()
  coord_sf(ylim = c(-2387082,1654989),crs="+proj=laea +lat_0=40 +lon_0=104")+ #设置坐标系
  annotation_scale(location = "bl") + #比例尺(不适用于默认坐标系)
  annotation_north_arrow(location = "tl", #指北针
                         height = unit(1, "cm"),width = unit(1, "cm"),
                         which_north = "false", style = north_arrow_fancy_orienteering)+
  ggtitle('points in Beibu Gulf')+
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = .5))

# 画九段线图
SouthChinaSea <- ggplot() +
  geom_sf(data = China, colour = "#525252",fill="white")+
  coord_sf(xlim = c(117131.4,2115095),ylim = c(-4028017,-1877844),
           crs="+proj=laea +lat_0=40 +lon_0=104")+
  theme(aspect.ratio = 1.25,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="#525252",fill = NA),
        plot.margin=unit(c(0,0,0,0),"mm"))
# 把大地图和小地图拼起来
ggdraw() +
  draw_plot(ChinaMap) +
  draw_plot(SouthChinaSea, 
            x = 0.86, y = 0.00, #控制小图的位置
            width = 0.1, height = 0.3) #调整小图的长宽比例

dev.off()

# !循环开始-------------------------------------------------

library(dismo)
library(ENMeval)
library(blockCV)
library(spocc)
library(spThin)
library(raster)
library(sp)
library(sf)
library(corrplot)
library(openxlsx)
library(dplyr)
library(readr) 
library(ggplot2)
library(patchwork)
library(beepr)

# 有效样本数量在5个以上的
for(i in 1:4){
  
  # 1 导入相关数据-------------------------------------------
  load('clim2_cors2_cors3_bg.Rdata') #这个要放到循环里面
  # clim2 变量保存了该流程需要的所有环境因子图层
  
  #提取接下来将要批量操作的物种分布数据所在路径
  fiels_list <- list.files('data/more_than5/') 
  # head(fiels_list) #检查格式
  
  # i=31 #用来设置具体样本
  species_name <- gsub('.csv','',fiels_list[i]) 
  file_name <- paste0('data/more_than5/',species_name,'.csv')
  occs <- read.csv(file_name,sep = " ")
  # 挑选出经纬度数据
  occs <- subset(occs, select = c("lon", "lat"))
  
  # 2 分布记录的稀疏化-----------------------------------------------
  # （保证每个栅格当中只有一条记录）
  # 必须要用x,y命名，以保证和backgrounds的数据框名称相符
  colnames(occs) <- c("x", "y") 
  
  # 定义稀疏化的函数
  thin <- function(x) {
    cells <- cellFromXY(clim2[[1]], x)
    dups <- duplicated(cells)
    res <- x[!dups,]
    cat("Number of occs before thinning:", nrow(x), "\n")
    cat("Number of occs after thinning:", nrow(res), "\n")
    return(res)
  }
  
  occs_thinned <- thin(occs) 
  #下面所有需要使用全部物种分布数据的时候，都使用occs_thinned
  
  # 由于有效分布点低于5，Maxent无法运行，所以需要判断一下
  # 有效坐标低于5就跳过，直接做下一个物种
  if(nrow(occs_thinned)<5){next}else{
    
    # 设置种子便于重现结果
    set.seed(12345)
    
    # 3 生成背景点------------------------
    # 在环境图层的范围内生成背景点（背景数据不是试图猜测缺失的位置，
    # 而是用来描述研究区域的环境特征，不要和伪不分布点混淆）
    # bg <- dismo::randomPoints(clim2[[1]], n = 10000) %>% as.data.frame()
    # 把bg和之前的clim2_cors2_cors3存到一起，方便调用
    # save(clim2,cors2,cors3,bg,file='clim2_cors2_cors3_bg.Rdata')
    # beep(8)
    
    # 4 划分测试集和训练集----------------
    # 4.1方法1:blockCV--------------
    
    occs_thinned$Species <- 1 #给稀疏化后的分布点添加一列，注明1表示存在
    bg$Species <- 0 #给背景点添加一列，注明0表示不存在
    PB <- rbind(occs_thinned,bg)
    # 把分布点和背景点转换为具有和clim2一直空间坐标系的点
    pb_data <- st_as_sf(PB, coords = c("x", "y"), crs = crs(clim2))
    
    
    # 空间区块使用systematic assignment
    sb <- spatialBlock(speciesData = pb_data, # presence-background data
                        species = "Species",
                        rasterLayer = clim2,
                        rows = 5,
                        cols = 6,
                        k = 4,
                        selection = "systematic",
                        biomod2Format = TRUE)
    
    
    # 画出背景点图
    # sb$plots + geom_sf(data = pb_data, alpha = 0.5,size=0.7)

    
    # 下面3行是出背景点分区的图
    occs.grp <- sb$foldID[1:nrow(occs_thinned)]
    bg.grp <- sb$foldID[(nrow(occs_thinned)+1):length(sb$foldID)]
    # evalplot.grps(pts = bg, pts.grp = bg.grp, envs = clim2)
    
    
    # 所属不同背景点区域的分布点（这边可以看出分布点和背景点的对应关系）
    # evalplot.grps(pts = occs_thinned, pts.grp = occs.grp, envs = clim2) + 
    #   ggplot2::ggtitle("Spatial block partitions: occurrences")
    
    # # 拼图
    # p1 <- sb$plots #分块图
    # p2 <- sb$plots + geom_sf(data = pb_data, alpha = 0.5,size=0.7) #背景点图
    # p3 <- evalplot.grps(pts = bg, pts.grp = bg.grp, envs = clim2) +
    #   ggplot2::ggtitle("Spatial block partitions: background points") #背景点分区块图
    # p4 <- evalplot.grps(pts = occs_thinned, pts.grp = occs.grp, envs = clim2) + 
    #   ggplot2::ggtitle("Spatial block partitions: occurrences") #分布点图
    # 
    #  # 保存拼图
    # pdf('blockCV.pdf',width =12 ,height = 8)
    # (p1|p3)/(p2|p4) + plot_annotation(tag_levels = 'A') #摆放图的顺序
    # dev.off()
    
    # 现在要从1-4折中选择一份出来作为测试集，
    # 结合所有样本分布点达到地图上，发现1折的分布点比较密集，
    # 不适合留作测试集，因为如果留作测试集就会牺牲训练的精度
    # 最终在2-4里面选了2，因为海南岛一半是1一半是2
    
    # 结合blockCV分区挑选合适的训练集和测试集
    occs.plot <- cbind(occs_thinned, partition = factor(occs.grp)) #得到分布点所在分区
    bg.plot <- cbind(bg, partition = factor(bg.grp)) #得到背景点所在分区
    all.plot <- rbind(occs.plot,bg.plot)
    
    # 挑选出测试集，然后拆分为分布和背景
    test.ste <- filter(all.plot,partition==2) #后面留着测试模型选阈值
    test.ste.occs <- filter(test.ste,Species==1) #测试集里面的分布点
    test.ste.bg <- filter(test.ste,Species==0) #测试集里面的背景点
    
    # 挑选出训练集，然后拆分为分布和背景
    train.set <- filter(all.plot,partition!=2) 
    train.ste.occs <- filter(train.set,Species==1) #训练集里面的分布点，放到ENMeval里面
    train.ste.bg <- filter(train.set,Species==0) #训练集里面的背景点
    
    
    # 无用的中间变量这时候可以删除,节省内存
    rm(occs)
    rm(cors2)
    rm(cors3)
    gc()
    
    # 查看R运行内存
    memory.limit()
    # 由于下面模型预测的矢量结果会很大，所以得调高R的运行内存
    memory.limit(25000)
    
    # 4 用ENMeval包驱动Maxent筛选模型--------------
    # ENMevaluate()的使用参考：
    # https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html
    
    eval_occs_all <- ENMevaluate(
      occs = train.ste.occs[,1:2], #我觉得应该使用稀疏化后所有的分布数据
      envs = clim2,  #环境数据
      bg = bg[,1:2],       #背景点
      partitions = 'checkerboard2', #分区块的模型选择
      tune.args = list(fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), #这边是该参数所有的类型,一般会全测
                       rm = c(0.1, seq(0.5, 6, 0.5))), 
      algorithm = 'maxent.jar',
      parallel = TRUE  #并行处理
    )
    beep(8)
     
    #创建该物种的文件夹，存放该物种所有数据
    dir.create(paste0("EMN/",species_name)) 
    #保存模型结果，便于后期调整阈值
    file_name <- paste0("EMN/",species_name,"/",species_name,".RData")
    save(eval_occs_all,file = file_name) #保存
    
    
    # 5 Delta AICc最小的模型---------------------
    # 关于模型选择有几点注意事项：
    # AUC值作为选择最优模型一直被诟病，原因在于模型选择的研究区域越大，AUC值就越大
    # Lobo et al.(2008)等人指出，验证AUC作为存在背景enm的绝对性能度量是不合适的，
    # 但它对于使用相同数据构建的模型的相对比较是有效的。
    # 这边选用delta.AICc值来选择最优模型
    # 在实践中，delta AICc分数小于2的模型通常被认为在统计上是等价的。
    
    # 绘制AICc曲线
    delta_AICc <- evalplot.stats(e = eval_occs_all,
                                 stats = c("delta.AICc"),
                                 color = "fc",
                                 x.var = "rm",
                                 error.bars = FALSE)
    
    name_p1 <- paste0('EMN/',species_name,'/',species_name,'_delta_AICc.pdf')
    
    pdf(name_p1,width = 5,height = 4)
    print(delta_AICc) #循环中的作图一定要有print()，不然生成的PDF文件会打不开
    dev.off()
    
    # 如果有多个相同的AICc，则只选取第一个，最好是AICc=0
    #第一个delta.AICc=0的模型
    bestmodel_index <- which(eval_occs_all@results$delta.AICc == 0)[1]
    #所有delta.AICc
    bestmodel_index2 <- which(eval_occs_all@results$delta.AICc == 0)
    # 查看所有delta.AICc=0的Model，并保存
    Model <- eval_occs_all@results[bestmodel_index2, ] 
    name_t1 <- paste0('EMN/',species_name,'/',species_name,'_model.csv')
    write.csv(Model,file = name_t1)
    
    # 6 反应曲线可视化----------
    # 方便了解每个物种适宜生存的范围
    
    res <- eval.results(eval_occs_all)
    opt.aicc <- res %>% filter(delta.AICc == 0)
    mod.seq <- eval.models(eval_occs_all)[[opt.aicc[1,]$tune.args]]
    
    # 显示每个因子的贡献率，和每个因子的重要性比较了一下，是一样的
    name <- paste0('EMN/',species_name,'/',species_name,'_contribution.pdf')
    pdf(name,width = 5,height = 4)
    plot(mod.seq, type = "cloglog",
         main=paste0(species_name,' variable contribution')) 
    dev.off()
    
    # 显示反应曲线
    # 没有贡献率的环境因子就是一条水平线
    name <- paste0('EMN/',species_name,'/',species_name,'_response.pdf')
    pdf(name,width = 8,height = 4)
    dismo::response(eval.models(eval_occs_all)[[opt.aicc[1,]$tune.args]])
    dev.off()
    
    
    # 7 模型的评估-----------
    ### Get the threshold of the binary output 获取二进制输出的阈值
    
    e <- evaluate(test.ste.occs[,1:2], test.ste.bg[,1:2], 
                  eval_occs_all@models[[bestmodel_index]], 
                  clim2)
    
    
    # 8 绘制适宜分布区图-----------------------------
    ### converet to binary output based on spec_sens，这边也可以是其他的阈值
    
    tr_all <- threshold(e) #不指定具体的阈值，会显示所有的阈值
    
    # 会出来几个指标，结合自己的生物学问题，选择合适的参数作为下一步显示的标准
    #            kappa     spec_sens no_omission prevalence equal_sens_spec sensitivity
    # thresholds 0.3156031 0.3156031 0.3156035   0.00096493 0.3156031       0.3367208
    
    tr_spec_sens <- threshold(e, 'spec_sens') #指定'spec_sens'
    
    name_t2 <- paste0('EMN/',species_name,'/',species_name,'_thresholds.csv')
    write.table(tr_all,file = name_t2)
    
    # 显示模型的预测结果(初步看看)
    eval_occs_all_predict <- eval.predictions(eval_occs_all)[[opt.aicc[1,]$tune.args]]
    raster::plot(eval_occs_all_predict,
         legend = FALSE, 
         main = paste0(species_name,' prediction'))
    # Error in plot.new() : figure margins too large 奇怪报错
    # 试着把绘图区域拉大一点，然后清空一下绘图区，就OK了
    
    # #选出阈值中最高的阈值作为后续判断是否分布的标准
    # bigest <- as.matrix(which.max(occs_all_th1[1,])) 
    # occs_all_binary_bigest <-
    #   eval_occs_all_cloglog >= occs_all_th1[,rownames(bigest)]
    # # 阈值的选择有时候会出现当选择最大的阈值时，栅格图层里面没有大于这个阈值的值，
    # # 不知道为什么会出现这种情况，
    # # 这边设置一个条件：如果出现occs_all_binary_bigest=0，那么就选用第二大的阈值，以此类推。
    # # 上面这个变量很重要，后续可以 raster2comm
    
    # 选出大于阈值的栅格
    occs_binary <- eval_occs_all_predict >= tr_spec_sens
    
    
    # 得先保存为.tif文件
    name_p2 <- paste0('EMN/',species_name,'/',species_name,'.tif')
    writeRaster(occs_binary,name_p2) #输出栅格备用
    # 画出来的区域就是可能分布的区域
    
    # 8 适生区面积----------------------------------
    
    # 根据栅格大小，粗略估计分布区面积（单位为平方千米）
    get_area_presence <- function(x) {
      cell_size = raster::area(x)@data@values
      return(sum(na.omit(cell_size[x@data@values == 1])))
    }
    
    area_occs_binary <- get_area_presence(occs_binary)
    
    # print("The suitable area for Castanopsis fargesii is:")
    # print(paste0(round(area_occs_all_binary_bigest, 2), "km^2"))
    # 5_Hippotion boerhaviae "157986.11km^2"
    area <- paste0(round(area_occs_binary, 0), "km^2")
    name_t3 <- paste0('EMN/',species_name,'/',species_name,'_area.csv')
    write.table(area,name_t3)
    
    # 8 各因子的相对重要性
    aic.opt <-
      eval_occs_all@models[[which(eval_occs_all@results$delta.AICc == 0)[1]]]
    
    # 定义一个函数来计算变量重要性
    # 这里有一点需要注意，这个重要性是通过AIC值去计算的，
    # AIC是越小越好，但是这边重要性的值是越大越直观
    var.importance <- function(x) {
      temp <-
        x@results[grepl('permutation.importance', rownames(x@results)),]
      names(temp) <-
        gsub("\\.permutation\\.importance", "", names(temp))
      return(temp)
    }
    
    df <- var.importance(aic.opt) #数值越大越重要
    # bio13 bio18 bio19  bio2  bio5 
    # 0     0     0      0     0 
    # Bio13 Bio18 Bio19  Bio2  Bio5  Bio7  Bio8 （bg从1万调整为1千之后）
    # 0     0     0      0     0     100   0 
    # 可能是分布数据太少结果很差
    
    name_t4 <- paste0('EMN/',species_name,'/',species_name,'_aic_opt.csv')
    write.csv(df,name_t4)
    
    par(mar = c(8, 4.1, 4.1, 2.1)) #设置图的位置，下左上右
    
    # 可视化比较每个环境因子的AIC值
    name_p3 <- paste0('EMN/',species_name,'/',species_name,'_aic_opt.pdf')
    pdf(name_p3,width = 5,height = 4)
    barplot(
      df,
      names.arg = row.names(df),
      las = 2,
      ylab = "Permutation Importance",
      main = species_name
    )
    dev.off()
    
    # 10 地图绘制-------------------------------------------------------
    
    library(sf)
    library(tmap) #加载比较慢
    library(raster)
    
    world <- read_sf("./mapfiles/world20200121_polygon.shp")
    china <- read_sf("./mapfiles/bou1_4l.shp")
    provinces <- read_sf("./mapfiles/province_polygon.shp")
    
    # 10.1 带梯度的适宜分布区图---------------------------------------
    
    # 基于tmap写的一个函数，展示适宜梯度
    show_map <-
      function(r,
               title = "",
               palette = NULL,
               legend.show = TRUE,
               main_title) {
        
        res <- tm_shape(r) +
          tm_raster(title = title,
                    legend.show = legend.show,
                    palette = palette) +
          tm_shape(provinces) +  #加上省图层
          tm_borders(col = "grey", lwd = 0.4, lty = 2) + #省界线的属性
          tm_shape(china) +      #加上中国图层
          tm_lines(col = "black") + #国界线颜色
          tm_scale_bar(position = c(.78, .045),width=0.15) + # 比例尺的位置
          tm_compass(type = "arrow", position = c("left", "top")) + #指北针的款式和位置
          tm_layout(main.title = main_title,
                    main.title.position = "center", #主标题居中
                    # inner.margins=c(0.03,0.03,0.03,0.03), # 底部、左侧、顶部和右侧空白
                    legend.position = c(.78, .1)) #图例的位置
        return(res)
      }
    
    # 出完整的全国预测图
    tmap_color <-
      show_map(
        eval_occs_all_predict,
        "Suitability",
        palette = "-Spectral",
        main_title = paste0("The suitability of ",species_name),
        legend.show = TRUE)
    
    # tmap_color
    
    # 保存
    name_p4 <- paste0('EMN/',species_name,'/',species_name,'_suitability_color1.pdf')
    pdf(name_p4,width = 6,height = 5)
    print(tmap_color)
    dev.off()
    
    
    # 10.2 两色适宜分布区图 -------------------------------------
    # 基于tmap写的一个函数，（绿色代表适宜，橙色代表不适宜）
    
    show_map2 <-
      function(r,
               title = "",
               palette = NULL,
               legend.show = TRUE,
               main_title) {
        r0 <- r
        r[r0 > 0] <- "Suitable"
        r[!r0 > 0] <- "Not suitable"
        
        res <- tm_shape(r) +
          tm_raster(title = title,
                    legend.show = legend.show,
                    palette = palette) +
          tm_shape(provinces) +
          tm_borders(col = "grey", lwd = 0.6, lty = 4) +
          tm_shape(china) +
          tm_lines(col = "black") + #国界线颜色
          tm_scale_bar(position = c(.78, .045),width=0.15) + # 比例尺的位置
          tm_compass(type = "arrow", position = c("left", "top")) + #指北针的款式和位置
          tm_layout(main.title = main_title,
                    main.title.position = "center", #主标题居中
                    # inner.margins=c(0.03,0.03,0.03,0.03), # 底部、左侧、顶部和右侧空白
                    legend.position = c(.78, .1)) #图例的位置
        return(res)
      }
    
    tmap2 <-
      show_map2(
        occs_binary,
        palette = "-Pastel2",
        main_title = paste0("The suitability of ",species_name),
        legend.show = TRUE
      )
    
    # tmap2
    
    # 保存
    name_p5 <- paste0('EMN/',species_name,'/',species_name,'_suitability_color2.pdf')
    pdf(name_p5,width = 6,height = 5)
    print(tmap2)
    dev.off()
    
    
    rm(list = ls()) # 清除这一步的环境变量
    gc()
    
  }
  
}
# 及时检查完成的物种出图的效果
# 如果有出图比较奇怪的，要及时调整代码
beep(8)
