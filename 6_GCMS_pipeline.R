
# GC MS 的水稻根系分泌物#---------

library(ggClusterNet)
ps = readRDS("./data/MS.data/b561GCMS/ps.GCMS.rds")

ps1 = ps %>% subset_samples.wt("group","WT")

map = ps1 %>% sample_data()

map
map$Group = map$treated
sample_data(ps) = map



# GC MS 的水稻根系分泌物#---------

library(ggClusterNet)
ps = readRDS("./data/MS.data/b561LCMS/ps_GC.rds")

ps1 = ps %>% subset_samples.wt("group","WT")

map = ps1 %>% sample_data()

map
map$Group = map$treated
sample_data(ps) = map





#--代谢组学数据
#-非靶向代谢组学数据类似环境数据
#但是由于代谢物数量较多，所以需要单独分析

library(tidyverse)
library(phyloseq)

#---0 构建phyloseq对象#------

dat = readxl::read_excel("./data/GCMSdata/GCMS.xlsx",sheet = 1) %>% as.data.frame()
head(dat)
dat = dat %>% distinct(ID, .keep_all = TRUE)  
row.names(dat) = dat$ID
dat$ID = NULL
dat = as.matrix(dat)
dat[is.na(dat)] = 0



map = read.csv("data/GCMSdata/map.csv")
head(map)
row.names(map) = map$ID

ps = phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(dat),taxa_are_rows = TRUE),
  phyloseq::sample_data(map)
  
)

saveRDS(ps,"./data/ps_LC.rds")
#---开始分析#--------

ps = readRDS("./data/ps_GC.rds")
ps = readRDS("./data/ps_LC.rds")


#---代谢物过滤：去除QC样本和非代谢物样本
# ps <- subset_samples(ps,Group %in% c("QC));ps


#---主题颜色设置#-------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//total_amplicon.R")
#---扩增子环境布置
ps0 = ps
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

#--提取有多少个分组
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
# 设定排序顺序
axis_order = phyloseq::sample_data(ps)$Group %>% unique()


#-设定结果保存路径
repath = "./result_and_plot/GCMS_result_and_plot2/"
fs::dir_create(repath)


#--这部分代码正在调试--尚不可使用-这部分问题先不解决
#-本来设计的使用爬虫实时调用，但是网络问题很多
# -所以就使用下载数据库来进行，现在正在权衡中

#--1 代谢物注释HMDB和KEGG数据库#------
id = ps %>% ggClusterNet::vegan_otu() %>% t() %>%
  as.data.frame() %>% row.names()
#-HMDB数据库注释
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//ann.HMDB.R")
repath = "E:/Shared_Folder/Function_local/R_function/micro/"
tax1 = ann.HMDB (id = id,repath  = repath )
colnames(tax1)

tax1 = tax1 %>% distinct(id,.keep_all = TRUE) %>%
  column_to_rownames("id")
head(tax1)
# #-2 代谢物注释KEGG数据库-模糊注释#----
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro//ann.kegg.compounds.R")
# tax2 = ann.kegg(id,repath = "E:/Shared_Folder/Function_local/R_function/micro/")
# head(tax2)
# #--注释kegg数据库算法2
# tax2 = ann.kegg2(id,repath = "E:/Shared_Folder/Function_local/R_function/micro/")
# head(tax2)
# tax= cbind(tax1,tax2)

# tax0 = ps %>% vegan_tax() %>% as.data.frame() #%>% rownames_to_column("ID")
# head(tax0)
# tax2 = tax0 %>% left_join(tax,by = "ID") %>% column_to_rownames("ID")
# head(tax2)

# tax =  readxl::read_excel("./data/GCMSdata/GCMS.xlsx",sheet = 2) %>% 
#   as.data.frame()
# head(tax)
# tax = tax %>% distinct(ID, .keep_all = TRUE)  
# row.names(tax) = tax$ID
phyloseq::tax_table(ps) = as.matrix(tax1)

ps
# tax_table(ps)
saveRDS(ps,"./data/GCMSdata/ps_GC_upper.rds")


ps

#--3 分类堆叠柱状图#--------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot_GC.R")
barpath = paste(repath,"/Microbial_composition/",sep = "")
dir.create(barpath)

phyloseq::rank_names(ps)
j = "Class"

strbar = c("Super_class","Class" , "Sub_class"  )
# strbar = c("Superclass","Class"  )

for (j in strbar) {
  result = barMainplot(ps = ps,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 12)
  p4_1 <- result[[1]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset3) +
    # scale_x_discrete(limits = axis_order) +
    mytheme1
  p4_1
  p4_2  <- result[[3]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset3) +
    # scale_x_discrete(limits = axis_order) + 
    mytheme1
  p4_2
  
  databar <- result[[2]] %>% 
    dplyr::group_by(Group,aa) %>%
    dplyr::summarise(sum(Abundance)) %>% as.data.frame()
  head(databar)
  colnames(databar) = c("Group",j,"Abundance(%)")
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
  ggsave(FileName1, p4_2, width = (5+ gnum), height =8 )
  FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
  ggsave(FileName2, p4_2, width = (5+ gnum), height =8 )
  
  FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
  ggsave(FileName1, p4_1, width = (5+ gnum), height =8 )
  FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
  ggsave(FileName2, p4_1, width = (5+ gnum), height =8 )
  
  FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
  write.csv(databar,FileName)
}



#--5-分类化合物分组差异#-------
library(EasyStat)
library(ggClusterNet)
barpath = paste(repath,"/Different_Class_EasyStat/",sep = "")
dir.create(barpath)

map = sample_data(ps)
head(map)
map = map[,1:2]
sample_data(ps) = map
for (j in strbar) {
  dat <- ps %>% scale_micro(method = "rela") %>%
    tax_glom_wt(ranks = j) %>%
    vegan_otu() %>% 
    as.data.frame()
  head(dat)
  
  dat$id = row.names(dat)
  
  dat2 = dat %>% 
    dplyr::left_join(as.tibble(sample_data(ps)),by = c("id" = "ID")) %>%
    # dplyr::filter(Group != "qiao") %>%
    dplyr::rename(group = Group) %>%
    select(id,group,everything())
  # dat2 %>%
  #   dim()
  
  dat2$group = as.factor(dat2$group)
  head(dat2)
  
  result = MuiKwWlx2(data = dat2,num = c(3:dim(dat2)[2]))
  
  FileName <- paste(barpath,"/",j,"_classification_different_label.csv", sep = "")
  write.csv(result,FileName,sep = "")
  FileName <- paste(barpath,"/",j,"_classification_data.csv", sep = "")
  write.csv(dat2,FileName,sep = "")
  
  result1 = EasyStat::FacetMuiPlotresultBox(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                                            ncol = 4 )
  p1_1 = result1[[1]] + 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset1)
  p1_1
  
  res = FacetMuiPlotresultBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                              ncol = 4)
  p1_2 = res[[1]]+
    # scale_x_discrete(limits = axis_order) + 
    guides(color = FALSE) +
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2
  
  res = FacetMuiPlotReBoxBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",ncol = 4)
  p1_3 = res[[1]]+ 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3
  
  h = dim(dat2)[2]%/%4
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_box", ".pdf", sep = "")
  ggsave(FileName, p1_1, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_bar", ".pdf", sep = "")
  ggsave(FileName, p1_2, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_boxbar", ".pdf", sep = "")
  ggsave(FileName, p1_3, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_box", ".jpg", sep = "")
  ggsave(FileName, p1_1, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_bar", ".jpg", sep = "")
  ggsave(FileName, p1_2, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_boxbar", ".jpg", sep = "")
  ggsave(FileName, p1_3, width = 12, height =3*h,limitsize = FALSE)
  
  
}


#---4-分组热图#-----------
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\GC_ggheatmap_buplot.R")
rank.names(ps)

barpath = paste(repath,"/Different_Class_heatmap_EasyStat/",sep = "")
dir.create(barpath)



for (j in strbar) {
  ps_rela <- ps %>% scale_micro(method = "rela") %>%
    tax_glom_wt(ranks = "Class")
  
  result <- GCheatmap (ps_rela,
                       label =  F,
                       col_cluster = F,
                       row_cluster = F)
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  
  h = taxa_names(ps_rela) %>% length()
  w = sample_names(ps_rela)  %>% length()
  
  
  filename = paste(barpath,"/",j,"_classification_","ggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggbubble.pdf",sep = "")
  ggsave(filename,p2,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggheatmap.png",sep = "")
  ggsave(filename,p1,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggbubble.png",sep = "")
  ggsave(filename,p2,width = w/1.7,height = h/3)
  
}



#---6 排序分析PCA等#-------------
betapath = paste(repath,"/beta_orda_total/",sep = "")
dir.create(betapath)


# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")


methodlist = c("t-sne","LDA", "PCA")
for (method in methodlist) {
  result = BetaDiv(ps = ps, group = "Group", dist = "bray",
                   method = method, Micromet = "anosim",
                   pvalue.cutoff = 0.05,pair = F)
  p3_1 = result[[1]] + 
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) +
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_1
  #带标签图形出图
  p3_2 = result[[3]] +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 12)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 12)
  
  # 提取出图数据
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  #---------排序-精修图
  plotdata =result[[2]]
  head(plotdata)
  # 求均值
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  # 合并到样本坐标数据中
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  
  # p2$layers[[2]] = NULL
  # library(ggcor)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 8)
  
}

map
#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_anosim.csv", sep = "")
write.csv(TResult,FileName)

#--换用adonis差异分析

# title1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
# title1
# FileName <- paste(betapath,"Total_adonis.csv", sep = "")
# write.csv(title1,FileName)
# pairResult = pairMicroTest(ps = ps, Micromet = "adonis", dist = "bray")
# FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
# write.csv(pair,FileName)


# 6.2 PLS-DA排序#---------
count = vegan_otu(ps)
map = as.data.frame(sample_data(ps))
map$Group
#PLS-DA分析，这里也是选取2个主成分
plsda.datatm <-plsda(count, map$Group, ncomp = 2)
#PLS-DA without centroids
mi=c("#1B9E77" ,"#D95F02")
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2' )
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2',ellipse = TRUE, 
          size.xlabel = 20, size.ylabel = 20, size.axis = 25, pch = 15, cex = 5)
#----提取数据作图
a = unclass(plsda.datatm)
#--提取坐标值
plotdata = as.data.frame(a$variates$X)
plotdata$SampleType = map$Group
#-提取解释度
eig = a$explained_variance$X
eig[1]
library(ggalt)
library(BiocManager)
# install("ggalt")
p = ggplot(data = plotdata,aes(x=comp1,y=comp2,group=SampleType,color=SampleType))+geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  geom_encircle(s_shape=1, expand=0) +
  labs(x=paste("X-variate 1 (", format(100 * eig[1]), "%)", sep=""),
       y=paste("X-variate 2 (", format(100 * eig[2] ), "%)", sep=""))+
  labs(title = "PLS-DA") 
p
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = mi)+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p

method = "pls-da"
FileName <- paste(betapath,"/a2_",method,"_plot.pdf", sep = "")
ggsave(FileName, p3_3, width = 8, height = 8)
FileName1 <- paste(betapath,"/a2_",method,"_plot.jpg", sep = "")
ggsave(FileName1 , p3_3, width = 8, height = 8)


#---7 层次聚类#--------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/cluster_plot.R")
clupath = paste(repath,"/cluster_plot/",sep = "")
dir.create(clupath)
res = cluster_plot (ps= ps,hcluter_method = "complete",
                    dist = "bray",cuttree = gnum,row_cluster = T,col_cluster =  T)

p0 = res[[1]]
p0



FileName <- paste(clupath,"cluster", ".jpg", sep = "")
ggsave(FileName, p0, width = 6, height =8,limitsize = FALSE)
FileName <- paste(clupath,"cluster", ".pdf", sep = "")
ggsave(FileName, p0, width = 6 , height = 8,limitsize = FALSE)

p1 = res[[2]]
p2 = res[[3]]

FileName <- paste(clupath,"heatmap_cluster", ".jpg", sep = "")
ggsave(FileName, p1, width = 8, height =8,limitsize = FALSE)
FileName <- paste(clupath,"heatap_cluster", ".pdf", sep = "")
ggsave(FileName, p1, width = 8 , height = 8,limitsize = FALSE)

FileName <- paste(clupath,"bubble_cluster", ".jpg", sep = "")
ggsave(FileName, p2, width = 8, height =8,limitsize = FALSE)
FileName <- paste(clupath,"bubble_cluster", ".pdf", sep = "")
ggsave(FileName, p2, width = 8 , height = 8,limitsize = FALSE)

dat = res[4]
FileName <- paste(clupath,"clu_data.csv", sep = "")
write.csv(dat,FileName)


#-----差异代谢物#----------
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
alppath = paste(repath,"/All_different_metabolites/",sep = "")
dir.create(alppath)

#--非参数检验
result = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "wilcox")
head(result)
FileName <- paste(alppath,"/data_wlx_all.compounds.csv", sep = "")
write_csv(result,FileName)
#--t检验检验--建议四个重复以上
result = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "ttext")
head(result)
FileName <- paste(alppath,"/data_ttest_all.compounds.csv", sep = "")
write_excel_csv(result,FileName)


#---8 单变量统计分析-箱线图等可视化#--------

#--提取差异代谢物标签
# head(result)

alppath = paste(repath,"/summary_stat_plot/",sep = "")
dir.create(alppath)

dat = ps %>% 
  ggClusterNet::vegan_otu() %>% 
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"
num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25 
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {
  
 if (j == 1) {
   A[[j]] = num[1:25]
 } else if(j != n.fac2){
   x = (25*(j - 1) + 1)
   y = 25*j
   A[[j]] = num[x:y]
 }else if (j == n.fac2){
   x = (25*(j - 1) + 1)
   y = 25*j
   A[[j]] = num[x:length(num)]
   
 }

}



for (i in 1:n.fac2) {
  result = EasyStat::MuiaovMcomper2(data = data,num = A[[i]])
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = A[[i]],
                                            result = result,
                                            sig_show ="abc",ncol = 5 )
  p1_1 = result1[[1]] + 
    ggplot2::scale_x_discrete(limits = axis_order) + 
    mytheme2 +
    ggplot2::guides(fill = guide_legend(title = NULL)) +
    ggplot2::scale_fill_manual(values = colset1)
  p1_1
  
  res = EasyStat::FacetMuiPlotresultBar(data = data,num = A[[i]],
                                        result = result,sig_show ="abc",ncol = 5)
  p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
    mytheme2+ 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2
  
  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = A[[i]],
                                       result = result,sig_show ="abc",ncol = 5)
  p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_box", ".pdf", sep = "")
  ggsave(FileName, p1_1, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_bar", ".pdf", sep = "")
  ggsave(FileName, p1_2, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_boxbar", ".pdf", sep = "")
  ggsave(FileName, p1_3, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_box", ".jpg", sep = "")
  ggsave(FileName, p1_1, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_bar", ".jpg", sep = "")
  ggsave(FileName, p1_2, width = 18, height =16,limitsize = FALSE)
  
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_boxbar", ".jpg", sep = "")
  ggsave(FileName, p1_3, width = 18, height =16,limitsize = FALSE)
  # result
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"_data_aov_abc.csv", sep = "")
  write.csv(result,FileName,sep = "")
  FileName <- paste(alppath,paste("part_",i,sep = ""),"_data.csv", sep = "")
  write.csv(data,FileName,sep = "")
  
  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    result = result,
    sample = TRUE,
    scale = TRUE
  )
  p1 = res[[1]]
  p1
  h = sample_names(ps) %>% length()
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatmap", ".jpg", sep = "")
  ggsave(FileName, p1, width = 12, height =h/2,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatap", ".pdf", sep = "")
  ggsave(FileName, p1, width = 12, height =h/2,limitsize = FALSE)
  
  p2 = res[[2]]
  p2
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble", ".jpg", sep = "")
  ggsave(FileName, p2, width = 12, height =h/2,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble", ".pdf", sep = "")
  ggsave(FileName, p2, width = 12 , height = h/2,limitsize = FALSE)
  
  
  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    result = result,
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    sample = FALSE,
    scale = TRUE
    
    
  )
  
  p1 = res[[1]]
  p1
  
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatmap_group", ".jpg", sep = "")
  ggsave(FileName, p1, width = gnum*1.5, height =8,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatap_group", ".pdf", sep = "")
  ggsave(FileName, p1, width = gnum*1.5, height =8,limitsize = FALSE)
  
  
  p2 = res[[2]]
  p2
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble_group", ".jpg", sep = "")
  ggsave(FileName, p2, width = gnum*1.5, height =6,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble_group", ".pdf", sep = "")
  ggsave(FileName, p2, width = gnum*1.5, height =6,limitsize = FALSE)
  
  
  res = EasyStat::value_stackBar(
    data = data,
    i =A[[i]],
    result = result,
    add_abc = TRUE)
  
  
  p1 = res[[1]]
  p1
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"sample_relative_abundacne", ".jpg", sep = "")
  ggsave(FileName, p1, width = 6, height =5,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"sample_relative_abundacne", ".pdf", sep = "")
  ggsave(FileName, p1, width = 6, height = 5,limitsize = FALSE)
  
}

#---9 载荷矩阵挑选重要代谢物#------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/loadingPCA.R")
pcapath = paste(repath,"/loadingPCA/",sep = "")
dir.create(pcapath)
res = loadingPCA(ps = ps,Top = 20)

otu = ps %>% vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)
p = res[[1]]
p
dat = res[[2]]

filemane = paste(pcapath,"/PCALoading.pdf",sep = "")
ggsave(filemane, p, width = 8, height = 6)
filemane = paste(pcapath,"/PCALoading.jpg",sep = "")
ggsave(filemane, p, width = 8, height = 6)
FileName <- paste(pcapath,"/Loadsing_pca.csv", sep = "")
write.csv(dat,FileName,sep = "")


#--10 机器学习#-------
matpath = paste(repath,"/Machine_learing/",sep = "")
dir.create(matpath )
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")
ROC  = F
rfcv = F
# library(randomForest)
# library(caret)
# library(ROCR) ##用于计算ROC
# library(e1071)

if (ROC ) {
  #--三种机器学习方法评测
  result = MicroRoc( ps = ps,group  = "Group")
  #--提取roc曲线
  p <- result[[1]] + 
    mytheme1
  p
  #提取AUC值
  data <- result[[2]]
  
  filename = paste(matpath,"/three_method_AUCvalue.csv",sep = "")
  write.csv(data,filename,quote = F)
  
  data <- result[[3]]
  filename = paste(matpath,"/three_method_AUCdata.csv",sep = "")
  write.csv(data,filename,quote = F)
  
  filename = paste(matpath,"/three_method_AUC_plot.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 8)
  filename = paste(matpath,"/three_method_AUC_plot.jpg",sep = "")
  ggsave(filename,p,width = 8,height = 8)
  
}


mapping = as.data.frame(phyloseq::sample_data(ps))

#--随机森林全套-如果圈图尚未显示前面几个，就设定max大一点
result = MicroRF_GC(ps = ps,group  = "Group",optimal = 40,
                 rfcv =F,nrfcvnum = 5,
                 min = -1,max = 5)
#火柴图展示前二十个重要的OTU
p <- result[[1]] + 
  mytheme1
p

filename = paste(matpath,"/randonforest_loading.pdf",sep = "")
ggsave(filename,p,width = 8,height = optimal/3)
filename = paste(matpath,"/randonforest_loading.jpg",sep = "")
ggsave(filename,p,width = 8,height = optimal/3)
# 圈图展示
p <- result[[2]]
p
filename = paste(matpath,"/randonforest_loading_circle.pdf",sep = "")
ggsave(filename,p,width = 8,height = 10)
filename = paste(matpath,"/randonforest_loading_circle.jpg",sep = "")
ggsave(filename,p,width = 8,height = 10)

p <- result[[6]]
p
filename = paste(matpath,"/Show_model.pdf",sep = "")
ggsave(filename,p,width = 8,height = 4)
filename = paste(matpath,"/Show_model.jpg",sep = "")
ggsave(filename,p,width = 8,height = 4)


if (rfcv) {
  # 展示交叉验证结果
  p <- result[[3]]
  filename = paste(matpath,"/randonforest_cross_check.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 12)
  data <- result[[4]]
  filename = paste(matpath,"/randomforest_cross_data.csv",sep = "")
  write.csv(data,filename,quote = F)
}

data <- result[[5]]
filename = paste(matpath,"/randomforest_data.csv",sep = "")
write.csv(data,filename,quote = F)




#----11代谢网络#-------
netpath = paste(repath,"/network_metabolites/",sep = "")
dir.create(netpath)
library(ggClusterNet)
library(igraph)
library(sna)

source("E:\\Shared_Folder\\Function_local\\R_function/micro/lizi_network.R")
result = network(ps = ps,
                 N = 250,
                 big = TRUE,
                 select_layout = F,
                 layout_net = "model_maptree",
                 r.threshold=0.9,
                 p.threshold=0.05,
                 label = FALSE,
                 ncol = gnum,
                 path = netpath,
                 fill = "#8DD3C7",
                 zipi = FALSE)
# 全部样本的网络比对
p = result[[1]]
p
# 全部样本网络参数比对
data = result[[2]]

plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 16*gnum,height = 16,limitsize = FALSE)

tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)


#-12通路富集分析#-------
library(MetaboAnalystR)

enrichpath = paste(repath,"/enrich/",sep = "")
dir.create(enrichpath)
tempath0 = getwd()
setwd(enrichpath)
# PID of current job: 55420
mSet<-InitDataObjects("conc", "pathora", FALSE)
tax = ps %>% ggClusterNet::vegan_tax() %>%as.data.frame()
head(tax)
cmpd.vec <- tax$KEGG[!is.na(tax$KEGG)]

mSet <- Setup.MapData(mSet, cmpd.vec);
mSet <- CrossReferencing(mSet, "kegg");
mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, "ath", "current")
mSet <- SetMetabolomeFilter(mSet, F);
# mSet$api
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")
mSet <- PlotPathSummary(mSet, F, "path_view_1_", "pdf", 72, width=NA)
mSet<-SaveTransformedData(mSet)
mSet<-PlotPathSummary(mSet, T, "path_view_1_", "png", 72, width=NA, NA, NA )

#--通路匹配的物质
tab = mSet$dataSet$map.table %>% as.data.frame()
head(tab)
# 富集的通路表格
tabpath = mSet$analSet$ora.mat %>% as.data.frame()


# for (i in 1:nrow(tabpath)) {
#   library(KEGGREST)
#   # listDatabases()
#   query <- keggGet(c(row.names(tabpath)[i]))
#   id = query[[1]]$NAME
#   com = query[[1]]$COMPOUND
#   datfil = tab %>% filter(KEGG  %in% names(com))
#   tempath = paste("./",id,sep = "")
#   fs::dir_create(tempath)
#   write.csv(datfil,paste(tempath,"/conpound_reaction.csv",sep = ""),quote = FALSE)
# }


# 富集的通路表格
dt = mSet$analSet$ora.mat %>% as.data.frame()
head(dt)

tem = mSet$analSet$ora.hits
A = c()
B = c()
C = c()
for (i in 1:length(names(tem))) {
  tem2 = tem[[names(tem)[i]]] %>% length()
  if (tem2 != 0) {
    A[i] = names(tem)[i]
    B[i] = tem[[names(tem)[i]]] %>% names() %>% str_c( collapse = "|")
    C[i] = tem[[names(tem)[i]]]  %>% str_c( collapse = "|")
  }
}

tem3 = data.frame(ID = A,name = B,kegg = C) %>% filter(!is.na(kegg))
head(tem3)
write_csv(tem3,paste0("./pathwat_contain.conpounds.csv"))


p1 <- ggplot(dt, aes(x =Impact, y = `-log(p)`))  +
  geom_point(pch = 21,aes(fill=`-log(p)`,size =`-log(p)`))  + 
  ggrepel::geom_text_repel(aes(x = Impact, y = `-log(p)`,label= row.names(dt))) +
  scale_fill_gradientn(colours =colorRampPalette(c("#F7F4F9","#FFFF33","#EE0000FF","#EE0000FF","#EE0000FF","#EE0000FF"))(60)) +
  theme_classic()
p1

ggsave(paste0("./compound.to.kegg.koplot.pdf"),p1,width = 8,height = 7)
ggsave(paste0("./compound.to.kegg.ko.png"),p1,width = 8,height = 7)


setwd(tempath0)


#---本人文章中一套分析#---------

library(tidyverse)
library(ggbiplot)
library(readxl)
library(phyloseq)
library(ggClusterNet)
# BiocManager::install("plsda")
# devtools::install_github("Dampeel/plsda")
heatmap_GSVA <- function(heat = sub_heat,
                         map = sub_map,
                         col_cluster =  TRUE,
                         row_cluster =  TRUE,
                         label =  TRUE){
  
  heat = as.data.frame(heat)
  heat$id = row.names(heat)
  head(heat)
  data = heat %>% 
    dplyr::select(id,everything())
  # data[data > 0.3]<-0.3
  mat <- data[,-1] #drop gene column as now in rows
  
  if (col_cluster) {
    clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }
  
  if (label ) {
    
    map$ID = row.names(map)
    labels = ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      scale_fill_brewer(palette = 'Set1',name="Cell Type") +
      theme_void()
  }
  
  pcm = reshape2::melt(data, id = c("id"))
  
  p1 = ggplot(pcm, aes(y = id, x = variable)) + 
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(12,"Spectral"))(60))+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
      
    )
  
  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p2 = ggplot(pcm, aes(y = id, x = variable)) + 
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(12,"Spectral"))(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
      
    )
  
  if (col_cluster) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
  }
  
  if (row_cluster) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02) 
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02) 
  }
  
  if (label) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }
  
  return(list(p1,p2))
  
  
}

repath = "./result_neg_CC/"
dir.create(repath)

ps = readRDS("./ps_GC.neg.rds")


ps = subset_samples(ps,Group %in% c("CC_10","CC_20"));ps
map = sample_data(ps)
# map$Group = map$Group2
# sample_data(ps) = map

count =  ps %>% vegan_otu() %>% t() %>%
  as.data.frame()
map = ps %>% sample_data()


# 进行PCA分析#----------------
pcapath <- paste(repath,"/PCA",sep = "")
dir.create(pcapath)

otu.pca <- prcomp(t(count), scale. = TRUE)
p=ggbiplot(otu.pca, obs.scale = 1, var.scale = 1,
           groups = map$Group, ellipse = TRUE,var.axes = F,circle = TRUE)+
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
p
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p <- p+theme_bw()+scale_colour_manual(values = mi)+labs(x=paste("PC 1 (", 61.7, "%)", sep=""),
                                                        y=paste("PC 2 (", 14.5, "%)", sep=""))+
  labs(title = "PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)

filename = paste(pcapath,"./pca_plot.pdf",sep = "")
ggsave(filename,p,width = 6,height = 5)
# 现在我们来提取作图坐标首先是样品排序轴坐标使用下面这条命令提取
yangpin<-otu.pca$x
yangpin=as.data.frame(yangpin)
yangpin$Group=map$Group
#提取荷载坐标
bianliang<-otu.pca$rotation
bianliang=as.data.frame(bianliang)
#提取特征根,这里提供的并不是特征值而是标准差，需要求其平方才是特征值
eig=otu.pca$sdev
eig=eig*eig
#在这里我设定了随机种子，方便两种形式图形比较
set.seed(10)
p=ggplot(data =yangpin,aes(x=yangpin$PC1,y=yangpin$PC2,group=Group,color=Group))+ 
  geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  labs(x=paste("PC 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PC 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  labs(title = "PCA") 

mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = mi)+labs(x=paste("PC 1 (", 61.7, "%)", sep=""),
                                                     y=paste("PC 2 (", 14.5, "%)", sep=""))+
  labs(title = "PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p

filename = paste(pcapath,"./pca_plot_2.pdf",sep = "")
ggsave(filename,p,width = 6,height = 5)

#转换成矩阵
# library(plsda)
library(mixOmics)
XXt <- as.matrix(count) %>% t()
#PLS-DA分析，这里也是选取2个主成分
plsda.datatm <-plsda(XXt, map$Group, ncomp = 2)

#PLS-DA without centroids
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")

pdf(paste(pcapath,"/a2_PLSDA图.pdf",sep = ""), width = 8, height = 6)
plotIndiv(plsda.datatm , comp = c(1,2),
          group = map$Group, style = 'ggplot2', ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, 
          size.xlabel = 20, size.ylabel = 20, size.axis = 25, pch = 15, cex = 5)

dev.off()

#------差异和热图可视化#-------
diffpath = paste(repath,"/difference",sep = "")
dir.create(diffpath)

library("gplots")
library("RColorBrewer")
library("ggplot2")
library("vegan")
ps

#########下面来分类统计分泌物的信息
##############全部峰值数据来做统计#########去除磷酸之后进行分析
# 读入mapping文件
design = map
# 读取OTU表，这里我选择的是整个otu表格，但是一般没有必要全部做差异的啊，相对丰度高的做做就可以了
otu_table = ps %>% vegan_otu() %>%
  t() %>%
  as.data.frame()
# idx = rownames(design) %in% colnames(otu_table) 
# sub_design = design[idx,]
# count = otu_table[, rownames(sub_design)]
# head(count)
# 转换原始数据为百分比，
norm = t(t(otu_table)/colSums(otu_table,na=T)) * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)

a=norm
head(a)
#预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
Pvalue<-c(rep(0,nrow(a)))
fdr<-c(rep(0,nrow(a)))
log2_FC<-c(rep(0,nrow(a)))
head(a)
a[is.na(a)] <- 0
###########开始运行脚本
for(i in 1:nrow(a)){
  if(sd(a[i,1:3])==0&&sd(a[i,4:6])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(a[i,1:3]),as.numeric(a[i,4:6]))
    Pvalue[i]<-y$p.value
    log2_FC[i]<-log2((mean(as.numeric(a[i,1:3]))+0.001)/(mean(as.numeric(a[i,4:6]))+0.001)) 
    fdr[i]=p.adjust(Pvalue[i], "BH") 
  }
}
# 在原文件后面加入log2FC，p value和FDR,共3列；
out<-cbind(a,log2_FC,Pvalue,fdr) %>% as.data.frame()
head(out)
tax = ps %>% vegan_tax() %>%
  as.data.frame()
out$tax= tax$Name
WT <-subset(out,fdr < 0.05 )

# WT1<- arrange(WT, desc(tax))
# row.names(WT)=row.names(WT)
# head(WT1)

filename = paste(diffpath,"/化合物t检验结果.csv",sep = "")
write.csv(WT,filename)

WT$ID = row.names(WT)

WT1 = WT %>% mutate(ord = log2_FC^2) %>%
  arrange(desc(ord)) %>%
  head(50)
wt=WT1 %>% as.data.frame()
row.names(wt) = wt$ID
wt$ID = NULL
row.names(wt) = paste(wt$ID,wt$tax,sep = "")

wt = wt[,1:6]


library(pheatmap)

#设置颜色梯度"#1B9E77"
color = colorRampPalette(c( "white", "#FFFFE5","#67001F"))(60)
wt2<--log10(wt+0.000001)
wt2<-sqrt(wt)

wt2<-sqrt(wt)
wt2[wt2>0.2]<-0.2


#我们支持分组,这里我随便造一个分组，作为纵向分组
annotation_row = data.frame(tax=WT1$tax)  
rownames(annotation_row) = rownames(wt2)
#我再造一个横向分组
annotation_col = data.frame(design$Group)  
rownames(annotation_col) = colnames(wt)
#install.packages("pheatmap")
#我们开始做分组


p=pheatmap(wt2,
           fontsize=6,
           cellwidth = 8,
           cellheight =6,
           cluster_rows = FALSE,
           cluster_cols = F,
           color = color,
           border_color = "white",
           #cutree_col = 2,
           #gaps_row = c(7,14,21),
           annotation_col = annotation_col,
           # annotation_row = annotation_row,
           labels_row = NULL,labels_col = NULL )

filename = paste(diffpath,"/热图差异.pdf",sep = "")

ggsave(filename, p, width = 15, height = 25,limitsize = FALSE)



result = heatmap_GSVA(heat = wt2,
                      map = sample_data(ps) )

p = result[[1]]
p

filename = paste(diffpath,"/热图差异.pdf",sep = "")
ggsave(filename, p, width = 10, height = 15,limitsize = FALSE)
p = result[[2]]
p

filename = paste(diffpath,"/气泡图差异.pdf",sep = "")
ggsave(filename, p, width = 10, height = 15,limitsize = FALSE)




#-------PCA载荷挑选#---------

count = otu_table

norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
# 进行PCA分析
otu.pca <- prcomp(t(norm), scale. = TRUE)
#现在我们来提取作图坐标#首先是样品排序轴坐标
predict(otu.pca) 
#也可以使用下面这条命令提取
yangpin<-otu.pca$x
yangpin=as.data.frame(yangpin)
yangpin$SampleType=sample_data(ps)$Group
#提取荷载坐标
bianliang<-otu.pca$rotation
bianliang=as.data.frame(bianliang)
head(bianliang)
dim(norm)
index = merge(norm ,bianliang, by="row.names",all=F)
head(index)
row.names(index)=index$Row.names
index$Row.names=NULL

index = merge(index ,tax, by="row.names",all=F)
head(index)
row.names(index)=index$Row.names
index$Row.names=NULL

##手动选择10个最终要的变量 PCA载荷矩阵挑选37个成分提取差异.txt
index$PCone = index$PC1^2
top = index %>% arrange(desc(PCone)) %>% 
  head(20)
top$com = paste(top$Name,sep = "")
#######开始出图，做火柴图########, "#7570B3","#E7298A")
library("ggplot2") 
p=ggplot(top, aes(x = PCone, y = reorder(com,PCone)))  +
  geom_segment(aes(yend=com),xend=0,size=3,colour = "#1B9E77" )+
  geom_point(size=4,pch=20, colour = "#1B9E77")
p 

p=  p+theme_bw()+theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
                       axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
p

filemane = paste(pcapath,"/PCA变量重要性火柴棒图.pdf",sep = "")
ggsave(filemane, p, width = 15, height = 8)

#----------随机森林重要代谢物挑选#-----------

rfpath = paste(repath,"/randomforest/",sep = "")
dir.create(rfpath)
# 加载随机森林包
library(randomForest)
#######使用随机森林做分类
set.seed(315)
t(norm) %>%dim()

iris.rf = randomForest(t(norm), as.factor(design$Group), importance=TRUE, proximity=TRUE,ntree=1000)
print(iris.rf)
a=round(importance(iris.rf), 2)
head(a)
index11 = merge(norm ,a, by="row.names",all=F)
head(index11)
dim(index11)

filename = paste(rfpath,"随机森林分类重要变量输出_norm.txt",sep = "")
write.table(index11,filename,quote = FALSE,row.names = T,
            col.names = T,sep = "\t")


varImpPlot(iris.rf)  

########这里做变量权重图的表格
colnames(index11)[1] = "ID"
tem = index11 %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(40)

# tax$ID = Row.namestax
tem2 <- tem %>% inner_join(tax)
tem2$Com = paste(tem2$Name,sep = "")

mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
library("ggplot2") 
p=ggplot(tem2, aes(x = MeanDecreaseAccuracy, y = reorder(Com,MeanDecreaseAccuracy)))  +
  geom_segment(aes(yend=Com),xend=0,size=3,colour = "#1B9E77" )
#geom_point(size=4,pch=20, colour = "#1B9E77")
p 

p = p+theme_bw()+theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
                       axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
p
filemane = paste(rfpath,"/a7_随机森林40个变量重要性火柴棒图.pdf",sep = "")
# ggsave(filemane, p, width = 15, height = 8)

ggsave(filemane, p, width = 12, height = 8)

#------代谢物通路富集分析#----------------
library(readxl)
library(GO.db)
library(DOSE)
library(GO.db)
library(GSEABase)
library(clusterProfiler)
library("GSVA")

id = read.delim("kegg.neg.txt")
head(id)
colnames(id)
dif = read.csv(paste0(diffpath,"/化合物t检验结果.csv"))
head(dif)
colnames(dif)[1] = "ID"
tem = dif %>% filter( fdr < 0.05, abs(log2_FC) > 0) %>% inner_join(id,by = "ID")
id.tem = tem$Kegg_ID %>% strsplit("[:]") %>%
  sapply( `[`, 2)


library(MetaboAnalystR)

enrichpath = paste(repath,"/enrich",sep = "")
dir.create(enrichpath)
tempath0 = getwd()
setwd(enrichpath)
# PID of current job: 55420
mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<- id.tem

mSet <- Setup.MapData(mSet, cmpd.vec);
mSet <- CrossReferencing(mSet, "kegg");
mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, "ath", "current")
mSet <- SetMetabolomeFilter(mSet, F);
# mSet$api
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")
mSet <- PlotPathSummary(mSet, F, "path_view_1_", "pdf", 72, width=NA)
mSet<-SaveTransformedData(mSet)
mSet<-PlotPathSummary(mSet, T, "path_view_1_", "png", 72, width=NA, NA, NA )

#--通路匹配的物质
tab = mSet$dataSet$map.table %>% as.data.frame()
head(tab)
# 富集的通路表格
tabpath = mSet$analSet$ora.mat %>% as.data.frame()


for (i in 1:nrow(tabpath)) {
  library(KEGGREST)
  # listDatabases()
  query <- keggGet(c(row.names(tabpath)[i]))
  id = query[[1]]$NAME
  com = query[[1]]$COMPOUND
  datfil = tab %>% filter(KEGG  %in% names(com))
  tempath = paste("./",id,sep = "")
  fs::dir_create(tempath)
  write.csv(datfil,paste(tempath,"/conpound_reaction.csv",sep = ""),quote = FALSE)
}


# 富集的通路表格
dt = mSet$analSet$ora.mat %>% as.data.frame()
head(dt)

tem = mSet$analSet$ora.hits
A = c()
B = c()
C = c()
for (i in 1:length(names(tem))) {
  tem2 = tem[[names(tem)[i]]] %>% length()
  if (tem2 != 0) {
    A[i] = names(tem)[i]
    B[i] = tem[[names(tem)[i]]] %>% names() %>% str_c( collapse = "|")
    C[i] = tem[[names(tem)[i]]]  %>% str_c( collapse = "|")
  }
}

tem3 = data.frame(ID = A,name = B,kegg = C) %>% filter(!is.na(kegg))
head(tem3)
write_csv(tem3,paste0("./pathwat_contain.conpounds.csv"))


p1 <- ggplot(dt, aes(x =Impact, y = `-log(p)`))  +
  geom_point(pch = 21,aes(fill=`-log(p)`,size =`-log(p)`))  + 
  ggrepel::geom_text_repel(aes(x = Impact, y = `-log(p)`,label= row.names(dt))) +
  scale_fill_gradientn(colours =colorRampPalette(c("#F7F4F9","#FFFF33","#EE0000FF","#EE0000FF","#EE0000FF","#EE0000FF"))(60)) +
  theme_classic()
p1

ggsave(paste0("./compound.to.kegg.koplot.pdf"),p1,width = 8,height = 7)
ggsave(paste0("./compound.to.kegg.ko.png"),p1,width = 8,height = 7)

setwd(tempath0)





