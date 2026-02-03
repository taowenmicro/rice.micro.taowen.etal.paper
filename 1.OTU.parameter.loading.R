# 参数文件; Tao Wen, 2023.8.28

#--输入和设定文件
# ps0 = base::readRDS("./data/dataNEW/ps.rds")
ps0 = readRDS("./231106/16s_ps.rds")

ps0
# ps0 = base::readRDS("./ps_16s.rds")

library(tidyverse)
library(phyloseq)
library(ggClusterNet)


#  整理map文件：核心是分组信息#---------
map = sample_data(ps0)
head(map)
map0 = readxl::read_xlsx("./original names and new names.xlsx") %>% as.data.frame()

head(map0)
map1 = map %>% 
  as.tibble() %>%
  left_join(map0,by = c("ID"="Original name")) %>% as.data.frame()
head(map1)
map1$X = NULL
row.names(map1) = map1$ID
colnames(map1)[4] = "newid"
sample_data(ps0) = map1

write_delim(map1,"map_deeping.txt")
#  整理结束后进行导入
map = read.delim("./map_deeping.txt")
head(map)
row.names(map) = map$ID
sample_data(ps0) = map
saveRDS(ps0,"ps_deeping.rds")

# 修改样本名称使用一下函数，newid的话就是将新的样本名字放在map里面作为一列，那就行啦；
# ps01 = changeSamplenames(ps = ps0,newid = "newid" )
# sample_data(ps01)
changeSamplenames = function(
    ps = ps0,
    newid = "newid"  #in sample data ,was one of the colname
){
  if (!is.null(ps@otu_table)) {
    otu = ps %>% vegan_otu() %>% t() %>%
      as.data.frame()
  }
  
  if (!is.null(ps@sam_data)) {
    map = ps %>% sample_data() %>% as.tibble() %>%column_to_rownames(newid ) %>%
      as.data.frame()
  }
  
  colnames(otu) = row.names(map)
  if (!is.null(ps@tax_table)) {
    psout = phyloseq(sample_data(map),
                     otu_table(as.matrix(otu), taxa_are_rows=TRUE),
                     tax_table(ps)
    )
  } else{
    psout = phyloseq(sample_data(map),
                     otu_table(as.matrix(otu), taxa_are_rows=TRUE))
  }
  return(psout)
}




## 用于筛选#-----
# ps0 %>% scale_micro() %>%
#   filter_taxa(function(x) sum(x ) > 0.01 , TRUE)
# #低丰度微生物挑选
# ps0 %>% scale_micro() %>%
# filter_taxa(function(x) sum(x ) < 0.0001 , TRUE)
# #--如何筛选样本:去除sample1
# ps_sub <- subset_samples.wt(ps0,"ID",c("sample1"),T);ps_sub

# #--如何筛选微生物
# ps0 <- ps0 %>% subset_taxa.wt("Kingdom", id) ;ps0
# # 是否需要仅仅注释为细菌或者真菌的做进一步分析
# ps0 <- ps0 %>% subset_taxa.wt("Kingdom", id) ;ps0


## 用于修改R包安装路径#-------
# ps0 = base::readRDS("./Error/221121/ps_its.rds")
# ps0
# .libPaths()
# .libPaths(new="C:/Program Files/R/R-4.1.1/library")


# 检查七个等级注释名字
# change.rank.name(ps0)


#1.1定义分组扩增子环境布置#-------
#  定义分组
# map = ps0 %>% sample_data()
# head(map)

# map$Group = as.factor(map$Group) %>% as.numeric()
# map$Group = gsub("1","Before",map$Group)
# map$Group = gsub("2","After",map$Group)

# map$Group %>% unique()

# #去除分组列中的特殊字符
# map$Group =  gsub("[-]",".", map$Group)
# map$species = map$Group %>%strsplit( "[.]") %>% sapply(`[`, 1) 

# sample_data(ps0) = map


#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my(ps0)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
# 构建保存结果文件夹
result<- dir.amp(ps0 = ps0,smart = TRUE)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps0 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps0)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL

# 热图展示的OTU数量
heatnum　=　30


#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40 # 选择多少个重要变量

#--功能预测
if (is.null(ps0@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps0@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps0 %>% ggClusterNet::filter_OTU_ps(1000)
if (Tax4Fun2) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}



#-----选择性功能

#设置CK，用于双向柱状图绘制-目前不绘制
CK = unique(phyloseq::sample_data(ps)$Group)[1]

# 用python做lefse
lefse.py = T
if (lefse.py) {
  lefsenum = 0
  ps_lefse <- ps %>%
    phyloseq::subset_taxa(
      # Kingdom == "Fungi"
      Kingdom == id
      # Genus  == "Genus1"
      # Species %in%c("species1") 
      # row.names(tax_table(ps0))%in%c("OTU1")
    )
  
  ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)
  
}


#1.2定义分组"1"扩增子环境布置#-------
#  定义分组
map = ps0 %>% sample_data()
head(map)
map$Group = map$group
map$Group %>% unique()
#去除分组列中的特殊字符
# map$Group =  gsub("[-]",".", map$Group)
map$species = map$Group %>%strsplit( "[_]") %>% sapply(`[`, 1) 
sample_data(ps0) = map

ps1 = ps0 %>% subset_samples.wt("species","one")
map = ps1 %>% sample_data()
head(map)
#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my(ps1)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
# 构建保存结果文件夹
result<- dir.amp(ps0 = ps0,smart = TRUE)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/Sub_OTU_",a,"/",sep = "");otupath
dir.create(otupath)


#--最终确定的phyloseq对象定义为ps
ps = ps1 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("Group2","Group3","Group1")

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
jj = j = "Phylum"

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL

# 热图展示的OTU数量
heatnum　=　30


#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40# 选择多少个重要变量

#--功能预测

if (is.null(ps1@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps1@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps1 %>% ggClusterNet::filter_OTU_ps(1000)
if (Tax4Fun2) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}


#-----选择性功能

#设置CK，用于双向柱状图绘制-目前不绘制
CK = unique(phyloseq::sample_data(ps)$Group)[1]

# 用python做lefse
lefse.py = T
if (lefse.py) {
  lefsenum = 0
  ps_lefse <- ps %>%
    phyloseq::subset_taxa(
      # Kingdom == "Fungi"
      Kingdom == id
      # Genus  == "Genus1"
      # Species %in%c("species1") 
      # row.names(tax_table(ps0))%in%c("OTU1")
    )
  
  ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)
  
}


#1.3定义分组"2"扩增子环境布置#-------
ps0 = readRDS("./231106/16s_ps.rds")

ps0
#  定义分组
map = ps0 %>% sample_data()
head(map)

map$Group %>% unique()
#去除分组列中的特殊字符
# map$Group =  gsub("[-]",".", map$Group)
map$species = map$Group %>%strsplit( "[_]") %>% sapply(`[`, 1) 
sample_data(ps0) = map

ps1 = ps0 %>% subset_samples.wt("species","two")
map = ps1 %>% sample_data()
head(map)
#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my(ps1)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
# 构建保存结果文件夹
result<- dir.amp(ps0 = ps0,smart = TRUE)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/Sub_OTU_two",a,"/",sep = "");otupath
dir.create(otupath)


#--最终确定的phyloseq对象定义为ps
ps = ps1 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("Group2","Group3","Group1")

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
jj = j = "Phylum"

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL

# 热图展示的OTU数量
heatnum　=　30


#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40# 选择多少个重要变量

#--功能预测

if (is.null(ps1@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps1@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps1 %>% ggClusterNet::filter_OTU_ps(1000)
if (Tax4Fun2) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}


#1.4定义分组"2"扩增子环境布置#-------
ps0 = readRDS("./231106/16s_ps.rds")

ps0
#  定义分组
map = ps0 %>% sample_data()
head(map)

map$Group %>% unique()
#去除分组列中的特殊字符
# map$Group =  gsub("[-]",".", map$Group)
map$species = map$Group %>%strsplit( "[_]") %>% sapply(`[`, 1) 
sample_data(ps0) = map

ps1 = ps0 %>% subset_samples.wt("species","three")
map = ps1 %>% sample_data()
head(map)
#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my(ps1)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
# 构建保存结果文件夹
result<- dir.amp(ps0 = ps0,smart = TRUE)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/Sub_OTU_three",a,"/",sep = "");otupath
dir.create(otupath)


#--最终确定的phyloseq对象定义为ps
ps = ps1 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("Group2","Group3","Group1")

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
jj = j = "Phylum"

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL

# 热图展示的OTU数量
heatnum　=　30


#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40# 选择多少个重要变量

#--功能预测

if (is.null(ps1@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps1@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps1 %>% ggClusterNet::filter_OTU_ps(1000)
if (Tax4Fun2) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}


#1.5定义分组"2"扩增子环境布置#-------
ps0 = readRDS("./231106/16s_ps.rds")

ps0
#  定义分组
map = ps0 %>% sample_data()
head(map)

map$Group %>% unique()
#去除分组列中的特殊字符
# map$Group =  gsub("[-]",".", map$Group)
map$species = map$Group %>%strsplit( "[_]") %>% sapply(`[`, 2) 
map$Group = map$species
sample_data(ps0) = map

ps1 = ps0 #%>% subset_samples.wt("species","three")
map = ps1 %>% sample_data()
head(map)
#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my(ps1)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
# 构建保存结果文件夹
result<- dir.amp(ps0 = ps0,smart = TRUE)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/Sub_OTU_allenv",a,"/",sep = "");otupath
dir.create(otupath)


#--最终确定的phyloseq对象定义为ps
ps = ps1 %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("Group2","Group3","Group1")

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
jj = j = "Phylum"

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL

# 热图展示的OTU数量
heatnum　=　30


#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40# 选择多少个重要变量

#--功能预测

if (is.null(ps1@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps1@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps1 %>% ggClusterNet::filter_OTU_ps(1000)
if (Tax4Fun2) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     
  
}