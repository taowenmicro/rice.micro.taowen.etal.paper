
library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果
ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1","No3"),T)
sample_data(ps02) %>% head()

map = sample_data(ps02)

map$Group = paste0(map$group,"_",map$treated,"_",map$species)
map$Group  %>% unique()
sample_data(ps02) = map
ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
     )
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_No2",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c("WT_NA_ZH11","WT_M.oryzae_ZH11","WT_NA_NIP","WT_M.oryzae_NIP" ,
    "Osmpk6_NA_ZH11","Osmpk6_M.oryzae_ZH11" ,"CA.OsMPK6_NA_NIP",  "CA.OsMPK6_M.oryzae_NIP",
   "WT_XOO_ZH11","Osmpk6_XOO_ZH11", "WT_XOO_NIP", "CA.OsMPK6_XOO_NIP")   
 
              
              



# axis_order = c(
#   "WT_CK","osb561.2_CK","osb561.1.2_CK","OE.OsB561.1_CK","OE.OsB561.2_CK",
#   "WT_M.oryzae","osb561.2_M.oryzae","osb561.1.2_M.oryzae","OE.OsB561.1_M.oryzae",
#   "OE.OsB561.2_M.oryzae", "WT_XOO","osb561.2_XOO" ,"osb561.1.2_XOO" ,"OE.OsB561.1_XOO",
#   "OE.OsB561.2_XOO" 
# ),
# axis_order =c("WT_CK_NA"   ,             "osb561.2_CK_4"      ,     "osb561.1.2_CK_9"    ,    
#   "osb561.1.2_CK_38"    ,    "OE.OsB561.1_CK_8"     ,   "OE.OsB561.1_CK_11"     , 
#   "OE.OsB561.2_CK_11"    ,   "OE.OsB561.2_CK_13" , "WT_M.oryzae_NA"   ,      
#   "osb561.2_M.oryzae_4"    , "osb561.1.2_M.oryzae_9" ,  "osb561.1.2_M.oryzae_38" ,
#   "OE.OsB561.1_M.oryzae_8" , "OE.OsB561.1_M.oryzae_11" ,"OE.OsB561.2_M.oryzae_11",
#   "OE.OsB561.2_M.oryzae_13" ,"WT_XOO_NA"            ,   "osb561.2_XOO_4"       ,  
#   "osb561.1.2_XOO_9"      ,  "osb561.1.2_XOO_38"    ,   "OE.OsB561.1_XOO_8"    ,  
#   "OE.OsB561.1_XOO_11"    ,  "OE.OsB561.2_XOO_11"    ,  "OE.OsB561.2_XOO_13"  )  
#     

# axis_order = c("WT_NA", "WT_CK" ,"Kitaake_CK",  "Osmpk6_NA", "CA.OsMPK6_NA" ,"oswrky24.53.70_CK" ,
#   "OsCPS4.OE_CK","oscps2oscps4_CK","OsCPS2.OE_CK","oscyp76m7oscyp76m8_CK",
#   "WT_XOO"  ,"Osmpk6_XOO","CA.OsMPK6_XOO",                    
#   "WT_M.oryzae","Kitaake_M.oryzae","Osmpk6_M.oryzae", "CA.OsMPK6_M.oryzae" ,
#   "oswrky24.53.70_M.oryzae","OsCPS2.OE_M.oryzae","oscps2oscps4_M.oryzae","OsCPS4.OE_M.oryzae",
#   "oscyp76m7oscyp76m8_M.oryzae")               
#                                             
                    
                        





axis_order.s = phyloseq::sample_data(ps) %>% row.names()
                      
                             
                               
              

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
if (is.null(ps@refseq)) {
  Tax4Fun2 = FALSE
} else if(!is.null(ps0@refseq)){
  Tax4Fun2 = TRUE
}

ps.t = ps %>% ggClusterNet::filter_OTU_ps(1000)
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




# No3的细菌群落特征#-------------



library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1","No2"),T)
sample_data(ps02) %>% head()

map = sample_data(ps02)

map$Group = paste0(map$group,"_",map$treated,"_",map$species)
map$Group  %>% unique()
sample_data(ps02) = map
ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_No3",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c("WT_CK_NIP", "WT_M.oryzae_NIP",
    "Kitaake_CK_Kitaake","Kitaake_M.oryzae_Kitaake",              
    "oswrky24.53.70_CK_NIP", "oswrky24.53.70_M.oryzae_NIP",                                    
    "oscps2oscps4_CK_Kitaake", "oscps2oscps4_M.oryzae_Kitaake",                          
    "OsCPS2.OE_CK_Kitaake","OsCPS2.OE_M.oryzae_Kitaake",
    "OsCPS4.OE_CK_Kitaake",                       
    "OsCPS4.OE_M.oryzae_Kitaake","oscyp76m7oscyp76m8_CK_Kitaake", 
    "oscyp76m7oscyp76m8_M.oryzae_Kitaake"              
  )   


axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"




# WT和接菌对照看看的细菌群落特征#-------------



library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 #%>% subset_samples.wt("experiment2",c("MAPK.Singl"))

sample_data(ps01)$group

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1"),T)
ps02 = ps02 %>% subset_samples.wt("group",c("WT","Kitaake"))
sample_data(ps02) %>% head()

map = sample_data(ps02)


map$Group = paste0(map$group,"_",map$treated,"_",map$species,"..",
                   map$experiment2,map$experiment1)
map$Group  %>% unique()
sample_data(ps02) = map
ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_WT.all",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c( "WT_NA_ZH11..MAPK.SinglNo2",
     "WT_M.oryzae_ZH11..MAPK.SinglNo2", 
     "WT_NA_NIP..MAPK.SinglNo2",
     "WT_M.oryzae_NIP..MAPK.SinglNo2",   
     
     "WT_CK_NIP..MAPK.SinglNo3",
     "WT_M.oryzae_NIP..MAPK.SinglNo3",   
     
     "WT_CK_NIP..B561No3",  
     "WT_M.oryzae_NIP..B561No3",
     
     "WT_XOO_ZH11..MAPK.SinglNo2",
     "WT_XOO_NIP..MAPK.SinglNo2",
     "WT_XOO_NIP..B561No3",
     "Kitaake_CK_Kitaake..MAPK.SinglNo3"    ,
     "Kitaake_M.oryzae_Kitaake..MAPK.SinglNo3"
  )   
           



axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"



# No3的细菌群落特征,使用No2对照#-------------



library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1"),T)

sample_data(ps02) %>% head()
map = sample_data(ps02)
map$Group = paste0(map$group,"_",map$experiment1)
map$Group  %>% unique()
sample_data(ps02) = map

ps02 = ps02 %>% subset_samples.wt("Group",c("Osmpk6_No2", "CA.OsMPK6_No2",
                                            "WT_No3"
                                            
                                            ),T)
map = sample_data(ps02)

map$Group = paste0(map$group,"_",map$treated,"_",map$species)
map$Group  %>% unique()
sample_data(ps02) = map


ps02 = ps02 %>% subset_samples.wt("Group",c("WT_NA_ZH11" ,
                                            "WT_XOO_ZH11","WT_XOO_NIP",
                                            "WT_M.oryzae_ZH11"
                                            
),T)

# map = sample_data(ps02)
# head(map)
# map$Group = paste0(map$group,"_",map$treated,"_",map$species,map$experiment1)
# map$Group  %>% unique()
# sample_data(ps02) = map


ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_No3_use_No2_CK",a,"/",sep = "");otupath
dir.create(otupath)
print(a)


#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)


#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c( "WT_NA_NIP","WT_M.oryzae_NIP",              
     "Kitaake_CK_Kitaake","Kitaake_M.oryzae_Kitaake",   
     "oswrky24.53.70_CK_NIP", "oswrky24.53.70_M.oryzae_NIP",
     "oscps2oscps4_CK_Kitaake","oscps2oscps4_M.oryzae_Kitaake",                            
     "OsCPS2.OE_CK_Kitaake","OsCPS2.OE_M.oryzae_Kitaake",
     "OsCPS4.OE_CK_Kitaake","OsCPS4.OE_M.oryzae_Kitaake",                        
     "oscyp76m7oscyp76m8_CK_Kitaake","oscyp76m7oscyp76m8_M.oryzae_Kitaake"          
  )   


axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"



# No1的细菌群落特征,使用No2对照#-------------



library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1"),F)

sample_data(ps02) %>% head()
map = sample_data(ps02)
map$Group = paste0(map$group,"_",map$experiment1)
map$Group  %>% unique()
sample_data(ps02) = map


map = sample_data(ps02)

map$Group = paste0(map$group,"_",map$treated,"_",map$species)
map$Group  %>% unique()
sample_data(ps02) = map



ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_No1",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c( "WT_NA_NIP",
     
     "Osmkk4_NA_NIP","OsMKK4WT_NA_NIP","OsMKK4DD_NA_NIP", "OsMKK4DD_DEX_NIP",
     
     "WT_NA_ZH11" ,"Osmpk6_NA_ZH11","CA.OsMPK6_NA_NIP"       
  )   


   
     

axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"

# No2微生物群落分析#---------

library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1","No3"),T)
sample_data(ps02) %>% head()

map = sample_data(ps02)

map$Group = paste0(map$group,"_",map$treated,"_",map$species)
map$Group  %>% unique()
sample_data(ps02) = map
ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_No2",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c("WT_NA_ZH11","WT_M.oryzae_ZH11","WT_NA_NIP","WT_M.oryzae_NIP",
    "Osmpk6_NA_ZH11" ,"Osmpk6_M.oryzae_ZH11" ,
    "CA.OsMPK6_NA_NIP" ,"CA.OsMPK6_M.oryzae_NIP",
    "WT_XOO_ZH11","Osmpk6_XOO_ZH11","WT_XOO_NIP","CA.OsMPK6_XOO_NIP"             
  )   




axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"

# B561结果#-------

library(ggClusterNet)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()
# 包含了四个批次的实验，第一个纯突变株种植；第二个是接病突变株种植；第三个是WARK接病结果；
# 第四个是植保素接病实验
# 但是测序的时候将第三批和第四批同时命名为NO3





#  同时还检测了一批B561突变株的结果，接种了细菌病原菌和真菌病原菌，测定的细菌群落
#  由于这批数据较为简单，首先挖掘被B561的结果

ps01 = ps00 %>% subset_samples.wt("experiment2",c("B561"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
sample_data(ps02) %>% head()

map = sample_data(ps02)
map$Group = paste0(map$group,"_",map$treated,"_",map$Line)
map$Group %>% unique()

sample_data(ps02) = map
ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_line_B561",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c(
#   "WT_CK","osb561.2_CK","osb561.1.2_CK","OE.OsB561.1_CK","OE.OsB561.2_CK",
#   "WT_M.oryzae","osb561.2_M.oryzae","osb561.1.2_M.oryzae","OE.OsB561.1_M.oryzae",
#   "OE.OsB561.2_M.oryzae", "WT_XOO","osb561.2_XOO" ,"osb561.1.2_XOO" ,"OE.OsB561.1_XOO",
#   "OE.OsB561.2_XOO" 
# ),
axis_order =c("WT_CK_NA"   ,    "WT_M.oryzae_NA" ,"WT_XOO_NA"       ,
              
              "osb561.2_CK_4"      ,  "osb561.2_M.oryzae_4"    , "osb561.2_XOO_4",
              "osb561.1.2_CK_9"    ,"osb561.1.2_M.oryzae_9" ,  "osb561.1.2_XOO_9"      ,
              "osb561.1.2_CK_38"    , "osb561.1.2_M.oryzae_38" , "osb561.1.2_XOO_38"    ,
                  
                 "OE.OsB561.1_CK_8"     ,  "OE.OsB561.1_M.oryzae_8" ,"OE.OsB561.1_XOO_8"    ,
              
              "OE.OsB561.1_CK_11"     ,  "OE.OsB561.1_M.oryzae_11" ,"OE.OsB561.1_XOO_11"    , 
              
              
              "OE.OsB561.2_CK_11" , "OE.OsB561.2_M.oryzae_11","OE.OsB561.2_XOO_11"    ,
              "OE.OsB561.2_CK_13" , "OE.OsB561.2_M.oryzae_13"   , "OE.OsB561.2_XOO_13"  )  



axis_order.s = phyloseq::sample_data(ps) %>% row.names()





# axis_order = c("KO","WT","OE")
# axis_order = c("Before","After")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件

#--物种分类树展示的OTU数量
Top_micro = 150

#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12

#  全部通路K6-WARY-CPS-X-一起分析#-------


library(ggClusterNet)
library(tidyverse)
library(phyloseq)
#  增碳MAPK幸好通路的一整个通路伴随着微生物群落的变化

ps00 = readRDS("./ps.16s.mapall.rds")
ps00

sample_data(ps00) %>% head()


ps01 = ps00 %>% subset_samples.wt("experiment2",c("MAPK.Singl"))
sample_data(ps01) %>% head()

#  首先分析根际的结果

ps02 = ps01 %>% subset_samples.wt("zone",c("rhi"))
ps02 = ps02 %>% subset_samples.wt("experiment1",c("No1"),T)
sample_data(ps02) %>% head()
map = sample_data(ps02)
map$Group = paste0(map$group,"_",map$treated,"_",map$species,map$experiment1)
map$Group  %>% unique()
sample_data(ps02) = map

ps02 = ps02 %>% subset_samples.wt("Group",c("WT_XOO_NIPNo2",
                                            "WT_XOO_ZH11No2",
                                            "CA.OsMPK6_XOO_NIPNo2",
                                            "Osmpk6_XOO_ZH11No2" ,
                                            "WT_CK_NIPNo3",
                                            "WT_M.oryzae_NIPNo3"),T)


map = sample_data(ps02)
map$Group = paste0(map$group,"_",map$treated,"_",map$species,map$experiment1,"_",map$Line)
map$Group  %>% unique()

sample_data(ps02) = map


ps = ps02



#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];
colset2 = res[[4]];
colset3 = res[[5]];
colset4 = res[[6]]

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A",  "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD",
     "#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C"
)
colset2 = mi

# 构建保存结果文件夹
result<- dir.amp(ps0 = ps,smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path

#-系统会判断这个数据是什么类型：细菌？真菌？
id = result[[2]];id



#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/OTU_MAPK_all_",a,"/",sep = "");otupath
dir.create(otupath)
print(a)

#--最终确定的phyloseq对象定义为ps
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)

#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order

axis_order =  
  c("WT_NA_ZH11No2_NA","WT_M.oryzae_ZH11No2_NA",
    "WT_NA_NIPNo2_NA",  "WT_M.oryzae_NIPNo2_NA",
    "Kitaake_CK_KitaakeNo3_NA","Kitaake_M.oryzae_KitaakeNo3_NA",
    "Osmpk6_NA_ZH11No2_NA", "Osmpk6_M.oryzae_ZH11No2_NA",                                          
    "CA.OsMPK6_NA_NIPNo2_3", "CA.OsMPK6_M.oryzae_NIPNo2_3",                                  
    "CA.OsMPK6_NA_NIPNo2_11","CA.OsMPK6_M.oryzae_NIPNo2_11",                                 
    "oswrky24.53.70_CK_NIPNo3_48", "oswrky24.53.70_M.oryzae_NIPNo3_48",                                   
    "oscps2oscps4_CK_KitaakeNo3_NA","oscps2oscps4_M.oryzae_KitaakeNo3_NA",             
    "OsCPS2.OE_CK_KitaakeNo3_15","OsCPS2.OE_M.oryzae_KitaakeNo3_15",               
    "OsCPS4.OE_CK_KitaakeNo3_15","OsCPS4.OE_M.oryzae_KitaakeNo3_15",                       
    "oscyp76m7oscyp76m8_CK_KitaakeNo3_11","oscyp76m7oscyp76m8_M.oryzae_KitaakeNo3_11" ,    
    "oscyp76m7oscyp76m8_CK_KitaakeNo3_15",
    "oscyp76m7oscyp76m8_M.oryzae_KitaakeNo3_15"  )   


  

