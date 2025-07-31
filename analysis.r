###########  ENV  ###########
library("dplyr")
library("tidyr")
library("ggplot2")
library(Rmisc)
library(agricolae)
library(lattice)
library(reshape)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library("patchwork")
library("ggplotify")
library("ggsci")
library(tidyverse)
library(ggthemes)
library(devtools)
library(tidyquant)
library(broom)
library(export)
library("Hmisc")
library("grafify")
library("job")
library(readxl)
library(ggpmisc)
library("scales")
library("paletteer")
library(ggalluvial)

theme_ykm<- function(..., bg='white'){
  theme_few() +
    theme(aspect.ratio =1,
          axis.text=element_text(size=7),
          axis.title=element_text(colour="black",size=7),
          legend.text=element_text(size=7),
          legend.title=element_text(size=7,colour="black"),
          plot.title=element_text(hjust = 0.5,size=7,colour="black"),
          plot.subtitle=element_text(size=7),
          plot.caption=element_text(size=7),
          plot.tag=element_text(size=7),
          strip.text=element_text(size=7))
  
}

dataframe2matrix  <- function(input){
  dataframe = as.data.frame(input)
  rownames(dataframe) = dataframe[,1] 
  return(dataframe[,-1])
}
format_percent <- function(input){
  out = input
  for (x in 1:length(input)){
    out[x] = input[x]/sum(input)
  }
  return(out)
}
filter_by_sumpercent <- function(input,cutoff){
  data.new3 = mutate(input,
                     percent = rowSums(dataframe2matrix(input))/sum(dataframe2matrix(input)))
  data.new4 = subset(data.new3,percent > cutoff)[,-ncol(data.new3)]
  return(data.new4)
}
filter_by_rank <- function(input,cutoff){
  data.new = mutate(input,
                    sum = apply(dataframe2matrix(input),1,sum))
  data.new2 = data.new[order(data.new$sum,decreasing=TRUE)[1:cutoff],-ncol(data.new)]
  return(data.new2)
}
my.t.test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj)
}
cor_matrix2matrix = function(matrix_a,matrix_b,cor_method){
  phage_num = as.numeric(nrow(matrix_a))
  species_num = as.numeric(nrow(matrix_b))
  phage_attribute = as.data.frame(matrix_a)[,1]
  species_attribute = as.data.frame(matrix_b)[,1]
  phage_attribute2 = rep(phage_attribute,each = species_num)
  cor_result = matrix(ncol = 2,nrow = phage_num*species_num)
  colnames(cor_result) = c("p","r")
  m_num = seq(1:phage_num)
  n_num = seq(1:species_num)
  h.phage.otu = matrix_a[,2:ncol(matrix_a)]
  h.species.otu = matrix_b[,2:ncol(matrix_b)]
  for (m in m_num){
    for (n in n_num){
      cor = cor.test(as.numeric(h.phage.otu[m,]),as.numeric(h.species.otu[n,]),method= cor_method )
      cor_result[(m-1)*species_num + n,1] = as.numeric(cor[3])
      cor_result[(m-1)*species_num + n,2] = as.numeric(cor[4])
    }
  }
  a_names = rep(phage_attribute,each = species_num)
  b_names = rep(species_attribute,phage_num)
  result = data.frame(a_names,b_names,cor_result)
  return(result)
}

table_order <- function(x) {
  table(x)[order(table(x))]
}


###########  map  ###########
library(ggspatial)
library(sf)

my_data <- subset(shp_data, NAME_1 == "Nei Mongol") |> 
  merge(show_data,by.y = "City",by.x="NAME_2",all=T)

sample_map.plot=ggplot(data = my_data) + geom_sf(aes( geometry = `geometry`)) + 
  geom_sf_text(aes(label = NAME_2,geometry = `geometry`), color = 'Black',size=2)+
  xlab("Longtitude (°E)") + ylab("Latitude (°N)") + 
  scale_fill_grafify()+
  theme(panel.background = element_rect(fill = "white",color = "black"),
        panel.grid = element_line(color = "lightgrey"))+
  ##添加指北针，“style”参数可以更改样式
  annotation_north_arrow(location='tl', which_north='false',
                         style=north_arrow_orienteering())+
  geom_point(aes(x=longtitude,y=latitude,shape=group,color=location),size=2)+
  scale_color_jco()+
  theme(legend.position="none")








###########  arg+mge+tax  ###########
library(ggtree)
library(aplot)
library(ggridges)
library(corrplot)

arg.gather =gather(arg.select,-1,key="sample",value ="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  unite("field", location, field, sep = "|", remove = F)
input = arg.gather[,c(1,3,5,7)]
ttest.result = data.frame(field="y",Type="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$Type)){
    data.pro = subset(data.sub,Type == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,Type=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("arg", field, Type, sep = "-", remove = FALSE)
arg.summary <- summarySE(data=arg.gather, measurevar="value",groupvars=c("Type","field","location","group"))
arg.spread = spread(arg.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("arg", field, Type, sep = "-", remove = FALSE) %>% 
  merge(ttest.result2,"arg") %>% 
  mutate(updown="none")
arg.spread$updown <- case_when(  
  arg.spread$fc < 0 & arg.spread$p < 0.05 ~ 'grass'   ,  
  arg.spread$fc > 0 & arg.spread$p < 0.05 ~ 'farm',
  .default = "none")
arg.spread2 = spread(arg.spread[,c(2,3,7)],key = "field.x",value = "fc")
arg.clust.plot = hclust(dist((dataframe2matrix(arg.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 1)+
  hexpand(2) 
clust.data = subset(arg.clust.plot$data,is.na(label)==F)
arg.spread$Type.x = factor(arg.spread$Type.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
arg.plot = ggplot(arg.spread,aes(x=field.x,y=Type.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3,5), range = c(1, 8))+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(24,25,21)) 

mge.gather =gather(mge.select,-1,key="sample",value ="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  unite("field", location, field, sep = "|", remove = F)
input = mge.gather[,c(1,3,5,7)]
ttest.result = data.frame(field="y",Type="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$Type)){
    data.pro = subset(data.sub,Type == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,Type=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("mge", field, Type, sep = "-", remove = FALSE)
mge.summary <- summarySE(data=mge.gather, measurevar="value",groupvars=c("Type","field","location","group"))
mge.spread = spread(mge.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("mge", field, Type, sep = "-", remove = FALSE) %>% 
  merge(ttest.result2,"mge") %>% 
  mutate(updown="none")
mge.spread$updown <- case_when(  
  mge.spread$fc < 0 & mge.spread$p < 0.05 ~ 'grass'   ,  
  mge.spread$fc > 0 & mge.spread$p < 0.05 ~ 'farm',
  .default = "none")
mge.spread2 = spread(mge.spread[,c(2,3,7)],key = "field.x",value = "fc")
mge.clust.plot = hclust(dist((dataframe2matrix(mge.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 0.2)+
  hexpand(2) 
clust.data = subset(mge.clust.plot$data,is.na(label)==F)
mge.spread$Type.x = factor(mge.spread$Type.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
mge.plot = ggplot(mge.spread,aes(x=field.x,y=Type.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+theme(aspect.ratio = 0.2)+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3, 5), range = c(1, 8))+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(24,25,21))

arg.plot %>% insert_left(arg.clust.plot)   / mge.plot
(arg.clust.plot /mge.clust.plot) | (arg.plot /mge.plot)

arg_group.plot=ggboxplot(arg.gather,x="group",y="value",color = "group",
                         add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+
  facet_wrap(~Type,ncol=6,scales = "free")

arg.gather.polymyxin = subset(arg.gather,Type == "polymyxin")
arg.gather.multidrug = subset(arg.gather,Type == "multidrug")
arg_polymyxin.plot=ggboxplot(arg.gather.polymyxin,x="group",y="value",color = "group",
                             add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

ggboxplot(arg.gather.polymyxin,x="group",y="value",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.format",vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

arg_multidrug.plot=ggboxplot(arg.gather.multidrug,x="group",y="value",color = "group",
                             add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

resistome_shannon_group.plot | arg_multidrug.plot | arg_polymyxin.plot

multidrug.group.ave=summarySE(data=arg.gather.multidrug, measurevar="value",groupvars="group")
polymyxin.group.ave=summarySE(data=arg.gather.polymyxin, measurevar="value",groupvars="group")

### sum

arg.sum = data.frame(sample=colnames(arg)[-1],sum =apply(arg[,-1],2,sum)) %>% 
  merge(meta,1) %>% 
  mutate(resistome="arg")
ggplot(arg.sum, aes(x = sum, y = field, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  scale_fill_jco()+
  theme_ykm()
mge.sum = data.frame(sample=colnames(mge)[-1],sum =apply(mge[,-1],2,sum)) %>% 
  merge(meta,1)%>% 
  mutate(resistome="mge")
ggplot(mge.sum, aes(x = sum, y = field, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  scale_fill_jco()+
  theme_ykm()

resistome = rbind.data.frame(arg.sum,mge.sum)
resistome_sum.plot=ggplot(resistome,aes(x=field,y=sum,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")+
  facet_wrap(~resistome,ncol=1,scales = "free")

sample_map.plot / resistome_sum.plot

input = resistome[,c(1:4,20)]
ttest.result = data.frame(resistome="x",field="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$resistome)){
  data.sub = subset(input,resistome == i )
  for (j in unique(data.sub$field)){
    data.pro = subset(data.sub,field == j )
    result1 <- my.t.test(sum~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(resistome=i,field=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(resistome=i,field=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
resistome.summary <- summarySE(data=resistome, measurevar="sum",groupvars=c("field","group","resistome"))
resistome.ave = spread(resistome.summary[,c(1,2,3,5)],key = "group",value = "sum") %>% 
  unite("class", field, resistome, sep = "-", remove = FALSE)
resistome.ttest = read.table("/home/shared/project/2025/ecotone90/analysis/result/arg_mge_sum.ttest.txt", sep = "\t",head = T,quote = "",check.names = F)%>% 
  unite("class", field, resistome, sep = "-", remove = FALSE)
resistome.merge=merge(resistome.ave,resistome.ttest,"class")
#write.table (resistome.merge, file = "/home/shared/project/2025/ecotone90/analysis/result/resistome.sum.merge.txt", sep = "\t", row.names = F ,col.names = TRUE, quote = FALSE)

resistome_group.plot=ggboxplot(resistome,x="group",y="sum",color = "group",
                               add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",paired = T)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+
  facet_wrap(~resistome,ncol=2,scales = "free")

#### correlation
resistome.line = spread(resistome[,c(1:4,16,20)],key = "resistome",value = "sum")
ggscatter(resistome.line, x = "arg", y = "mge", add = "reg.line",cor.coef = T,color="group") +
  stat_poly_eq(formula = y ~ x,aes(color=group,label=paste( ..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_grafify(palette = "kelly")+
  theme(legend.position=c(0.9,0.9))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))+
  xlab("arg") + ylab("mge")+
  facet_wrap(~location,scales = "free")

arg_mge.cor = cor_matrix2matrix(arg,mge,"spearman")
cor.p = dataframe2matrix(spread(arg_mge.cor[,-4],key = "b_names",value = "p"))
cor.r = dataframe2matrix(spread(arg_mge.cor[,-3],key = "b_names",value = "r"))
corrplot(
  t(cor.r), 
  p.mat = t(cor.p), 
  sig.level = 0.05,
  insig = "blank",
  method = "circle",
  is.corr = TRUE,
  col.lim = c(-1, 1),
  tl.col = "black",
  addgrid.col = "black",  # 关键参数：设置边框颜色
  tl.cex = 0.5
)

arg.fa = subset(meta,group=="FA")
arg.gr = subset(meta,group=="GR")
arg.fa.data = arg[,colnames(arg) %in% c("Type",arg.fa$sample)]
arg.gr.data = arg[,colnames(arg) %in% c("Type",arg.gr$sample)]
fa.cor = cor_matrix2matrix(arg,mge,"spearman")
mge.fa = subset(meta,group=="FA")
mge.gr = subset(meta,group=="GR")
mge.fa.data = mge[-4,colnames(mge) %in% c("Type",mge.fa$sample)]
mge.gr.data = mge[-4,colnames(mge) %in% c("Type",mge.gr$sample)]
fa.cor = cor_matrix2matrix(arg.fa.data,mge.fa.data,"spearman")
fa.cor.p = dataframe2matrix(spread(fa.cor[,-4],key = "b_names",value = "p"))
fa.cor.r = dataframe2matrix(spread(fa.cor[,-3],key = "b_names",value = "r"))
corrplot(
  t(fa.cor.r), 
  p.mat = t(fa.cor.p), 
  sig.level = 0.05,
  insig = "blank",
  method = "circle",
  is.corr = TRUE,
  col.lim = c(-1, 1),
  tl.col = "black",
  addgrid.col = "black",  # 关键参数：设置边框颜色
  tl.cex = 0.5
)
gr.cor = cor_matrix2matrix(arg.gr.data,mge.gr.data,"spearman")
gr.cor.p = dataframe2matrix(spread(gr.cor[,-4],key = "b_names",value = "p"))
gr.cor.r = dataframe2matrix(spread(gr.cor[,-3],key = "b_names",value = "r"))
corrplot(
  t(gr.cor.r), 
  p.mat = t(gr.cor.p), 
  sig.level = 0.05,
  insig = "blank",
  method = "circle",
  is.corr = TRUE,
  col.lim = c(-1, 1),
  tl.col = "black",
  addgrid.col = "black",  # 关键参数：设置边框颜色
  tl.cex = 0.5
)

######  tax phylum

phylum.sum = data.frame(phylum = names(rowSums(dataframe2matrix(phylum))),value=rowSums(dataframe2matrix(phylum)))
top20 = phylum.sum[order(phylum.sum[,"value"],decreasing=T)[1:20],] 


phylum.select = subset(phylum,phylum %in% top20$phylum)
phylum.gather =gather(phylum.select,-1,key="sample",value ="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  unite("field", location, field, sep = "|", remove = F)
input = phylum.gather[,c(1,3,5,7)]
ttest.result = data.frame(field="y",phylum="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$phylum)){
    data.pro = subset(data.sub,phylum == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,phylum=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,phylum=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("class", field, phylum, sep = "-", remove = FALSE)
phylum.summary <- summarySE(data=phylum.gather, measurevar="value",groupvars=c("phylum","field","location","group"))
phylum.spread = spread(phylum.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("class", field, phylum, sep = "-", remove = FALSE) %>% 
  merge(ttest.result2,"class") %>% 
  mutate(updown="none")
phylum.spread$updown <- case_when(  
  phylum.spread$fc < 0 & phylum.spread$p < 0.05 ~ 'grass'   ,  
  phylum.spread$fc > 0 & phylum.spread$p < 0.05 ~ 'farm',
  .default = "none")
phylum.spread2 = spread(phylum.spread[,c(2,3,7)],key = "field.x",value = "fc")
phylum.clust.plot = hclust(dist((dataframe2matrix(phylum.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 1)+
  hexpand(2) 
clust.data = subset(phylum.clust.plot$data,is.na(label)==F)
phylum.spread$phylum.x = factor(phylum.spread$phylum.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
phylum.plot = ggplot(phylum.spread,aes(x=field.x,y=phylum.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3,5), range = c(1, 8))+
  #theme(legend.position = "none") +
  scale_shape_manual(values = c(24,25,21)) 

phylum.clust.plot | phylum.plot


#############  tax shannon richness  ############

eazy_alpha <- function(input){
  bacteria.m = as.matrix(t(dataframe2matrix(input) ))
  Shannon=diversity(bacteria.m,index='shannon')
  Simpson=diversity(bacteria.m,"simpson")
  invers_simpson=diversity(bacteria.m,"inv")
  Species_richness = specnumber(bacteria.m)
  Pielou_evenness = Shannon/log(Species_richness)
  sample = rownames(bacteria.m)
  alpha = data.frame(sample,Shannon,Simpson,invers_simpson,Species_richness,Pielou_evenness)
  return(alpha)
}
my.t.test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj)
}


tax.alpha=eazy_alpha(tax) %>% 
  merge(meta,1)
tax.shannon.plot=ggplot(tax.alpha,aes(x=field,y=Shannon,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(aes( group=group),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")
tax.richness.plot=ggplot(tax.alpha,aes(x=field,y=Species_richness,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(aes( group=group),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")

input = tax.alpha[,c(2,7,8)]
ttest.result = data.frame(field="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.pro = subset(input,field == i )
  result1 <- my.t.test(Shannon~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(field=i,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(field=i,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  
  ttest.result = rbind.data.frame(ttest.result,ttest)
}
ttest.result = ttest.result[-1,]  
shannon.summary <- summarySE(data=tax.alpha, measurevar="Shannon",groupvars=c("field","group"))
shannon.ave = spread(shannon.summary[,c(1,2,3,4)],key = "group",value = "Shannon") 
shannon.merge=merge(shannon.ave,ttest.result,1)

input = tax.alpha[,c(5,7,8)]
ttest.result = data.frame(field="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.pro = subset(input,field == i )
  result1 <- my.t.test(Species_richness~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(field=i,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(field=i,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  
  ttest.result = rbind.data.frame(ttest.result,ttest)
}
ttest.result = ttest.result[-1,]  
richness.summary <- summarySE(data=tax.alpha, measurevar="Species_richness",groupvars=c("field","group"))
richness.ave = spread(richness.summary[,c(1,2,3,4)],key = "group",value = "Species_richness") 
richness.merge=merge(richness.ave,ttest.result,1)

tax.richness.plot / tax.shannon.plot

tax_shannon_group.plot=ggboxplot(tax.alpha,x="group",y="Shannon",color = "group",
                                 add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")
ggboxplot(tax.alpha,x="group",y="Shannon",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")

tax_richness_group.plot=ggboxplot(tax.alpha,x="group",y="Species_richness",color = "group",
                                  add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")

(tax_shannon_group.plot | tax_richness_group.plot) / (tax_microbial_CN_group.plot| tax_avd_group.plot)

############## tax resistome alpha cor  ################################## 
tax_resistome.alpha = merge(tax.alpha[,c(1,2,5,8)],resistome.alpha[,c(1,2,5,8)],1)
ggscatter(tax_resistome.alpha, x = "Shannon.x", y = "Shannon.y", add = "reg.line",cor.coef = F,color = "group.x",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group.x,label=paste( ..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position=c(0.9,0.9))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))

ggscatter(tax_resistome.alpha, x = "Species_richness.x", y = "Species_richness.y", add = "reg.line",cor.coef = F,color = "group.x",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group.x,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position=c(0.9,0.9))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))

###########  resistome sum + shannon  ###########
eazy_alpha <- function(input){
  bacteria.m = as.matrix(t(dataframe2matrix(input) ))
  Shannon=vegan::diversity(bacteria.m,index='shannon')
  Simpson=vegan::diversity(bacteria.m,"simpson")
  invers_simpson=vegan::diversity(bacteria.m,"inv")
  Species_richness = specnumber(bacteria.m)
  Pielou_evenness = Shannon/log(Species_richness)
  sample = rownames(bacteria.m)
  alpha = data.frame(sample,Shannon,Simpson,invers_simpson,Species_richness,Pielou_evenness)
  return(alpha)
}
my.t.test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj)
}


arg_mge.gene = unique.data.frame(rbind.data.frame(arg.gene[,1:91],mge.gene[,1:91]))

arg.alpha=eazy_alpha(arg.gene[,1:91]) %>% 
  merge(meta,1)
mge.alpha=eazy_alpha(mge.gene[,1:91]) %>% 
  merge(meta,1)
ggplot(arg.alpha,aes(x=field,y=Shannon,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(aes(group=group),method = "t.test",vjust=0.1)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")

resistome.alpha=eazy_alpha(arg_mge.gene) %>% 
  merge(meta,1)
resistome.shannon.plot=ggplot(resistome.alpha,aes(x=field,y=Shannon,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(comparisons = list( c("FA", "GR") ),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")
resistome.richness.plot=ggplot(resistome.alpha,aes(x=field,y=Species_richness,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(comparisons = list( c("FA", "GR") ),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")

input = resistome.alpha[,c(2,7,8)]
ttest.result = data.frame(field="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.pro = subset(input,field == i )
  result1 <- my.t.test(Shannon~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(field=i,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(field=i,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  
  ttest.result = rbind.data.frame(ttest.result,ttest)
}
ttest.result = ttest.result[-1,]  
shannon.summary <- summarySE(data=resistome.alpha, measurevar="Shannon",groupvars=c("field","group"))
shannon.ave = spread(shannon.summary[,c(1,2,3,4)],key = "group",value = "Shannon") 
shannon.ttest = read.table("/home/shared/project/2025/ecotone90/analysis/result/resistome_shannon.ttest.txt", sep = "\t",head = T,quote = "",check.names = F)
shannon.merge=merge(shannon.ave,shannon.ttest,1)

resistome.sum = data.frame(sample=colnames(arg_mge.gene)[-1],sum =apply(arg_mge.gene[,-1],2,sum)) %>% 
  merge(meta,1)
input = resistome.sum[,c(1:4)]
ttest.result = data.frame(field="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.pro = subset(input,field == i )
  result1 <- my.t.test(sum~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(field=i,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(field=i,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  
  ttest.result = rbind.data.frame(ttest.result,ttest)
}
ttest.result = ttest.result[-1,]
resistome.summary <- summarySE(data=resistome.sum, measurevar="sum",groupvars=c("field","group"))
resistome.ave = spread(resistome.summary[,c(1,2,3,4)],key = "group",value = "sum") 
resistome.ttest = read.table("/home/shared/project/2025/ecotone90/analysis/result/resistome_sum.ttest.txt", sep = "\t",head = T,quote = "",check.names = F)
resistome.merge=merge(resistome.ave,resistome.ttest,1)
resistome.sum.plot=ggplot(resistome.sum,aes(x=field,y=sum,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(comparisons = list( c("FA", "GR") ),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")


resistome.sum.plot / resistome.shannon.plot 

resistome_shannon_group.plot=ggboxplot(resistome.alpha,x="group",y="Shannon",color = "group",
                                       add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

arg.alpha.plot=ggboxplot(arg.alpha,x="group",y="Shannon",color = "group",
                         add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test")+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

mge.alpha.plot=ggboxplot(mge.alpha,x="group",y="Shannon",color = "group",
                         add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test")+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")

arg.alpha.plot | mge.alpha.plot
###########  resistome +tax structure  ###########
library("ape")
library(Rtsne)

arg_mge.gene = unique.data.frame(rbind.data.frame(arg.gene[,1:91],mge.gene[,1:91]))
arg_mge.t = t(dataframe2matrix(arg_mge.gene))
distance.bray<-vegdist(arg_mge.t,method = 'bray')

tsne <- Rtsne(arg_mge.t,pca=T,dims=2,perplexity=12,theta=0)
tsne.resis = data.frame(sample = rownames(arg_mge.t),
                        tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])%>% 
  merge(meta,1)
plot.resis=ggplot(data=tsne.resis,aes(tsne1,tsne2)) +
  geom_point(aes(color=field,shape = group),size=2) +
  theme_ykm()+
  scale_color_grafify(palette = "kelly")

adonis.group = dataframe2matrix(meta) %>% 
  unite("landuse",group,plant,sep = "-",remove = F)
adonis.resis = adonis2(arg_mge.t~group+field+plant,data = adonis.group,permutations = 999,method="bray")

tax.t = t(dataframe2matrix(tax.data))
tsne <- Rtsne(tax.t,pca=T,dims=2,perplexity=12,theta=0)
tsne.tax = data.frame(sample = rownames(tax.t),
                      tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])%>% 
  merge(meta,1)
plot.tax=ggplot(data=tsne.tax,aes(tsne1,tsne2)) +
  geom_point(aes(color=field,shape = group),size=2) +
  theme_ykm()+
  scale_color_grafify(palette = "kelly")

adonis.tax = adonis2(tax.t~group+field+plant,data = adonis.group,permutations = 999,method="bray")

plot.resis | plot.tax

###########  resistome composition  ###########

arg.combine=arg
mge.combine=mge
arg.5 = filter_by_rank(arg,5)
mge.2 = filter_by_rank(mge,2)

arg.combine$Type[ arg.combine$Type %nin% arg.5$Type ] <- "Other_ARGs"
mge.combine$Type[ mge.combine$Type %nin% mge.2$Type ] <- "Other_MGEs"

arg.combine.sum = apply(arg.combine[,-1],2,aggregate,list(arg.combine$Type),sum)
arg.combine.df = data.frame(arg.combine.sum)[,c(1,seq(2,length(arg.combine.sum)*2,2))]
colnames(arg.combine.df)=colnames(arg.combine)
mge.combine.sum = apply(mge.combine[,-1],2,aggregate,list(mge.combine$Type),sum)
mge.combine.df = data.frame(mge.combine.sum)[,c(1,seq(2,length(mge.combine.sum)*2,2))]
colnames(mge.combine.df)=colnames(mge.combine)

arg_mge.combine = rbind.data.frame(arg.combine.df,mge.combine.df)
percent.arg = apply(arg_mge.combine[,-1],2,format_percent)
percent.arg2 = data.frame(Type =arg_mge.combine$Type,percent.arg)
percent.arg.gather =gather(percent.arg2,-1,key="sample",value="value") %>% 
  separate(sample,into = c("location","field","group"),sep = "\\.") %>% 
  unite("field", location, field, sep = "|", remove = FALSE)
percent.arg.sum = summarySE(data=percent.arg.gather, measurevar="value",groupvars=c("Type","group")) %>% 
  mutate(x="x")
percent.arg.sum.plot=ggplot(percent.arg.sum, aes(x = x, y = value*100,fill=Type))+   
  geom_col(width = 0.3,color="black",position = "stack")+  
  labs(x="", y="Proportion (%)")+  
  scale_fill_grafify(palette = "okabe_ito")+
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))+
  facet_wrap(~group,ncol=1)

###########  tax composition  ###########

tax.sep = separate(tax,Tax,into = c("kingdom","species"),sep = ";c__")[,-2]
tax.sum = apply(tax.sep[,-1],2,aggregate,list(tax.sep$kingdom),sum)
tax.sum.df = data.frame(tax.sum)[,c(1,seq(2,length(tax.sum)*2,2))]
colnames(tax.sum.df)=colnames(tax.sep)
Bacteria = tax.sum.df[grep(pattern="d__Bacteria",tax.sum.df$kingdom),]
Fungi = tax.sum.df[grep(pattern="d__Eukaryota;k__Fungi",tax.sum.df$kingdom),]
Other_eukaryota = rbind.data.frame(tax.sum.df[grep(pattern="d__Eukaryota;k__norank",tax.sum.df$kingdom),],
                                   tax.sum.df[grep(pattern="d__Eukaryota;k__Metazoa",tax.sum.df$kingdom),],
                                   tax.sum.df[grep(pattern="d__Eukaryota;k__Viridiplantae",tax.sum.df$kingdom),])
Archaea = tax.sum.df[grep(pattern="d__Archaea",tax.sum.df$kingdom),]
Viruses = tax.sum.df[grep(pattern="d__Viruses",tax.sum.df$kingdom),]
bacteria.5 = filter_by_rank(Bacteria,5)
Bacteria$kingdom[ Bacteria$kingdom %nin% bacteria.5$kingdom ] <- "Other bacteria"
Fungi$kingdom<- "Fungi"
Other_eukaryota$kingdom<- "Other eukaryota"
Archaea$kingdom<- "Archaea"
Viruses$kingdom<- "Viruses"

phylum.combine = rbind.data.frame(Bacteria,Fungi,Other_eukaryota,Archaea,Viruses)
phylum.sum = apply(phylum.combine[,-1],2,aggregate,list(phylum.combine$kingdom),sum)
phylum.sum.df = data.frame(phylum.sum)[,c(1,seq(2,length(phylum.sum)*2,2))]
colnames(phylum.sum.df)=colnames(phylum.combine)

percent.phylum = apply(phylum.sum.df[,-1],2,format_percent)
percent.phylum2 = data.frame(Type =phylum.sum.df$kingdom,percent.phylum)
percent.phylum.gather =gather(percent.phylum2,-1,key="sample",value="value") %>% 
  separate(sample,into = c("location","field","group"),sep = "\\.") %>% 
  unite("field", location, field, sep = "|", remove = FALSE)
percent.phylum.sum = summarySE(data=percent.phylum.gather, measurevar="value",groupvars=c("Type","group")) %>% 
  mutate(x="x")

percent.phylum.sum.plot=ggplot(percent.phylum.sum, aes(x = x, y = value*100,fill=Type))+   
  geom_col(width = 0.3,color="black",position = "stack")+  
  labs(x="", y="Proportion (%)")+  
  scale_fill_grafify(palette = "kelly")+
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))+
  facet_wrap(~group,ncol=1)

ggboxplot(percent.phylum.gather,x="group",y="value",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")+
  facet_wrap(~Type,scales = "free",ncol=3)


###########  mag  ###########
library(treeio)

mag.arg = mag.arg%>% 
  separate(GeneID,into = c("strain","gene"),sep = "\\|") %>% 
  mutate(count=1) 
mag.arg$Type[ mag.arg$Type %nin% c("multidrug","rifamycin","vancomycin") ] <- "Others"
mag.arg.num = stats::aggregate(count ~ strain + Type , data = mag.arg, FUN = base::sum) %>% 
  spread(key = "Type",value = "count")
mag.arg.num[is.na(mag.arg.num)]<- 0 

mag.mge = mag.mge %>% 
  separate(GeneID,into = c("strain","gene"),sep = "\\|") %>% 
  mutate(count=1)
mag.mge.num = stats::aggregate(count ~ strain + Type , data = mag.mge, FUN = base::sum) %>% 
  spread(key = "Type",value = "count")
mag.mge.num[is.na(mag.mge.num)]<- 0 
mag.all.num = merge(mag.arg.num,mag.mge.num,1,all=T)
mag.all.num[is.na(mag.all.num)]<- 0 
mag.all.num$strain = gsub(pattern=".Bin", replacement=".", mag.all.num$strain)
mag.tax = read.table("/home/shared/project/2025/ecotone90/11.binning/gnm_taxon.txt", sep = "\t",head = T,quote = "",check.names = F) 
mag.tax.num = merge(mag.tax[-166,c(1,2)],mag.all.num[-166,],1) %>% 
  separate(classification,into = c("domain","phylum","class","order","family","genus","species"),sep = ";",remove = F)


reserve_tip<-mag.tax.num$user_genome
to_drop<-tree$tip.label[-match(reserve_tip,tree$tip.label)]
tree_reduced<-treeio::drop.tip(tree,to_drop)
ggtree(tree_reduced)+
  geom_text2(aes(label=node))+
  geom_tiplab()

mag.tpm.gather = gather(mag.tpm,-1,key="sample",value=value) %>% 
  separate(MAGs,into = c("sample0","id"),sep = ".Bin",remove = F) %>% 
  subset(sample0 == sample)
mag.tpm.gather$MAGs = gsub(pattern=".Bin", replacement=".", mag.tpm.gather$MAGs)
mag.tpm.arg = merge(mag.tpm.gather,mag.tax.num,1)[,-c(2:3)]
mag.tpm.sep = separate(mag.tpm.arg,sample,into  = c("location","field","group","rep"),sep = "-",remove = F) 

mag.tpm.meta=merge(mag.tpm.sep,meta,by="sample") %>% 
  mutate(plant.y=plant) %>% 
  mutate(sum=multidrug+Others+rifamycin+vancomycin+insertion_sequence+integrase+ist+transposase) %>% 
  mutate(arg=multidrug+Others+rifamycin+vancomycin) %>% 
  mutate(mge=insertion_sequence+integrase+ist+transposase) 


mag.tpm.meta.out=mag.tpm.meta[,c(1,2,36,37,3,25,5,27,41,24,7:23,42:44)]

cols <- c("FA"="#0073C2","GR"="#EFC000")

ggscatter(mag.tpm.meta.out, x = "arg", y = "mge", add = "reg.line",cor.coef = F,color = "group.y") +
  stat_poly_eq(formula = y ~ x,aes(color = group.y,label=paste( ..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position=c(0.9,0.9))+guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))

ggplot(mag.tpm.meta.out, aes(x=arg, y=mge)) +
  geom_point(aes(size = value,color = phylum)) +
  theme_ykm()+
  xlab("Time (weeks)")+
  ylab("Sample")+
  scale_colour_nejm()+scale_fill_nejm()

mag.tpm.meta.arg = subset(mag.tpm.meta.out,arg > 0) %>% 
  mutate(count=1)
ggboxplot(mag.tpm.meta.arg,x="family",y="arg",color = "group.y",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(aes(group = group.y),method = "wilcox.test",label = "p.signif")+
  theme_ykm()+
  scale_colour_manual(values = cols)+
  theme(aspect.ratio = 0.5)
mag.arg.count = aggregate(count ~ group.y + phylum+class+family, data = mag.tpm.meta.arg, FUN = base::sum)




mag.tpm.meta.mge = subset(mag.tpm.meta.out,mge > 0) %>% 
  mutate(count=1)
ggboxplot(mag.tpm.meta.mge,x="family",y="mge",color = "group.y",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(aes(group = group.y),method = "wilcox.test",label = "p.signif")+
  theme_ykm()+
  scale_colour_manual(values = cols)+
  theme(aspect.ratio = 0.5)
mag.mge.count = aggregate(count ~ group.y + phylum+class+family, data = mag.tpm.meta.mge, FUN = base::sum)


mag.arg.mge = unique.data.frame(merge(mag.arg[,c(1,2,6,7)],mag.mge[,c(1,2,7,8)],1)) %>% 
  separate(strain,into = c("location","field","group","rep"),sep = "-",remove = F)
mag.arg.mge$strain = gsub(pattern=".Bin", replacement=".", mag.arg.mge$strain)
mag.arg.mge.tax = merge(mag.arg.mge,mag.tax[,1:2],1)

mag.arg.mge.tax.fa = subset(mag.arg.mge.tax,group=="FA")
mag.arg.mge.tax.gr = subset(mag.arg.mge.tax,group=="GR")
mag.arg.mge.tax.select = mag.arg.mge.tax[c(63,66),]





##### chord diagram
library(glue)
library(circlize)

mag.arg = mag.arg %>% 
  separate(GeneID,into = c("strain","gene"),sep = "\\|") %>% 
  mutate(count=1) 
mag.arg.num = stats::aggregate(count ~ strain + Type , data = mag.arg, FUN = base::sum) %>% 
  spread(key = "Type",value = "count")
mag.arg.num[is.na(mag.arg.num)]<- 0 

mag.mge = mag.mge%>% 
  separate(GeneID,into = c("strain","gene"),sep = "\\|") %>% 
  mutate(count=1)
mag.mge$Type[ mag.mge$Type %nin% c("insertion_sequence","transposase") ] <- "Other MGEs"
mag.mge.num = stats::aggregate(count ~ strain + Type , data = mag.mge, FUN = base::sum) %>% 
  spread(key = "Type",value = "count")
mag.mge.num[is.na(mag.mge.num)]<- 0 
mag.all.num = merge(mag.arg.num,mag.mge.num,1,all=T)
mag.all.num[is.na(mag.all.num)]<- 0 
mag.all.num$strain = gsub(pattern=".Bin", replacement=".", mag.all.num$strain)
mag.tax.num = merge(mag.tax[-166,c(1,2)],mag.all.num[-166,],1) %>% 
  separate(classification,into = c("domain","phylum","class","order","family","genus","species"),sep = ";",remove = F) %>% 
  separate(user_genome,into  = c("location","field","group","rep"),sep = "-",remove = F) 
mag.tax.num$phylum[ mag.tax.num$phylum %nin% c("p__Pseudomonadota","p__Actinomycetota","p__Acidobacteriota","p__Myxococcota") ] <- "Other Phylum"
mag.tax.num[is.na(mag.tax.num)]<- 0 

mag.fa = subset(mag.tax.num,group=="FA")
mag.fa.select = apply(mag.fa[,14:22],2,aggregate,list(mag.fa$phylum),sum)
mag.fa.df = data.frame(mag.fa.select)[,c(1,seq(2,length(mag.fa.select)*2,2))]
colnames(mag.fa.df)=colnames(mag.fa)[c(8,14:22)]
mag.fa.df2 =as.matrix(dataframe2matrix(mag.fa.df))
mag.fa.df3 = data.frame(from = rep(rownames(mag.fa.df2), times = ncol(mag.fa.df2)), 
                        to = rep(colnames(mag.fa.df2), each = nrow(mag.fa.df2)),
                        value = as.vector(mag.fa.df2), 
                        stringsAsFactors = FALSE)
opar<-par(no.readonly=TRUE)
par(mai=c(1.8,1.8,1.8,1.8),cex=1.2)

col_fun <- function(value) { 
  sapply(seq_along(value), function(i) {  
    base_color <- grid.col[mag.fa.df3$from[i]]
    adjustcolor(base_color, alpha.f = ifelse(value[i] <1000,0.5, attr(base_color,"alpha") %||%1)) 
  })}

chordDiagram( mag.fa.df3, big.gap =30, col = col_fun(mag.fa.df3$value),annotationTrack = c("grid"), grid.col = grid.col,   
              preAllocateTracks =0.1, directional =1, direction.type = c("diffHeight","arrows"),link.arr.type ="big.arrow")

circos.trackPlotRegion(track.index =1, panel.fun= function(x, y) { 
  xlim = get.cell.meta.data("xlim") 
  ylim = get.cell.meta.data("ylim") 
  sector.name =get.cell.meta.data("sector.index") 
  circos.text(  x = mean(xlim), y =1.1, labels = paste0(' ', sector.name,' '), 
                facing ="clockwise", niceFacing = TRUE, adj = c(0,0.5),  cex =.75, xpd = TRUE)},  bg.border = NA)


mag.gr = subset(mag.tax.num,group=="GR")
mag.gr.select = apply(mag.gr[,14:22],2,aggregate,list(mag.gr$phylum),sum)
mag.gr.df = data.frame(mag.gr.select)[,c(1,seq(2,length(mag.gr.select)*2,2))]
colnames(mag.gr.df)=colnames(mag.gr)[c(8,14:22)]
mag.gr.df2 =as.matrix(dataframe2matrix(mag.gr.df))
mag.gr.df3 = data.frame(from = rep(rownames(mag.gr.df2), times = ncol(mag.gr.df2)), 
                        to = rep(colnames(mag.gr.df2), each = nrow(mag.gr.df2)),
                        value = as.vector(mag.gr.df2), 
                        stringsAsFactors = FALSE)
opar<-par(no.readonly=TRUE)
par(mai=c(1.8,1.8,1.8,1.8),cex=1.2)
col_fun <- function(value) { 
  sapply(seq_along(value), function(i) {  
    base_color <- grid.col[mag.gr.df3$from[i]]
    adjustcolor(base_color, alpha.f = ifelse(value[i] <1000,0.5, attr(base_color,"alpha") %||%1)) 
  })}

chordDiagram( mag.gr.df3, big.gap =30, col = col_fun(mag.gr.df3$value),annotationTrack = c("grid"), grid.col = grid.col,   
              preAllocateTracks =0.1, directional =1, direction.type = c("diffHeight","arrows"),link.arr.type ="big.arrow")

circos.trackPlotRegion(track.index =1, panel.fun= function(x, y) { 
  xlim = get.cell.meta.data("xlim") 
  ylim = get.cell.meta.data("ylim") 
  sector.name =get.cell.meta.data("sector.index") 
  circos.text(  x = mean(xlim), y =1.1, labels = paste0(' ', sector.name,' '), 
                facing ="clockwise", niceFacing = TRUE, adj = c(0,0.5),  cex =.75, xpd = TRUE)},  bg.border = NA)







##### net
mag.tax.num = merge(mag.tax[-166,c(1,2)],mag.all.num[-166,],1) %>% 
  separate(classification,into = c("domain","phylum","class","order","family","genus","species"),sep = ";",remove = F)

mag.sarg.select = separate(mag.sarg[,c(1,6)],GeneID,into = c("sample","id"),sep = ".Bin",remove = F) %>% 
  merge(meta,by="sample") 
mag.sarg.fa = subset(mag.sarg.select,group=="FA")[,c(2,4)]
mag.sarg.gr = subset(mag.sarg.select,group=="GR")[,c(2,4)]

mag.mge.select = separate(mag.mge[,c(1,7)],GeneID,into = c("sample","id"),sep = ".Bin",remove = F) %>% 
  merge(meta,by="sample") 
mag.mge.fa = subset(mag.mge.select,group=="FA")[,c(2,4)]
mag.mge.gr = subset(mag.mge.select,group=="GR")[,c(2,4)]

edges.sarg.fa <- data.frame()
for (DBid in unique(mag.sarg.fa$DBid)) {
  # 获取具有相同 Subtype 的 GeneID
  same_subtype_genes <- mag.sarg.fa[mag.sarg.fa$DBid == DBid, "GeneID"]
  if (length(same_subtype_genes) > 1) {
    combinations <- t(combn(same_subtype_genes, 2))
    temp_edges <- data.frame(
      from = combinations[, 1],
      to = combinations[, 2],
      Subtype = DBid,
      group="arg"
    )
    edges.sarg.fa <- rbind(edges.sarg.fa, temp_edges)
  }
}
edges.sarg.gr <- data.frame()
for (DBid in unique(mag.sarg.gr$DBid)) {
  # 获取具有相同 Subtype 的 GeneID
  same_subtype_genes <- mag.sarg.gr[mag.sarg.gr$DBid == DBid, "GeneID"]
  if (length(same_subtype_genes) > 1) {
    combinations <- t(combn(same_subtype_genes, 2))
    temp_edges <- data.frame(
      from = combinations[, 1],
      to = combinations[, 2],
      Subtype = DBid,
      group="arg"
    )
    edges.sarg.gr <- rbind(edges.sarg.gr, temp_edges)
  }
}

edges.mge.fa <- data.frame()
for (ID in unique(mag.mge.fa$ID)) {
  # 获取具有相同 Subtype 的 GeneID
  same_subtype_genes <- mag.mge.fa[mag.mge.fa$ID == ID, "GeneID"]
  if (length(same_subtype_genes) > 1) {
    combinations <- t(combn(same_subtype_genes, 2))
    temp_edges <- data.frame(
      from = combinations[, 1],
      to = combinations[, 2],
      Subtype = ID,
      group="mge"
    )
    edges.mge.fa <- rbind(edges.mge.fa, temp_edges)
  }
}
edges.mge.gr <- data.frame()
for (ID in unique(mag.mge.gr$ID)) {
  # 获取具有相同 Subtype 的 GeneID
  same_subtype_genes <- mag.mge.gr[mag.mge.gr$ID == ID, "GeneID"]
  if (length(same_subtype_genes) > 1) {
    combinations <- t(combn(same_subtype_genes, 2))
    temp_edges <- data.frame(
      from = combinations[, 1],
      to = combinations[, 2],
      Subtype = ID,
      group="mge"
    )
    edges.mge.gr <- rbind(edges.mge.gr, temp_edges)
  }
}

edges.all.fa = rbind.data.frame(edges.sarg.fa,edges.mge.fa)
edges.all.gr = rbind.data.frame(edges.sarg.gr,edges.mge.gr)
node.all.fa = data.frame(ID = unique(c(edges.all.fa$from,edges.all.fa$to))) 
node.all.fa = mutate(node.all.fa,user_genome = gsub(pattern=".Bin", replacement=".", node.all.fa$ID)) %>% 
  merge(mag.tax.num,"user_genome") %>% 
  separate(classification,into = c("domain2genus","species"),sep = ";s__",remove = F)
node.all.gr = data.frame(ID = unique(c(edges.all.gr$from,edges.all.gr$to))) 
node.all.gr = mutate(node.all.gr,user_genome = gsub(pattern=".Bin", replacement=".", node.all.gr$ID)) %>% 
  merge(mag.tax.num,"user_genome") %>% 
  separate(classification,into = c("domain2genus","species"),sep = ";s__",remove = F)
edges.fa.tax = merge(edges.all.fa,node.all.fa[,c(2,7,4)],by.x = 1,by.y = 1) %>% 
  merge(node.all.fa[,c(2,7,4)],by.x = 2,by.y = 1) 
edges.fa.tax$class <- case_when(  
  edges.fa.tax$domain2genus.x == edges.fa.tax$domain2genus.y ~ 'intragenus'   ,  
  edges.fa.tax$phylum.x != edges.fa.tax$phylum.y ~ 'interphylum',
  .default = "intergenus")
edges.fa.tax$Weight <- case_when(  
  edges.fa.tax$domain2genus.x == edges.fa.tax$domain2genus.y ~ '1'   ,  
  edges.fa.tax$phylum.x != edges.fa.tax$phylum.y ~ '5',
  .default = "3")
edges.gr.tax = merge(edges.all.gr,node.all.gr[,c(2,7,4)],by.x = 1,by.y = 1) %>% 
  merge(node.all.gr[,c(2,7,4)],by.x = 2,by.y = 1) 
edges.gr.tax$class <- case_when(  
  edges.gr.tax$domain2genus.x == edges.gr.tax$domain2genus.y ~ 'intragenus'   ,  
  edges.gr.tax$phylum.x != edges.gr.tax$phylum.y ~ 'interphylum',
  .default = "intergenus")
edges.gr.tax$Weight <- case_when(  
  edges.gr.tax$domain2genus.x == edges.gr.tax$domain2genus.y ~ '1'   ,  
  edges.gr.tax$phylum.x != edges.gr.tax$phylum.y ~ '5',
  .default = "3")
fa.noself=edges.fa.tax$from != edges.fa.tax$to
edges.fa.tax2 = edges.fa.tax[fa.noself,]
gr.noself=edges.gr.tax$from != edges.gr.tax$to
edges.gr.tax2 = edges.gr.tax[gr.noself,]

colnames(edges.fa.tax2)[1:2]=c("source","target")
colnames(edges.gr.tax2)[1:2]=c("source","target")

edge.type= data.frame(class=c("intergenus" ,"interphylum"  ,"intragenus"  ,"intergenus" ,"interphylum"  ,"intragenus"),
                      count=c(table(edges.fa.tax2$class),table(edges.gr.tax2$class)),
                      group=c(rep("FA",3),rep("GR",3)),
                      node = c(rep(nrow(node.all.fa),3),rep(nrow(node.all.gr),3)))
colnames(mag.mge)[c(1,7:9)] <- colnames(mag.sarg)[c(1,6:8)]
arg_mge.anno = unique.data.frame(rbind.data.frame(mag.sarg[,c(6:8)],mag.mge[,c(7:9)]))

colnames(edges.fa.tax2)[3]="DBid"
edges.fa.tax3 = merge(edges.fa.tax2[,3:9],arg_mge.anno,1,all.x=T)
interphylum.fa = subset(edges.fa.tax3,class=="interphylum") %>% 
  unite("name", Type, Subtype, sep = "-", remove = FALSE) %>% 
  mutate(count=1)
interphylum.fa2 = as.data.frame(table_order(interphylum.fa$name)) %>% 
  separate(x,into = c("type","subtype"),sep = "-",remove = F)
intergenus.fa = subset(edges.fa.tax3,class=="intergenus") %>% 
  unite("name", Type, Subtype, sep = "-", remove = FALSE) %>% 
  mutate(count=1)
intergenus.fa2 = as.data.frame(table_order(intergenus.fa$name)) %>% 
  separate(x,into = c("type","subtype"),sep = "-",remove = F)

interphylum.fa.plot=ggplot(interphylum.fa2,aes(x=subtype,y=Freq,fill=type))+
  geom_bar(position =position_dodge(0.6),width = 0.5,stat = "identity",color="black")+
  xlab("")+
  ylab("")+
  theme_ykm()+theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev(interphylum.fa2$subtype))+
  scale_fill_manual(values = grid.col)+
  ylim(0,25)


colnames(edges.gr.tax2)[3]="DBid"
edges.gr.tax3 = merge(edges.gr.tax2[,3:9],arg_mge.anno,1,all.x=T)
interphylum.gr = subset(edges.gr.tax3,class=="interphylum") %>% 
  unite("name", Type, Subtype, sep = "-", remove = FALSE) %>% 
  mutate(count=1)
interphylum.gr2 = as.data.frame(table_order(interphylum.gr$name)) %>% 
  separate(x,into = c("type","subtype"),sep = "-",remove = F)
intergenus.gr = subset(edges.gr.tax3,class=="intergenus") %>% 
  unite("name", Type, Subtype, sep = "-", remove = FALSE) %>% 
  mutate(count=1)
intergenus.gr2 = as.data.frame(table_order(intergenus.gr$name)) %>% 
  separate(x,into = c("type","subtype"),sep = "-",remove = F)

interphylum.gr.plot=ggplot(interphylum.gr2,aes(x=subtype,y=Freq,fill=type))+
  geom_bar(position =position_dodge(0.6),width = 0.5,stat = "identity",color="black")+
  xlab("")+
  ylab("")+
  theme_ykm()+theme(aspect.ratio = 1)+
  scale_x_discrete(limits=rev(interphylum.gr2$subtype))+
  scale_fill_manual(values = grid.col)+
  ylim(0,25)

interphylum.fa.plot | interphylum.gr.plot



###########  protest  ###########
matrix2dataframe <- function(mat){
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  row <- rep(rownames(mat),ncol)
  col <- rep(colnames(mat), each=nrow)
  frame <- data.frame(row,col,value =as.numeric(mat))
  return(frame)
}

resistome.data = unique.data.frame(rbind.data.frame(arg.gene[,1:91],mge.gene[,1:91]))
s.dist <- vegdist(t(dataframe2matrix(tax.data)),method = "bray")
r.dist <- vegdist(t(dataframe2matrix(resistome.data)),method = "bray")
mantel(s.dist,r.dist,method = "pearson")
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)
Y <- cbind(data.frame(pro.s.r$X), data.frame(pro.s.r$Yrot))
X <- data.frame(pro.s.r$rotation)
Y$sample <- rownames(Y)
pro.all = merge(Y,meta,by="sample")
plot.protest=ggplot(pro.all) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "gray90", size = 0.5) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "gray90", size = 0.5) +
  geom_point(aes(MDS1, MDS2, fill = group), size = 2, shape = 21) +
  geom_point(aes(X1, X2, fill = group), size = 2, shape = 22) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  annotate('text', label = 'Mantel test: R2 = 0.83***\nProcrustes correlation: 0.76***',
           x = 0.5, y = 1.2, size = 4,hjust = 0) +
  theme_ykm()+scale_fill_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position = "none") 

tax.fa = tax.data[,colnames(tax.data) %in% c("Tax",subset(meta,group=="FA")[,1])]
resistome.fa = resistome.data[,colnames(resistome.data) %in% c("GeneID",subset(meta,group=="FA")[,1])]
s.dist <- vegdist(t(dataframe2matrix(tax.fa)),method = "bray")
r.dist <- vegdist(t(dataframe2matrix(resistome.fa)),method = "bray")
mantel(s.dist,r.dist,method = "pearson")
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)
Y <- cbind(data.frame(pro.s.r$X), data.frame(pro.s.r$Yrot))
X <- data.frame(pro.s.r$rotation)
Y$sample <- rownames(Y)
pro.fa = merge(Y,meta,by="sample")
plot.fa=ggplot(pro.fa) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#bcbcbc", size = 0.5) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#bcbcbc", size = 0.5) +
  geom_point(aes(MDS1, MDS2, fill = group), size = 2, shape = 21) +
  geom_point(aes(X1, X2, fill = group), size = 2, shape = 22) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Chahar") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  annotate('text', label = 'Mantel test: R2 = 0.91***\nProcrustes correlation: 0.83***',
           x = -1, y = 1.2, size = 4,hjust = 0) +
  theme_ykm()+scale_fill_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position = "none") 

tax.gr = tax.data[,colnames(tax.data) %in% c("Tax",subset(meta,group=="GR")[,1])]
resistome.gr = resistome.data[,colnames(resistome.data) %in% c("GeneID",subset(meta,group=="GR")[,1])]
s.dist <- vegdist(t(dataframe2matrix(tax.gr)),method = "bray")
r.dist <- vegdist(t(dataframe2matrix(resistome.gr)),method = "bray")
mantel(s.dist,r.dist,method = "pearson")
mds.s <- monoMDS(s.dist)
mds.r <- monoMDS(r.dist)
pro.s.r <- procrustes(mds.s,mds.r)
protest(mds.s,mds.r)
Y <- cbind(data.frame(pro.s.r$X), data.frame(pro.s.r$Yrot))
X <- data.frame(pro.s.r$rotation)
Y$sample <- rownames(Y)
pro.gr = merge(Y,meta,by="sample")
plot.gr=ggplot(pro.gr) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#bcbcbc", size = 0.5) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#bcbcbc", size = 0.5) +
  geom_point(aes(MDS1, MDS2, fill = group), size = 2, shape = 21) +
  geom_point(aes(X1, X2, fill = group), size = 2, shape = 22) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Chahar") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  annotate('text', label = 'Mantel test: R2 = 0.85***\nProcrustes correlation: 0.80***',
           x = -1, y = 1.2, size = 4,hjust = 0) +
  theme_ykm()+scale_fill_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position = "none") 

distance.tax = as.matrix(s.dist)
distance.tax.data = matrix2dataframe (distance.tax) 
distance.resis = as.matrix(r.dist)
distance.resis.data = matrix2dataframe (distance.resis) 

distance.tax.sep = separate(distance.tax.data,row,c("location.x","field.x","group.x","rep.x"),sep = "-",remove = F) %>% 
  separate(col,c("location.y","field.y","group.y","rep.y"),sep = "-",remove = F)


###########  linkET  ###########
library(linkET)
tax.fa = tax.resistome.data[,colnames(tax.resistome.data) %in% c("GeneID",subset(meta,group=="FA")[,1])]
tax.gr = tax.resistome.data[,colnames(tax.resistome.data) %in% c("GeneID",subset(meta,group=="GR")[,1])]
env.fa = subset(meta,group=="FA")[,c(1,18:20,8:15,6:7)]
env.gr = subset(meta,group=="GR")[,c(1,18:20,8:15,6:7)]

mantel.fa <- mantel_test(spec = t(dataframe2matrix(tax.fa)), env = dataframe2matrix(env.fa), spec_select = list(Microbiome = 1:64634, Resistome = 64635:68517), mantel_fun = 'mantel') %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c('< 0.2', '0.2 - 0.4', '>= 0.4')),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c('< 0.01', '0.01 - 0.05', '>= 0.05')))
linket.fa=qcorrplot(correlate(dataframe2matrix(env.fa), method = 'spearman'), type = 'upper', diag = FALSE) +  #环境变量矩阵计算 Spearman 相关系数
  geom_square() +  #绘制 Spearman 相关系数热图
  geom_mark(sep = '\n', size = 2.5, sig.thres = 0.05,only_mark=T) +  #显示 Spearman 相关系数和显著性
  geom_couple(aes(color = pd, size = rd), data = mantel.fa, curvature = nice_curvature()) +  #环境和微生物的相关性展示为上述 Mantel 相关
  scale_fill_gradientn(colors = c('#67001F', '#F7B394', 'white','#68A8CF', '#053061' ), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  scale_size_manual(values = c(0, 1, 2)) +  #根据 Mantel 相关指定线条粗细
  scale_color_manual(values = c('#D95F02', '#1B9E77', 'white')) +  #根据 Mantel 相关 p 值指定线条颜色
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +
  theme(legend.key = element_blank())

mantel.gr <- mantel_test(spec = t(dataframe2matrix(tax.gr)), env = dataframe2matrix(env.gr), spec_select = list(Microbiome = 1:64634, Resistome = 64635:68517), mantel_fun = 'mantel') %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c('< 0.2', '0.2 - 0.4', '>= 0.4')),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c('< 0.01', '0.01 - 0.05', '>= 0.05')))
#write.table (mantel.gr, file = "/home/shared/project/2025/ecotone90/analysis/result/mantel.gr.txt", sep = "\t", row.names = F ,col.names = TRUE, quote = FALSE)
linket.gr=qcorrplot(correlate(dataframe2matrix(env.gr), method = 'spearman'), type = 'upper', diag = FALSE) +  #环境变量矩阵计算 Spearman 相关系数
  geom_square() +  #绘制 Spearman 相关系数热图
  geom_mark(sep = '\n', size = 2.5, sig.thres = 0.05,only_mark=T) +  #显示 Spearman 相关系数和显著性
  geom_couple(aes(color = pd, size = rd), data = mantel.gr, curvature = nice_curvature()) +  #环境和微生物的相关性展示为上述 Mantel 相关
  scale_fill_gradientn(colors = c('#67001F', '#F7B394', 'white','#68A8CF', '#053061'), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  scale_size_manual(values = c(0, 1, 2)) +  #根据 Mantel 相关指定线条粗细
  scale_color_manual(values = c('#D95F02', '#1B9E77', 'white')) +  #根据 Mantel 相关 p 值指定线条颜色
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +
  theme(legend.key = element_blank())

linket.fa | linket.gr

###########  NST  ###########
library(ape)
library(NST) 


tax.CHA1 = t(dataframe2matrix(tax.data[,colnames(tax.data) %in% c("Tax",subset(meta,field=="CHA|1")[,1])]))
group.CHA1 = data.frame(group=subset(meta,field=="CHA|1")[,3])
rownames(group.CHA1)=rownames(tax.CHA1)
tNST.CHA1 <- tNST(comm=tax.CHA1, group=group.CHA1,
                  dist.method="bray", abundance.weighted=TRUE, rand=999,
                  output.rand=TRUE, nworker=16, LB=FALSE, null.model="PF",
                  between.group=FALSE, SES=TRUE, RC=TRUE)
tNST.CHA1.data=tNST.CHA1[[3]] %>% 
  gather(2:3,key="name",value ="sample") %>% 
  summarySE( measurevar="NST.ij.bray",groupvars=c("group","sample")) %>% 
  mutate(field="CHA|1")
tNST.CHA1.boot <- nst.boot(nst.result=tNST.CHA1,group=group.CHA1,rand=999,trace=TRUE,
                           two.tail=T,out.detail=F,between.group=T,nworker=16)
tNST.CHA1.bootdata=tNST.CHA1.boot[[2]] %>% 
  mutate(field="CHA|1")
tNST.CHA1.aov <- nst.panova(nst.result=tNST.CHA1,group = group.CHA1, rand = 999, trace = TRUE)%>% 
  mutate(field="CHA|1")

tNST.data = tNST.CHA1.data[1,]
tNST.boot = tNST.CHA1.bootdata[1,]
tNST.aov = tNST.CHA1.aov[1,]

for (i in unique(meta$field)) {
  tax.CHA1 = t(dataframe2matrix(tax.data[,colnames(tax.data) %in% c("Tax",subset(meta,field==i)[,1])]))
  group.CHA1 = data.frame(group=subset(meta,field==i)[,3])
  rownames(group.CHA1)=rownames(tax.CHA1)
  tNST.CHA1 <- tNST(comm=tax.CHA1, group=group.CHA1,
                    dist.method="bray", abundance.weighted=TRUE, rand=999,
                    output.rand=TRUE, nworker=64, LB=FALSE, null.model="PF",
                    between.group=T, SES=TRUE, RC=TRUE)
  tNST.CHA1.data=tNST.CHA1[[3]] %>% 
    gather(2:3,key="name",value ="sample") %>% 
    summarySE( measurevar="NST.ij.bray",groupvars=c("group","sample")) %>% 
    mutate(field=i)
  tNST.CHA1.boot <- nst.boot(nst.result=tNST.CHA1,group=group.CHA1,rand=999,trace=TRUE,
                             two.tail=T,out.detail=F,between.group=T,nworker=16)
  tNST.CHA1.bootdata=tNST.CHA1.boot[[2]] %>% 
    mutate(field=i)
  tNST.CHA1.aov <- nst.panova(nst.result=tNST.CHA1,group = group.CHA1, rand = 999, trace = TRUE)%>% 
    mutate(field=i)
  tNST.data = rbind.data.frame(tNST.data,tNST.CHA1.data)
  tNST.boot = rbind.data.frame(tNST.boot,tNST.CHA1.bootdata)
  tNST.aov = rbind.data.frame(tNST.aov,tNST.CHA1.aov)
}
tNST.data = tNST.data[-1,]
tNST.boot = tNST.boot[-1,]
tNST.aov = tNST.aov[-1,]
tNST.aov.p =subset(tNST.aov,index=="NST"& group1=="FA"&P.anova < 0.05)

tNST.summary = summarySE(data=tNST.data, measurevar="NST.ij.bray",groupvars=c("group","field")) %>% 
  mutate(deterministic=1-NST.ij.bray)
colnames(tNST.summary)[4] = "stochastic"
tNST.summary.gather = gather(tNST.summary[,c(1,2,4,8)],3:4,key="process",value ="value") %>% 
  unite("sample", field, group, sep = "-", remove = FALSE) %>% 
  separate(field,into = c("city","location"),sep = "\\|", remove = FALSE)

ggplot(tNST.summary.gather, aes( x = sample,y=value,fill = process, stratum = process, alluvium = process))+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=0,xmax=1.5,fill="#0073C2",alpha=0.1 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=1.5,xmax=2.5,fill="#EFC000",alpha=0.1 )+
  geom_stratum(width = 0.5, color='black')+
  geom_alluvium(alpha = 0.5,width = 0.5,curve_type = "linear") +
  scale_fill_manual(values = c("stochastic"="#FEDFB2","deterministic"="#EE7E77"))+
  facet_wrap(~field, scales = "free_x",nrow=1) +
  theme_ykm()+theme(aspect.ratio = 2)+
  xlab("")+ylab("")


tax.all = t(dataframe2matrix(tax.data))
group.all = data.frame(group=meta$group)
rownames(group.all)=rownames(tax.all)

tNST.all <- tNST(comm=tax.all, group=group.all,
                 dist.method="bray", abundance.weighted=TRUE, rand=999,
                 output.rand=TRUE, nworker=64, LB=FALSE, null.model="PF",
                 between.group=T, SES=TRUE, RC=TRUE)
tNST.all.data=tNST.all[[3]] %>% 
  separate(name1,into = c("location1","field1","group1","rep1"),sep = "-", remove = FALSE)%>% 
  separate(name2,into = c("location2","field2","group2","rep2"),sep = "-", remove = FALSE) %>% 
  unite("sample1", location1, field1, group1,sep = "-", remove = FALSE) %>% 
  unite("sample2", location2, field2, group2,sep = "-", remove = FALSE) %>% 
  subset(sample1 == sample2) %>% 
  gather(c(2,8),key="name",value ="sample") %>% 
  summarySE( measurevar="NST.ij.bray",groupvars=c("group","sample"))
tNST.all.boot <- nst.boot(nst.result=tNST.all,group=group.all,rand=999,trace=TRUE,
                          two.tail=T,out.detail=TRUE,between.group=T,nworker=16)
tNST.boot.data=tNST.all.boot[[2]]
tNST.all.aov <- nst.panova(nst.result=tNST.all,group = group.all, rand = 999, trace = TRUE)


###########  AVD  ###########

ai <- abs(dataframe2matrix(tax.data)-apply(dataframe2matrix(tax.data), 1, mean))/apply(dataframe2matrix(tax.data), 1, sd)
AVD <- apply(ai,2,sum)/(1*nrow(dataframe2matrix(tax.data)))
avd.data=data.frame(sample=names(AVD),avd=AVD)%>% 
  merge(meta,1)
tax.avd.plot=ggplot(avd.data,aes(x=field,y=avd,fill=group))+
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",width = .2,position = position_dodge(0.5)) +
  stat_summary(fun = "mean", geom = "point",color="black",shape = 21,position = position_dodge(0.5),size=2) +
  xlab("")+
  ylab("")+
  theme_ykm()+
  stat_compare_means(aes( group=group),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***","**","*","n.s.")),vjust=0.6)+
  scale_fill_jco()+
  theme(aspect.ratio = 0.5,legend.position = "none")

tax.shannon.plot / tax.avd.plot

tax_avd_group.plot=ggboxplot(avd.data,x="group",y="avd",color = "group",
                             add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")

############ soil compare ####################
meta.soil = gather(meta[,c(3,6:15)],-1,key="soil",value=value)

ggboxplot(meta.soil,x="group",y="value",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")+
  facet_wrap(~soil,ncol=4,scales = "free")

meta.soil.select = subset(meta.soil,soil %in% c("TP(g/100g)","AP(mg/kg)","pH","EC","Water(%)","MBN(mg/kg)"))
meta.soil.select$soil = factor(meta.soil.select$soil, levels=as.factor(c("TP(g/100g)","AP(mg/kg)","pH","EC","Water(%)","MBN(mg/kg)")) )
soil_group.plot=ggboxplot(meta.soil.select,x="group",y="value",color = "group",
                          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")+
  facet_wrap(~soil,ncol=3,scales = "free")

ggboxplot(meta.soil.select,x="group",y="value",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(aspect.ratio = 1.5,legend.position = "none")+
  facet_wrap(~soil,ncol=3,scales = "free")

soil.gather =gather(meta[,-c(16:20)],6:15,key="index",value ="value") 
input = soil.gather[,c(2,3,6,7)]
ttest.result = data.frame(field="y",index="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$index)){
    data.pro = subset(data.sub,index == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,index=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,index=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("soil", field, index, sep = "-", remove = FALSE)
soil.summary <- summarySE(data=soil.gather, measurevar="value",groupvars=c("index","field","group")) %>% 
  unite("soil", field, index, sep = "-", remove = FALSE)
soil.spread = spread(soil.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  merge(ttest.result2,"soil") %>% 
  mutate(updown="none")
soil.spread$updown <- case_when(  
  soil.spread$fc < 0 & soil.spread$p < 0.05 ~ 'grass'   ,  
  soil.spread$fc > 0 & soil.spread$p < 0.05 ~ 'farm',
  .default = "none")

soil.spread2 = spread(soil.spread[,c(2,3,6)],key = "field.x",value = "fc")
soil.clust.plot = hclust(dist((dataframe2matrix(soil.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 1)+
  hexpand(2) 
clust.data = subset(soil.clust.plot$data,is.na(label)==F)
soil.spread$index.x = factor(soil.spread$index.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
soil.plot = ggplot(soil.spread,aes(x=field.x,y=index.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3,5), range = c(1, 8))+
  #theme(legend.position = "none") +
  scale_shape_manual(values = c(24,25,21)) 

soil.clust.plot | soil.plot

###########  data combine  ###########


soil = scale(dataframe2matrix(meta[,c(1,6,7,9:13)]))
Zcore <- as.data.frame((soil-mean(soil))/sd(soil))
var = colnames(Zcore)
Zcore$emf <- apply(select(Zcore,var),1,mean)

sample.data = data.frame(meta,
                         resistome_shannon = resistome.alpha$Shannon,
                         resistome_richness = resistome.alpha$Species_richness,
                         resistome_abundance = resistome.sum$sum,
                         tax_shannon = tax.alpha$Shannon,
                         tax_richness = tax.alpha$Species_richness,
                         resistome_tsne = tsne.resis$tsne1,
                         tax_tsne = tsne.tax$tsne1,
                         stochastic = tNST.data$NST.ij.bray,
                         deterministic = 1-tNST.data$NST.ij.bray,
                         emf = Zcore$emf,
                         avd = avd.data$avd,
                         cogQ = COG.Q$value
) %>% 
  mutate(microbial_CN = `MBC.mg.kg.`/`MBN.mg.kg.`)

sample.data = read.table("/home/shared/project/2025/ecotone90/analysis/result/sample.data.txt", sep = "\t",head = T,quote = "",check.names = F)
arg_shannon.group.ave =  summarySE(data=sample.data, measurevar="arg_shannon",groupvars=c("group"))
mge_shannon.group.ave =  summarySE(data=sample.data, measurevar="mge_shannon",groupvars=c("group"))
resistone_shannon.group.ave =  summarySE(data=sample.data, measurevar="resistome_shannon",groupvars=c("group"))
tax_shannon.group.ave =  summarySE(data=sample.data, measurevar="tax_shannon",groupvars=c("group"))
tax_richness.group.ave =  summarySE(data=sample.data, measurevar="tax_richness",groupvars=c("group"))
microbial_CN.group.ave =  summarySE(data=sample.data, measurevar="microbial_CN",groupvars=c("group"))
avd.group.ave =  summarySE(data=sample.data, measurevar="avd",groupvars=c("group"))
tp.group.ave =  summarySE(data=sample.data, measurevar="TP.g.100g.",groupvars=c("group"))
ap.group.ave =  summarySE(data=sample.data, measurevar="AP.mg.kg.",groupvars=c("group"))
ec.group.ave =  summarySE(data=sample.data, measurevar="EC",groupvars=c("group"))
ph.group.ave =  summarySE(data=sample.data, measurevar="pH",groupvars=c("group"))
water.group.ave =  summarySE(data=sample.data, measurevar="Water...",groupvars=c("group"))
mbn.group.ave =  summarySE(data=sample.data, measurevar="MBN.mg.kg.",groupvars=c("group"))


cor.shannon.plot=ggscatter(sample.data, x = "tax_shannon", y = "resistome_shannon", add = "reg.line",cor.coef = F,color = "group",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position="none")
cor.richness.plot=ggscatter(sample.data, x = "tax_richness", y = "resistome_shannon", add = "reg.line",cor.coef = F,color = "group",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position="none")
cor.avd.plot=ggscatter(sample.data, x = "avd", y = "resistome_shannon", add = "reg.line",cor.coef = F,color = "group",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position="none")
cor.CN.plot=ggscatter(sample.data, x = "microbial_CN", y = "resistome_shannon", add = "reg.line",cor.coef = F,color = "group",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position="none")

(cor.shannon.plot | cor.richness.plot)/(cor.avd.plot|cor.CN.plot)

ggboxplot(sample.data,x="group",y="microbial_CN",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")

cor.long.plot=ggscatter(sample.data, x = "latitude", y = "resistome_shannon", add = "reg.line",cor.coef = F,color = "group",fullrange  =T) +
  stat_poly_eq(formula = y ~ x,aes(color = group,label=paste(..adj.rr.label..,..p.value.label..,sep = "~~~~")), parse = TRUE)+
  theme_ykm()+
  scale_color_manual(values = c("FA"="#0073C2","GR"="#EFC000"))+
  theme(legend.position="none")



###########  sem  ###########
library(plspm)
library(piecewiseSEM)
library(semPlot)
library(qgraph)
library(ggchicklet)

sample.data.fa = subset(sample.data,group=="FA")
sample.data.gr = subset(sample.data,group=="GR")

fix_blocks = list(c(18:20),c(6:7,9,12:13),c(24:25),c(27,29,31),21)
fix_modes = rep("A", 5)

fix_pls.fa = plspm(sample.data.fa, path, fix_blocks, modes = fix_modes,boot.val = TRUE,tol =0.00001)
innerplot(fix_pls.fa, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray20', box.lwd = 0)
outerplot(fix_pls.fa, colpos = 'red', colneg = 'blue')
summary(fix_pls.fa)
fix_pls.fa$inner_model
fix_pls.fa$outer_model
fix_pls.fa$effects
fix_pls.fa$gof
effects.fa=as.data.frame(fix_pls.fa$effects)[c(4,7,9,10),-4] %>% 
  gather(-1,key="type",value="value")

effects.fa.plot=ggplot(effects.fa,aes(x=relationships,y=value))+
  geom_chicklet(aes(fill=type),
                width=0.8,
                radius=grid::unit(3,"pt"))+
  theme_ykm()+
  scale_fill_manual(values = c("#01847F","#F9D2E4"))+
  geom_hline(aes(yintercept=0))+theme(legend.position = "none")+
  ylim(-1,1)

fix_pls.gr = plspm(sample.data.gr, path, fix_blocks, modes = fix_modes,boot.val = TRUE,tol =0.00001)
innerplot(fix_pls.gr, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray20', box.lwd = 0)
summary(fix_pls.gr)
fix_pls.gr$inner_model
fix_pls.gr$outer_model
fix_pls.gr$effects
fix_pls.gr$gof
effects.gr=as.data.frame(fix_pls.gr$effects)[c(4,7,9,10),-4] %>% 
  gather(-1,key="type",value="value")
effects.gr.plot=ggplot(effects.gr,aes(x=relationships,y=value))+
  geom_chicklet(aes(fill=type),
                width=0.8,
                radius=grid::unit(3,"pt"))+
  theme_ykm()+
  scale_fill_manual(values = c("#01847F","#F9D2E4"))+
  geom_hline(aes(yintercept=0))+theme(legend.position = "none")+
  ylim(-1,1)

effects.fa.plot | effects.gr.plot



###########  variable importance  ########### 
library(Boruta)
library(MuMIn)
library(rdacca.hp)


### shannon
shannon.variable = dataframe2matrix(sample.data[,-c(2,4:5,14:15,19:20,22:23,26,30)])
variable.group = data.frame(index = colnames(shannon.variable)[c(2:3,14:21,4:9,10:12)],group=c(rep("Microbial",10),rep("Soil",6),rep("Geography",3))) 

fa.shannon.variable = subset(shannon.variable,group=="FA")[,-1]
fa.shannon.boruta <- Boruta(resistome_shannon ~ ., data=fa.shannon.variable, pValue = 0.05, mcAdj = TRUE, doTrace = 2)
fa.shannon.boruta.importance <- fa.shannon.boruta$ImpHistory
plot(fa.shannon.boruta.importance, las = 2, xlab = '', main = 'Variable Importance')
fa.shannon.boruta.result <- data.frame(index = rownames(attStats(fa.shannon.boruta)),attStats(fa.shannon.boruta))
gr.shannon.variable = subset(shannon.variable,group=="GR")[,-1]
gr.shannon.boruta <- Boruta(resistome_shannon ~ ., data=gr.shannon.variable, pValue = 0.05, mcAdj = TRUE, doTrace = 2)
gr.shannon.boruta.importance <- gr.shannon.boruta$ImpHistory
plot(gr.shannon.boruta.importance, las = 2, xlab = '', main = 'Variable Importance')
gr.shannon.boruta.result <- data.frame(index = rownames(attStats(gr.shannon.boruta)),attStats(gr.shannon.boruta))

fa.shannon.variable.scale<-data.frame(scale(fa.shannon.variable))
fa.resistome_shannon.model<-lm(resistome_shannon~., data=fa.shannon.variable.scale,na.action = "na.fail" )
fa.resistome_shannon.model_all<-dredge(fa.resistome_shannon.model)  
fa.resistome_shannon.model_best<-summary(get.models(fa.resistome_shannon.model_all, 1)[[1]])

fa.resistome_shannon.result<-data.frame(index=rownames(fa.resistome_shannon.model_best$coefficients),fa.resistome_shannon.model_best$coefficients) %>% 
  merge(variable.group,1)
fa.resistome_shannon.effect=rdacca.hp(fa.shannon.variable.scale$resistome_shannon, list(Soil=fa.shannon.variable.scale[,3:8], Microbial=fa.shannon.variable.scale[,c(1:2,13:20)],Geography= fa.shannon.variable.scale[,9:11]), method = 'RDA', type = 'adjR2')
fa.resistome_shannon.effect_data<-data.frame(index = c(rownames(fa.resistome_shannon.effect$Hier.part),"Unexplained"),value=c(fa.resistome_shannon.effect$Hier.part[,3],1-fa.resistome_shannon.effect$Total_explained_variation))%>% 
  mutate(x="Index type")
fa.resistome_shannon.effect_data$index = factor(fa.resistome_shannon.effect_data$index, levels=as.factor(c("Microbial","Soil","Geography","Unexplained") ))

fa.resistome_shannon.result$sig[fa.resistome_shannon.result$Pr...t.. < 0.05]<-"*"
fa.resistome_shannon.result$sig[fa.resistome_shannon.result$Pr...t.. < 0.01]<-"**"
fa.resistome_shannon.result$sig[fa.resistome_shannon.result$Pr...t.. < 0.001]<-"***"

fa.shannon.p1=ggplot(fa.resistome_shannon.result,aes(index,Estimate))+  
  geom_rect(ymin=-Inf,ymax =Inf,xmin=6.5,xmax=Inf,fill="#4DBBD5",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=1.5,xmax=6.5,fill="#00A087",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=-Inf,xmax=2.5,fill="#E64B35",alpha=0.05 )+
  geom_point(size=3,shape=22,color="black",aes(fill=group))+  
  geom_errorbar(aes(ymin = Estimate - `Std..Error`, ymax = Estimate + `Std..Error`,color=group), width = 0,position = position_dodge(width = 0.7),cex=0.9)+  
  labs(x="", y="Parameter estimates")+  
  geom_text(aes(label = sig,y= Estimate, x=index),color="black",vjust=-0.05,size =5)+  
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_npg() +  scale_color_npg() +  
  scale_x_discrete(limits=fa.resistome_shannon.result$index[order(match(fa.resistome_shannon.result$index, variable.group$index),decreasing = T)])+
  theme_ykm()+coord_flip()+theme(legend.position = "none")

fa.shannon.p2=ggplot(fa.resistome_shannon.effect_data, aes(x = x, y = value*100,fill=index))+   
  geom_col(width = 0.6,color="black",position = "stack")+  
  labs(x="", y="Relative effect of esimates (%)")+  
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) +  scale_color_npg() +  
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,legend.position = "none",
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))

gr.shannon.variable.scale<-data.frame(scale(gr.shannon.variable))
gr.resistome_shannon.model<-lm(resistome_shannon~., data=gr.shannon.variable.scale,na.action = "na.fail" )
gr.resistome_shannon.model_all<-dredge(gr.resistome_shannon.model)  ###这一步很慢
gr.resistome_shannon.model_best<-summary(get.models(gr.resistome_shannon.model_all, 1)[[1]])

gr.resistome_shannon.result<-data.frame(index=rownames(gr.resistome_shannon.model_best$coefficients),gr.resistome_shannon.model_best$coefficients) %>% 
  merge(variable.group,1)
gr.resistome_shannon.effect=rdacca.hp(gr.shannon.variable.scale$resistome_shannon, list(Soil=fa.shannon.variable.scale[,3:8], Microbial=fa.shannon.variable.scale[,c(1:2,13:20)],Geography= fa.shannon.variable.scale[,9:11]), method = 'RDA', type = 'adjR2')
gr.resistome_shannon.effect_data<-data.frame(index = c(rownames(gr.resistome_shannon.effect$Hier.part),"Unexplained"),value=c(gr.resistome_shannon.effect$Hier.part[,3],1-gr.resistome_shannon.effect$Total_explained_variation))%>% 
  mutate(x="Index type")
gr.resistome_shannon.effect_data$index = factor(gr.resistome_shannon.effect_data$index, levels=as.factor(c("Microbial","Soil","Geography","Unexplained") ))

gr.resistome_shannon.result$sig[gr.resistome_shannon.result$Pr...t.. < 0.05]<-"*"
gr.resistome_shannon.result$sig[gr.resistome_shannon.result$Pr...t.. < 0.01]<-"**"
gr.resistome_shannon.result$sig[gr.resistome_shannon.result$Pr...t.. < 0.001]<-"***"

gr.shannon.p1=ggplot(gr.resistome_shannon.result,aes(index,Estimate))+  
  geom_rect(ymin=-Inf,ymax =Inf,xmin=3.5,xmax=Inf,fill="#4DBBD5",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=1.5,xmax=3.5,fill="#00A087",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=-Inf,xmax=1.5,fill="#E64B35",alpha=0.05 )+
  geom_point(size=3,shape=22,color="black",aes(fill=group))+  
  geom_errorbar(aes(ymin = Estimate - `Std..Error`, ymax = Estimate + `Std..Error`,color=group), width = 0,position = position_dodge(width = 0.7),cex=0.9)+  
  labs(x="", y="Parameter estimates")+  
  geom_text(aes(label = sig,y= Estimate, x=index),color="black",vjust=-0.05,size =5)+  
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_npg() +  scale_color_npg() +  
  scale_x_discrete(limits=gr.resistome_shannon.result$index[order(match(gr.resistome_shannon.result$index, variable.group$index),decreasing = T)])+
  theme_ykm()+coord_flip()+theme(legend.position = "none")

gr.shannon.p2=ggplot(gr.resistome_shannon.effect_data, aes(x = x, y = value*100,fill=index))+   
  geom_col(width = 0.6,color="black",position = "stack")+  
  labs(x="", y="Relative effect of esimates (%)")+  
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) +  scale_color_npg() +  
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,legend.position = "none",
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))

fa.shannon.p1 | fa.shannon.p2 | gr.shannon.p1 | gr.shannon.p2


### abundance
abundance.variable = dataframe2matrix(sample.data[,-c(4,8,16,17,22,21,26)])
fa.abundance.variable = subset(abundance.variable,group=="FA")[,-2]
fa.abundance.boruta <- Boruta(resistome_abundance ~ ., data=fa.abundance.variable, pValue = 0.05, mcAdj = TRUE, doTrace = 2)
fa.abundance.boruta.importance <- fa.abundance.boruta$ImpHistory
plot(fa.abundance.boruta.importance, las = 2, xlab = '', main = 'Variable Importance')
fa.abundance.boruta.result <- data.frame(index = rownames(attStats(fa.abundance.boruta)),attStats(fa.abundance.boruta))
gr.abundance.variable = subset(abundance.variable,group=="GR")[,-2]
gr.abundance.boruta <- Boruta(resistome_abundance ~ ., data=gr.abundance.variable, pValue = 0.05, mcAdj = TRUE, doTrace = 2)
gr.abundance.boruta.importance <- gr.abundance.boruta$ImpHistory
plot(gr.abundance.boruta.importance, las = 2, xlab = '', main = 'Variable Importance')
gr.abundance.boruta.result <- data.frame(index = rownames(attStats(gr.abundance.boruta)),attStats(gr.abundance.boruta))

fa.abundance.variable.scale<-data.frame(scale(fa.abundance.variable[,-c(1:2,21)]))
fa.resistome_abundance.model<-lm(resistome_abundance~., data=fa.abundance.variable.scale,na.action = "na.fail" )
fa.resistome_abundance.model_all<-dredge(fa.resistome_abundance.model)  ###这一步很慢
fa.resistome_abundance.model_best<-summary(get.models(fa.resistome_abundance.model_all, 1)[[1]])
fa.resistome_abundance.result<-data.frame(index=rownames(fa.resistome_abundance.model_best$coefficients),fa.resistome_abundance.model_best$coefficients) %>% 
  merge(variable.group,1)
fa.resistome_abundance.effect=rdacca.hp(fa.abundance.variable.scale$resistome_abundance, list(Soil=fa.abundance.variable.scale[,1:9], Microbial=fa.abundance.variable.scale[,c(14:20)],Geography= fa.abundance.variable.scale[,10:12]), method = 'RDA', type = 'adjR2')
fa.resistome_abundance.effect_data<-data.frame(index = c(rownames(fa.resistome_abundance.effect$Hier.part),"Unexplained"),value=c(fa.resistome_abundance.effect$Hier.part[,3],1-fa.resistome_abundance.effect$Total_explained_variation))%>% 
  mutate(x="Index type")
fa.resistome_abundance.effect_data$index = factor(fa.resistome_abundance.effect_data$index, levels=as.factor(c("Microbial","Soil","Geography","Unexplained") ))

fa.resistome_abundance.result$sig[fa.resistome_abundance.result$Pr...t.. < 0.05]<-"*"
fa.resistome_abundance.result$sig[fa.resistome_abundance.result$Pr...t.. < 0.01]<-"**"
fa.resistome_abundance.result$sig[fa.resistome_abundance.result$Pr...t.. < 0.001]<-"***"

fa.abundance.p1=ggplot(fa.resistome_abundance.result,aes(index,Estimate))+  
  geom_rect(ymin=-Inf,ymax =Inf,xmin=5.5,xmax=Inf,fill="#4DBBD5",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=-Inf,xmax=5.5,fill="#00A087",alpha=0.05 )+
  #geom_rect(ymin=-Inf,ymax =Inf,xmin=-Inf,xmax=1.5,fill="#E64B35",alpha=0.05 )+
  geom_point(size=3,shape=22,color="black",aes(fill=group))+  
  geom_errorbar(aes(ymin = Estimate - `Std..Error`, ymax = Estimate + `Std..Error`,color=group), width = 0,position = position_dodge(width = 0.7),cex=0.9)+  
  labs(x="", y="Parameter estimates")+  
  geom_text(aes(label = sig,y= Estimate, x=index),color="black",vjust=-0.05,size =5)+  
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) + 
  scale_color_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) + 
  scale_x_discrete(limits=fa.resistome_abundance.result$index[order(match(fa.resistome_abundance.result$index, variable.group$index),decreasing = T)])+
  theme_ykm()+coord_flip()+theme(legend.position = "none")

fa.abundance.p2=ggplot(fa.resistome_abundance.effect_data, aes(x = x, y = value*100,fill=index))+   
  geom_col(width = 0.6,color="black",position = "stack")+  
  labs(x="", y="Relative effect of esimates (%)")+  
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) +  
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,legend.position = "none",
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))

gr.abundance.variable.scale<-data.frame(scale(gr.abundance.variable[,-c(1:2,21)]))
gr.resistome_abundance.model<-lm(resistome_abundance~., data=gr.abundance.variable.scale,na.action = "na.fail" )
gr.resistome_abundance.model_all<-dredge(gr.resistome_abundance.model)  ###这一步很慢
gr.resistome_abundance.model_best<-summary(get.models(gr.resistome_abundance.model_all, 1)[[1]])
gr.resistome_abundance.result<-data.frame(index=rownames(gr.resistome_abundance.model_best$coefficients),gr.resistome_abundance.model_best$coefficients) %>% 
  merge(variable.group,1)
gr.resistome_abundance.effect=rdacca.hp(gr.abundance.variable.scale$resistome_abundance, list(Soil=gr.abundance.variable.scale[,1:9], Microbial=gr.abundance.variable.scale[,c(14:20)],Geography= gr.abundance.variable.scale[,10:12]), method = 'RDA', type = 'adjR2')
gr.resistome_abundance.effect_data<-data.frame(index = c(rownames(gr.resistome_abundance.effect$Hier.part),"Unexplained"),value=c(gr.resistome_abundance.effect$Hier.part[,3],1-gr.resistome_abundance.effect$Total_explained_variation))%>% 
  mutate(x="Index type")
gr.resistome_abundance.effect_data$index = factor(gr.resistome_abundance.effect_data$index, levels=as.factor(c("Microbial","Soil","Geography","Unexplained") ))

gr.resistome_abundance.result$sig[gr.resistome_abundance.result$Pr...t.. < 0.05]<-"*"
gr.resistome_abundance.result$sig[gr.resistome_abundance.result$Pr...t.. < 0.01]<-"**"
gr.resistome_abundance.result$sig[gr.resistome_abundance.result$Pr...t.. < 0.001]<-"***"

gr.abundance.p1=ggplot(gr.resistome_abundance.result,aes(index,Estimate))+  
  geom_rect(ymin=-Inf,ymax =Inf,xmin=5.5,xmax=Inf,fill="#4DBBD5",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=-1.5,xmax=5.5,fill="#00A087",alpha=0.05 )+
  geom_rect(ymin=-Inf,ymax =Inf,xmin=-Inf,xmax=1.5,fill="#E64B35",alpha=0.05 )+
  geom_point(size=3,shape=22,color="black",aes(fill=group))+  
  geom_errorbar(aes(ymin = Estimate - `Std..Error`, ymax = Estimate + `Std..Error`,color=group), width = 0,position = position_dodge(width = 0.7),cex=0.9)+  
  labs(x="", y="Parameter estimates")+  
  geom_text(aes(label = sig,y= Estimate, x=index),color="black",vjust=-0.05,size =5)+  
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) + 
  scale_color_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) + 
  scale_x_discrete(limits=gr.resistome_abundance.result$index[order(match(gr.resistome_abundance.result$index, variable.group$index),decreasing = T)])+
  theme_ykm()+coord_flip()+theme(legend.position = "none")

gr.abundance.p2=ggplot(gr.resistome_abundance.effect_data, aes(x = x, y = value*100,fill=index))+   
  geom_col(width = 0.6,color="black",position = "stack")+  
  labs(x="", y="Relative effect of esimates (%)")+  
  scale_fill_manual(values = c("Soil"="#00A087","Microbial"="#4DBBD5","Geography"="#E64B35","Unexplained"="#3C5488")) +  
  scale_y_continuous(expand = c(0,0.001)) +scale_x_discrete(expand = c(0,0.001)) +
  theme_minimal()+
  theme(aspect.ratio = 6,legend.position = "none",
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))

fa.abundance.p1 | fa.abundance.p2 | gr.abundance.p1 | gr.abundance.p2

(fa.abundance.p1 | fa.abundance.p2 | gr.abundance.p1 | gr.abundance.p2)/(fa.shannon.p1 | fa.shannon.p2 | gr.shannon.p1 | gr.shannon.p2)




###########  resistome tax cor ###########

cor_matrix2matrix = function(matrix_a,matrix_b,cor_method){
  phage_num = as.numeric(nrow(matrix_a))
  species_num = as.numeric(nrow(matrix_b))
  phage_attribute = as.data.frame(matrix_a)[,1]
  species_attribute = as.data.frame(matrix_b)[,1]
  phage_attribute2 = rep(phage_attribute,each = species_num)
  cor_result = matrix(ncol = 2,nrow = phage_num*species_num)
  colnames(cor_result) = c("p","r")
  m_num = seq(1:phage_num)
  n_num = seq(1:species_num)
  h.phage.otu = matrix_a[,2:ncol(matrix_a)]
  h.species.otu = matrix_b[,2:ncol(matrix_b)]
  for (m in m_num){
    for (n in n_num){
      cor = cor.test(as.numeric(h.phage.otu[m,]),as.numeric(h.species.otu[n,]),method= cor_method )
      cor_result[(m-1)*species_num + n,1] = as.numeric(cor[3])
      cor_result[(m-1)*species_num + n,2] = as.numeric(cor[4])
    }
  }
  a_names = rep(phage_attribute,each = species_num)
  b_names = rep(species_attribute,phage_num)
  result = data.frame(a_names,b_names,cor_result)
  return(result)
}

resistome.fa=subset(sample.data,group=="FA")[,c(1,23,22,21)]
resistome.gr=subset(sample.data,group=="GR")[,c(1,23,22,21)]
tax.fa=subset(sample.data,group=="FA")[,c(1,24,25,28,29,31,32)]
tax.gr=subset(sample.data,group=="GR")[,c(1,24,25,28,29,31,32)]
resistome.fa.data = data.frame(index=colnames(resistome.fa)[-1],t(dataframe2matrix(resistome.fa)))
resistome.gr.data = data.frame(index=colnames(resistome.gr)[-1],t(dataframe2matrix(resistome.gr)))
tax.fa.data = data.frame(index=colnames(tax.fa)[-1],t(dataframe2matrix(tax.fa)))
tax.gr.data = data.frame(index=colnames(tax.gr)[-1],t(dataframe2matrix(tax.gr)))

resistome_tax.cor.fa = cor_matrix2matrix(resistome.fa.data,tax.fa.data,"pearson")
cor.p = dataframe2matrix(spread(resistome_tax.cor.fa[,-4],key = "b_names",value = "p"))
cor.r = dataframe2matrix(spread(resistome_tax.cor.fa[,-3],key = "b_names",value = "r"))
corrplot(
  t(cor.r), 
  p.mat = t(cor.p), 
  sig.level = 0.05,
  insig = "blank",
  method = "circle",
  is.corr = TRUE,
  col.lim = c(-1, 1),
  tl.col = "black",
  addgrid.col = "black",  # 关键参数：设置边框颜色
  tl.cex = 0.5
)

resistome_tax.cor.gr = cor_matrix2matrix(resistome.gr.data,tax.gr.data,"pearson")
cor.p = dataframe2matrix(spread(resistome_tax.cor.gr[,-4],key = "b_names",value = "p"))
cor.r = dataframe2matrix(spread(resistome_tax.cor.gr[,-3],key = "b_names",value = "r"))
corrplot(
  t(cor.r), 
  p.mat = t(cor.p), 
  sig.level = 0.05,
  insig = "blank",
  method = "circle",
  is.corr = TRUE,
  col.lim = c(-1, 1),
  tl.col = "black",
  addgrid.col = "black",  # 关键参数：设置边框颜色
  tl.cex = 0.5
)


###########  host  ###########
library(ggsankey)
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}

arg_mge.anno = rbind.data.frame(arg.anno[,c(1,6:8)],mge.anno[,c(1,7:9)])
arg_mge.nr.anno = merge(arg_mge.nr[,c(1,97)],arg_mge.anno,1)
arg_mge.anno.sep = separate(arg_mge.nr.anno,Tax,into = c("phylum","species"),sep = ";c__") %>% 
  separate(species,into = c("class","species"),sep = ";o__")%>% 
  separate(species,into = c("order","species"),sep = ";f__")

arg_mge.anno.sep$phylum = gsub(pattern="Myxococcota_A", replacement="Myxococcota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$phylum = gsub(pattern="Myxococcota_B", replacement="Myxococcota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$phylum = gsub(pattern="Nitrospinota_B", replacement="Nitrospinota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$phylum = gsub(pattern="Nitrospirota_A", replacement="Nitrospirota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$phylum = gsub(pattern="Nitrospirota_B", replacement="Nitrospirota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$Subtype = gsub(pattern="Streptomyces rishiriensis parY mutant conferring resistance to aminocoumarin", replacement="parY", arg_mge.anno.sep$Subtype)
arg_mge.anno.sep$phylum[arg_mge.anno.sep$phylum %nin% c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota") ] <- "AAA_Other phylum"
arg_mge.anno.sep$class[arg_mge.anno.sep$class %nin% c("Actinomycetia","Alphaproteobacteria","Gammaproteobacteria","Thermoleophilia","Anaerolineae") ] <- "AAA_Other class"
arg_mge.anno.sep$order[arg_mge.anno.sep$order %nin% c("Rhizobiales","Mycobacteriales","Burkholderiales","Xanthomonadales","Streptomycetales","Propionibacteriales","Steroidobacterales","Enterobacterales","Gaiellales","Actinomycetales","Streptosporangiales","Sphingomonadales","Ardenticatenales") ] <- "AAA_Other order"
arg_mge.anno.sep$Type[arg_mge.anno.sep$Type %in% c("ist","integrase","plasmid") ] <- "Other_MGEs"
arg_mge.anno.sep$Type[arg_mge.anno.sep$Type %nin% c("vancomycin","multidrug","transposase","rifamycin","insertion_sequence","aminocoumarin","quinolone","Other_MGEs") ] <- "Other_ARGs"
arg_mge.anno.sep$Subtype[arg_mge.anno.sep$Subtype %nin% c("tnpA","IS91","vanRO","mtrA","vanSO","MexF","QepA4","Bado_rpoB_RIF","rsmA","novA","ceoB","MuxB","rphA","parY") ] <- "AAA_Other genes"

arg_mge.anno.sep$class = gsub(pattern="Thermoleophilia", replacement="Ad_Thermoleophilia", arg_mge.anno.sep$class)
arg_mge.anno.sep$class = gsub(pattern="Anaerolineae", replacement="AAB_Anaerolineae", arg_mge.anno.sep$class)
arg_mge.anno.sep$order = gsub(pattern="Streptosporangiales", replacement="P_Streptosporangiales", arg_mge.anno.sep$order)
arg_mge.anno.sep$order = gsub(pattern="Streptomycetales", replacement="P_Streptomycetales", arg_mge.anno.sep$order)
arg_mge.anno.sep$order = gsub(pattern="Enterobacterales", replacement="X_Enterobacterales", arg_mge.anno.sep$order)
arg_mge.anno.sep$order = gsub(pattern="Burkholderiales", replacement="X_Burkholderiales", arg_mge.anno.sep$order)
arg_mge.anno.sep$phylum = gsub(pattern="d__Bacteria;p__Gemmatimonadota", replacement="d__Bacteria;p__A_Gemmatimonadota", arg_mge.anno.sep$phylum)
arg_mge.anno.sep$phylum = gsub(pattern="d__Bacteria;p__Chloroflexota", replacement="d__Bacteria;p__Acj_Chloroflexota", arg_mge.anno.sep$phylum)


arg_mge.anno.count = arg_mge.anno.sep[,-c(1,5)] %>% 
  mutate(count=1)
arg_mge.anno.sum = aggregate(count ~ phylum + class + order +Subtype+Type+Mechanism, data = arg_mge.anno.count, FUN = base::sum)
arg_mge.sankey = make_long(arg_mge.anno.sum[,c(6,5,4,3,2,1,7)],Mechanism,Type,Subtype,order,class,phylum,value = count)

cols <- c("#FCF7D5FF","#F5F3C1FF","#EAF0B5FF","#DDECBFFF","#D0E7CAFF","#C2E3D2FF","#B5DDD8FF","#A8D8DCFF","#E8ECFBFF","#DDD8EFFF","#D1C1E1FF","#C3A8D1FF","#B58FC2FF","#A778B4FF","#9B62A7FF","#8C4E99FF","#6F4C9BFF","#6059A9FF","#5568B8FF","#4E79C5FF","#4D8AC6FF","#4E96BCFF","#549EB3FF","#59A5A9FF","#60AB9EFF","#69B190FF","#77B77DFF","#8CBC68FF","#A6BE54FF","#BEBC48FF","#D1B541FF","#DDAA3CFF","#E49C39FF","#E78C35FF","#E67932FF","#E4632DFF","#DF4828FF","#DA2222FF","#B8221EFF","#95211BFF","#721E17FF","#521A13FF")
ggplot(arg_mge.sankey,        
       aes(x = x,next_x = next_x,node = node,next_node = next_node,fill = factor(node))) +  
  geom_sankey(width = 0.15,smooth = 5,space = 15,na.rm = TRUE,position = "identity",flow.alpha = 0.4,node.color = "transparent") +  
  geom_sankey_text(aes(label = node),width = 0.2,space = 15,position = "identity",size = 3,color = "black",hjust = 0.1) +  
  theme_sankey(base_size = 7) + paletteer::scale_fill_paletteer_d("palettesForR::Web")+  
  theme(plot.title = element_text(hjust = 0.08), plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"),axis.title.x = element_blank(),axis.text.x = element_text(color = "black", size =7), legend.position = "none")

ggplot(arg_mge.sankey,        
       aes(x = x,next_x = next_x,node = node,next_node = next_node,fill = factor(node))) +  
  geom_sankey(width = 0.15,smooth = 5,space = 15,na.rm = TRUE,position = "identity",flow.alpha = 0.4,node.color = "transparent") +  
  geom_sankey_text(aes(label = node),width = 0.2,space = 15,position = "identity",size = 3,color = "black",hjust = 0.1) +  
  theme_sankey(base_size = 7) + scale_fill_viridis_d(drop=FALSE) +  
  theme(plot.title = element_text(hjust = 0.08), plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"),axis.title.x = element_blank(),axis.text.x = element_text(color = "black", size =7), legend.position = "none")



arg_mge.tax = unique.data.frame(merge(arg_mge.nr[,97],tax,1))
arg_mge.tax.sep = separate(arg_mge.tax,1,into = c("phylum","species"),sep = ";c__", remove = F) %>% 
  separate(species,into = c("class","species"),sep = ";o__")%>% 
  separate(species,into = c("order","species"),sep = ";f__")
arg_mge.tax.sep$phylum = gsub(pattern="Myxococcota_A", replacement="Myxococcota", arg_mge.tax.sep$phylum)
arg_mge.tax.sep$phylum = gsub(pattern="Myxococcota_B", replacement="Myxococcota", arg_mge.tax.sep$phylum)
arg_mge.tax.sep$phylum = gsub(pattern="Nitrospinota_B", replacement="Nitrospinota", arg_mge.tax.sep$phylum)
arg_mge.tax.sep$phylum = gsub(pattern="Nitrospirota_A", replacement="Nitrospirota", arg_mge.tax.sep$phylum)
arg_mge.tax.sep$phylum = gsub(pattern="Nitrospirota_B", replacement="Nitrospirota", arg_mge.tax.sep$phylum)

host.order = apply(arg_mge.tax.sep[,6:95],2,aggregate,list(arg_mge.tax.sep$order),sum)
host.order.df = data.frame(host.order)[,c(1,seq(2,length(host.order)*2,2))]
colnames(host.order.df)=c("Type",colnames(arg_mge.tax.sep[,6:95]))

order.gather =gather(host.order.df,-1,key="sample",value ="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  unite("field", location, field, sep = "|", remove = F)
input = order.gather[,c(1,3,5,7)]
ttest.result = data.frame(field="y",Type="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$Type)){
    data.pro = subset(data.sub,Type == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,Type=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("order", field, Type, sep = "-", remove = FALSE)
ttest.sig = subset(ttest.result2,sig=="sig")
ttest.result3 = subset(ttest.result2,Type %in% unique(ttest.sig$Type))
order.summary <- summarySE(data=order.gather, measurevar="value",groupvars=c("Type","field","location","group"))
order.spread = spread(order.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("order", field, Type, sep = "-", remove = FALSE) %>% 
  merge(ttest.result3,"order") %>% 
  mutate(updown="none")
order.spread$updown <- case_when(  
  order.spread$fc < 0 & order.spread$p < 0.05 ~ 'grass'   ,  
  order.spread$fc > 0 & order.spread$p < 0.05 ~ 'farm',
  .default = "none")
order.spread.select = subset(order.spread,Type.x %in% c("Rhizobiales","Mycobacteriales","Burkholderiales","Xanthomonadales","Streptomycetales","Propionibacteriales","Steroidobacterales","Enterobacterales","Gaiellales","Actinomycetales","Streptosporangiales","Sphingomonadales","Ardenticatenales"))

order.spread2 = spread(order.spread.select[,c(2,3,7)],key = "field.x",value = "fc")
order.clust.plot = hclust(dist((dataframe2matrix(order.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 1)+
  hexpand(2) 
clust.data = subset(order.clust.plot$data,is.na(label)==F)
order.spread.select$Type.x = factor(order.spread.select$Type.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
order.plot = ggplot(order.spread.select,aes(x=field.x,y=Type.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3,5), range = c(1, 8))+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(24,25,21)) 

order.clust.plot | order.plot

host.phylum = apply(arg_mge.tax.sep[,6:95],2,aggregate,list(arg_mge.tax.sep$phylum),sum)
host.phylum.df = data.frame(host.phylum)[,c(1,seq(2,length(host.phylum)*2,2))]
colnames(host.phylum.df)=c("Type",colnames(arg_mge.tax.sep[,6:95]))

phylum.gather =gather(host.phylum.df,-1,key="sample",value ="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  unite("field", location, field, sep = "|", remove = F)
input = phylum.gather[,c(1,3,5,7)]
ttest.result = data.frame(field="y",Type="x",df=1,t=1,p=1,stderr=1)
for (i in unique(input$field)){
  data.sub = subset(input,field == i )
  for (j in unique(data.sub$Type)){
    data.pro = subset(data.sub,Type == j )
    result1 <- my.t.test(value~group,data=data.pro)
    if(is.na(result1)[1]==T){ttest = data.frame(field=i,Type=j,df=NA,t=NA,p=NA,stderr=NA)}
    else{ttest = data.frame(field=i,Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
    ttest.result = rbind.data.frame(ttest.result,ttest)
  }
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("phylum", field, Type, sep = "-", remove = FALSE)
ttest.sig = subset(ttest.result2,sig=="sig")
ttest.result3 = subset(ttest.result2,Type %in% unique(ttest.sig$Type))
phylum.summary <- summarySE(data=phylum.gather, measurevar="value",groupvars=c("Type","field","location","group"))
phylum.spread = spread(phylum.summary[,c(1,2,3,4,6)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("phylum", field, Type, sep = "-", remove = FALSE) %>% 
  merge(ttest.result3,"phylum") %>% 
  mutate(updown="none")
phylum.spread$updown <- case_when(  
  phylum.spread$fc < 0 & phylum.spread$p < 0.05 ~ 'grass'   ,  
  phylum.spread$fc > 0 & phylum.spread$p < 0.05 ~ 'farm',
  .default = "none")
phylum.spread.select = subset(phylum.spread,Type.x %in% c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota"))

phylum.spread2 = spread(phylum.spread.select[,c(2,3,7)],key = "field.x",value = "fc")
phylum.clust.plot = hclust(dist((dataframe2matrix(phylum.spread2)))) %>% 
  ggtree(layout="rectangular",branch.length="none")+
  theme(legend.position = "none")+theme(aspect.ratio = 0.3)+
  hexpand(2) 
clust.data = subset(phylum.clust.plot$data,is.na(label)==F)
phylum.spread.select$Type.x = factor(phylum.spread.select$Type.x , levels=as.factor(clust.data$label[order(clust.data$y,decreasing=F)]) )
phylum.plot = ggplot(phylum.spread.select,aes(x=field.x,y=Type.x)) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Orange-White-Blue", limits = c(-7, 7))+
  theme_ykm()+
  geom_point(aes(size=abs(fc),shape= updown, fill=fc))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=c(4.5,9.5,11.5,13.5),size=.4)+
  scale_size(breaks = c(1, 3,5), range = c(1, 8))+
  theme(legend.position = "none") +theme(aspect.ratio = 0.3)+
  scale_shape_manual(values = c(24,25,21)) 

(phylum.clust.plot | phylum.plot)/(order.clust.plot | order.plot)

ggboxplot(phylum.gather,x="group",y="value",color = "group",
          add = "jitter", outlier.shape = NA,error.plot="errorbar")+
  stat_compare_means(comparisons = list( c("FA", "GR")),method = "t.test",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***","**","*","n.s.")),vjust=0.6)+
  theme_ykm()+
  xlab("")+
  ylab("")+
  scale_colour_jco()+scale_fill_jco()+theme(legend.position = "none")+
  facet_wrap(~Type,ncol=6,scales = "free")


input = phylum.gather[,c(1,5,7)]
ttest.result = data.frame(Type="x",df=1,t=1,p=1,stderr=1)
for (j in unique(input$Type)){
  data.pro = subset(input,Type == j )
  result1 <- my.t.test(value~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(Type=j,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  ttest.result = rbind.data.frame(ttest.result,ttest)
  
}
ttest.result = ttest.result[-1,]
ttest.result2 = mutate(ttest.result,sig=ifelse(ttest.result$p<0.05,"sig","nonsig")) %>% 
  unite("phylum", Type, sep = "-", remove = FALSE)
ttest.sig = subset(ttest.result2,sig=="sig")
ttest.result3 = subset(ttest.result2,Type %in% unique(ttest.sig$Type))
phylum_group.summary <- summarySE(data=phylum.gather, measurevar="value",groupvars=c("Type","group"))
phylum_group.spread = spread(phylum_group.summary[,c(1,2,4)],key = "group",value = "value") %>% 
  mutate(fc = log2(FA/GR))%>% 
  unite("phylum",  Type, sep = "-", remove = FALSE) %>% 
  merge(ttest.result3,"phylum") %>% 
  mutate(updown="none")
phylum_group.spread$updown <- case_when(  
  phylum_group.spread$fc < 0 & phylum_group.spread$p < 0.05 ~ 'grass'   ,  
  phylum_group.spread$fc > 0 & phylum_group.spread$p < 0.05 ~ 'farm',
  .default = "none")
#write.table (phylum_group.spread, file = "/home/shared/project/2025/ecotone90/analysis/result/host.phylum_group.ttest.txt", sep = "\t", row.names = F ,col.names = TRUE, quote = FALSE)

phylum_group.spread$phylum = gsub(pattern="d__Bacteria;p__", replacement="", phylum_group.spread$phylum)
phylum_group.dotchart=ggdotchart(phylum_group.spread, x="phylum", y="fc", color = "updown", 
                                 palette = c("#0073C2","#EFC000"), 
                                 sorting = "descending", add = "segment", 
                                 add.params = list(color="lightgray", size=2), 
                                 group = "updown", dot.size = 4,
                                 ggtheme = theme_ykm())+ geom_hline(aes(yintercept=0))+
  theme(aspect.ratio = 0.5,legend.position = "none",axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  ylim(-3.2,3.2)


phylum_group.spread2 =spread(phylum_group.summary[,c(1,2,4)],key = "group",value = "value") 
phylum_group.6 = filter_by_rank(phylum_group.spread2,6)
phylum_group.spread2$Type[ phylum_group.spread2$Type %nin% phylum_group.6$Type ] <- "Others"
phylum_group.sum = apply(phylum_group.spread2[,-1],2,aggregate,list(phylum_group.spread2$Type),sum)
phylum_group.sum.df = data.frame(phylum_group.sum)[,c(1,seq(2,length(phylum_group.sum)*2,2))]
colnames(phylum_group.sum.df)=c("phylum","FA","GR")

percent.phylum = apply(phylum_group.sum.df[,-1],2,format_percent)
percent.phylum2 = data.frame(Type =phylum_group.sum.df$phylum,percent.phylum)
percent.phylum_group.gather =gather(percent.phylum2,-1,key="sample",value="value")
percent.phylum_group.gather$Type = factor(percent.phylum_group.gather$Type, levels=as.factor(c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota","d__Bacteria;p__Desulfobacterota_B","Others") ) )

percent.phylum_group.plot=ggplot(percent.phylum_group.gather, aes(x = sample, y=value*100, fill = Type,
                                                                  stratum = Type, alluvium = Type)) +
  scale_fill_grafify(palette = "kelly") +
  scale_y_continuous(expand = c(0,0)) +
  geom_col(width = 0.6,
           color= NA, size = 0.5) +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0, #直线连线效果
            color= 'white', size = 0.5) +
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0,
                fill= NA, color = 'white', size = 0.5)+
  theme_minimal()+
  theme(aspect.ratio = 1.5,
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title.y = element_text(hjust=0.5,size=7),axis.text.x = element_text(angle = 0,hjust = 0.5,size=7))+
  theme(legend.position = "none") 

percent.phylum_group.plot+phylum_group.dotchart

arg.nr=merge(arg,nr,1)
arg.tax = unique.data.frame(merge(arg.nr[,102],tax,1))
arg.tax.gather = gather(arg.tax,-1,key="sample",value="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  separate(x,into = c("phylum","species"),sep = ";c__",remove = F) %>% 
  separate(species,into = c("class","species"),sep = ";o__") %>% 
  separate(species,into = c("order","species"),sep = ";f__")

arg.tax.group.ave = summarySE(data=arg.tax.gather, measurevar="value",groupvars=c("x","phylum","group"))
arg.tax.group.spread = spread(arg.tax.group.ave[,c(1,2,3,5)],key = "group",value = "value") %>% 
  mutate(logfa=log10(FA)) %>% 
  mutate(loggr=log10(GR))
arg.tax.group.spread$phylum[ arg.tax.group.spread$phylum %nin% c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota","d__Bacteria;p__Desulfobacterota_B") ] <- "Others"
arg.tax.group.spread$phylum = factor(arg.tax.group.spread$phylum, levels=as.factor(c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota","Others","d__Bacteria;p__Desulfobacterota_B") ) )

input = arg.tax.gather[,c(1,9,11)]
colnames(input)[1]="Type"
ttest.result = data.frame(Type="x",df=1,t=1,p=1,stderr=1)
for (j in unique(input$Type)){
  data.pro = subset(input,Type == j )
  result1 <- my.t.test(value~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(Type=j,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  ttest.result = rbind.data.frame(ttest.result,ttest)
  
}
ttest.result = ttest.result[-1,]
arg.tax.group.ttest = merge(arg.tax.group.spread,ttest.result,1) %>% 
  mutate(sig=ifelse(arg.tax.group.ttest$p<0.05,"sig","nonsig")) %>% 
  mutate(alpha=ifelse(arg.tax.group.ttest$p<0.05,1,0.25)) %>% 
  mutate(FA01=normalization(arg.tax.group.ttest$FA))%>% 
  mutate(GR01=normalization(arg.tax.group.ttest$GR))


arg.tax.fagr.plot=ggplot(arg.tax.group.ttest, aes(x=log10(FA), y=log10(GR))) +
  geom_point(aes(color = phylum),alpha=arg.tax.group.ttest$alpha) +
  theme_ykm()+
  xlab("Farmland")+
  ylab("Pastureland")+
  scale_fill_grafify(palette = "kelly") +scale_color_grafify(palette = "kelly") +
  xlim(-4.2,4.2)+ylim(-4.2,4.2)+geom_abline(intercept = 0, slope = 1)+
  theme(legend.position = "none") 
arg.tax.fagr.plot=ggplot(arg.tax.group.ttest, aes(x=FA01, y=GR01)) +
  geom_point(aes(color = phylum),alpha=arg.tax.group.ttest$alpha) +
  theme_ykm()+
  xlab("Farmland")+
  ylab("Pastureland")+
  scale_fill_grafify(palette = "kelly") +scale_color_grafify(palette = "kelly") +
  xlim(0,1)+ylim(0,1)+geom_abline(intercept = 0, slope = 1)+
  theme(legend.position = "none") 
mge.nr=merge(mge,nr,1)
mge.tax = unique.data.frame(merge(mge.nr[,100],tax,1))
mge.tax.gather = gather(mge.tax,-1,key="sample",value="value") %>% 
  separate(sample,into = c("location","field","group","rep"),sep = "-",remove = F) %>% 
  separate(x,into = c("phylum","species"),sep = ";c__",remove = F) %>% 
  separate(species,into = c("class","species"),sep = ";o__") %>% 
  separate(species,into = c("order","species"),sep = ";f__")
mge.tax.gather$phylum = gsub(pattern="Myxococcota_A", replacement="Myxococcota", mge.tax.gather$phylum)
mge.tax.gather$phylum = gsub(pattern="Myxococcota_B", replacement="Myxococcota", mge.tax.gather$phylum)
mge.tax.gather$phylum = gsub(pattern="Nitrospinota_B", replacement="Nitrospinota", mge.tax.gather$phylum)
mge.tax.gather$phylum = gsub(pattern="Nitrospirota_A", replacement="Nitrospirota", mge.tax.gather$phylum)
mge.tax.gather$phylum = gsub(pattern="Nitrospirota_B", replacement="Nitrospirota", mge.tax.gather$phylum)

mge.tax.group.ave = summarySE(data=mge.tax.gather, measurevar="value",groupvars=c("x","phylum","group"))
mge.tax.group.spread = spread(mge.tax.group.ave[,c(1,2,3,5)],key = "group",value = "value") %>% 
  mutate(logfa=log10(FA)) %>% 
  mutate(loggr=log10(GR))
mge.tax.group.spread$phylum[ mge.tax.group.spread$phylum %nin% c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota","d__Bacteria;p__Desulfobacterota_B") ] <- "Others"
mge.tax.group.spread$phylum = factor(mge.tax.group.spread$phylum, levels=as.factor(c("d__Bacteria;p__Actinomycetota","d__Bacteria;p__Pseudomonadota","d__Bacteria;p__Acidobacteriota","d__Bacteria;p__Chloroflexota","d__Bacteria;p__Gemmatimonadota","Others","d__Bacteria;p__Desulfobacterota_B") ) )

input = mge.tax.gather[,c(1,9,11)]
colnames(input)[1]="Type"
ttest.result = data.frame(Type="x",df=1,t=1,p=1,stderr=1)
for (j in unique(input$Type)){
  data.pro = subset(input,Type == j )
  result1 <- my.t.test(value~group,data=data.pro)
  if(is.na(result1)[1]==T){ttest = data.frame(Type=j,df=NA,t=NA,p=NA,stderr=NA)}
  else{ttest = data.frame(Type=j,df=result1$parameter,t=result1$statistic,p=result1$p.value,stderr=result1$stderr)}
  ttest.result = rbind.data.frame(ttest.result,ttest)
  
}
ttest.result = ttest.result[-1,]
mge.tax.group.ttest = merge(mge.tax.group.spread,ttest.result,1) %>% 
  mutate(sig=ifelse(mge.tax.group.ttest$p<0.05,"sig","nonsig")) %>% 
  mutate(alpha=ifelse(mge.tax.group.ttest$p<0.05,1,0.25)) 
#write.table (mge.tax.group.ttest, file = "/home/shared/project/2025/ecotone90/analysis/result/mge.tax.group.ttest.txt", sep = "\t", row.names = F ,col.names = TRUE, quote = FALSE)

mge.tax.fagr.plot=ggplot(mge.tax.group.ttest, aes(x=log10(FA), y=log10(GR))) +
  geom_point(aes(color = phylum),alpha=mge.tax.group.ttest$alpha) +
  theme_ykm()+
  xlab("Farmland")+
  ylab("Pastureland")+
  scale_fill_grafify(palette = "kelly") +scale_color_grafify(palette = "kelly") +
  xlim(-3.5,3.5)+ylim(-3.5,3.5)+geom_abline(intercept = 0, slope = 1)+
  theme(legend.position = "none") 

arg.tax.fagr.plot | mge.tax.fagr.plot

########### arg_mge cornet   ###################
library(Hmisc)
library(ggnetwork)
library(igraph)
library(sna)
library(intergraph)

table_order <- function(x) {
  table(x)[order(table(x))]
}

arg_mge.fa =arg_mge.subtype[,colnames(arg_mge.subtype) %in% c("Subtype",subset(meta,group=="FA")[,1])]
arg_mge.gr =arg_mge.subtype[,colnames(arg_mge.subtype) %in% c("Subtype",subset(meta,group=="GR")[,1])]

res = rcorr(t(dataframe2matrix(arg_mge.fa)),type = "spearman")
res.p = res$P
res.r = res$r
mol.names<-rownames(res.p[,-1])
myindex <- which(lower.tri(res.p[,-1]) == TRUE, arr.ind = TRUE)
target <- mol.names[myindex[, 1]]
source <- mol.names[myindex[, 2]]
p_value<-res.p[,-1][lower.tri(res.p[,-1])]
cov_value<-res.r[,-1][lower.tri(res.r[,-1])]
out.fa<-data.frame(source,target,cov_value,p_value)
sig.fa = subset(out.fa, p_value < 0.05 & cov_value > 0.5 ) %>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=1,by.y=2)%>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=2,by.y=2) 
sig.fa$group=ifelse(sig.fa$group.x != sig.fa$group.y,"ARG_MGE","Others")
#sig.fa =subset(sig.fa,Type.x != Type.y)
node.fa = data.frame(ID=unique(c(sig.fa$target,sig.fa$source))) %>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=1,by.y=2)
node.fa$class = node.fa$Type
node.fa$class[ node.fa$class %in% c("transposase","insertion_sequence","integrase","ist","plasmid") ] <- "MGEs"
node.fa$class[ node.fa$class %nin% c("multidrug","beta_lactam","aminoglycoside","MLS","vancomycin","tetracycline","MGEs") ] <- "Other_ARGs"

igraph.fa <- graph_from_edgelist(as.matrix(sig.fa[,1:2]), directed = F)
igraph.node.fa = data.frame(ID = names(V(igraph.fa)),degree = igraph::degree(igraph.fa,mode="in")) %>% 
  join(node.fa,by = "ID")
E(igraph.fa)$group <- sig.fa$group
V(igraph.fa)$class <- igraph.node.fa$class
V(igraph.fa)$degree <- igraph.node.fa$degree
ggnetwork_net <- ggnetwork(x = igraph.fa,layout=igraph::layout.gem(igraph.fa))
ggplot(data = ggnetwork_net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color=group),alpha=0.5) +
  geom_nodes(aes(color=class,size=degree)) +
  theme_blank()+
  scale_size(breaks = c(0.5, 1.5, 3), range = c(0.5, 5))

res = rcorr(t(dataframe2matrix(arg_mge.gr)),type = "spearman")
res.p = res$P
res.r = res$r
mol.names<-rownames(res.p[,-1])
myindex <- which(lower.tri(res.p[,-1]) == TRUE, arr.ind = TRUE)
target <- mol.names[myindex[, 1]]
source <- mol.names[myindex[, 2]]
p_value<-res.p[,-1][lower.tri(res.p[,-1])]
cov_value<-res.r[,-1][lower.tri(res.r[,-1])]
out.gr<-data.frame(source,target,cov_value,p_value)
sig.gr = subset(out.gr, p_value < 0.05 & cov_value > 0.5 ) %>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=1,by.y=2)%>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=2,by.y=2) 
sig.gr$group=ifelse(sig.gr$group.x != sig.gr$group.y,"ARG_MGE","Others")
#sig.gr =subset(sig.gr,Type.x != Type.y)
node.gr = data.frame(ID=unique(c(sig.gr$target,sig.gr$source))) %>% 
  merge(arg_mge.subtype[,c(1:2,93)],by.x=1,by.y=2)
node.gr$class = node.gr$Type
node.gr$class[ node.gr$class %in% c("transposase","insertion_sequence","integrase","ist","plasmid") ] <- "MGEs"
node.gr$class[ node.gr$class %nin% c("multidrug","beta_lactam","aminoglycoside","MLS","vancomycin","tetracycline","MGEs") ] <- "Other_ARGs"

igraph.gr <- graph_from_edgelist(as.matrix(sig.gr[,1:2]), directed = F)
igraph.node.gr = data.frame(ID = names(V(igraph.gr)),degree = igraph::degree(igraph.gr,mode="in")) %>% 
  join(node.gr,by = "ID")
E(igraph.gr)$group <- sig.gr$group
V(igraph.gr)$class <- igraph.node.gr$class
V(igraph.gr)$degree <- igraph.node.gr$degree
ggnetwork_net <- ggnetwork(x = igraph.gr,layout=igraph::layout.gem(igraph.gr))
ggplot(data = ggnetwork_net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color=group),alpha=0.5) +
  geom_nodes(aes(color=class,size=degree)) +
  theme_blank()+
  scale_size(breaks = c(0.5, 1.5, 3), range = c(0.5, 5))

igraph.node.fa$fagr = "FA"
igraph.node.gr$fagr = "GR"
igraph.node.all = rbind.data.frame(igraph.node.fa,igraph.node.gr)
igraph.node.spread = spread(igraph.node.all,key = "fagr",value = "degree")
igraph.node.spread[is.na(igraph.node.spread)]<- 0
igraph.node.spread2 =mutate(igraph.node.spread,sum=FA+GR) %>% 
  mutate(difference=abs(FA-GR))
igraph.node.sum10=igraph.node.spread2[order(igraph.node.spread2[,"sum"],decreasing=T),][1:11,] 
igraph.node.diff10=igraph.node.spread2[order(igraph.node.spread2[,"difference"],decreasing=T),][1:11,] 
igraph.node.20 = rbind.data.frame(igraph.node.sum10,igraph.node.diff10)
igraph.node.20[10,1] = "parY"
igraph.node.20.unite=unite(igraph.node.20,"name", Type, ID, sep = "-", remove = FALSE)
igraph.node.20.gather = gather(igraph.node.20.unite[,1:7],6:7,key="fagr",value ="value")

ggplot(igraph.node.20.gather,aes(x=name,y=value,fill=fagr))+
  geom_bar(position =position_dodge(0.6),width = 0.5,stat = "identity",color="black")+
  theme_ykm()+theme(aspect.ratio = 1.5,legend.position = "none")+
  scale_fill_jco()+
  coord_flip()


edge.fa.cross = subset(edge.fa,group=="ARG_MGE")
node.fa.cross = subset(node.fa[,-2],ID %in% unique(c(edge.fa.cross$source,edge.fa.cross$target)))
edge.gr.cross = subset(edge.gr,group=="ARG_MGE")
node.gr.cross = subset(node.gr[,-2],ID %in% unique(c(edge.gr.cross$source,edge.gr.cross$target)))

mge.fa = as.data.frame(table_order(edge.fa.cross[,1])) %>% 
  merge(unique.data.frame(edge.fa.cross[,c(1,7)]),1) %>% 
  mutate(group="FA")
mge.gr = as.data.frame(table_order(edge.gr.cross[,1])) %>% 
  merge(unique.data.frame(edge.gr.cross[,c(1,7)]),1) %>% 
  mutate(group="GR")
mge.fa.gr = rbind.data.frame(mge.fa,mge.gr)
ggplot(mge.fa.gr,aes(x=x,y=Freq,fill=group))+
  geom_bar(position =position_dodge(0.6),width = 0.5,stat = "identity",color="black")+
  theme_ykm()+theme(aspect.ratio = 1.5,legend.position = "none")+
  scale_fill_jco()+
  coord_flip()




