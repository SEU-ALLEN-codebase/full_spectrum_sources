library(circlize)
library(dendextend)
library(openxlsx)
library(stringr)
library(devtools)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)
library(dendextend)
library(dendsort)



axonF<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/axon_feature.txt")
colnames(axonF)<-c('name','arbor_id','region','max_density','num_nodes','total_path_length','volume','branch','d_to_soma','d_to_soma2','hub','variance_ratio')

table<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/1891_type.txt",header=TRUE)
ssp_layer<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/SSp_layer.txt")
colnames(ssp_layer)<-c('name','L1','L2/3','L4','L5','L6a','L6b','Non-cortical')
ssp_arbor_layer<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/SSp_arbors_layer.txt")
colnames(ssp_arbor_layer)<-c('name','L1','L2/3','L4','L5','L6a','L6b','Non-cortical')

LD_order<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/category.txt",header=TRUE,sep=",")

ssp_table<-matrix(0,nrow = nrow(ssp_arbor_layer),ncol=8)
for(i in 1:nrow(ssp_arbor_layer)){
  n<-nchar(ssp_arbor_layer$name[i])
  ssp_table[i,1]<-substr(ssp_arbor_layer$name[i],1,n-7)
  ssp_table[i,2]<-substr(ssp_arbor_layer$name[i],n,n)
  if(ssp_arbor_layer$L1[i]>1000){
    ssp_table[i,3]<-ssp_arbor_layer$L1[i]
  }
  else{
    ssp_table[i,3]<-0
  }
  if(ssp_arbor_layer$`L2/3`[i]>1000){
    ssp_table[i,4]<-ssp_arbor_layer$`L2/3`[i]
  }
  else{
    ssp_table[i,4]<-0
  }
  if(ssp_arbor_layer$L4[i]>1000){
    ssp_table[i,5]<-ssp_arbor_layer$L4[i]
  }
  else{
    ssp_table[i,5]<-0
  }
  if(ssp_arbor_layer$L5[i]>1000){
    ssp_table[i,6]<-ssp_arbor_layer$L5[i]
  }
  else{
    ssp_table[i,6]<-0
  }
  if(ssp_arbor_layer$L6a[i]>1000){
    ssp_table[i,7]<-ssp_arbor_layer$L6a[i]
  }
  else{
    ssp_table[i,7]<-0
  }
  if(ssp_arbor_layer$L6b[i]>1000){
    ssp_table[i,8]<-ssp_arbor_layer$L6b[i]
  }
  else{
    ssp_table[i,8]<-0
  }
}

ssp_table2<-as.data.frame(ssp_table)
colnames(ssp_table2)<-c('name','arbor_id','L1','L2/3','L4','L5','L6a','L6b')
ssp_table2$L1<-as.numeric(ssp_table2$L1)
ssp_table2$`L2/3`<-as.numeric(ssp_table2$`L2/3`)
ssp_table2$L4<-as.numeric(ssp_table2$L4)
ssp_table2$L5<-as.numeric(ssp_table2$L5)
ssp_table2$L6a<-as.numeric(ssp_table2$L6a)
ssp_table2$L6b<-as.numeric(ssp_table2$L6b)
ssp_table2<-ssp_table2[which(!((ssp_table2$L1==0)&(ssp_table2$`L2/3`==0)&(ssp_table2$L4==0)&(ssp_table2$L5==0)&(ssp_table2$L6a==0)&(ssp_table2$L6b==0))),]

ssp_ctype<-list()
ssp_DL<-list()
ssp_ptype<-list()
ssp_layer<-list()
for(i in 1:nrow(ssp_table2)){
  name<-ssp_table2$name[i]
  line<-table[which(table$Name==name),]
  ssp_ctype[i]<-line$Soma_region[1]
  ssp_ptype[i]<-line$Subclass_or_type[1]
  ssp_layer[i]<-line$Cortical_layer[1]
  arborID<-ssp_table2$arbor_id[i]
  line2<-LD_order[which((LD_order$name==name)&(LD_order$arbor_id==arborID)),]
  ssp_DL[i]<-line2$LD[1]
}
## 统计subtype
ssp<-c("SSp-bfd","SSp-ll","SSp-m","SSp-n","SSp-ul","SSp-un")
L<-c("2/3","4","5","6a","6b")
ssp_subtype_table<-matrix(0,nrow = length(ssp)*length(L),ncol=6)
N<-1
del_lines<-list()
for(i in ssp){
  for (j in L) {
    lines<-ssp_table2[which((ssp_ctype==i)&(ssp_layer==j)),]
    if(nrow(lines)<1){
      for(k in 1:6){
        ssp_subtype_table[N,k]<-0
      }
      del_lines<-c(del_lines,N)
      N<-N+1
      next
    }
    m<-colMeans(lines[3:8])
    ssp_subtype_table[N,]<-m
    s<-sum(m)
    if(s==0){
      del_lines<-c(del_lines,N)
    }
    N<-N+1
  }
}

rev_del<-rev(del_lines)
ssp_subtype_table2<-ssp_subtype_table
for(i in rev_del){
  ssp_subtype_table2<-ssp_subtype_table2[-i,]
}

## standardize
ssp_subtype_table3<-ssp_subtype_table2
for(i in 1:17){
  ssp_subtype_table3[i,]<-(ssp_subtype_table3[i,]-min(ssp_subtype_table3[i,]))/(max(ssp_subtype_table3[i,])-min(ssp_subtype_table3[i,]))
}
ssp_subtype_table3<-ssp_subtype_table3[,-6]
colnames(ssp_subtype_table3)<-c("L2/3","L4","L5","L6a","L6b")
rownames(ssp_subtype_table3)<-c("SSp-bfd-L2/3","SSp-bfd-L4","SSp-bfd-L5","SSp-bfd-L6a",
                                "SSp-ll-L2/3","SSp-ll-L5","SSp-ll-L6a",
                                "SSp-m-L2/3","SSp-m-L4","SSp-m-L5",
                                "SSp-n-L2/3","SSp-n-L4","SSp-n-L5",
                                "SSp-ul-L2/3","SSp-ul-L4","SSp-ul-L5",
                                "SSp-un-L5")

ssp_subtype_R<-c("SSp-un","SSp-bfd","SSp-n","SSp-bfd",
                 "SSp-m","SSp-n","SSp-m",
                 "SSp-ll","SSp-ul","SSp-ul","SSp-m","SSp-bfd",
                 "SSp-n","SSp-bfd","SSp-ul",
                 "SSp-ll","SSp-n")
ssp_subtype_R<-rev(ssp_subtype_R)
ssp_subtype_R<-data.frame(ssp_subtype_R)
rownames(ssp_subtype_R)<-c("SSp-n\n-L2/3","SSp-ll\n-L2/3","SSp-ul\n-L4","SSp-bfd\n-L4",
                           "SSp-n\n-L4","SSp-bfd\n-L2/3","SSp-m\n-L4",
                           "SSp-ul\n-L2/3","SSp-ul\n-L5","SSp-ll\n-L5",
                           "SSp-m\n-L5","SSp-ll\n-L6a","SSp-m\n-L2/3",
                           "SSp-bfd\n-L5","SSp-n\n-L5","SSp-bfd\n-L6a","SSp-un\n-L5")


## figure
circos.clear() 
circos.par(gap.after=c(50))
col_r=structure(brewer.pal(6,"Set1"),names=ssp)
circos.heatmap(ssp_subtype_R,col = col_r,track.height=0.01,
               rownames.side="outside",rownames.col="black",rownames.cex=1,rownames.font=2,
               )
col_p=colorRamp2(c(-0.1,0.5,1.1),c("lightseagreen","white", "violet"))
circos.heatmap(ssp_subtype_table3,col = col_p,track.height = 0.4,
               dend.side="inside", cluster = TRUE,dend.track.height=0.2,
               #rownames.side="outside",rownames.col="black",rownames.cex=0.5,rownames.font=2,
               dend.callback=function(dend,m,si){color_branches(dend,k=8,col=1:8)})

lg=Legend(title="Projection Strength",col_fun=col_p,direction = c("horizontal"),
          grid_width = unit(4, "cm"),legend_height = unit(2, "cm"))
circle_size= unit(0.07,"snpc")
grid.draw(lg)
