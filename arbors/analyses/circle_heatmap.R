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

### read data
axonF<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/axon_feature2.txt")
colnames(axonF)<-c('name','arbor_id','region','max_density','num_nodes','total_path_length','volume','branch','d_to_soma','d_to_soma2','hub','variance_ratio')

table<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/1891_type.txt",header=TRUE)
ssp_layer<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/SSp_layer.txt")
colnames(ssp_layer)<-c('name','L1','L2/3','L4','L5','L6a','L6b','Non-cortical')
ssp_arbor_layer_initial<-read.table("/Users/krystal/Desktop/project/cloud_paper/2023/SSp_arbors_layer2.txt")
colnames(ssp_arbor_layer_initial)<-c('name','L1','L2/3','L4','L5','L6a','L6b','Non-cortical')
filter_name<-read.table("/Users/krystal/Desktop/project/cloud_paper/0511/Final_names.txt")
colnames(filter_name)<-c('name')

ssp_arbor_layer<-ssp_arbor_layer_initial
keep<--1
for(i in 1:nrow(ssp_arbor_layer_initial)){
  n<-nchar(ssp_arbor_layer_initial$name[i])
  name<-substr(ssp_arbor_layer_initial$name[i],1,n-7)
  r<-filter_name[which(filter_name$name==name),"name"]
  if (length(r)!=0){
    keep<-c(keep,i)
  }
}
keep<-keep[-1]
# remove arbors with few projections
ssp_arbor_layer<-ssp_arbor_layer[keep,]


LD_order<-read.table("/Users/krystal/Desktop/project/cloud_paper/0511/category.txt",header=TRUE,sep=",")

# combine all necessary information into a united table
ssp_table<-matrix(0,nrow = nrow(ssp_arbor_layer),ncol=12)
count<-1
for(i in 1:nrow(ssp_arbor_layer)){
  if(sum(ssp_arbor_layer[i,2:7])>10){
    n<-nchar(ssp_arbor_layer$name[i])
    ssp_table[count,1]<-substr(ssp_arbor_layer$name[i],1,n-7) # swc name
    ssp_table[count,2]<-substr(ssp_arbor_layer$name[i],n,n)   # arbor id
    ssp_table[count,3]<-ssp_arbor_layer$L1[i] # axonal length in cortical laminar L1
    ssp_table[count,4]<-ssp_arbor_layer$`L2/3`[i]
    ssp_table[count,5]<-ssp_arbor_layer$L4[i]
    ssp_table[count,6]<-ssp_arbor_layer$L5[i]
    ssp_table[count,7]<-ssp_arbor_layer$L6a[i]
    ssp_table[count,8]<-ssp_arbor_layer$L6b[i]
    ssp_table[count,9]<-ssp_arbor_layer$`Non-cortical`[i]
    tmp<-LD_order[which((LD_order$name==ssp_table[count,1])&(LD_order$arbor_id)==as.numeric(ssp_table[count,2])),]
    ssp_table[count,10]<-tmp$LD[1]
    tmp2<-table[which(table$Name==ssp_table[count,1]),]
    ssp_table[count,11]<-tmp2$Soma_region[1]
    ssp_table[count,12]<-tmp2$Cortical_layer[1]
    count<-count+1
  }
}
ssp_table<-ssp_table[1:199,]

# Extract the information from previous table
ssp_table2<-as.data.frame(ssp_table)
colnames(ssp_table2)<-c('name','arbor_id','L1','L2/3','L4','L5','L6a','L6b','Non_cortical','LD','S','CL')
ssp_table2$L1<-as.numeric(ssp_table2$L1)
ssp_table2$`L2/3`<-as.numeric(ssp_table2$`L2/3`)
ssp_table2$L4<-as.numeric(ssp_table2$L4)
ssp_table2$L5<-as.numeric(ssp_table2$L5)
ssp_table2$L6a<-as.numeric(ssp_table2$L6a)
ssp_table2$L6b<-as.numeric(ssp_table2$L6b)
ssp_table2$`Non_cortical`<-as.numeric(ssp_table2$`Non_cortical`)
# Remove non-cortical layer projection items
#ssp_table2<-ssp_table2[which(!((ssp_table2$L1==0)&(ssp_table2$`L2/3`==0)&(ssp_table2$L4==0)&(ssp_table2$L5==0)&(ssp_table2$L6a==0)&(ssp_table2$L6b==0))),]



ssp_proj<-ssp_table2

ssp<-c('SSp-un','SSp-ul','SSp-n','SSp-m','SSp-ll','SSp-bbfd')
### based on stype
ssp_stype_projD<-matrix(0,nrow = 6,ncol = 14)
N<-1
for (i in ssp) {
  lines<-ssp_proj[which(ssp_proj$S==i),]
  local<-lines[which(lines$LD==0),]
  distal<-lines[which(lines$LD==1),]
  m1<-colMeans(local[3:9])
  m2<-colMeans(distal[3:9])
  ssp_stype_projD[N,]<-c(m1,m2)
  N<-N+1
}

ssp_stype_projD2<-ssp_stype_projD
maxp<-max(ssp_stype_projD2)
minp<-min(ssp_stype_projD2)
for(i in 1:6){
  ssp_stype_projD2[i,]<-(ssp_stype_projD[i,]-minp)/(maxp-minp)
}

ssp_layer_projD<-matrix(0,nrow = 4,ncol = 14)
N<-1  
L2<-c("2/3","4","5")
for (i in L2) {
  lines<-ssp_proj[which(ssp_proj$CL==i),]
  local<-lines[which(lines$LD==0),]
  distal<-lines[which(lines$LD==1),]
  m1<-colMeans(local[3:9])
  m2<-colMeans(distal[3:9])
  ssp_layer_projD[N,]<-c(m1,m2)
  N<-N+1
}
lines<-ssp_proj[which((ssp_proj$CL=="6a")|(ssp_proj$CL=="6b")|(ssp_proj$CL=="6")),]
local<-lines[which(lines$LD==0),]
distal<-lines[which(lines$LD==1),]
m1<-colMeans(local[3:9])
m2<-colMeans(distal[3:9])
ssp_layer_projD[4,]<-c(m1,m2)

#ssp_layer_projD<-ssp_layer_projD[-c(1,6),]#remove empty types stype$L1 and stype$L6b
ssp_layer_projD2<-ssp_layer_projD
minpl<-min(ssp_layer_projD2)
maxpl<-max(ssp_layer_projD2)
for(i in 1:4){
  ssp_layer_projD2[i,]<-(ssp_layer_projD[i,]-minpl)/(maxpl-minpl)
}

write.table(ssp_stype_projD2,file = "/Users/krystal/Desktop/project/cloud_paper/0511/ssp_stype_proj.txt",col.names = FALSE,row.names = FALSE)
write.table(ssp_layer_projD2,file = "/Users/krystal/Desktop/project/cloud_paper/0511/ssp_layer_proj.txt",col.names = FALSE,row.names = FALSE)


## 统计subtype
ssp<-c("SSp-bfd","SSp-ll","SSp-m","SSp-n","SSp-ul","SSp-un")
L<-c("1","2/3","4","5","6a","6b")
ssp_subtype_table<-matrix(0,nrow = length(ssp)*length(L),ncol=12)
N<-1
del_lines<-list()
for(i in ssp){
  for (j in L) {
    lines<-ssp_table2[which((ssp_table2$S==i)&(ssp_table2$CL==j)),]
    if(nrow(lines)<1){
      for(k in 1:12){
        ssp_subtype_table[N,k]<-0
      }
      del_lines<-c(del_lines,N)
      N<-N+1
      next
    }
    local<-lines[lines$LD==0,]
    distal<-lines[lines$LD==1,]
    if(nrow(local)==0){
      m1<-c(0,0,0,0,0,0)
    }
    else{
      m1<-colMeans(local[3:8])
    }
    if(nrow(distal)==0){
      m2<-c(0,0,0,0,0,0)
    }
    else{
      m2<-colMeans(distal[3:8])
    }
    ssp_subtype_table[N,]<-c(m1,m2)
    s1<-sum(m1)
    s2<-sum(m2)
    if((s1+s2)==0){
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
#max3<-max(ssp_subtype_table3)
#min3<-min(ssp_subtype_table3)
for(i in 1:16){
  max3<-max(ssp_subtype_table3[i,])
  min3<-min(ssp_subtype_table3[i,])
  ssp_subtype_table3[i,]<-(ssp_subtype_table3[i,]-min3)/(max3-min3)
}

#ssp_subtype_table3<-ssp_subtype_table3[-c(17,4),]
ssp_subtype_table3<-ssp_subtype_table3[-c(16),]
colnames(ssp_subtype_table3)<-c("L1","L2/3","L4","L5","L6a","L6b","L1","L2/3","L4","L5","L6a","L6b")
rownames(ssp_subtype_table3)<-c("SSp-bfd-L2/3","SSp-bfd-L4","SSp-bfd-L5",
                                "SSp-ll-L2/3","SSp-ll-L5","SSp-ll-L6a",
                                "SSp-m-L2/3","SSp-m-L4","SSp-m-L5",
                                "SSp-n-L2/3","SSp-n-L4","SSp-n-L5",
                                "SSp-ul-L2/3","SSp-ul-L4","SSp-ul-L5",
                                "SSp-un-L5")

ssp_subtype_R<-c("SSp-ul","SSp-bfd","SSp-un","SSp-m",
                 "SSp-bfd","SSp-n","SSp-ul",
                 "SSp-ll","SSp-ul","SSp-m","SSp-ll",
                 "SSp-n","SSp-ll","SSp-m","SSp-bfd",
                 "SSp-n")
ssp_subtype_R<-data.frame(ssp_subtype_R)
rownames(ssp_subtype_R)<-c("SSp-ul\n-L4","SSp-ll\n-L2/3","SSp-n\n-L4","SSp-m\n-L5",
                           "SSp-bfd\n-L5","SSp-n\n-L5","SSp-ul\n-L5",
                           "SSp-ll\n-L5","SSp-ul\n-L2/3",
                           "SSp-m\n-L2/3","SSp-ll\n-L6a","SSp-n\n-L4",
                           "SSp-ll\n-L2/3","SSp-m\n-L4","SSp-bfd\n-L4","SSp-n\n-L2/3")


## figure
circos.clear() 
circos.par(gap.after=c(50))
col_r=structure(brewer.pal(6,"Set1"),names=ssp)
circos.heatmap(ssp_subtype_R,col = col_r,track.height=0.03,
               rownames.side="outside",rownames.col="black",rownames.cex=0.8,rownames.font=2,
               )
col_p=colorRamp2(c(-0.1,0.5,1.1),c("lightseagreen","white", "violet"))
circos.heatmap(ssp_subtype_table3,col = col_p,track.height = 0.65,
               dend.side="inside", cluster = TRUE,dend.track.height=0.2,
               rownames.side="outside",rownames.col="black",rownames.cex=0.4,rownames.font=2,
               dend.callback=function(dend,m,si){color_branches(dend,k=8,col=1:8)})

circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 cn = colnames(ssp_subtype_table3)
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                             1:n /1.77+2, cn, 
                             cex = 1, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)

lg=Legend(title="Projection Strength",col_fun=col_p,direction = c("horizontal"),
          grid_width = unit(4, "cm"),legend_height = unit(2, "cm"))
circle_size= unit(0.07,"snpc")
grid.draw(lg)

library("lessR")
subtype_names<-c("","SSp-ul-L4","SSp-ll-L2/3","SSp-n-L4","SSp-bfd-L4",
                 "SSp-ul-L5","SSp-ll-L5","SSp-ul-L2/3",
                 "SSp-m-L2/3","SSp-ll-L6a",
                 "SSp-un-L5","SSp-bfd-L2/3","SSp-m-L5",
                 "SSp-bfd-L5","SSp-n-L5","SSp-n-L2/3","SSp-m-L4")
#subtype_names_df<-data.frame(variable=tmp)
#PieChart(variable,data = subtype_names_df,values_size = 1,labels_cex = 1.2)
PieChart(subtype_names,c(50,310/16,310/16,310/16,310/16,310/16,310/16,310/16,310/16,310/16,310/16,310/16,310/16,
                         310/16,310/16,310/16,310/16),main=NULL,labels_cex = 1.3)


### count numbers
namelist<-unique(ssp_table2$name)
count_type<-c(0,0,0,0,0,0)
count_layer<-c(0,0,0,0,0,0)
count_n<-0
for (i in namelist){
  lines<-ssp_table2[which(ssp_table2$name==i),]
  stype<-lines$S[1]
  cl<-lines$CL[1]
  ID1<-which(ssp==stype)
  count_type[ID1]<-count_type[ID1]+1
  if(stype %in% ssp){
    if (cl=="n/a"){
      count_n<-count_n+1
    }
    else if (cl=="6"){
      count_n<-count_n+1
      print(1)
    }
    else{
      ID2<-which(L==cl)
      count_layer[ID2]<-count_layer[ID2]+1
    }
  }
}

for(i in ssp){
  for(j in L){
    r<-ssp_table2[which((ssp_table2$S==i)&(ssp_table2$CL==j)),]
    if (nrow(r)==0){
      print(nrow(r))
    }
    else{
      n<-length(unique(r$name))
      print(n)
    }
  }
}

write.table(ssp_table2,file = "/Users/krystal/Desktop/project/cloud_paper/0511/ssp_table.txt",col.names = TRUE,row.names = FALSE)
