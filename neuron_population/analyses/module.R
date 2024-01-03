#set working  dir
setwd("E:\\platform-paper\\fig2modules_figs3neurites_1224\\\\")
print(getwd())

#import df
df = read.csv("middle_files\\fig2.csv")

#create corr
b_num = dim(df)[[1]]
r_num = dim(df)[[2]]-2
corr = data.frame(matrix(ncol = r_num, nrow = r_num))
colnames(corr) = colnames(df)[1:r_num+1]
rownames(corr) = colnames(df)[1:r_num+1]

#clac corr with w
library(wCorr)
w = c(df[[dim(df)[[2]]]])
for (ii in seq(1,r_num)){
  x = c(df[[ii+1]])
  for (jj in seq(1,r_num)){
    y = c(df[[jj+1]])
    
    #replace 0 with nan in x and y
    x0 = list()
    y0 = list()
    w0 = list()
    for (ib in seq(1,b_num)){
      if (x[ib]!=0 & y[ib]!=0){
        x0 = append(x0,x[ib])
        y0 = append(y0,y[ib])
        w0 = append(w0,w[ib])
      }
    }
    
    len0 = length(x0)
    if (len0>=20){
      c = weightedCorr(x0,y0,
                       method = "Spearman", #c("Pearson", "Spearman", "Polyserial", "Polychoric"),
                       weights = w0, #rep(1, length(x)),ML = FALSE,fast = TRUE
      )      
      corr[ii,jj] = c
    }
    
    }
}
write.table(corr,'middle_files\\corr_314region.csv',sep=',')

corr = read.csv('middle_files\\corr_314region.csv')


#calc hclust
new_corr = corr
rownames(new_corr) = 1:nrow(new_corr)

hc = hclust(dist(new_corr))
merge = hc[['merge']]
write.table(merge,'./merge.csv',sep=',')
order = hc[['order']]
write.table(order,'./order.csv',sep=',')

#draw initial intermodule
structure_list = fromJSON(file="middle_files\\structure_list.json")
structure_list = data.matrix(structure_list)
rownames(structure_list) = rownames(corr)
#structure_list = rbind(structure_list,structure_list)
#colnames(structure_list) = colnames(corr)
#corr = rbind(corr,structure_list)
#corr = t(corr)
#print(dim(corr))

#install.packages("dendextend","circlize","openxlsx")
library(circlize)
library(dendextend)
library(openxlsx)

#install.packages("devtools")
#library(usethis)
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(grid)
library(gridBase)
library(ComplexHeatmap)



new_corr = new_corr[hc[["order"]],]
dendlist = as.numeric(rownames(new_corr))
new_structure_list = structure_list[dendlist]
new_structure_list = data.matrix(new_structure_list)
rownames(new_structure_list) = rownames(structure_list)[dendlist]
outnames = rep(" ",length(new_structure_list))
innames = rep(" ",length(new_structure_list))
for (i in seq(1,length(new_structure_list),by=2)){
    outnames[i] = rownames(new_structure_list)[i]
    innames[i+1] = rownames(new_structure_list)[i+1]
}
  
plot.new()
circle_size = unit(1, "snpc")
circos.clear()
circos.par(gap.after=c(90),"start.degree" = 315)
col_fun0 = colorRamp2(c(313,315,354,477,519,528,549,698,703,771,803,1089,1097),
                      c("darkred","darkgreen","purple3","cyan","khaki1","orange1","lightpink",
                                 "darkseagreen1","darkslategray","violet","blue","green","red"))#"darkolivegreen3"
circos.heatmap(new_structure_list,col=col_fun0,track.height = 0.05,
                #rownames.side = "outside",
                #rownames.col = 'black',#label color
                #rownames.cex = 0.3,#label size
                cluster = FALSE
 )
                                 
plot.new()
circle_size = unit(1, "snpc")
circos.clear()
circos.par(gap.after=c(90),"start.degree" = 315)
col_fun0 = colorRamp2(c(313,315,354,477,519,528,549,698,703,771,803,1089,1097),
                      c("darkred","darkgreen","purple3","cyan","khaki1","orange1","lightpink",
                        "darkseagreen1","darkslategray","violet","blue","green","red"))#"darkolivegreen3"
rownames(new_structure_list) = outnames
circos.heatmap(new_structure_list,col=col_fun0,track.height = 0.01,
               rownames.side = "outside",
               rownames.col = 'black',#label color
               rownames.cex = 0.25,#label size
               cluster = FALSE
)

plot.new()
circle_size = unit(1, "snpc")
circos.clear()
circos.par(gap.after=c(90),"start.degree" = 315)
col_fun0 = colorRamp2(c(313,315,354,477,519,528,549,698,703,771,803,1089,1097),
                      c("darkred","darkgreen","purple3","cyan","khaki1","orange1","lightpink",
                                 "darkseagreen1","darkslategray","violet","blue","green","red"))
rownames(new_structure_list) = innames
circos.heatmap(new_structure_list,col=col_fun0,track.height = 0.01,
               rownames.side = "inside",
               rownames.col = 'black',#label color
               rownames.cex = 0.25,#label size
               rownames.font = 2,#font width size
               cluster = FALSE
)

#col_fun0 = colorRamp2(c(313,315,354,477,519,528,549,698,703,771,803,1089,1097),
#                      c("goldenrod1","lightpink","lightpink","lightpink","darkgreen","cyan","darkgreen",
#                                "goldenrod1","darkgreen","lightpink","cyan","darkgreen","lightpink"))#"darkolivegreen3"
#circos.heatmap(new_structure_list,col=col_fun0,track.height = 0.05,
#               cluster = FALSE
#)


plot.new()
circle_size = unit(3, "snpc")
circos.clear()
circos.par(gap.after=c(90),"start.degree" = 315)
#col_fun = colorRamp2(c(-1,0,1),
#                     c("darkmagenta","white","darkgreen"))
col_fun = colorRamp2(c(-1,0,1),
                     c("purple3","white","darkgreen"))
circos.heatmap(corr, col = col_fun, track.height = 0.5,
#               rownames.side = "outside",
#               rownames.col = 'black',#label color
#               rownames.cex = 0.5,#label size
               cluster = TRUE,
               dend.side = "inside",
               dend.track.height = 0.3,
               dend.callback = function(dend, m, si){
              color_branches(dend, k = 16, col = c(2,4,7,2,4,7,2,4,7,2,4,7,2,4,7,2))#
               }
)

plot.new()
circle_size = unit(3, "snpc")
circos.clear()
lg=Legend(col_fun=col_fun,direction = c("vertical"),at=c(-1,0,1))
draw(lg, x = unit(0.86, "snpc"), y = unit(0.64, "snpc"), just = "topright")

