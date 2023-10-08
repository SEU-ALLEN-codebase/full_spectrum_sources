setwd("D:\\temp_need\\paper_copy\\fig3distance\\")
print(getwd())

#install.packages("rjson")
print(any(grepl("rjson",installed.packages())))
library(rjson)
color_list = fromJSON(file="color_list.json")
print(color_list)

corr = read.csv("corr.csv")
corr = data.matrix(corr)
#corr = apply(corr,as.numeric)
print(corr)

rownames(corr) = colnames(corr)
print(dim(corr))

#corr = corr[,1:10]
  
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
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))

circos.clear()
circos.par(gap.after=c(30))
col_fun = colorRamp2(c(-1,0,1),c("lightseagreen","white", "violet"))
circos.heatmap(corr, col = col_fun, track.height = 0.5, rownames.side = "outside",
             rownames.col = color_list,#label color
              rownames.cex = 0.3,#label size
#               rownames.font = 1:nrow(mat1) %% 4 + 1#label font type
               cluster = TRUE,
               dend.side = "inside",
               dend.track.height = 0.4,
               dend.callback = function(dend, m, si) {
                 color_branches(dend, k = 6, col = 2:7)
               }
)

#col_r=structure(brewer.pal(6,"Set1"),names=ssp)
#circos.heatmap(ssp_subtype_R,col = col_r,track.height=0.01,
#               rownames.side="outside",rownames.col="black",rownames.cex=1,rownames.font=2,
#)

lg=Legend(col_fun=col_fun,direction = c("vertical"),at=c(1,0,-1),
          grid_width = unit(0.1, "snpc"),legend_height = unit(0.1, "snpc"))
draw(lg, x = unit(0.99, "snpc"), y = unit(0.6, "snpc"), just = "topright")

#circos.heatmap(corr,split=color_list,col=col_fun)
