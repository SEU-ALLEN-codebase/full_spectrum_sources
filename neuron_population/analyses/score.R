setwd("E:\\platform-paper\\fig2modules_figs3neurites_1224\\")
print(getwd())

#import df
df = read.csv("middle_files\\df_183brain_316region_density.csv")

#create corr
b_num = dim(df)[[1]]
r_num = dim(df)[[2]]-1
corr = data.frame(matrix(ncol = b_num, nrow = b_num))
colnames(corr) = colnames(df)[1:r_num+1]
rownames(corr) = colnames(df)[1:r_num+1]

#clac corr with w
library(wCorr)
for (ii in seq(1,b_num)){
  x <- c(df[ii,])
  for (jj in seq(1,b_num)){
    y <- c(df[jj,])
    
    #replace 0 with nan in x and y
    x0 = list()
    y0 = list()
    for (ib in seq(1,b_num)){
      if (x[ib+1]!=0 & y[ib+1]!=0){
        x0 = append(x0,x[ib+1])
        y0 = append(y0,y[ib+1])
      }
    }

    len0 = length(x0)
    w0 <- rep(1/len0,len0)
    if (len0>=20){
      c = weightedCorr(x0,y0,
                       method = "Spearman", #c("Pearson", "Spearman", "Polyserial", "Polychoric"),
                       weights = w0, #rep(1, length(x)),ML = FALSE,fast = TRUE
      )      
      corr[ii,jj] = c
    }
    
  }
}
write.table(corr,'middle_files\\corr_183brain_316region_density.csv',sep=',')
