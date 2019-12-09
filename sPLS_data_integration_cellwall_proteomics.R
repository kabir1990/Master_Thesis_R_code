setwd("D:/LIST_DOCS/Single_omics_data_integration_for_thesis/Images for sPLS data integration")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("mixOmics")

library(mixOmics)
#library(ggplot2)
#install.packages("extrafont")
#library(extrafont)

data<- read.csv("Transcriptomedata.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = 1)
X<-t(data)

data2<- read.csv("EGTA.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",")

row.names(data2)<-make.names(data2[,1],TRUE)
data3<-data2[,-c(1,6)]

Y<-t(data3)

######sPLS analysis Transcirptome vs Metabolomics ##########

head(cbind(rownames(X), rownames(Y)))

cratero.spls <- spls(X, Y, ncomp =3, keepX = c(200,200,200),
                     keepY= c(20,20,20), mode = "regression")

##### plots###

#cormetaplot=plotVar(cratero.spls, comp =1:2, cex = c(1.5, 3))
# 
# tiff(filename = "cormetaplot.tif",
#      width = 6.7, height = 5.9, type = c("windows"), units = "in",
#      bg = "white", res = 600, family = "sans", restoreConsole = TRUE)

plotVar(cratero.spls, comp =1:2, cex = c(1.5, 3))

dev.off()

###### sPLS analysis Transcirptome vs CaCl2 fracitons ##########

head(cbind(rownames(X), rownames(Y)))

cratero.spls <- spls(X, Y, ncomp =3, keepX = c(200,200,200),
                   keepY= c(20,20,20), mode = "regression")

##################################################################
#### sPLS analysis Transcirptome vs LiCl fracitons ###########

data4<- read.csv("Transcriptomedata.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = 1)
A<-t(data)

data5<- read.csv("LiCl - spls.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",")

row.names(data2)<-make.names(data2[,1],TRUE)
data6<-data2[,-c(1,6)]

B<-t(data4)

## variableplots
head(cbind(rownames(A), rownames(B)))

cratero.spls <- spls(A, B, ncomp =3, keepX = c(200,200,200),
                     keepY= c(50,50,50), mode = "regression")


### Variable plots
plotVar(cratero.spls, comp =1:2, cex = c(1.5, 3))

##################################################################
#### sPLS analysis Transcirptome vs EGTA fracitons ###############

data7<- read.csv("Transcriptomedata.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = 1)
C<-t(data)

data8<- read.csv("EGTA - SPLS.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",")

row.names(data2)<-make.names(data2[,1],TRUE)
data9<-data2[,-c(1,6)]

D<-t(data4)

## variableplots
head(cbind(rownames(C), rownames(D)))

cratero.spls <- spls(C, D, ncomp =3, keepX = c(200,200,200),
                     keepY= c(50,50,50), mode = "regression")

### Variable plots
plotVar(cratero.spls, comp =1:2, cex = c(1, 3))






