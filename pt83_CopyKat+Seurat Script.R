library(Seurat)
library(ggplot2)
library(gplots)
data_dir <- '/Users/Pancho/STEMREM Project/83_RAW'
list.files(data_dir)
library(Seurat)
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix)
library(Seurat)
raw83 <- Read10X(data.dir = '/Users/Pancho/STEMREM Project/83_RAW')
raw83 <- CreateSeuratObject(counts = raw83, project = "copykat83", min.cells = 0, min.features = 0)
exp.raw83data <- as.matrix(raw83@assays$RNA@counts)
# Save the matrix for future usage

write.table(exp.raw83data, file="exp.rawdata83.txt", sep="\t", quote = FALSE, row.names = TRUE)

library(copykat)
copykat83 <- copykat(rawmat=exp.raw83data, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")

pred.test <- data.frame(copykat83$prediction)
CNA.test <- data.frame(copykat83$CNAmat)
# you should get a copykat_prediction.txt file that has what we need to input as a column in Seurat

library(dplyr)
library(Seurat)
library(patchwork)
library(Cairo)

#load patients data, returns UMI count matrix
data_dir <- '/Users/Pancho/STEMREM Project/83_RAW'
tnbc83.data <- Read10X(data.dir = '/Users/Pancho/STEMREM Project/83_RAW')
#make a seurat object with raw data, though i think this has already been
#normalized by authors
tnbc83 <- CreateSeuratObject(counts = tnbc83.data, project = "tnbcpt83", min.cells = 3, min.features = 200)
#tnbc83 = An object of class Seurat
#33538 features (genes) across 1065 samples (individual cells) within 1 assay 
#Active assay: RNA (33528 features, 0 variable features)

# upload prediction data from copykat
copykat83_prediction <- read.delim("/Users/Pancho/pt83_copykat_prediction.txt")
str(copykat83_prediction)

 
head(tnbc83$nFeature_RNA)

# tmp <- copykat83_prediction$copykat.pred
# names(tmp) <- copykat83_prediction$cell.names

preds2 <- sapply(names(tnbc83$orig.ident),
                 function(ID){
                   pred <- copykat83_prediction$copykat.pred[copykat83_prediction$cell.names == ID]
                   if (length(pred) != 0){
                     return(pred)
                   } else {
                     return(NA)
                   }
                 })


tnbc83_ck <- AddMetaData(object = tnbc83, metadata = preds2, col.name = "Prediction")
head(tnbc_ck)
tnbc83_a <- subset(tnbc83_ck, subset = Prediction == 'aneuploid')
head(tnbc83_a)

# this is all i have up to - i am consistently getting an error in trying to add the prediction column
# i was thinking that it was due to the preds2 chunk not working on mine but unsure why 
# the erroe seems to try to only replace 1 thing at a time, but not sure why 

