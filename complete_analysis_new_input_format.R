#########################################################################################
# Copyright (C) 2018 {Bulat Zagidullin} {bulat.zagidullin@helsinki.fi}; {Jing Tang} {jing.tang@helsinki.fi}
# MIT License
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#  
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Version: 0.2a
#
# This file is part of project Visualization of Flow-cytometry Drug Screening Data.
# Network Pharmacology for Precision Medicine Group
# https://www.fimm.fi/en/research/groups/tang
# Institute for Molecular Medicine Finland, University of Helsinki
#
# Short Description: this Software is used to generate and draw heatmaps of Phosphoflow Flow Cytometry experiments.
# 
# Manual:
# Input is one FlowJo-generated excel (.xls) file.
# Input files should contain an already applied gating strategy. 
# Gating should be consistent for all WellIDs.
# "Name" column should contain gating hierarchy in form of, e.g. "IntelliCyt_iQue3.fcs/Cells/Singlets/Q2: BL1-A+ , VL6-A+/CD110+/pAKT".
# "Name" column should contain mean fluorescent intensity values in form of, e.g. "Mean : VL1-A = 92428.3515625".
# "Statistic" column should contain numerical values of the mean fluorescent intensity or "n/a", if not available.
# "IC$WellID" column should contain a valid WellID value in *EACH* row.
# Normalization of mean fluorescent intensity values is done using Log10 method, such that fold change is log10(x) ??? log10(control) 
# as per https://my.vanderbilt.edu/irishlab/protocols/scales-and-transformation/ visited on 12/07/2018.
#
# 1. Before using the script it is necessary to input correct setwd() and filename
# 2. To modify the plate layout or Cell_Types identified in an experiment, please modify "Cell_type", "keys" and "values" variables 
# 3. To modify any vizualization settings please refer to heatmap.2 function and its help page
# 4. Plots and csv with fold change numerical values are saved in your getwd() directory
#
# ToDo's:
# i dunno, can't be fucked right now
# would be nice to have a correct extension of the file when saving
# would be nice to have a function that runs the data transform
# would be nice to use ggplot for plotting, not this weird heatmap2 contraption. Urgh
# would be freaking amazing to not have manually input the keys and values
######################################################################################### 
library('readxl')
library('dplyr')
library('reshape2')
library('gplots')
setwd('/Users/zagidull/Documents/fimm_files/fcs files/heatmaps_phosphoflow')
filename <- 'Sample 3 Phosphoflow on MKs(EP,TPO, and drugs-2.xls'

keeps <- c("IC$WellID","Cell_Type","Statistic")
Cell_Type <- c('MK+','MK+/CD110+','IMK+','IMK+/CD110+','MK-') # this has to be consistent with the gating hierarchy in the excel file
keys <- c('A01','B01','C01','D01','E01','F01','G01','H01',
          'A02','B02','C02','D02','E02','F02','G02','H02',
          'A03','B03','C03','D03','E03','F03','G03','H03',
          'A04','B04','C04','D04','E04','F04','G04','H04',
          'A05','B05','C05','D05','E05','F05','G05','H05')
values <-c('Control','TPO 100ng','EPAG 1.8uM','TPO+EPAG','Control','TPO 100ng','EPAG 1.8uM','TPO+EPAG',
           'Veneto 20nM','Veneto 40nM','Veneto 20nM + EPAG','Veneto 40nM + EPAG','Veneto 20nM','Veneto 40nM','Veneto 20nM + EPAG','Veneto 40nM + EPAG',
           'ETOPO 100nM','ETOPO 200nM','ETOPO 100nM + EPAG','ETOPO 200nM + EPAG','ETOPO 100nM','ETOPO 200nM','ETOPO 100nM + EPAG','ETOPO 200nM + EPAG',
           'Mido 20nM','Mido 40nM','Mido 20nM + EPAG', 'Mido 40nM + EPAG', 'Mido 20nM','Mido 40nM','Mido 20nM + EPAG', 'Mido 40nM + EPAG',
           'Navito 20nM','Navito 40nM','Navito 20nM + EPAG','Navito 40nM + EPAG','Navito 20nM','Navito 40nM','Navito 20nM + EPAG','Navito 40nM + EPAG'
           )

#creates a list of key:value pairs. Where is key is WellID and value is Condition
mylist <- list()
for (i in seq_along(keys)) {
  mylist[keys[i]] <- values[i]
}

# greps the data from the excel output
mod_testF <- function(x) {
  if (x == 'pAKT') {
    pat <- 'Median\\s:\\sVL1'
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) # get "Mean : VL1" anywhere in a string or "/pAKT" in the end of a string in "Name" column 
    return(out)
  } else if (x == 'pERK') {
    pat <- 'Median\\s:\\sBL2'
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) 
    return(out)  
  } else if (x == 'pSTAT3') {
    pat <- 'Median\\s:\\sBL3'
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) 
    return(out)
  } else if (x == 'pSTAT5') {
    pat <- 'Median\\s:\\sRL1'
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) 
    return(out)      
  } else {
    print('Wrong downstream signal chosen! Please check again')
  }
  return(out)
}

# prepares conditions based on WellID
monst <- function(x) {
  out_final = character()
  un <- unique(x$`IC$WellID`)
  for (i in seq_along(un)) {
    wellID = un[i]
    data_subset <- x[x[,1] == wellID,]
    out <- rep(as.character(mylist[wellID]),nrow(data_subset))
    out_final <- c(out_final,out)
  }
  return(out_final)
}

new_PF <- read_xls(filename)

# testing
#pAKT
mod_testF('pAKT') -> temp.pAKT_raw
temp.pAKT_raw$Cell_Type <- Cell_Type
temp.pAKT_raw[keeps] -> temp.pAKT_hm
as.data.frame(temp.pAKT_hm) -> temp.pAKT_hm
temp.pAKT_hm[,3] <- as.numeric(temp.pAKT_hm[,3])
monst(temp.pAKT_hm) -> temp.out_final
temp.pAKT_hm$`IC$WellID` <- temp.out_final
colnames(temp.pAKT_hm)[1] <- 'Condition' 
        #relevant when we need to remove something
        #rem <- with(pAKT_hm, which(pAKT_hm$Condition == 'Unstained')) #prep to remove unstained
        #rem1 <- with(pAKT_hm, which(is.na(pAKT_hm$Condition)))
        #pAKT_hm <- pAKT_hm[-c(rem,rem1),]
dcast(temp.pAKT_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> temp.pAKT_ready
temp.mat <- matrix(,nrow=length(rownames(temp.pAKT_ready)),ncol=length(colnames(temp.pAKT_ready))-1) #create matrix of the same dimension as the one from the input, minus one column, as it is cell-type
colnames(temp.mat) <- colnames(temp.pAKT_ready)[-1] #first colnames is cell-type we dont need taht, hence we use [-1] index
temp.mat <- temp.pAKT_ready[,-1] #copy over the values into newly created mat, minus the cell type column 
rownames(temp.mat) <- temp.pAKT_ready[,1] #apply first matrix as a column
temp.mat <- log10(temp.mat) - log10(temp.mat[,1]) 
pAKT <- temp.mat
rm(list=ls(pattern='^temp'))
#pERK
mod_testF('pERK') -> temp.pERK_raw
temp.pERK_raw$Cell_Type <- Cell_Type
temp.pERK_raw[keeps] -> temp.pERK_hm
as.data.frame(temp.pERK_hm) -> temp.pERK_hm
temp.pERK_hm[,3] <- as.numeric(temp.pERK_hm[,3])
monst(temp.pERK_hm) -> temp.out_final
temp.pERK_hm$`IC$WellID` <- temp.out_final
colnames(temp.pERK_hm)[1] <- 'Condition' 
#relevant when we need to remove something
#rem <- with(pAKT_hm, which(pAKT_hm$Condition == 'Unstained')) #prep to remove unstained
#rem1 <- with(pAKT_hm, which(is.na(pAKT_hm$Condition)))
#pAKT_hm <- pAKT_hm[-c(rem,rem1),]
dcast(temp.pERK_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> temp.pERK_ready
temp.mat <- matrix(,nrow=length(rownames(temp.pERK_ready)),ncol=length(colnames(temp.pERK_ready))-1) #create matrix of the same dimension as the one from the input, minus one column, as it is cell-type
colnames(temp.mat) <- colnames(temp.pERK_ready)[-1] #first colnames is cell-type we dont need taht, hence we use [-1] index
temp.mat <- temp.pERK_ready[,-1] #copy over the values into newly created mat, minus the cell type column 
rownames(temp.mat) <- temp.pERK_ready[,1] #apply first matrix as a column
temp.mat <- log10(temp.mat) - log10(temp.mat[,1]) 
pERK <- temp.mat
rm(list=ls(pattern='^temp'))
#pSTAT3
mod_testF('pSTAT3') -> temp.pSTAT3_raw
temp.pSTAT3_raw$Cell_Type <- Cell_Type
temp.pSTAT3_raw[keeps] -> temp.pSTAT3_hm
as.data.frame(temp.pSTAT3_hm) -> temp.pSTAT3_hm
temp.pSTAT3_hm[,3] <- as.numeric(temp.pSTAT3_hm[,3])
monst(temp.pSTAT3_hm) -> temp.out_final
temp.pSTAT3_hm$`IC$WellID` <- temp.out_final
colnames(temp.pSTAT3_hm)[1] <- 'Condition' 
#relevant when we need to remove something
#rem <- with(pAKT_hm, which(pAKT_hm$Condition == 'Unstained')) #prep to remove unstained
#rem1 <- with(pAKT_hm, which(is.na(pAKT_hm$Condition)))
#pAKT_hm <- pAKT_hm[-c(rem,rem1),]
dcast(temp.pSTAT3_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> temp.pSTAT3_ready
temp.mat <- matrix(,nrow=length(rownames(temp.pSTAT3_ready)),ncol=length(colnames(temp.pSTAT3_ready))-1) #create matrix of the same dimension as the one from the input, minus one column, as it is cell-type
colnames(temp.mat) <- colnames(temp.pSTAT3_ready)[-1] #first colnames is cell-type we dont need taht, hence we use [-1] index
temp.mat <- temp.pSTAT3_ready[,-1] #copy over the values into newly created mat, minus the cell type column 
rownames(temp.mat) <- temp.pSTAT3_ready[,1] #apply first matrix as a column
temp.mat <- log10(temp.mat) - log10(temp.mat[,1]) 
pSTAT3 <- temp.mat
rm(list=ls(pattern='^temp'))
#STAT5
mod_testF('pSTAT5') -> temp.pSTAT5_raw
temp.pSTAT5_raw$Cell_Type <- Cell_Type
temp.pSTAT5_raw[keeps] -> temp.pSTAT5_hm
as.data.frame(temp.pSTAT5_hm) -> temp.pSTAT5_hm
temp.pSTAT5_hm[,3] <- as.numeric(temp.pSTAT5_hm[,3])
monst(temp.pSTAT5_hm) -> temp.out_final
temp.pSTAT5_hm$`IC$WellID` <- temp.out_final
colnames(temp.pSTAT5_hm)[1] <- 'Condition' 
#relevant when we need to remove something
#rem <- with(pAKT_hm, which(pAKT_hm$Condition == 'Unstained')) #prep to remove unstained
#rem1 <- with(pAKT_hm, which(is.na(pAKT_hm$Condition)))
#pAKT_hm <- pAKT_hm[-c(rem,rem1),]
dcast(temp.pSTAT5_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> temp.pSTAT5_ready
temp.mat <- matrix(,nrow=length(rownames(temp.pSTAT5_ready)),ncol=length(colnames(temp.pSTAT5_ready))-1) #create matrix of the same dimension as the one from the input, minus one column, as it is cell-type
colnames(temp.mat) <- colnames(temp.pSTAT5_ready)[-1] #first colnames is cell-type we dont need taht, hence we use [-1] index
temp.mat <- temp.pSTAT5_ready[,-1] #copy over the values into newly created mat, minus the cell type column 
rownames(temp.mat) <- temp.pSTAT5_ready[,1] #apply first matrix as a column
temp.mat <- log10(temp.mat) - log10(temp.mat[,1]) 
pSTAT5 <- temp.mat
rm(list=ls(pattern='^temp'))
# end testing

#colour settings. To change colour scheme please modify colour names in colorRampPalette() only
breaks = c(seq(-1,-0.08,length=50),seq(-0.07,0.07,length=50),seq(0.08,1,length=50))
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = length(breaks)-1)


#png(file = "heatmap2.png")
heatmap.2(as.matrix(pAKT),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
          margins = c(15,15), col = my_palette, breaks = breaks, main='pAKT fold change vs SFEMII',na.color = 'magenta', srtCol=60)
heatmap.2(as.matrix(pERK),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
          margins = c(15,15), col = my_palette, breaks = breaks, main='pERK fold change vs SFEMII',na.color = 'magenta', srtCol=60)
heatmap.2(as.matrix(pSTAT3),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
          margins = c(15,15), col = my_palette, breaks = breaks, main='pSTAT3 fold change vs SFEMII',na.color = 'magenta', srtCol=60)
heatmap.2(as.matrix(pSTAT5),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
          margins = c(15,15), col = my_palette, breaks = breaks, main='pSTAT5 fold change vs SFEMII',na.color = 'magenta', srtCol=60)
#dev.off()

# saving file to your getwd() with the name original name, but with current time and date prepended, it is comma-delimited. When opening with excel it might give some security warning - ignore. So, use text-to-columns in excel to view it correctly
mat_out <- list()
mat_out[["pAKT"]] <- pAKT
mat_out[["pERK"]] <- pERK
mat_out[["pSTAT3"]] <- pSTAT3
mat_out[["pSTAT5"]] <- pSTAT5

current.Time <- Sys.time() 
csv.filename <- paste(current.Time,filename,sep=" ") 

out_file <- file(csv.filename, open="a")
for (i in seq_along(mat_out)){
  write.table(names(mat_out)[i], file=out_file, sep=",", dec=".", col.names=FALSE, row.names=FALSE)
  write.table(mat_out[[i]], file=out_file, sep=",", dec=".", quote=FALSE, col.names=NA, row.names=TRUE)  
}
close(out_file)
