
library('readxl')
library('dplyr')
library('reshape2')
library('gplots')
library('grid')


setwd('/Users/zagidull/Documents/fimm_files/fcs files/heatmaps_phosphoflow')
filename <- 'Phosphoflow on MKs New -3_for_abstract.xls'

keeps <- c("IC$WellID","Cell_Type","Statistic")
Cell_Type <- c('MK+/CD110+','MK+','IMK+/CD110+','IMK+') # this has to be consistent with the gating hierarchy in the excel file
keys <- c('A01','B01','C01','D01',
          'A02','B02',
          'E02')
values <-c('SFEM','SFEM','SFEM+TPO','SFEM+EPAG',
           'SFEM','SFEM',
           'Unstained'
)

#creates a list of key:value pairs. Where is key is WellID and value is Condition
mylist <- list()
for (i in seq_along(keys)) {
  mylist[keys[i]] <- values[i]
}

# greps the data from the excel output
testF <- function(x) {
  if (x == 'pAKT') {
    pat = paste0('Mean\\s:\\sVL1|\\/',x) 
    pat = paste0(pat,'$')
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) # get "Mean : VL1" anywhere in a string or "/pAKT" in the end of a string in "Name" column 
    return(out)
  } else if (x == 'pERK') {
    pat = paste0('Mean\\s:\\sBL2|\\/',x)
    pat = paste0(pat,'$')
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) 
    return(out)  
  } else if (x == 'pSTAT3') {
    pat = paste0('Mean\\s:\\sBL3|\\/',x)
    pat = paste0(pat,'$')
    out <- filter(new_PF, grepl(pattern = pat, new_PF$Name)) 
    return(out)
  } else if (x == 'pSTAT5') {
    pat = paste0('Mean\\s:\\sRL1|\\/',x)
    pat = paste0(pat,'$')
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

testF('pAKT') -> pAKT_raw #extract mean intesity values and laser channels
pAKT_raw <- pAKT_raw[seq(2, nrow(pAKT_raw), 2), ] #leave only intensity values
pAKT_raw$Cell_Type <- Cell_Type #apply cell types
pAKT_raw[keeps] -> pAKT_hm #keep only desired columns
as.data.frame(pAKT_hm) -> pAKT_hm
pAKT_hm[,3] <- as.numeric(pAKT_hm[,3]) #since intensity values where chr
monst(pAKT_hm) -> out_final #prepare conditions according to plate layout
pAKT_hm$`IC$WellID` <- out_final #apply conditions
colnames(pAKT_hm)[1] <- 'Condition' 
rem <- with(pAKT_hm, which(pAKT_hm$Condition == 'Unstained')) #prep to remove unstained
rem1 <- with(pAKT_hm, which(is.na(pAKT_hm$Condition))) #prep to remove when NA
pAKT_hm <- pAKT_hm[-c(rem,rem1),] #remove unstained and NA
dcast(pAKT_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> pAKT_ready #get mean values of wells with matching conditions
mat <- matrix(,nrow=4,ncol=15)
colnames(mat) <- colnames(pAKT_ready)[-1]
mat <- pAKT_ready[,-1]
rownames(mat) <- pAKT_ready[,1]
mat[,c(13,14,15,1,2,3,4,5,6,7,8,9,10,11,12)] -> mat #rearrange so that SFEMII and SFEMII-EPO and SFEMII-TPO are in front
mat <- log10(mat) - log10(mat[,1]) # log transform to get fold change
pAKT <- mat[-7,] # get rid of -MK

testF('pERK') -> pERK_raw
pERK_raw <- pERK_raw[seq(2, nrow(pERK_raw), 2), ]
pERK_raw$Cell_Type <- Cell_Type
pERK_raw[keeps] -> pERK_hm
as.data.frame(pERK_hm) -> pERK_hm
pERK_hm[,3] <- as.numeric(pERK_hm[,3]) #since intensity values where chr
monst(pERK_hm) -> out_final1 #prepare conditions according to plate layout
pERK_hm$`IC$WellID` <- out_final1 #apply conditions
colnames(pERK_hm)[1] <- 'Condition' 
rem2 <- with(pERK_hm, which(pERK_hm$Condition == 'Unstained')) #prep to remove unstained
rem3 <- with(pERK_hm, which(is.na(pERK_hm$Condition))) #prep to remove when NA
pERK_hm <- pERK_hm[-c(rem2,rem3),] #remove unstained and na
dcast(pERK_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> pERK_ready #get mean values of matching wells
mat1 <- matrix(,nrow=4,ncol=15)
colnames(mat1) <- colnames(pERK_ready)[-1]
mat1 <- pERK_ready[,-1]
rownames(mat1) <- pERK_ready[,1]
mat1[,c(13,14,15,1,2,3,4,5,6,7,8,9,10,11,12)] -> mat1 #rearrange so that SFEMII and SFEMII-EPO and SFEMII-TPO are in front
mat1 <- log10(mat1) - log10(mat1[,1]) # log transform to get fold change
pERK <- mat1[-7,] # get rid of -MK

testF('pSTAT3') -> pSTAT3_raw
pSTAT3_raw <- pSTAT3_raw[seq(2, nrow(pSTAT3_raw), 2), ]
pSTAT3_raw$Cell_Type <- Cell_Type
pSTAT3_raw[keeps] -> pSTAT3_hm
as.data.frame(pSTAT3_hm) -> pSTAT3_hm
pSTAT3_hm[,3] <- as.numeric(pSTAT3_hm[,3]) #since intensity values where chr
monst(pSTAT3_hm) -> out_final2 #prepare conditions according to plate layout
pSTAT3_hm$`IC$WellID` <- out_final2 #apply conditions
colnames(pSTAT3_hm)[1] <- 'Condition' 
rem4 <- with(pSTAT3_hm, which(pSTAT3_hm$Condition == 'Unstained')) #prep to remove unstained
rem5 <- with(pSTAT3_hm, which(is.na(pSTAT3_hm$Condition))) #prep to remove when NA
pSTAT3_hm <- pSTAT3_hm[-c(rem4,rem5),] #remove unstained and na
dcast(pSTAT3_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> pSTAT3_ready #get mean values of matching wells
mat2 <- matrix(,nrow=4,ncol=15)
colnames(mat2) <- colnames(pSTAT3_ready)[-1]
mat2 <- pSTAT3_ready[,-1]
rownames(mat2) <- pSTAT3_ready[,1]
mat2[,c(13,14,15,1,2,3,4,5,6,7,8,9,10,11,12)] -> mat2 #rearrange so that SFEMII and SFEMII-EPO and SFEMII-TPO are in front
mat2 <- log10(mat2) - log10(mat2[,1]) # log transform to get fold change
pSTAT3 <- mat2[-7,] # get rid of -MK

testF('pSTAT5') -> pSTAT5_raw
pSTAT5_raw <- pSTAT5_raw[seq(2, nrow(pSTAT5_raw), 2), ]
pSTAT5_raw$Cell_Type <- Cell_Type
pSTAT5_raw[keeps] -> pSTAT5_hm
as.data.frame(pSTAT5_hm) -> pSTAT5_hm
pSTAT5_hm[,3] <- as.numeric(pSTAT5_hm[,3]) #since intensity values where chr
monst(pSTAT5_hm) -> out_final3 #prepare conditions according to plate layout
pSTAT5_hm$`IC$WellID` <- out_final3 #apply conditions
colnames(pSTAT5_hm)[1] <- 'Condition' 
rem6 <- with(pSTAT5_hm, which(pSTAT5_hm$Condition == 'Unstained')) #prep to remove unstained
rem7 <- with(pSTAT5_hm, which(is.na(pSTAT5_hm$Condition))) #prep to remove when NA
pSTAT5_hm <- pSTAT5_hm[-c(rem6,rem7),] #remove unstained and na
dcast(pSTAT5_hm, Cell_Type ~ Condition, fun.aggregate = mean) -> pSTAT5_ready #get mean values of matching wells
mat3 <- matrix(,nrow=4,ncol=15)
colnames(mat3) <- colnames(pSTAT5_ready)[-1]
mat3 <- pSTAT5_ready[,-1]
rownames(mat3) <- pSTAT5_ready[,1]
mat3[,c(13,14,15,1,2,3,4,5,6,7,8,9,10,11,12)] -> mat3 #rearrange so that SFEMII and SFEMII-EPO and SFEMII-TPO are in front
mat3 <- log10(mat3) - log10(mat3[,1]) # log transform to get fold change
pSTAT5 <- mat3[-7,] # get rid of -MK

#colour settings. To change colour scheme please modify colour names in colorRampPalette() only
breaks = c(seq(-1,-0.1,length=30),seq(-0.09,0.09,length=30),seq(0.10,1,length=30))
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = length(breaks)-1)

#library(gridGraphics)
#grab_grob <- function(){
#  grid.echo()
#  grid.grab()
#}


#grid.newpage()

# library(gridExtra)
# grid.arrange(g,g, ncol=2, clip=TRUE)

#lay <- grid.layout(nrow = 2, ncol=2)
#pushViewport(viewport(layout = lay))

#png(file = "heatmap2.png")
#heatmap.2(as.matrix(pAKT),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
#          margins = c(15,15), col = my_palette, breaks = breaks, main='pAKT fold change vs SFEMII',na.color = 'magenta', srtCol=60,cexRow=0.75, cexCol = 0.75)
#g <- grab_grob()
#grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1, 
#                                  layout.pos.col = 1, clip=TRUE)))
#heatmap.2(as.matrix(pERK),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
#          margins = c(15,15), col = my_palette, breaks = breaks, main='pERK fold change vs SFEMII',na.color = 'magenta', srtCol=60,cexRow=0.75, cexCol = 0.75)
#g <- grab_grob()
#grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1, 
#                                  layout.pos.col = 2, clip=TRUE)))
#heatmap.2(as.matrix(pSTAT3),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
#          margins = c(15,15), col = my_palette, breaks = breaks, main='pSTAT3 fold change vs SFEMII',na.color = 'magenta', srtCol=60,cexRow=0.75, cexCol = 0.75)
#g <- grab_grob()
#grid.draw(editGrob(g, vp=viewport(layout.pos.row = 2, 
#                                  layout.pos.col = 1, clip=TRUE)))
#heatmap.2(as.matrix(pSTAT5),density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', 
#          margins = c(15,15), col = my_palette, breaks = breaks, main='pSTAT5 fold change vs SFEMII',na.color = 'magenta', srtCol=60,cexRow=0.75, cexCol = 0.75)
#g <- grab_grob()
#grid.draw(editGrob(g, vp=viewport(layout.pos.row = 2, 
#                                  layout.pos.col = 2, clip=TRUE)))
#heatmap.2(as.matrix(mtcars))

#upViewport(1)


library(gridGraphics) 

grid.newpage() 


pushViewport(viewport(0, 0.5, 0.5, 0.5, just=c("left", "bottom"))) 
grid.echo(function() { 
  heatmap.2(as.matrix(pAKT),key.title = 'nothing"',key.xlab = "Fold Change",key.ylab = "blah", density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', margins = c(10,10), col = my_palette, breaks = breaks,na.color = 'magenta', main='pAKT',cexRow=0.85, cexCol = 0.85,srtCol=0)
}, newpage=FALSE) 
popViewport() 

pushViewport(viewport(0.5, 0.5, 0.5, 0.5, just=c("left", "bottom"))) 
grid.echo(function() { 
  heatmap.2(as.matrix(pERK), key = F,density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none',margins = c(10,10), col = my_palette, breaks = breaks,na.color = 'magenta', main='pERK',srtCol=0,cexRow=0.85, cexCol = 0.85)
}, newpage=FALSE) 
popViewport()

pushViewport(viewport(0, 0, 0.5, 0.5, just=c("left", "bottom"))) 
grid.echo(function() { 
  heatmap.2(as.matrix(pSTAT3),key = F,density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none',margins = c(10,10), col = my_palette, breaks = breaks, na.color = 'magenta',main='pSTAT3', srtCol=0,cexRow=0.85, cexCol = 0.85)
}, newpage=FALSE) 
popViewport()

pushViewport(viewport(0.5, 0, 0.5, 0.5, just=c("left", "bottom"))) 
grid.echo(function() { 
  heatmap.2(as.matrix(pSTAT5), key = F,density.info="none", trace="none",Colv="NA",Rowv = 'NA', dendrogram = 'none', margins = c(10,10), col = my_palette, breaks = breaks,na.color = 'magenta', main='pSTAT5', srtCol=0,cexRow=0.85, cexCol = 0.85)
}, newpage=FALSE) 
popViewport()