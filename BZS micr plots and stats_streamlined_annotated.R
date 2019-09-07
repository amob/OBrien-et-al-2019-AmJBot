#From O'Brien et al 2019 American Journal of Botany
#Analyze 16s rDNA sequencing ASV data from QIIME2 and produce manuscript products

###load required libraries
library(Polychrome) # for less overlapping colors in barplot
library(SDMTools) #legend.gradient
library(ade4)#for pca

#recall ASVs are amplicon sequence variants, not otus.

###read in data
bzstab.s <- read.csv("~/bzsfieldtab.csv",header=T,row.names=1) 
bzstax <- read.csv("~/bzstax.csv",header=T,row.names=1,stringsAsFactors=F) 
#these 2 files are output from qiime2 analysis discussed in script. mitochondria and streptophyta (vascular plant chloroplast) sequences have been removed
#they are organized exactly the same way, i.e. abundance for row 1409 in bzstab.s is linked to the taxa id for row 1409 in bzstax
#all  1409 rows


###calculations 
dim(bzstax)#1409
length(which(rowSums(sign(bzstab.s)) == 1))#1209
sapply(1:6, function(z) sum(sign(bzstab.s[,z])))
#shared master-field in pairs.
dim(bzstax[ which(sign(bzstab.s[,5])==1  & sign(bzstab.s[,6])==1),])
dim(bzstax[ which(sign(bzstab.s[,3])==1  & sign(bzstab.s[,4])==1),])
dim(bzstax[ which(sign(bzstab.s[,1])==1  & sign(bzstab.s[,2])==1),]) # 1 so dim is NA
# #master in any field
# bzstax[ which( (sign( bzstab.s[,5])==1  & sign(bzstab.s[,6])==1)  | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,6])==1) | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,6])==1)),]
# bzstax[ which( (sign( bzstab.s[,5])==1  & sign(bzstab.s[,4])==1)  | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,4])==1) | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,4])==1)),]
# bzstax[ which( (sign( bzstab.s[,5])==1  & sign(bzstab.s[,2])==1)  | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,2])==1) | (sign(bzstab.s[,3])==1  & sign(bzstab.s[,2])==1)),]
# #show taxa that are in inocula
# bzstax[which(bzstab.s[,2]>0 | bzstab.s[,4]>0 | bzstab.s[,6]>0),]
# 


###subsetting

bzsItax.tab <- bzstab.s[rowSums(bzstab.s[,c(2,4,6)]) >0 ,]
bzsItax.tax <- bzstax[rowSums(bzstab.s[,c(2,4,6)]) >0 ,]
famItax <-  sapply(1:nrow(bzsItax.tax), function(z) ifelse(nchar(bzsItax.tax[z,5])>3 & !is.na(bzsItax.tax[z,5]), strsplit(bzsItax.tax[z,5],split="__")[[1]][[2]],"" ))
famItax[which(famItax=="[Chromatiaceae]")] <- "Chromatiaceae"
famItax[which(famItax=="[Exiguobacteraceae]")] <- "Exiguobacteraceae"
famItax[which(famItax=="[Weeksellaceae]")] <- "Weeksellaceae"
genItax <-  sapply(1:nrow(bzsItax.tax), function(z) ifelse(nchar(bzsItax.tax[z,6])>3 & !is.na(bzsItax.tax[z,6]), strsplit(bzsItax.tax[z,6],split="__")[[1]][[2]],"" ))
famgenItax <- paste(genItax,famItax)


###further calculations

#inoc id'd to sp? 
table(bzsItax.tax[,7]) # 
#result shows only 7 id'd to sp. rest "s__" or NA
 #           s__ s__azotifigens     s__facilis   s__johnsonii     s__veronii s__viridiflava 
#             33              1              1              1              3              1 
# 54/61
#0.8852

sum(rowSums(bzsItax.tab[,c(1,3,5)]) > 0) #ASV of inocula in any field
# [1] 23
 dim(bzsItax.tab) #total ASVs in inocual
# [1] 61  6
# > 1/14 #shared bruce pairs over total bruce inoc ASV
# [1] 0.07142857
# > 5/17 #shared KSR pairs over total KSR inoc asv
# [1] 0.2941176
# > 13/38#shared MT pairs over total MT inoc asv
# [1] 0.3421053

##inoc-inoc sharing
#BH, KSR, MT in that order sharing num with either other inocula
dim(bzstax[ which(sign(bzstab.s[,2])==1  & ( sign(bzstab.s[,4])==1 | sign(bzstab.s[,6])==1  )),])
dim(bzstax[ which(sign(bzstab.s[,4])==1  & ( sign(bzstab.s[,2])==1 | sign(bzstab.s[,6])==1  )),])
dim(bzstax[ which(sign(bzstab.s[,6])==1  & ( sign(bzstab.s[,2])==1 | sign(bzstab.s[,4])==1  )),])
#results over number ASVs
 6/14
#[1] 0.4285714
 4/17
#[1] 0.2352941
 5/38
#[1] 0.1315789


###make summary table of abundance at class order

taxsum <- matrix(NA, ncol = ncol(bzstab.s),nrow = length(unique(bzstax[,3])) )
for(i in 1:length(unique(bzstax[,3]))) {
	rows <- which(bzstax[,3] == unique(bzstax[,3])[i])
	subtax <- bzstab.s[rows,]
	if(!is.null(nrow(subtax)) ){	
	taxsum[i,] <- colSums(subtax,na.rm=T)
	} else{
	taxsum[i,] <- subtax	
	}
} 

#family sums. ONLY FOR TAXA IN INOC
famsum <- matrix(NA, ncol = ncol(bzsItax.tab),nrow = length(unique(bzsItax.tax [,5])) )
for(i in 1:length(unique(bzsItax.tax [,5]))) {
	rows <- which(bzsItax.tax[,5] == unique(bzsItax.tax [,5])[i])
	subtax <- bzsItax.tab[rows,]
	if(!is.null(nrow(subtax)) ){	
	famsum[i,] <- colSums(subtax,na.rm=T)
	} else{
	famsum[i,] <- subtax	
	}
}

#gen sums. ONLY FOR TAXA IN INOC
gensum <- matrix(NA, ncol = ncol(bzsItax.tab),nrow = length(unique(bzsItax.tax [,6])) )
for(i in 1:length(unique(bzsItax.tax [,6]))) {
	rows <- which(bzsItax.tax[,6] == unique(bzsItax.tax [,6])[i])
	subtax <- bzsItax.tab[rows,]
	if(!is.null(nrow(subtax)) ){	
	gensum[i,] <- colSums(subtax,na.rm=T)
	} else{
	gensum[i,] <- subtax	
	}
}

prettyclassnames <- sapply(1:length(unique(bzstax[,3])), function(z) ifelse(nchar(unique(bzstax[,3])[z])>3 & !is.na(unique(bzstax[,3])[z]), strsplit(unique(bzstax[,3])[z],split="__")[[1]][[2]],"" ))
#prettyfamnames <- sapply(1:length(unique(classininoc.tax[,5])), function(z) ifelse(nchar(unique(classininoc.tax[,5])[z])>3 & !is.na(unique(classininoc.tax[,5])[z]), strsplit(unique(classininoc.tax[,5])[z],split="__")[[1]][[2]],"" ))
prettyclassnames[5] <- "Saprospirae"
prettyclassnames[21] <- "Unidentified"
prettyclassnames[24] <- "Chloracidobacteria"
prettyclassnames[27] <- "Spartobacteria"
prettyclassnames[31] <- "Pedosphaerae"
prettyclassnames[35] <- "Unidentified"
prettyclassnames[39] <- "Methylacidiphilae"
prettyclassnames[41] <- "Fimbriimonadia"

#subbed tax by class table so less than X percent are added together
lowabundclass <-  !sapply(1:nrow(taxsum),function(z) any(taxsum[z,]>0.02))
hiabundsum <- taxsum[!lowabundclass,]
simplertaxsum <- rbind(hiabundsum[order(rowSums(hiabundsum),decreasing=T),],colSums(taxsum[lowabundclass,]) )
simplertaxsumnames <- c(prettyclassnames[!lowabundclass][order(rowSums(hiabundsum),decreasing=T)], "Sum of groups < 2%")

prettyfamnames <- sapply(1:length(unique(bzsItax.tax[,5])), function(z) ifelse(nchar(unique(bzsItax.tax[,5])[z])>3 & !is.na(unique(bzsItax.tax[,5])[z]), strsplit(unique(bzsItax.tax[,5])[z],split="__")[[1]][[2]],"" ))
prettyfamnames[2] <- "Chromatiaceae"
prettyfamnames[6] <- "Exiguobacteraceae"
prettyfamnames[14] <- "Weeksellaceae"

newpal <- createPalette(nrow(simplertaxsum), c(rgb(0.7,0,0),rgb(1,1,0),rgb(0,0,0.7)), M=1000)

fampca <- dudi.pca(t(famsum),scannf=F,nf=5)

#genus in any field?
#fam in any field?
gen.uglynames <-  cbind(as.data.frame(gensum),unique(bzsItax.tax[,6])) 
fam.uglynames <-  cbind(as.data.frame(famsum),unique(bzsItax.tax[,5])) 
taxgen.in.field <- as.numeric(sapply(1:nrow(bzsItax.tab) , function(z)  bzsItax.tax[z,6]%in% gen.uglynames[rowSums(gen.uglynames[,c(1,3,5)])>0,7]) )
taxgen.in.field[bzsItax.tax[,6]=="g__" | is.na(bzsItax.tax[,6])] <- 0.5
taxfam.in.field <- as.numeric(sapply(1:nrow(bzsItax.tab) , function(z)  bzsItax.tax[z,5]%in% fam.uglynames[rowSums(fam.uglynames[,c(1,3,5)])>0,7]) )
taxfam.in.field[bzsItax.tax[,5]=="f__" | is.na(bzsItax.tax[,5])] <- 0.5
#CAREFUL IF USING ANY ALTERED/UPDATED FILES. the above line *happens* to eliminate ones without genus names, so no extra code written. 

bw <- colorRampPalette( c( rgb(1,1,1),rgb(0,0,0) ) )

pdf("~/Figure1.pdf",height=3,width=7)
layout(matrix(c(1:3), ncol=3), widths=c(3.75,1.75,1.75))
par(oma=c(0,0,0,2))
par(mar=c(4,5,3,0))
	bg <- barplot(simplertaxsum[,], col= newpal,ylab="Proportion",xlab="",names.arg=c("Field","Inoculum","Field","Inoculum","Field","Inoculum") )
		axis(  at=c((bg[1]+bg[2])/2, (bg[3]+bg[4])/2 ,(bg[5]+bg[6])/2 ), side=1, labels = c("BH", "KSR","MT") ,lty=0,line=1.2)
		mtext("a.",side=3,adj=-0.2,line=0.5)
par(mar=c(4,0,0,2))
	plot(rep(1,times=3)~c(1:3), pch=NA,bty="n",xaxt="n",yaxt="n",ylab="",xlab="")
 		legend(x=1,y=1.3,legend=simplertaxsumnames[],fill=newpal,bty="n",ncol=1,cex=1,xpd=NA)
par(mar=c(6,4,3,0))
plot(fampca$li[,1]~c(1.5,1.5,3.5,3.5,5.5,5.5),ylab="PC Axis 1",xlab="", cex=2, ylim=c(-5,8),xlim=c(0,7),xaxt="n",
	col="black",pch=c(5,2,5,2,5,2)) #field is 16, bruce is red, ksr is black,mt is blue
	mtext("b.",side=3,adj=-0.2,line=0.5)
	legend(x=0,y=7,legend=c("Field","Inocula"),col="black",pch=c(18,17),bty="n", pt.cex=2)
	abline(v=c(2.5,4.5),lty=3)
	abline(h=0,lwd=1.2)
	axis(side=1,at=c(1.5,3.5,5.5),labels=c("BH", "KSR","MT"))
dev.off()

sitenames <- c("BH Field","BH Inoculum","KSR Field","KSR Inoculum","MT Field","MT Inoculum")

halfcolscale <- round(exp(log(3.568115e-05)/2),digits=3)
pdf("~/abundpanel_for_supp_text.pdf",height=4,width=10)
layout(matrix(c(1:8), ncol=2,byrow=F), heights=c(0.1,0.1,0.1,2) , widths=c(5,1) )
par(oma=c(2,0,2,2))
par(mar=c(0,8,0,2))
	image(as.matrix(  				    taxfam.in.field[order(rowSums(bzsItax.tab),decreasing=T)] ) , col= bw(3),yaxt="n",xaxt="n")
	 axis(side=2,at=c(0.5),labels="Family in Field",las=1)
	image(as.matrix(  				    taxgen.in.field[order(rowSums(bzsItax.tab),decreasing=T)] ) , col= bw(3),yaxt="n",xaxt="n")
	 axis(side=2,at=c(0.5),labels="Genus in Field",las=1)
	image(sign(as.matrix(rowSums(bzsItax.tab[,c(1,3,5)][order(rowSums(bzsItax.tab),decreasing=T),] )) ), col=bw(2),yaxt="n",xaxt="n")
	 axis(side=2,at=c(0.5),labels="ASV in Field",las=1)
par(mar=c(18,8,0,2))
	image(as.matrix(log(bzsItax.tab[order(rowSums(bzsItax.tab),decreasing=T),])),col=rev(heat.colors(20)),zlim=c(log(3.568115e-05),log(1)),yaxt="n",xaxt="n")
	 mtext("ASV taxonomy",side=1,line=16.25,cex=0.75)
	 axis(side=2,at=seq(from=0,to=1,length.out=6),labels=(sitenames),las=1)
	 axis(side=1,at=seq(from=0,to=1,length.out=nrow(bzsItax.tab)),labels=(famgenItax)[order(rowSums(bzsItax.tab),decreasing=T)],las=2)
par(mar=c(0,0,0,0))
	plot(1,1,bty="n",pch=NA,xaxt="n",yaxt="n",ylab="",xlab="",xlim=c(0,2),ylim=c(0,2))
	plot(1,1,bty="n",pch=NA,xaxt="n",yaxt="n",ylab="",xlab="")
	plot(1,1,bty="n",pch=NA,xaxt="n",yaxt="n",ylab="",xlab="")
par(mar=c(0,0,0,0))
	plot(0:6,0:6,bty="n",pch=NA,xaxt="n",yaxt="n",ylab="",xlab="")
legend.gradient(cbind(c(0,1,1,0),c(1,1,5,5)),cols=c(rgb(1,1,1),rev(heat.colors(20))),title="",limits=c("0","1"),cex=1.5,bty="o")
rect(0, 1, 1, 5)
     text(1.5,3, "ln(Proportion) of Inocula",adj=c(0.5,0),srt=270,cex=1.5)
legend(x=-1,y=6.75,fill=c(rgb(0,0,0),rgb(.5,.5,.5),rgb(1,1,1)),legend=c("Taxonomy present","Unknown","Taxonomy absent"),cex=1.25,bty="n",xpd=NA)
dev.off()
