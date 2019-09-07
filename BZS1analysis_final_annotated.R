
##LOAD REQUIRED LIBRARY

library(MCMCglmm)



####DEFINING UNCTIONS



bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}
	
range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

#
circ2perim <- function(A,C){ return(sqrt((A*4*pi)/C))}
std.error <- function(x){sd(x)/sqrt(length(x))}


#function for analyzing data in map and dat files
MapToWellsB <- function(dat,map,platecol,firstcol, sumcols, meancols){ #columns labeled "time"[not anymore] "roi", "image"
	n <- nrow(map)
	outdata <- matrix(NA,nrow=n, ncol=(length(sumcols)+length(meancols)+1))
	for(i in 1:n){
		rois <- map[i,firstcol:ncol(map)] #what happens to NA values
		p <- map$plate[i]
		#day <- map$time[i]
		welldat <- dat[dat$plate==p & dat$roi%in%rois,]#& dat$time==day 
		welldat.sums <- colSums(welldat[,sumcols])
		welldat.means <- colMeans(welldat[,meancols],na.rm=T)
		wellstats <- c(nrow(welldat),welldat.sums,welldat.means)
		outdata[i,] <- wellstats
		}
	mappeddata <- cbind(map[,1:(firstcol-1)],outdata)
	colnames(mappeddata) <- c(colnames(map)[1:(firstcol-1)],"particles",colnames(dat)[sumcols],colnames(dat)[meancols])
	return(mappeddata)
}


#This function converts OD to cells/ul to a desired number of cells in inocula, according to the formula at http://www.reric.org/wordpress/archives/169
# Build a function to calculate cells/ul for each sample
#we used this function to standardize inocula concentrations.
ODtoCells <- function(x){
	log10cellspml <- 9.9 - 4/(1 + (x/.14)^0.8 ) #log base 10 cellspermL
	cellpul <-  (10^log10cellspml)/1000 #per ul
	return(cellpul)
} #note that the minimum estimate is 794 cells / ul... (relevant for expressing OD as cell count, if of interest, below. not in O'Brien et al.)


#ratio calculations for pixel area from bottom or top images
#well diam (start) end, plate1 b3: 445
#well diam nov23 (mid): 565
#ratio to multiply bottom plates:
botrat <- 445/565




####READING DATA AND MAKING DATAFRAMES FOR ANALYSIS
##THIS ASSUMES FILES ARE IN "~/" directory and saves output to same


#
#read in start and end image analyses and roi to well map files
bzdatnov16start <- read.csv("~/StartDatBZS1Nov16.csv",header=T,stringsAsFactors=F)#
bzmapnov16start <- read.csv("~/StartMapBZS1Nov16.csv",header=T,stringsAsFactors=F)
bzdatnov20 <- read.csv("~/Nov20datBZS1.csv",header=T,stringsAsFactors=F)
bzmapnov20 <- read.csv("~/Nov20mapBZS1.csv",header=T,stringsAsFactors=F)
bzdatnov23 <- read.csv("~/Nov23datBZS1.csv",header=T,stringsAsFactors=F)
bzmapnov23 <- read.csv("~/Nov23mapBZS1.csv",header=T,stringsAsFactors=F)
bzdatnov26end <- read.csv("~/EndDatBZS1Nov26.csv",header=T,stringsAsFactors=F)
bzmapnov26end <- read.csv("~/EndMapBZS1Nov26.csv",header=T,stringsAsFactors=F)

#merge map and dat data, for columns of interest, not all files the same.
mapdatn16 <- MapToWellsB(bzdatnov16start,bzmapnov16start,firstcol=4,sumcols=c(4,6),meancols=c(5,7,9:14))#
mapdatn20 <- MapToWellsB(bzdatnov20,bzmapnov20,firstcol=4,sumcols=c(4,6),meancols=c(5,7,9:14))
mapdatn23 <- MapToWellsB(bzdatnov23,bzmapnov23,firstcol=4,sumcols=c(4,6),meancols=c(5,7,9:14))#
mapdatn26 <- MapToWellsB(bzdatnov26end,bzmapnov26end,firstcol=4,sumcols=c(4,6),meancols=c(5,7,9:14))#different measures

trts <- read.csv("~/SaltBTduckweeddesTRTS.csv") #data is already sorted the same design

ODdat <- read.csv("~/AO BZS.1 12mar2019_sorted.csv",header=T)
logOD600 <- log(round(ODdat$od600mean.blanked,digits=4) + 0.01) #CORRECTING FOR DETECTION LIMIT: a bit larger than the most negative value, and approximately one third of the minimum positive value: min(round(ODdat$od600mean.blanked,digits=4)[min(which(round(ODdat$od600mean.blanked,digits=4)>0))])
logOD750 <- log(round(ODdat$od750mean.blanked,digits=4) + 0.01)
flOD600 <- log(ODdat$od600.fancyblankmean + 0.01)
flOD750 <- log(ODdat$od750.fancyblankmean + 0.01)
#cell counts. use fancy blanked, in case of interest in cell count vs OD (we consider OD more accurate and report in paper)
mincells <- min(ODtoCells(ODdat$od600.fancyblankmean))
mindetect10p <- min(ODtoCells(ODdat$od600.fancyblankmean)[ODtoCells(ODdat$od600.fancyblankmean)>mincells])*0.1
Lcellfb6 <- log(ODtoCells(ODdat$od600.fancyblankmean) - mincells + mindetect10p) # minus the minimum detectable (approximately 793) + 10% of the minimum detected over 0

#NOTES SAY, Nov26End plate 13, C5 AND plate 24 D5, pass filters, but might (prb are) be dead – covered with black fungus., remove them here
rmrows <- c(which(trts$plate==13&trts$col==5 & trts$row=="C"),which(trts$plate==24&trts$col==5 & trts$row=="D"))
mapdatn26[rmrows,4:14] <- NA

bzs.dattrt <- data.frame(pixstart = mapdatn16$area,  pix20 =  mapdatn20$area*botrat, pix23 = mapdatn23$area*botrat, pixend = mapdatn26$area,
	alivestart = as.numeric(mapdatn16$area>0), alive20 =as.numeric(mapdatn20$area>0), alive23 = as.numeric(mapdatn23$area>0), aliveend = as.numeric(mapdatn26$area>0),
	deadstart = 1-as.numeric(mapdatn16$area>0), dead20 = 1- as.numeric(mapdatn20$area>0), dead23= 1-as.numeric(mapdatn23$area>0), deadend = 1- as.numeric(mapdatn26$area>0),
	 dpixse =  mapdatn26$area -  mapdatn16$area ,
	dpixs20 = mapdatn20$area*botrat - mapdatn16$area, dpix2023 =  (mapdatn23$area - mapdatn20$area)*botrat, dpixs23 =  mapdatn23$area*botrat - mapdatn16$area, dpix23e =  mapdatn26$area -  mapdatn23$area*botrat  , dpixse =  mapdatn26$area -  mapdatn16$area ,
	bzt = trts$bzt, nacl=trts$nacl, inoc=trts$inoc, plnt.micr=trts$plnt.micr, plate = trts$plate,
	lbzt = log(trts$bzt+.01), lnacl=log(trts$nacl+.08), setuperror= trts$SetupErrorPossible,
	lOD600 = logOD600, lOD750 = logOD750, flOD600 = flOD600, Lcellfb6 = Lcellfb6, cellfb6 = ODtoCells(ODdat$od600.fancyblankmean) - 792,
	red =mapdatn26$red / (3*mapdatn26$meangray), green = 1- mapdatn26$red / (3*mapdatn26$meangray) - mapdatn26$blue / (3*mapdatn26$meangray), blue =mapdatn26$blue / (3*mapdatn26$meangray), #there's an error in the RGB data. plate 3 has identical values for green and blue., given their range (below all other green values), they are clearly the blue values duplicated, #HOWEVER, these values add to 1, nearly perfectly, so just get green from the other two. check out  plot(bzs.dattrt$green ~ I(mapdatn26$green / (3*mapdatn26$meangray)) )
	aggrE = mapdatn26$area/mapdatn26$perim, daggrE = mapdatn26$area/mapdatn26$perim - mapdatn16$area/mapdatn16$perim) 
bzs.dattrt.c <- bzs.dattrt[bzs.dattrt$setuperror != "y",]


bzs.dattrt$x <- ( (ceiling(trts$plate/8) -1) * 127.63 ) + ( trts$column * (127/6) )
bzs.dattrt$y <- -1* ( ( (ifelse( (trts$plate %% 8)==0,8,(trts$plate %% 8)) -1) * 85.4 ) + ( as.numeric(trts$row) * (85/4) ) )# across plate + within plate y # %% gives the remainder
#in case interest in spatial effects

sizetime <- data.frame(pix = c(bzs.dattrt[,1],bzs.dattrt[,2],bzs.dattrt[,3],bzs.dattrt[,4]),
			days = rep(c(0,4,7,10), each = nrow(bzs.dattrt)),
			bzt = rep(bzs.dattrt$bzt,times=4), nacl = rep(bzs.dattrt$nacl,times=4), plate = rep(bzs.dattrt$plate, times=4),
			inoc = rep(bzs.dattrt$inoc,times=4), plnt.micr = rep(bzs.dattrt$plnt.micr,times=4), indiv = rep(1:576,times=4), 
			lnacl = rep(log(bzs.dattrt$nacl + .08),times=4), lbzt = rep(log(bzs.dattrt$bzt+.01),times=4))
sizetime$rpix <- round(sizetime$pix)
#
bzs.dattrtI <- bzs.dattrt[bzs.dattrt$inoc=="y",]

#BZT and NACL were chosen at power of 10 increasing levels. 
#note that OD and cell count data seem to have similar inflated 0 issues causing a similar non-normal trend. but "0" is not strictly 0, instead the detection limit. Log of cell count, where we subtract the minimum calculateable value, and then re-add 10% of it. -- see above




##LINEAR MODEL HYPOTHESIS TESTING, AND FIGURES


####RGB COLOR
#cor(cbind(bzs.dattrt$red,bzs.dattrt$green,bzs.dattrt$blue),use="complete.obs") # green and blue are essentially identical and opposite.
green.1anla.np <- MCMCglmm(green~nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc + inoc:bzt:nacl, random = ~  plnt.micr ,	 family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
green.1anla.np.tmp <- MCMCglmm(green~nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc , random = ~  plnt.micr ,	 family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
green.1anla.np.tmp2 <- MCMCglmm(green~nacl + bzt + inoc  , random = ~  plnt.micr ,	 family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
green.1anla.refit.np <- MCMCglmm(green~nacl , random = ~  plnt.micr ,	 family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
####AGGEGATION
aggrE.Lr.mcmcglmm4anla.np <-  MCMCglmm(aggrE~  nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc + inoc:bzt:nacl, random = ~  plnt.micr ,family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
aggrE.Lr.mcmcglmm4anla.np.tmp <-  MCMCglmm(aggrE~  nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc , random = ~  plnt.micr ,family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
aggrE.Lr.mcmcglmm4anla.np.tmp2 <-  MCMCglmm(aggrE~  nacl + bzt + inoc , random = ~  plnt.micr ,family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
aggrE.Lr.mcmcglmm4anla.refit.np <-  MCMCglmm(aggrE~  nacl, random = ~  plnt.micr ,family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
####OD
 OD600I.mcmcglmm4a.np <- MCMCglmm(flOD600~ nacl + bzt + bzt:nacl, random = ~  plnt.micr ,family="gaussian",data=bzs.dattrtI,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
 OD600I.mcmcglmm4a.np.tmp <- MCMCglmm(flOD600~ nacl + bzt , random = ~  plnt.micr ,family="gaussian",data=bzs.dattrtI,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
 OD600I.mcmcglmm4a.refit.np <- MCMCglmm(flOD600~ nacl , random = ~  plnt.micr ,family="gaussian",data=bzs.dattrtI,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
 OD600I.inoc.test.np <- MCMCglmm(flOD600~ inoc , random = ~  plnt.micr ,family="gaussian",data=bzs.dattrt,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
####INOC TEST
inocpost <- OD600I.inoc.test.np$Sol
plantinoc <- rowMeans(inocpost[,3:5]) 
#plateinoc <- rowMeans(inocpost[,6:29]) 
inoc.mn <- sapply(c(0,1), function(n) mean(inocpost[,1] + inocpost[,2]*n + plantinoc ))
inoc.CI <- sapply(c(0,1), function(n) HPDinterval(inocpost[,1] + inocpost[,2]*n + plantinoc , ,prob=0.95))

##TRAITS FIGURE
#posteriors and predictions
nacl.s <- sort(unique(bzs.dattrt$nacl))
#
greenpost <- green.1anla.refit.np$Sol
plantgr <- rowMeans(greenpost[,3:5]) 
green.mn <- sapply(nacl.s, function(n) mean(greenpost[,1] + greenpost[,2]*n + plantgr ))
green.CI <- sapply(nacl.s, function(n) HPDinterval(greenpost[,1] + greenpost[,2]*n + plantgr , ,prob=0.95))
#
agrpost <- aggrE.Lr.mcmcglmm4anla.refit.np$Sol
plantag <- rowMeans(agrpost[,3:5]) 
aggr.mn <- sapply(nacl.s, function(n) mean(agrpost[,1] + agrpost[,2]*n + plantag ))
aggr.CI <- sapply(nacl.s, function(n) HPDinterval(agrpost[,1] + agrpost[,2]*n + plantag , ,prob=0.95))
#
odpost <- OD600I.mcmcglmm4a.refit.np$Sol
plantod <- rowMeans(odpost[,3:5]) 
od.mn <- sapply(nacl.s, function(n) mean(odpost[,1] + odpost[,2]*n + plantod ))
od.CI <- sapply(nacl.s, function(n) HPDinterval(odpost[,1] + odpost[,2]*n + plantod , ,prob=0.95))
#data averages
grnmn <- tapply(bzs.dattrt$green, (bzs.dattrt$nacl), mean, na.rm=T)
grnse <- tapply(bzs.dattrt$green, (bzs.dattrt$nacl), sd, na.rm=T)/sqrt(table((bzs.dattrt$nacl)[!is.na(bzs.dattrt$green)]))
agmn <- tapply(bzs.dattrt$aggrE, (bzs.dattrt$nacl), mean, na.rm=T)
agse <- tapply(bzs.dattrt$aggrE, (bzs.dattrt$nacl), sd, na.rm=T)/sqrt(table((bzs.dattrt$nacl)[!is.na(bzs.dattrt$aggrE)]))
odmn <- tapply(bzs.dattrtI$flOD600, (bzs.dattrtI$nacl), mean, na.rm=T)
odse <- tapply(bzs.dattrtI$flOD600, (bzs.dattrtI$nacl), sd, na.rm=T)/sqrt(table((bzs.dattrtI$nacl)[!is.na(bzs.dattrtI$flOD600)]))
#plot
pdf("~/BZS1_Trait_results_plot.pdf",width=7,height=3)
par(mfrow=c(1,3))
par(mar=c(4,4,1,0))
par(oma=c(1,0,2,1))
plot(grnmn~c(1:3),pch=1,cex=3, ylim =c(0.4,0.48),
		xaxt="n",ylab = "",xlab="",xlim=c(0,4),cex.axis=1.5)
	mtext("Proportion green pixels",side=2,line=2.5)
	mtext("a.",side=3,adj=-0.3,line=1.5)
	arrows(1:3, y0 = grnmn -grnse, y1 =grnmn +grnse,length=0,lwd=3)
	axis(side = 1,at=1:3, labels=nacl.s,,cex.axis=1.5)
plot(agmn~c(1:3),pch=1,cex=3, ylim=c(4,10),
		xaxt="n",ylab = "",xlim=c(0,4),xlab="",cex.axis=1.5) 
	mtext("Aggregation (area:perimeter)",side=2,line=2.5)
	mtext("NaCl g/L",side=1,line=2.5)
	mtext("b.",side=3,adj=-0.3,line=1.5)
	arrows(1:3, y0 = agmn - agse, y1 =agmn + agse,length=0,lwd=3)
	axis(side = 1,at=1:3, labels=nacl.s,cex.axis=1.5)
plot(odmn~c(1:3),pch=1,cex=3, ylim =c(-4.25,-3.75),
		xaxt="n",xlab = "",ylab="",xlim=c(0,4),cex.axis=1.5) 
	mtext("ln Optical Density",side=2,line=2.5)
	mtext("c.",side=3,adj=-0.3,line=1.5)
	arrows(1:3, y0 = odmn - odse , y1 = odmn+odse,length=0,lwd=3)
	axis(side = 1,at=1:3, labels=nacl.s,cex.axis=1.5)
dev.off()


###########SURVIVAL AND GROWTH
#######


### TESTS FOR NORMALITY, FIGURES.
aliveatend <- which(bzs.dattrt$pixend>0)
sizetime.alive <- sizetime[ c(aliveatend, 576 + aliveatend, 576*2 + aliveatend, 576*3 + aliveatend),]
sizetime.aliveHS <- sizetime.alive[sizetime.alive$nacl>1,]
sizetime.aliveMS <- sizetime.alive[sizetime.alive$nacl==0.8,]#
sizetime.aliveNS <- sizetime.alive[sizetime.alive$nacl==0,]#
stestHS <- shapiro.test(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl >1])
stest0S <- shapiro.test(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0])
stestMS <- shapiro.test(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0.8])
pdf("~/qqPix_liveplantsAS.pdf",width=7,height=5)
par(mfrow=c(2,3))
par(mar=c(4,4,1,0))
par(oma=c(0,0,2,1))
hist(sizetime$pix[sizetime$nacl==0 & sizetime$days==10],ylab="Frequency, All plants",xlab="",main="",breaks=20)
mtext("NaCl 0 g/L",line=1)
hist(sizetime$pix[sizetime$nacl==0.8 & sizetime$days==10],ylab="",xlab="",main="",breaks=20)
mtext("NaCl 0.8 g/L",line=1)
mtext("Pixel area",line=2,side=1)
hist(sizetime$pix[sizetime$nacl==10 & sizetime$days==10],ylab="",xlab="",main="",breaks=20)
mtext("NaCl 10 g/L",line=1)
qqnorm(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0],main="",ylab="Sample Quantiles, Live plants",xlab="")
qqline(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0])
mtext(paste("Shapiro-Wilk W=",round(stest0S$statistic,digits=2),"*",sep=""),side=3)
qqnorm(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0.8],main="",ylab="",xlab="")
qqline(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl ==0.8])
mtext(paste("Shapiro-Wilk W=",round(stestMS$statistic,digits=2),"*",sep=""),side=3)
mtext("Theoretical Quantiles",line=2,side=1)
qqnorm(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl >1],main="",ylab="",xlab="")
qqline(sizetime.alive$pix[sizetime.alive$days==10 & sizetime.alive$nacl >1])
mtext(paste("Shapiro-Wilk W=",round(stestHS$statistic,digits=2),"*",sep=""),side=3)
dev.off()
##CONCLUSION IS TO SPLIT BY SURVIVED OR DIED, AND GROWTH AT DIFFERENT SALT LEVELS

#
######MODEL SURVIVAL
aliveend.np <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc + inoc:bzt:nacl,
								random = ~  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)#
aliveend.tmp <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc,
								random = ~  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)# 
aliveend.npF <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc + inoc:bzt:nacl +  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)#
aliveend.tmpF <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + inoc + bzt:nacl + bzt:inoc + nacl:inoc +  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)#
aliveend.np.rm.n <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + nacl:bzt ,
								random = ~  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)#
aliveend.np.rm.nF <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + nacl:bzt +  plnt.micr  ,
								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
								nitt=1000000, thin=100, burnin=2000)#
#RE rand vs fixed genos: do not change final selected model  but then nacl:bzt is insig when inoc interact terms rm'd. further, genoeffs but are sig for fixed and ns rand
# aliveend.np.rm.nx2samps <- MCMCglmm(cbind(aliveend,deadend)~ nacl + bzt + nacl:bzt ,
# 								random = ~  plnt.micr  ,
# 								family="multinomial2", data=bzs.dattrt,verbose=FALSE,pr=TRUE,
# 								nitt=2000000, thin=100, burnin=2000)#
logistic <- function(x){exp(x)/(1+ exp(x))}
DuckLivePost <- aliveend.np.rm.n$Sol
days.s <- seq(from=0,to=10,length.out=1000)
nacl.s <- sort(unique(sizetime$nacl))
bzt.s <- sort(unique(sizetime$bzt))
nacl.ss <- seq(from=0,to=10,length.out=1000)
bzt.ss <- seq(from=0,to=10,length.out=1000)
 lnacl.s <- sort(unique(sizetime$lnacl))
plntmnefB <- rowMeans(DuckLivePost[,5:7])
# 
alive.pred.mn <- lapply(bzt.s, function(b) sapply(nacl.ss, function(n)
				mean(logistic(DuckLivePost[,1] + DuckLivePost[,2]*n + DuckLivePost[,3]*b + DuckLivePost[,4]*n*b + plntmnefB ))
				))
alive.pred.CI <- lapply(bzt.s, function(b) sapply(nacl.ss, function(n)
				HPDinterval(as.mcmc(logistic(DuckLivePost[,1] + DuckLivePost[,2]*n + DuckLivePost[,3]*b + DuckLivePost[,4]*n*b + plntmnefB )),prob=0.9))
				)
alive.pred.mnB <- lapply(nacl.s, function(n) sapply(bzt.ss, function(b)
				mean(logistic(DuckLivePost[,1] + DuckLivePost[,2]*n + DuckLivePost[,3]*b + DuckLivePost[,4]*n*b + plntmnefB ))
				))
alive.pred.CIB <- lapply(nacl.s, function(n) sapply(bzt.ss, function(b)
				HPDinterval(as.mcmc(logistic(DuckLivePost[,1] + DuckLivePost[,2]*n + DuckLivePost[,3]*b + DuckLivePost[,4]*n*b + plntmnefB )),prob=0.9))
				)

nsurva.mn <- tapply(bzs.dattrt$aliveend,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), mean,na.rm=T)
nsurva.se <- tapply(bzs.dattrt$aliveend,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), sd,na.rm=T)/ sqrt(table(paste(bzs.dattrt$nacl, bzs.dattrt$bzt)[!is.na(bzs.dattrt$aliveend)]))
nsurva.bzt <- tapply(bzs.dattrt$bzt,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), mean,na.rm=T)
nsurva.lbzt <- tapply(bzs.dattrt$lbzt,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), mean,na.rm=T)
nsurva.lnacl <- tapply(bzs.dattrt$lnacl,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), mean,na.rm=T)
nsurva.nacl <- tapply(bzs.dattrt$nacl,paste(bzs.dattrt$nacl, bzs.dattrt$bzt), mean,na.rm=T)
pdf("~/saltsurvdatonly.pdf",height=2.25,width=3)
par(mar=c(3,3,1,5))
	plot(nsurva.mn ~ I(rep(1:4,times=3) + 0.05*as.numeric(as.factor(nsurva.nacl)))  , pch=21,xaxt="n", yaxt="n",
	  col = rgb(range01(nsurva.lnacl),0,0),bg = rgb(range01(nsurva.lnacl),0,0,alpha=0.5), cex=1.25,
	 ylim=c(0,1),xlim=c(0.75,4.25),ylab="",xlab="")
	axis(side=1,at=1:4, labels=c(0,0.1,1,10))
	axis(side=2,at=c(0,0.25,0.5,0.75,1), labels=c(0,"",0.5,"",1))
	mtext("Proportion Survived",side=2,line=2)
	mtext("Benzotriazole mg/L",side=1,line=2)
	arrows(x0 = I(rep(1:4,times=3) + 0.05*as.numeric(as.factor(nsurva.nacl))),
	 y0 =nsurva.mn -nsurva.se , y1=nsurva.mn + nsurva.se,
	 col = rgb(range01(nsurva.lnacl),0,0), length=0)
	text(5.5,0.8,"NaCl g/L",xpd=NA)
	legend(4.5,0.8, c("0","0.8","10"),pch=c(21,21,21), col = rgb(range01(unique(nsurva.lnacl)),0,0),pt.bg = rgb(range01(unique(nsurva.lnacl)),0,0,alpha=0.5),bty="n",xpd=NA)
dev.off()

###MODEL GROWTH.  
###SPLIT INTO 3 because of normality/multidistribution issues.
growth.time.nlHSP.np <- MCMCglmm(rpix ~ days + I(days^2) + days:bzt + days:inoc + days:bzt:inoc ,
								random = ~ indiv + plnt.micr  ,
								family="poisson", data=sizetime.aliveHS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlHSP2.np <- MCMCglmm(rpix ~ days + I(days^2) + I(days^2):bzt + I(days^2):inoc + I(days^2):bzt:inoc,
								random = ~ indiv + plnt.micr  ,
								family="poisson", data=sizetime.aliveHS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlHSP2.nptmp <- MCMCglmm(rpix ~ days + I(days^2) + I(days^2):bzt + I(days^2):inoc ,
								random = ~ indiv + plnt.micr  ,
								family="poisson", data=sizetime.aliveHS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlHSP.rm.np <- MCMCglmm(rpix ~ days + I(days^2) + I(days^2):inoc , ### days ^2 interacts fit better
								random = ~ indiv + plnt.micr  ,
								family="poisson", data=sizetime.aliveHS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlMSG1.np <- MCMCglmm(pix ~  days + I(days^2) + days:bzt + days:inoc + days:bzt:inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveMS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlMSG2.np <- MCMCglmm(pix ~  days + I(days^2) + I(days^2):bzt + I(days^2):inoc + I(days^2):bzt:inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveMS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#days squared ints better.
growth.time.nlMSG2.np.tmp <- MCMCglmm(pix ~  days + I(days^2) + I(days^2):bzt + I(days^2):inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveMS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlMSG2.np.rm <- MCMCglmm(pix ~  days + I(days^2) + I(days^2):inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveMS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#BEST.
growth.time.nlNSG1.np <-MCMCglmm(pix ~  days + I(days^2) + days:bzt + days:inoc + days:bzt:inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveNS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#
growth.time.nlNSG2.np <-  MCMCglmm(pix ~  days + I(days^2) + I(days^2):bzt + I(days^2):inoc + I(days^2):bzt:inoc,
								random = ~ indiv + plnt.micr ,
								family="gaussian", data=sizetime.aliveNS,verbose=FALSE,pr=TRUE,
								nitt=100000, thin=100, burnin=2000)#days squared ints better., this one needs no removal 

DuckGrowPostHS <- growth.time.nlHSP.rm.np$Sol
DuckGrowPostMS <- growth.time.nlMSG2.np.rm$Sol
DuckGrowPostNS <- growth.time.nlNSG2.np$Sol

indivmnefHS <- rowMeans(DuckGrowPostHS[,5:107]) 
plntmnefHS <- rowMeans(DuckGrowPostHS[,108:110]) #HPDinterval(DuckGrowPostHS[,108:110]) and similar for MS, NS, show no sig diff
indivmnefMS <- rowMeans(DuckGrowPostMS[,5:172]) 
plntmnefMS <- rowMeans(DuckGrowPostMS[,173:175])
indivmnefNS <- rowMeans(DuckGrowPostNS[,7:185]) 
plntmnefNS <- rowMeans(DuckGrowPostNS[,186:188])

##Growth FIGURE
NS10 <- sizetime.aliveNS[sizetime.aliveNS$days==10,]
MS10 <- sizetime.aliveMS[sizetime.aliveMS$days==10,]
HS10 <- sizetime.aliveHS[sizetime.aliveHS$days==10,]
hs10mn <- (tapply(log(HS10$rpix),HS10$inoc,mean,na.rm=T)) # fit on log scale, so take averages on log scale also.
hs10se <- (tapply(log(HS10$pix),HS10$inoc,sd,na.rm=T) / sqrt(table(HS10$inoc)) ) #filtered to survived.
ms10mn <- tapply(MS10$pix,MS10$inoc,mean,na.rm=T)
ms10se <- tapply(MS10$pix,MS10$inoc,sd,na.rm=T)/ sqrt(table(MS10$inoc))
ns10mn <- tapply(NS10$pix,paste(NS10$inoc,NS10$bzt),mean,na.rm=T)
ns10se <- tapply(NS10$pix,paste(NS10$inoc,NS10$bzt),sd,na.rm=T) /  sqrt(table(paste(NS10$inoc,NS10$bzt)))
ns10bzt <- tapply(NS10$bzt,paste(NS10$inoc,NS10$bzt),mean,na.rm=T)

pdf("~/summarygrowthres.pdf",height=2.5,width=4.5)
layout(matrix(1:3,ncol=3),widths=c(1.25,1,1))
par(mar=c(4,0,2,0))
par(oma=c(0,4,2,2))
plot(ns10mn ~ rep(c(1:4),times=2), cex=2,pch=21, xaxt = "n",
			col=rep(c(rgb(0.5,0.5,0.5),rgb(0,0,0)),each=4),
			bg=rep(c(rgb(0.5,0.5,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),each=4),
			ylim=c(0,6000),xlim=c(0.75,4.25),xlab="",ylab="")
	axis(at=c(1:4),labels=c(0,0.1,1,10),side=1)
	mtext("Pixel area", side=2,line=2.5)
	mtext("NaCl 0 g/L", side=3,line=0.75)
	mtext("a.", side=3,line=2.5,adj = 0)
 	mtext("BZT mg/L", side=1,line=2.5)
	arrows(x0= rep(c(1:4),times=2), y0=ns10mn-ns10se, y1= ns10mn+ns10se, length=0,col=rep(c(rgb(0.5,0.5,0.5,),rgb(0,0,0)),each=4))
plot(ms10mn~c(1,2), ylim=c(0,6000),xaxt="n" ,xlim=c(0.5,2.5),xlab="",ylab="",cex=2,pch=21,
		col=c(rgb(0.5,0.5,0.5,),rgb(0,0,0)),bg=c(rgb(0.5,0.5,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)) ,yaxt="n")
	mtext("b.", side=3,line=2.5,adj = 0)
	mtext("NaCl 0.8 g/L", side=3,line=0.75)
	arrows(x0= rep(1:2, times=2), y0=ms10mn-ms10se, y1= ms10mn+ms10se,length=0,col=c(rgb(0.5,0.5,0.5,),rgb(0,0,0)) )
plot(exp(hs10mn)~c(1,2),pch = 21,ylim=c(0,6000),xaxt="n" ,xlim=c(0.5,2.5),xlab="",ylab="",cex=2,
		col=c(rgb(0.5,0.5,0.5,),rgb(0,0,0)),bg=c(rgb(0.5,0.5,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),yaxt="n")
	mtext("c.", side=3,line=2.5,adj = 0)
	mtext("NaCl 10 g/L", side=3,line=0.75)
	arrows(x0= 1:2, y0=exp(hs10mn-hs10se), y1= exp(hs10mn+hs10se),length=0 ,col=c(rgb(0.5,0.5,0.5,),rgb(0,0,0)))
	legend(0.35,6000, c("Not inoculated","Inoculated"),
		border=c(rgb(0.5,0.5,0.5,),rgb(0,0,0)),fill=c(rgb(0.5,0.5,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),bty="n" )
dev.off()	

	
######SUPFIG, FULL TIMESERIES.
# (Intercept)      days  I(days^2) I(days^2):inocy
grow.pmnHS.inoc <- sapply(days.s, function(d) 
				exp(mean(DuckGrowPostHS[,1] + DuckGrowPostHS[,2]*d + DuckGrowPostHS[,3]*(d^2) + DuckGrowPostHS[,4]*d^2*1 + indivmnefHS + plntmnefHS ))
				)
grow.pmnHS.ninoc <- sapply(days.s, function(d) 
				exp(mean(DuckGrowPostHS[,1] + DuckGrowPostHS[,2]*d + DuckGrowPostHS[,3]*(d^2) + DuckGrowPostHS[,4]*d^2*0 + indivmnefHS + plntmnefHS ))
				)
grow.pCIHS.inoc <- sapply(days.s, function(d) 
				exp(HPDinterval(as.mcmc(DuckGrowPostHS[,1] + DuckGrowPostHS[,2]*d + DuckGrowPostHS[,3]*(d^2) + DuckGrowPostHS[,4]*d^2*1 + indivmnefHS + plntmnefHS ),prob=0.95))
				)
grow.pCIHS.ninoc <- sapply(days.s, function(d) 
				exp(HPDinterval(as.mcmc(DuckGrowPostHS[,1] + DuckGrowPostHS[,2]*d + DuckGrowPostHS[,3]*(d^2) + DuckGrowPostHS[,4]*d^2*0 + indivmnefHS + plntmnefHS ),prob=0.95))
				)
#"(Intercept)"     "days"            "I(days^2)"       "I(days^2):inocy"
grow.pmnMS.inoc <-  sapply(days.s, function(d) 
				mean(DuckGrowPostMS[,1] + DuckGrowPostMS[,2]*d + DuckGrowPostMS[,3]*d^2 + DuckGrowPostMS[,4]*1*d^2 + indivmnefMS + plntmnefMS )
				)
grow.pmnMS.ninoc <-  sapply(days.s, function(d) 
				mean(DuckGrowPostMS[,1] + DuckGrowPostMS[,2]*d + DuckGrowPostMS[,3]*d^2 + DuckGrowPostMS[,4]*0*d^2 + indivmnefMS + plntmnefMS )
				)
grow.pCIMS.inoc <-  sapply(days.s, function(d) 
				HPDinterval(as.mcmc(DuckGrowPostMS[,1] + DuckGrowPostMS[,2]*d + DuckGrowPostMS[,3]*d^2 + DuckGrowPostMS[,4]*1*d^2 + indivmnefMS + plntmnefMS ),prob=0.95)
				)
grow.pCIMS.ninoc <-  sapply(days.s, function(d) 
				HPDinterval(as.mcmc(DuckGrowPostMS[,1] + DuckGrowPostMS[,2]*d + DuckGrowPostMS[,3]*d^2 + DuckGrowPostMS[,4]*0*d^2 + indivmnefMS + plntmnefMS ),prob=0.95)
				)
# "(Intercept)"         "days"                "I(days^2)"           "I(days^2):bzt"       "I(days^2):inocy"     "I(days^2):bzt:inocy" 
grow.pmnNS.inoc <- lapply(bzt.s, function(b) sapply(days.s, function(d) 
				mean(DuckGrowPostNS[,1] + DuckGrowPostNS[,2]*d + DuckGrowPostNS[,3]*d^2 + DuckGrowPostNS[,4]*d^2*b + DuckGrowPostNS[,5]*d^2*1 + DuckGrowPostNS[,6]*d^2*b*1 + indivmnefNS + plntmnefNS ))
				)
grow.pmnNS.ninoc <- lapply(bzt.s, function(b) sapply(days.s, function(d) 
				mean(DuckGrowPostNS[,1] + DuckGrowPostNS[,2]*d + DuckGrowPostNS[,3]*d^2 + DuckGrowPostNS[,4]*d^2*b + DuckGrowPostNS[,5]*d^2*0 + DuckGrowPostNS[,6]*d^2*b*0 + indivmnefNS + plntmnefNS ))
				)
grow.pCINS.inoc <- lapply(bzt.s, function(b) sapply(days.s, function(d) 
				HPDinterval(as.mcmc(DuckGrowPostNS[,1] + DuckGrowPostNS[,2]*d + DuckGrowPostNS[,3]*d^2 + DuckGrowPostNS[,4]*d^2*b + DuckGrowPostNS[,5]*d^2*1 + DuckGrowPostNS[,6]*d^2*b*1 + indivmnefNS + plntmnefNS ),prob=0.95))
				)
grow.pCINS.ninoc <- lapply(bzt.s, function(b) sapply(days.s, function(d) 
				HPDinterval(as.mcmc(DuckGrowPostNS[,1] + DuckGrowPostNS[,2]*d + DuckGrowPostNS[,3]*d^2 + DuckGrowPostNS[,4]*d^2*b + DuckGrowPostNS[,5]*d^2*0 + DuckGrowPostNS[,6]*d^2*b*0 + indivmnefNS + plntmnefNS ),prob=0.95))
				)
lbzt.s <- sort(unique(sizetime$lbzt))
bzt.s <- sort(unique(sizetime$bzt))
 lnacl.s <- sort(unique(sizetime$lnacl))


trtgropixhs <-tapply(log(sizetime.aliveHS$pix + min(sizetime.aliveHS$pix[sizetime.aliveHS$pix>0] )  ),paste(sizetime.aliveHS$days, sizetime.aliveHS$inoc),mean,na.rm=T)
trtgropixhsse <- tapply(log(sizetime.aliveHS$pix + min(sizetime.aliveHS$pix[sizetime.aliveHS$pix>0]) ),paste(sizetime.aliveHS$days, sizetime.aliveHS$inoc),sd,na.rm=T) / sqrt( table( paste(paste(sizetime.aliveHS$days, sizetime.aliveHS$inoc)) ) )
subfromhspix <-  min(sizetime.aliveHS$pix[sizetime.aliveHS$pix>0] )
trtgropixms <-tapply(sizetime.aliveMS$pix,paste(sizetime.aliveMS$days, sizetime.aliveMS$inoc),mean,na.rm=T)
trtgropixmsse <- tapply(sizetime.aliveMS$pix,paste(sizetime.aliveMS$days, sizetime.aliveMS$inoc),sd,na.rm=T) / sqrt( table( paste(paste(sizetime.aliveMS$days, sizetime.aliveMS$inoc)) ) )
trtgrodaysxs <- c(0,0.5,10,10.5,4,4.5,7,7.5)
trtgropixns <-tapply(sizetime.aliveNS$pix[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],mean,na.rm=T)
trtgropixnsse <-tapply(sizetime.aliveNS$pix[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],sd,na.rm=T) / sqrt( table( paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0] ) )
trtgrodaysns <-tapply(sizetime.aliveNS$days[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],mean,na.rm=T)
trtgrolbztns <-tapply(sizetime.aliveNS$lbzt[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],mean,na.rm=T)
trtgrobztns <-tapply(sizetime.aliveNS$bzt[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],paste(sizetime.aliveNS$days, sizetime.aliveNS$lbzt, sizetime.aliveNS$inoc)[sizetime.aliveNS$bzt== 10 | sizetime.aliveNS$bzt==0],mean,na.rm=T)
jtrtgrodaysns <- trtgrodaysns + rep(c(-0.6,-0.2,0.2,0.6),times=4)
pdf("~/timeseriesPixFitAS.pdf",width=7,height=3)
layout(matrix(1:4,ncol=4),widths=c(1,1,1,0.5))
par(oma=c(3,1,1,0))
par(mar=c(1,3,3,0))
plot(trtgropixns ~ jtrtgrodaysns, ylab="",xlab="",ylim=c(0,5000),xlim=c(-1,11),main="NaCl 0 g/L",pch=NA)
	sapply(c(1,4), function(b) lines(grow.pmnNS.inoc[[b]]~days.s,	col=rgb(0,0,range01(lbzt.s)[b])) ) 
	sapply(c(1,4), function(b) lines(grow.pmnNS.ninoc[[b]]~days.s,	col=rgb(0,0,range01(lbzt.s)[b]) ,lty=2) )
	sapply(c(1,4), function(b) polygon(c(days.s,rev(days.s)),c(grow.pCINS.inoc[[b]][1,],rev(grow.pCINS.inoc[[b]][2,])),	col=rgb(0,0,range01(lbzt.s)[b],alpha=0.35),border=NA ) ) 
	sapply(c(1,4), function(b) polygon(c(days.s,rev(days.s)),c(grow.pCINS.ninoc[[b]][1,],rev(grow.pCINS.ninoc[[b]][2,])),	border=rgb(0,0,range01(lbzt.s)[b],alpha=1),col=NA ,lty=3  )) 
	points(trtgropixns ~ jtrtgrodaysns, pch=rep(c(1,16),times=length(trtgropixns)/2),cex=2,
					col=rgb(0,0,range01(trtgrolbztns),alpha=0.75))
	arrows(x0 = jtrtgrodaysns, y0= trtgropixns - trtgropixnsse, y1= trtgropixns + trtgropixnsse,length=0,col=rgb(0,0,range01(trtgrolbztns),alpha=0.75))
	mtext("a.", side=3,line=2.5,adj = 0)
	mtext("Pixel area",side=2,line=2.5)
plot(trtgropixms~ trtgrodaysxs, pch=NA, ylab="",xlab="",ylim=c(0,5000),xlim=c(-1,11),main="NaCl 0.8 g/L")
	lines(grow.pmnMS.inoc~days.s,	col=rgb(range01(lnacl.s)[2],0,0) ) 
	lines(grow.pmnMS.ninoc~days.s,	col=rgb(range01(lnacl.s)[2],0,0),lty=2) 
	polygon(c(days.s,rev(days.s)),c(grow.pCIMS.inoc[1,],rev(grow.pCIMS.inoc[2,])),	col=rgb(range01(lnacl.s)[2],0,0,alpha=0.35),border=NA ) 
	polygon(c(days.s,rev(days.s)),c(grow.pCIMS.ninoc[1,],rev(grow.pCIMS.ninoc[2,])),border=rgb(range01(lnacl.s)[2],0,0,alpha=1),col=NA,lty=3 ) 
	points(trtgropixms~ trtgrodaysxs, pch=rep(c(1,16),times=length(trtgropixms)/2),cex=2,
					col=rgb(range01(lnacl.s)[2],0,0,alpha=0.5))
	arrows(x0 = trtgrodaysxs, y0= trtgropixms - trtgropixmsse, y1= trtgropixms + trtgropixmsse,length=0,col=rgb(range01(lnacl.s)[2],0,0,alpha=0.75))
	mtext("b.", side=3,line=2.5,adj = 0)
	mtext("Days",side=1,line=2.5)
plot(exp(trtgropixhs)~ trtgrodaysxs, pch = NA, ylab="",xlab="",ylim=c(0,5000),xlim=c(-1,11),main="NaCl 0.8 g/L")
	lines(grow.pmnHS.inoc~days.s,col=rgb(1,0,0) ,lwd=1 )
	lines(grow.pmnHS.ninoc~days.s,col=rgb(1,0,0) ,lwd=1 ,lty=2)
	polygon(c(days.s,rev(days.s)),c(grow.pCIHS.inoc[1,],rev(grow.pCIHS.inoc[2,])),col=rgb(1,0,0,alpha=0.5),border=NA ) 
	polygon(c(days.s,rev(days.s)),c(grow.pCIHS.ninoc[1,],rev(grow.pCIHS.ninoc[2,])),col=NA,border= rgb(1,0,0,alpha=0.35),lty=3) 
	points((exp(trtgropixhs) -subfromhspix)~ trtgrodaysxs, pch=rep(c(1,16),times=length(trtgropixhs)/2),cex=2,
					col=rgb(range01(lnacl.s)[3],0,0,alpha=0.5))
	arrows(x0 = trtgrodaysxs, y0= exp(trtgropixhs - trtgropixhsse)-subfromhspix, y1= exp(trtgropixhs + trtgropixhsse)-subfromhspix,length=0,col=rgb(range01(lnacl.s)[3],0,0,alpha=0.75))
	mtext("c.", side=3,line=2.5,adj = 0)
	snames <- c("0","0","0.8","10")
	bnames <- c("0","10","all","all")
par(mar=c(1,0,3,0))
plot(c(6:3,6:3)~rep(c(1,2),each=4),ylab="",xlab="",yaxt="n",xaxt="n",pch=rep(c(1,19),each=4),bty="n", ylim=c(1,10), xlim=c(0,5), cex=1.2,
		col=  	c( rgb(range01(lnacl.s)[1],0,range01(lbzt.s),alpha=0.75)[c(1,4)], rgb(range01(lnacl.s)[2],0,0,,alpha=0.75), rgb(1,0,0,alpha=0.75) ) )
	arrows(x0=0.25,x1=1.5,y0=c(6:3)-0.25, lty=2, length=0,
		col=c( rgb(range01(lnacl.s)[1],0,range01(lbzt.s))[c(1,4)], rgb(range01(lnacl.s)[2],0,0), rgb(1,0,0) ))
	arrows(x0=1.75,x1=2.75,y0=c(6:3)-0.25, lty=1, length=0,
		col=c( rgb(range01(lnacl.s)[1],0,range01(lbzt.s))[c(1,4)], rgb(range01(lnacl.s)[2],0,0), rgb(1,0,0) ))
	text(x=3,y=6:3,snames,adj=c(0,0.5))
	text(x=4,y=6:3,bnames,adj=c(0,0.5))
	text(c("Not inoculated","Inoculated","NaCl g/L","BZT mg/L"),x=c(0.75,1.75,3,4),y=6.5,adj=c(0,1),srt=90)
dev.off()	

#to for discussion of relative sizes. Using averages calculated above for growth and suvrvival, but also the following:
#view HS mns
exp(trtgropixhs) -subfromhspix
 exp(trtgropixhs - trtgropixhsse)-subfromhspix
 exp(trtgropixhs + trtgropixhsse)-subfromhspix
(exp(trtgropixhs) -subfromhspix)/(exp(trtgropixhs) -subfromhspix)[c(1,2,1,2,1,2,1,2)]
( exp(trtgropixhs - trtgropixhsse)-subfromhspix)/(exp(trtgropixhs) -subfromhspix)[c(1,2,1,2,1,2,1,2)]
( exp(trtgropixhs + trtgropixhsse)-subfromhspix)/(exp(trtgropixhs) -subfromhspix)[c(1,2,1,2,1,2,1,2)]
	
#genotype effects	
tapply(bzs.dattrt$aliveend,bzs.dattrt$plnt.micr, mean,na.rm=T)
tapply(bzs.dattrt$pixend[bzs.dattrt$aliveend==1],bzs.dattrt$plnt.micr[bzs.dattrt$aliveend==1], mean,na.rm=T) #not correcting for starting differences.
	
	
	
####BZT DEGRADATAION analysis and figure
degrade <- read.csv("~/BZT in BZS1 for trts pooled within geno.csv")
whichwells <- read.csv("~/ORB Samples all genos.csv")
growthdegdat <- cbind(whichwells,bzs.dattrt[ paste(trts$row,trts$col,trts$plate)%in%paste(whichwells$row,whichwells$col,whichwells$plate),])
#order of growthdegdat and poolgrowth by trt/geno SHOULD be the same. good to always double check
poolgrowth <- data.frame(bzt = rep(rep(c(0.1,10),each=4),times=3), nacl = rep(rep(c(0,10),each=2),times=6), inoc=rep(c("n","y"),times=12),plnt.micr=rep(c("BruceH","KSR","MoccTr"),each=8),
			simmeas.bztstart = rep(rep(c(0.101580039208,10.1580039208),each=4),times=3) ,# these are measurement from the starting 10mgL concentration measured after storage.
			px.mn = tapply(growthdegdat$pixend, paste(growthdegdat$plnt.micr,growthdegdat$bzt,growthdegdat$nacl,growthdegdat$inoc),mean,na.rm=T),
			px.se = tapply(growthdegdat$pixend, paste(growthdegdat$plnt.micr,growthdegdat$bzt,growthdegdat$nacl,growthdegdat$inoc),sd,na.rm=T)/sqrt(3),
			od.mn = tapply(growthdegdat$flOD600, paste(growthdegdat$plnt.micr,growthdegdat$bzt,growthdegdat$nacl,growthdegdat$inoc),mean,na.rm=T),
			od.se = tapply(growthdegdat$flOD600, paste(growthdegdat$plnt.micr,growthdegdat$bzt,growthdegdat$nacl,growthdegdat$inoc),sd,na.rm=T)/sqrt(3))
poolgrowth$bzt.post <- degrade$BZT.pool
poolgrowth$bzt.propfinal <- poolgrowth$bzt.post/poolgrowth$simmeas.bztstart
write.csv(poolgrowth,"~/BZTinBZS1_trts_within_geno_plusgrowth.csv",row.names=F)


degmodfullprop.rand <- MCMCglmm(bzt.propfinal ~ inoc +bzt + nacl, random= ~plnt.micr,family="gaussian",data=poolgrowth,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)
degmodpix.prR <-   MCMCglmm(bzt.propfinal ~ px.mn, random= ~plnt.micr ,family="gaussian",data=poolgrowth,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)#
degmodod.prR <-   MCMCglmm(bzt.propfinal ~ od.mn, random= ~plnt.micr ,family="gaussian",data=poolgrowth,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)#
degmodpixod.prRAND <-   MCMCglmm(bzt.propfinal ~ od.mn + px.mn ,random=~plnt.micr ,family="gaussian",data=poolgrowth,verbose=FALSE,pr=TRUE,nitt=100000, thin=100, burnin=2000)#

bzprpx.mn = tapply(growthdegdat$pixend, paste(growthdegdat$nacl,growthdegdat$bzt,growthdegdat$inoc),mean,na.rm=T)
bzprpx.se = tapply(growthdegdat$pixend, paste(growthdegdat$nacl,growthdegdat$bzt,growthdegdat$inoc),sd,na.rm=T)/sqrt(9)
bzprod.mn = tapply(growthdegdat$flOD600, paste(growthdegdat$nacl,growthdegdat$bzt,growthdegdat$inoc),mean,na.rm=T)
bzprod.se = tapply(growthdegdat$flOD600, paste(growthdegdat$nacl,growthdegdat$bzt,growthdegdat$inoc),sd,na.rm=T)/sqrt(9)

bzprpmn <- tapply(poolgrowth$bzt.propfinal, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc),mean)
bzprpse <- tapply(poolgrowth$bzt.propfinal, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc),sd) / sqrt(table(paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc))) #no NAs
bzprpbzt <- tapply(poolgrowth$bzt, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc),mean)
bzprpnacl <- tapply(poolgrowth$nacl, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc),mean)
bzprpinoc <- tapply(as.numeric(poolgrowth$inoc)-1, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$inoc),mean)
bzprpmnpm <- tapply(poolgrowth$bzt.propfinal, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr),mean)
bzprpsepm <- tapply(poolgrowth$bzt.propfinal, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr),sd) / sqrt(table(paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr))) #no NAs
bzprpbztpm <- tapply(poolgrowth$bzt, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr),mean)
bzprpnaclpm <- tapply(poolgrowth$nacl, paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr),mean)
bzprppmpm <- tapply(as.numeric(poolgrowth$plnt.micr), paste(poolgrowth$nacl,poolgrowth$bzt,poolgrowth$plnt.micr),mean)
bzprpmnpmO <- tapply(poolgrowth$bzt.propfinal, poolgrowth$plnt.micr,mean)
bzprpsepmO <- tapply(poolgrowth$bzt.propfinal, poolgrowth$plnt.micr,sd) / sqrt(table(poolgrowth$plnt.micr)) #no NAs
bzprppmpmO <- tapply(as.numeric(poolgrowth$plnt.micr), poolgrowth$plnt.micr,mean)
pixodsol <- degmodpixod.prRAND$Sol
mnplntef <- rowMeans(pixodsol[,4:6])
pgpxmn <- mean(poolgrowth$px.mn)
pgodmn <- mean(poolgrowth$od.mn)
 odseq <- seq(from = min(poolgrowth$od.mn),to=max(poolgrowth$od.mn), length.out=1000)
 pxseq <- seq(from = min(poolgrowth$px.mn),to=max(poolgrowth$px.mn), length.out=1000)
trypxseqCI <-sapply(pxseq, function(z) HPDinterval(as.mcmc(pixodsol[,1] + pixodsol[,2]*pgodmn + pixodsol[,3]*z + mnplntef),prob=0.95))
trypxseqmn <-sapply(pxseq, function(z) mean(pixodsol[,1] + pixodsol[,2]*pgodmn + pixodsol[,3]*z + mnplntef) )
tryodseqCI <-sapply(odseq, function(z) HPDinterval(as.mcmc(pixodsol[,1] + pixodsol[,2]*z + pixodsol[,3]*pgpxmn + mnplntef),prob=0.95))
tryodseqmn <-sapply(odseq, function(z) mean(pixodsol[,1] + pixodsol[,2]*z + pixodsol[,3]*pgpxmn + mnplntef) )

pdf("~/treatmenteffsBZTloss.pdf",width=5.5,height=2.5)
layout(matrix(1:4,ncol=4,nrow=1,byrow=T),widths=c(2,2,2,1.5))
par(oma=c(0,1,1,1))
par(mar=c(5,5,1,0))
plot(bzprpmn*100~as.numeric(as.factor(paste(bzprpnacl, bzprpbzt))), pch=21, #bty="n",
			col = rgb(ifelse(bzprpinoc,0.75,1)*bzprpnacl/10,0, ifelse(bzprpinoc,0.75,1)*floor(bzprpbzt)/10), 
			bg = rgb(ifelse(bzprpinoc,0.75,1)*bzprpnacl/10,0,ifelse(bzprpinoc,0.75,1)*floor(bzprpbzt)/10,alpha=0.75*ifelse(bzprpinoc,1,0.5)),
			xaxt="n",ylab="% BZT remaining",xlab="",xlim=c(0.75,4.25),ylim=c(0,100),cex=2,cex.lab=1.5)#BZT mg/L
 mtext("a.",side=3, line = 0.5 ,adj=0)
 arrows(x0=as.numeric(as.factor(paste(bzprpnacl, bzprpbzt))),
		y0= (bzprpmn- bzprpse)*100, y1 =(bzprpmn + bzprpse)*100,length=0,
		col = rgb(ifelse(bzprpinoc,0.75,1)*bzprpnacl/10,0, ifelse(bzprpinoc,0.75,1)*floor(bzprpbzt)/10))
par(mar=c(5,0,1,0))
plot(poolgrowth$bzt.propfinal*100~poolgrowth$px.mn,ylab="",yaxt="n",xlab="Pixel area", ylim=c(0,100) ,cex=1,cex.lab=1.5,
	pch=21, col = rgb(ifelse(poolgrowth$inoc=="n",1,0.75)*poolgrowth$nacl/10,0,ifelse(poolgrowth$inoc=="n",1,0.75)*floor(poolgrowth$bzt)/10), 
	bg = rgb(ifelse(poolgrowth$inoc=="n",1,0.75)*poolgrowth$nacl/10,0,ifelse(poolgrowth$inoc=="n",1,0.75)*floor(poolgrowth$bzt)/10,alpha= 0.75*ifelse(poolgrowth$inoc=="y",1,0.5) ))
 lines(trypxseqmn*100~pxseq,lwd=1.5)
 polygon(x=c(pxseq, rev(pxseq)), y= c(trypxseqCI[1,]*100,rev(trypxseqCI[2,]*100)) ,col=rgb(0,0,0,alpha=0.1),border=NA)
 mtext("b.",side=3, line = 0.5 ,adj=0)
plot(poolgrowth$bzt.propfinal*100~poolgrowth$od.mn,ylab="",xlab="log(optical density)",ylim=c(0,100) ,yaxt="n",cex=1,cex.lab=1.5,
	pch=21, col = rgb(ifelse(poolgrowth$inoc=="n",1,0.75)*poolgrowth$nacl/10,0,ifelse(poolgrowth$inoc=="n",1,0.75)*floor(poolgrowth$bzt)/10), 
	bg = rgb(ifelse(poolgrowth$inoc=="n",1,0.75)*poolgrowth$nacl/10,0,ifelse(poolgrowth$inoc=="n",1,0.75)*floor(poolgrowth$bzt)/10,alpha= 0.75*ifelse(poolgrowth$inoc=="y",1,0.5) ))
 lines(tryodseqmn*100 ~ odseq,lwd=1.5)
 polygon(x=c(odseq, rev(odseq)), y= c(tryodseqCI[1,]*100,rev(tryodseqCI[2,]*100)) ,col=rgb(0,0,0,alpha=.1),border=NA)
 mtext("c.",side=3, line = 0.5 ,adj=0)
plot(rep(4:1,times=2)~rep(c(1.1,2.1),each=4), pch=21, yaxt="n",xaxt="n", xlim=c(-1,6),ylim=c(0,9), cex=1.5,xlab="",
		col = rgb(c(0,0,1,1,0,0,0.75,0.75),0,c(0,1,0,1,0,0.75,0,0.75)), 
		bg = rgb(c(0,0,1,1,0,0,0.75,0.75),0,c(0,1,0,1,0,0.75,0,0.75),alpha=rep(c(0.375,0.75),each=4)),  bty="n")
text(x=c(2.75),y=c(4:1),c("0.1", "10","0.1","10"),adj=c(0,0.5))
text(x=c(4),y=c(4:1),c("0","0","10","10"),adj=c(0,0.5))
text(x=c(1,2,3,4.25),y=4.5,c("Not Inoculated","Inoculated","BZT mg/L","NaCl g/L"),adj=c(0,0.5),srt=90)

dev.off()

##numbers for genotype differences
tapply(poolgrowth$bzt.propfinal,poolgrowth$plnt.micr,mean)



