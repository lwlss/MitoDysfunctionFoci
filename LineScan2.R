library(caTools)

# https://twitter.com/robustgar/status/798962929475457024
effectsize=function(x,y){
 m <- outer(x,y,FUN="-")
 qxly <- sum(m<0)/length(m)
 return(qxly)
}

# https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

hypertest=function(x,m,y,n) {
  1-phyper(x-1,y,n-y,m)
}

binomtest=function(x,m,w) {
  1-pbinom(x-1,m,w)
}

# https://stackoverflow.com/questions/5012516/count-how-many-consecutive-values-are-true
f7 <- function(x){ tmp<-cumsum(x);tmp-cummax((!x)*tmp)}

findRegions=function(dt,DAPI_cutoff,Ratio_cutoff,Diff_cutoff,focdef="Ratio",wind=25,minFlen=40,maxFlen=9999999999){
  dt$Perinucleus=dt$DAPI>DAPI_cutoff
  dt$PerinuclearRegion=runmax(dt$Perinucleus,wind)>0
  if(focdef=="Ratio") dt$MitoDysfunction=dt$Ratio<Ratio_cutoff
  if(focdef=="Diff") dt$MitoDysfunction=dt$Diff>Diff_cutoff
  dt$Focus=runmax(dt$MitoDysfunction,wind)>0
  dt$Run=f7(dt$Focus)
  runpeaks = which(diff(sign(diff(dt$Run)))==-2)+1
  peakheights = dt$Run[runpeaks]
  runstarts = runpeaks - peakheights + 1 
  for (i in seq_along(runstarts)){
    dt$Run[runstarts[i]:runpeaks[i]]=peakheights[i]
  } 
  dt$Focus[(dt$Run<minFlen)|(dt$Run>maxFlen)]=FALSE
  
  dt$Region="Subsarcolemmal"
  dt$Region[dt$PerinuclearRegion]="Perinuclear"
  return(dt)
}

makeTest=function(dat,var,wind){
 testCutoff = function(cutoff){
  pos = runmax(dat[[var]][dat$CellType=="Pos"]<cutoff,wind)
  neg = runmax(dat[[var]][dat$CellType=="Neg"]>cutoff,wind)
  num = ifelse(var=="Ratio",sum(1-pos)+sum(1-neg),sum(pos)+sum(neg))
  return(sum(num)/(2*length(pos)))
 }
return(testCutoff)
}

croot="LineScans"
dat1=read.delim(paste(croot,".txt",sep=""),sep="\t",stringsAsFactors=FALSE)
dat1$DAPI.1=-99
dat1$MTCOI.1=-99
dat1$SDHA.1=-99
dat1$FociDeficient=NA
dat1$Patient.no=NA
dat1$Foi.no=NA
dat1$Filename=croot
colnames(dat1)[2]="Distance"

#froot="Focix10"
froot="Combined_manual_foci"
if(froot=="Focix10"){
  focidat=read.delim(paste(froot,".txt",sep=""),sep="\t",stringsAsFactors=FALSE,skip=1)
  focidat$CellLabel=sprintf("Foci%02d",focidat$Foci.no)
  focidat$Foci.no=NULL
  focidat$FociDeficient=NA
  focidat$Patient.no=NA
  focidat$Foi.no=NA
  focidat$Filename=froot
}
if(froot=="Combined_manual_foci"){
  focidat=read.delim(paste(froot,".txt",sep=""),sep="\t",stringsAsFactors=FALSE,skip=0)
  focidat$CellLabel=sprintf("Foci%02d_%02d",focidat$Patient.no,focidat$Foi.no)
  focidat$DAPI.1=-99
  focidat$MTCOI.1=-99
  focidat$SDHA.1=-99
  focidat$Filename=froot
}

dat=rbind(dat1,focidat)
dat=dat[!is.na(dat$DAPI),]

focdef="Diff"

dat$CellType=sub("[^[:alpha:]]+", "", dat$CellLabel)
dat$CellLabel=gsub("_","S",gsub("Foci","FociP",dat$CellLabel))
nd=c()
for(ct in unique(dat$CellType)){
  labs=sort(unique(dat$CellLabel[dat$CellType==ct]))
  ndict=seq_along(labs)
  names(ndict)=labs
  nd=c(nd,ndict) 
}
dat$RepNo=nd[dat$CellLabel]

# Remove first round of Foci measurements (data from focidat should replace existing data, not add to it)
dat=dat[!((dat$CellType=="Foci")&(dat$Filename==croot)),]

# Subtract background from DAPI signal
# Note that need to filter saturated DAPI values (intensity of 4095) as this skews mode
DAPIbgrnds = aggregate(dat$DAPI,by=list(dat$CellLabel),FUN = function(x) Mode(x[x<4095]))
DAPIdict = DAPIbgrnds$x
names(DAPIdict)= DAPIbgrnds$Group.1
dat$DAPI = dat$DAPI-DAPIdict[dat$CellLabel]

dat$Diff = dat$SDHA-dat$MTCOI
dat$Ratio = dat$MTCOI/dat$SDHA
delta = dat$Distance[2]-dat$Distance[1]

DAPI_cutoff = quantile(dat$DAPI,0.85)
Diff_cutoff = quantile(dat$Diff[dat$CellType=="Pos"],0.998)
Ratio_cutoff = quantile(dat$Ratio[dat$CellType=="Pos"],0.001)
wind = 24

cutDiff = makeTest(dat,"Diff",wind)
cuts = seq(0,max(dat$Diff),length.out=5000)
res = sapply(cuts,cutDiff)

cutRatio = makeTest(dat,"Ratio",wind)
cuts = seq(0,max(dat$Ratio),length.out=5000)
res = sapply(cuts,cutRatio)

dat=findRegions(dat,DAPI_cutoff,Ratio_cutoff,Diff_cutoff,focdef=focdef)

boxplot(Diff~CellType,data=dat)

alpha=0.25
cols=c(rgb(1,0,0,alpha),rgb(0,0,1,alpha),rgb(0,1,0,alpha))
names(cols)=c("Pos","Neg","Foci")
dat$Colour=cols[dat$CellType]

# Examine distribution of manually classified focus lengths
fl=read.delim("FocusLengths.txt",sep="\t",stringsAsFactors=FALSE)
hist(fl$Focus,xlim=c(0,1.1*max(fl$Focus)))
delt = dat$Distance[2]-dat$Distance[1] # Distance between observations (um)
frange=quantile(fl$Focus,c(0.025,0.975))/delt
fmin=frange[1]
fmax=frange[2]

head(dat)

mkplt = function(mlab="",foctype="Focus",counter=1){
  dt = findRegions(dt,DAPI_cutoff,Ratio_cutoff,Diff_cutoff,focdef=focdef,wind=wind,minFlen=fmin,maxFlen=9999999999999999)
  Ncells=length(unique(dt$CellLabel))
  dt$SumDist = dt$Distance
  if(length(unique(dt$RepNo))>1){
   for(n in 2:max(dt$RepNo)){
     dt$SumDist[dt$RepNo==n]=dt$SumDist[dt$RepNo==n]+max(dt$SumDist[dt$RepNo==(n-1)])
   }
  }else{dt$SumDist=dt$Distance}
  starts=dt[(dt$Distance==0.0)&(dt$SumDist>0),]

  if(mlab=="") {
     mainlab=paste(ct," (N = ",Ncells,")",sep="")
   }else{
     P=sum(dt$PerinuclearRegion)
     S=sum(!dt$PerinuclearRegion)
     F=sum(dt[[foctype]])
     FnP=sum(dt[[foctype]]&dt$PerinuclearRegion)  
     focal_fraction=F/(P+S)
     nuclear_fraction=P/(P+S)
     overlap_fraction=FnP/(P+S)
     mainlab=paste(mlab,"Focal frac.:",formatC(focal_fraction,3),"Perinuclear frac.:",formatC(nuclear_fraction,3),"Obs. overlap:",formatC(overlap_fraction,3),"Pred. overlap:",formatC(focal_fraction*nuclear_fraction,3))
   }
  
  plot(dt$SumDist,dt$DAPI,type="n",main=mainlab,xlab="",ylab="",ylim=c(0,max(dat$DAPI)))
  abline(v=starts$SumDist,col="black",lty=2,lwd=1)
  polygon(c(dt$SumDist[1],dt$SumDist,dt$SumDist[length(dt$SumDist)],dt$SumDist[1]),99999999*c(0,dt$PerinuclearRegion,0,0),col=rgb(0,0,1,0.2),border=NA)
  polygon(c(dt$SumDist[1],dt$SumDist,dt$SumDist[length(dt$SumDist)],dt$SumDist[1]),99999999*c(0,dt[[foctype]],0,0),col=rgb(1,0,0,0.2),border=NA)
  points(dt$SumDist,dt$DAPI,type="l",col="blue")
  mtext(side = 2, line = 2.5, "DAPI")
  par(new=TRUE)
  plot(dt$SumDist,pmax(0,dt[[focdef]]),type="l",col="red",axes=FALSE,ann=FALSE,ylim=c(0,max(dat[[focdef]])))
  if(foctype=="Focus") abline(h = ifelse(focdef=="Ratio",Ratio_cutoff,Diff_cutoff),col="green",lwd=1)
  axis(side = 4)
  mtext(side = 4, line = 3, ifelse(focdef=="Ratio","MTCOI/SDHA","SDHA-MTCOI"))
  if(counter%%3==0) mtext(side = 1, line = 2.5, "Distance along fibre section perimeter (um)")
  legend("topright",lwd=1,legend=c("DAPI",ifelse(focdef=="Ratio","MTCOI/SDHA","SDHA-MTCOI")),col=c("blue","red"),bg="white")
}

#foctype="Focus" # Automatic focus classification
foctype="FociDeficient" # Manual focus classification

# Overview of all sections
pdf(paste(froot,ifelse(foctype=="Focus","_OVERVIEW_auto_class.pdf","_OVERVIEW_manual_class.pdf"),sep=""),width=16,height=8)

# Where should we segment Diff?
dneg=density(dat$Diff[dat$CellType=="Neg"])
dpos=density(dat$Diff[dat$CellType=="Pos"])

op=par(mfrow=c(3,1),mai=c(0.25,0.75,0.35,1),oma=c(3,0,0,0))
cts=unique(dat$CellType)
for(i in seq_along(cts)){
  ct=cts[i]
  dt = dat[dat$CellType==ct,]
  mkplt(foctype=foctype,counter=i)
}
par(op)
#mtext(side = 1, line = 1.5, "Distance along fibre section perimeter (um)")
dev.off()

# Plot each section separately
pdf(paste(froot,ifelse(foctype=="Focus","_SEPARATE_auto_class.pdf","_SEPARATE_manual_class.pdf"),sep=""),width=16,height=8)

# Predicted-observed test
df = dat[dat$CellType=="Foci",]
section_ids=sort(unique(df$CellLabel))
overlap_fractions=rep(0,length(section_ids))
focal_fractions=rep(0,length(section_ids))
nuclear_fractions=rep(0,length(section_ids))
for(i in seq_along(section_ids)){
  dt=df[df$CellLabel==section_ids[i],]
  P=sum(dt$PerinuclearRegion)
  S=sum(!dt$PerinuclearRegion)
  F=sum(dt[[foctype]])
  FnP=sum(dt[[foctype]]&dt$PerinuclearRegion)  
  focal_fractions[i]=F/(P+S)
  nuclear_fractions[i]=P/(P+S)
  overlap_fractions[i]=FnP/(P+S)
}

predicted_fractions=focal_fractions*nuclear_fractions
maxfrac=max(c(overlap_fractions,predicted_fractions))
patientno=as.numeric(substr(section_ids,6,7))

# Summary
pdf(paste(froot,ifelse(foctype=="Focus","_SUMMARY_auto_class.pdf","_SUMMARY_manual_class.pdf"),sep=""),width=10,height=5)
op=par(mfrow=c(1,2))
plot(NULL,xlab="Observed overlap between focal and perinuclear regions",ylab="Predicted overlap, assuming random locations",xlim=c(0,maxfrac),ylim=c(0,maxfrac),main="Focal fibre section perimeters")
abline(a=0,b=1,col="darkgrey",lwd=3,lty=2)
clist=c("black","red","blue")
points(overlap_fractions,predicted_fractions,pch=16,cex=1.25,col=clist[patientno])
#text(overlap_fractions,predicted_fractions,gsub("Foci","",section_ids),cex=0.35,pos=2,offset=0.25)
legend("bottomright",col=clist[unique(patientno)],legend=paste("Patient",unique(patientno)),pch=16,pt.cex=1.25)
res = t.test(predicted_fractions-overlap_fractions)
hist(predicted_fractions-overlap_fractions,breaks=15,xlab="Predicted - Observed",main=paste("p-value:",formatC(res$p.value,4)))
abline(v=0,col="darkgrey",lwd=3,lty=2)
par(op)
dev.off()

# Plot each section separately
pdf(paste(froot,ifelse(foctype=="Focus","_SEPARATE_auto_class.pdf","_SEPARATE_manual_class.pdf"),sep=""),width=16,height=8)

op=par(mfrow=c(3,1),mai=c(0.25,0.75,0.35,1),oma=c(3,0,0,0))
if(foctype=="Focus") clabs = sort(unique(dat$CellLabel))
if(foctype=="FociDeficient") clabs = sort(unique(dat$CellLabel[dat$CellType=="Foci"]))
for(i in seq_along(clabs)){
  cl=clabs[i]
  dt = dat[dat$CellLabel==cl,]
  mkplt(mlab=cl,foctype=foctype,counter=ifelse(i==length(clabs),3,i))
}
#legend("topright",lwd=1,legend=c("DAPI",ifelse(focdef=="Ratio","MTCOI/SDHA","SDHA-MTCOI")),col=c("blue","red"),bg="white")
par(op)
dev.off()

# Plot one section
pdf(paste(froot,ifelse(foctype=="Focus","_INDIVIDUAL_auto_class.pdf","_INDIVIDUAL_manual_class.pdf"),sep=""),width=8,height=8)
op=par(mai=c(0.25,0.75,0.35,1),oma=c(3,0,0,0))
  cl="FociP01S02"
  dt = dat[dat$CellLabel==cl,]
  mkplt(mlab=cl,foctype=foctype,counter=ifelse(i==length(clabs),3,i))
  legend("topright",lwd=1,legend=c("DAPI",ifelse(focdef=="Ratio","MTCOI/SDHA","SDHA-MTCOI")),col=c("blue","red"),bg="white")
par(op)
dev.off()

# Test for over-representation of foci in perinuclear region of muscle fibre perimeters marked "foci": early development of mito deficiency
dt = dat[dat$CellType=="Foci",]
dt = findRegions(dt,DAPI_cutoff,Ratio_cutoff,Diff_cutoff,focdef=focdef,wind=wind,minFlen=fmin,maxFlen=999999999999999)
P=sum(dt$PerinuclearRegion)
S=sum(!dt$PerinuclearRegion)
F=sum(dt[[foctype]])
FnP=sum(dt[[foctype]]&dt$PerinuclearRegion)

hypertest(FnP,F,P,P+S)
binomtest(FnP,F,P/(P+S))


