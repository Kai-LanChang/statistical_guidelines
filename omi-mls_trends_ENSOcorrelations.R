#This example is made for ozone trend analysis of 
#the OMI/MLS satellite data 2004-2021 (Ziemke et al., 2006; 2019)
#Dr Kai-Lan Chang (CIRES & NOAA CSL, 15 Mar 2023)
#setwd("/Users/kai-lanchang/Dropbox/data")
library(dichromat)
library(RColorBrewer)
library(viridis)
library(fields)
library(quantreg)
library(ncdf4)
library(rworldmap)
newworld =getMap(resolution = "li")
difcol=colorschemes$DarkRedtoBlue.18
difcol=c('blue 4',difcol,'red 4')
toar.col=c(rgb(0.2081, 0.1663, 0.5292),rgb(0.3961, 0.3176, 0.8000),
rgb(0.0123, 0.4213, 0.8802),rgb(0.4941, 0.7647, 0.8980),
rgb(0.1157, 0.7022, 0.6843),rgb(0.5216, 0.6980, 0.1725),
rgb(0.9968, 0.7513, 0.2325),rgb(1.0000, 0.4863, 0.0000),
rgb(0.8000, 0.3176, 0.3176),rgb(0.6980, 0.1725, 0.1725))

##moving block bootstrap function
mbfun=function(formula,data,tau){
 #data=rbind(data,data)
 n=nrow(data)
 b=ceiling(n^0.25) 
 nblocks=ceiling(n/b)
 blocks=lapply(seq_len(n-b+1), function(i) seq(i, i+b-1))
 bn=sample(1:length(blocks),nblocks,replace=T)
 sampdata=data[unlist(blocks[bn]), ]  
 mod=rq(formula, data=sampdata, tau=tau)
 coef(mod)
}

set.seed(2013)
#=========================================================
#read data
#=========================================================
#data can be found at https://psl.noaa.gov/enso/mei/
enso=read.table("meiv2.data",header=F, sep="", nrows=45, skip=1)
enso=data.frame(year=rep(seq(1979,1978+nrow(enso)),each=12), month=seq(1,12), enso=as.vector(t(enso)[-1,]))

#data can be found at https://acd-ext.gsfc.nasa.gov/Data_services/cloud_slice/.
sat=ncvar_get(nc_open("tco_omimls.nc"), "TropoO3")
ny=18
lat=seq(-57.5, by=5, length=24)
lon=seq(-177.5, by=5, length=72)

#=========================================================
#enso correlations and trends
#=========================================================
dt=sat[25,12,] #lon -57.5W and lat 2.5S
dt=data.frame(cbind(o3=dt, year=c(rep(2004,3),rep(seq(2005,2021),each=12)), month=c(seq(10,12), rep(1:12, ny-1)), ind=1:length(dt)))
seasonality=predict(lm(o3~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=dt), newdata=data.frame(month=1:12))
dt=merge(dt, data.frame(month=1:12, o3d=seasonality), by="month")
dt$o3d=dt$o3-dt$o3d
dt=merge(dt,enso,by=c("year","month"), sort=T)
dt=dt[order(dt$ind),]
#a cross correlation example
png("omi-mls_lag_ccf_2004-2021.png", height=5, width=5, units ="in", res =300)
par(mar=c(3.5, 3.5, 3.5, 0.5), mgp=c(2.4, 0.8, 0), las=1)
ccf(dt$enso, dt$o3d,main="Lag correlations between enso & ozone", ylab="correlation", xlab="temporal lag")
dev.off()

#simple linear trends
ops_slope=matrix(0,nrow=72,ncol=24)
ops_sigma=matrix(0,nrow=72,ncol=24)
for (j in 1:24){
    for (i in 1:72){
        dt=sat[i,j,]
        dt=data.frame(cbind(o3=dt, year=c(rep(2004,3),rep(seq(2005,2021),each=12)), month=c(seq(10,12), rep(1:12, ny-1)), ind=1:length(dt)))
        seasonality=predict(lm(o3~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=dt), newdata=data.frame(month=1:12))
        dt=merge(dt, data.frame(month=1:12, o3d=seasonality), by="month")
        dt$o3d=dt$o3-dt$o3d
        dt=merge(dt,enso,by=c("year","month"), sort=T)
        dt=dt[order(dt$ind),]
        mbb=t(replicate(1000, mbfun(formula=o3d~ind,data=dt,tau=0.5)))
        ops_slope[i,j]=coef(rq(o3d~ind, data=dt, tau=0.5))[2]*12
        ops_sigma[i,j]=t(apply(mbb, 2, sd, na.rm=TRUE))[2]*12 
    }
print(paste(100*j/24, "% complete", sep=""))
}

#zero-lag enso correlations & trends
op0_enso=matrix(0,nrow=72,ncol=24)
op0_slope=matrix(0,nrow=72,ncol=24)
op0_sigma=matrix(0,nrow=72,ncol=24)
for (j in 1:24){
    for (i in 1:72){
        dt=sat[i,j,]
        dt=data.frame(cbind(o3=dt, year=c(rep(2004,3),rep(seq(2005,2021),each=12)), month=c(seq(10,12), rep(1:12, ny-1)), ind=1:length(dt)))
        seasonality=predict(lm(o3~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=dt), newdata=data.frame(month=1:12))
        dt=merge(dt, data.frame(month=1:12, o3d=seasonality), by="month")
        dt$o3d=dt$o3-dt$o3d
        dt=merge(dt,enso,by=c("year","month"), sort=T)
        dt=dt[order(dt$ind),]
        op0_enso[i,j]=cor(dt$enso, dt$o3d)
        mbb=t(replicate(1000, mbfun(formula=o3d~ind+enso,data=dt,tau=0.5)))
        op0_slope[i,j]=coef(rq(o3d~ind+enso, data=dt, tau=0.5))[2]*12
        op0_sigma[i,j]=t(apply(mbb, 2, sd, na.rm=TRUE))[2]*12 
    }
print(paste(100*j/24, "% complete", sep=""))
}

#maximum enso correlations & trends
lag_fn = function(x, k=1, pad=NA){
  if(k == 0)
    return(x)
  nas <- rep(pad, min(length(x), abs(k)))
  if(k < 0)
    c(tail(x, k), nas) else c(nas, head(x, -k))
}

op_enso=matrix(0,nrow=72,ncol=24)
op_ensolag=matrix(0,nrow=72,ncol=24)
op_slope=matrix(0,nrow=72,ncol=24)
op_sigma=matrix(0,nrow=72,ncol=24)
for (j in 1:24){
    for (i in 1:72){
        dt=sat[i,j,]
        dt=data.frame(cbind(o3=dt, year=c(rep(2004,3),rep(seq(2005,2021),each=12)), month=c(seq(10,12), rep(1:12, ny-1)), ind=1:length(dt)))
        seasonality=predict(lm(o3~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=dt), newdata=data.frame(month=1:12))
        dt=merge(dt, data.frame(month=1:12, o3d=seasonality), by="month")
        dt$o3d=dt$o3-dt$o3d
        dt=merge(dt,enso,by=c("year","month"), sort=T)
        dt=dt[order(dt$ind),]
        pt=ccf(dt$enso, dt$o3d, lag.max=6, pl=F)
        op_enso[i,j]=pt$acf[which.max(abs(pt$acf))]
        op_ensolag[i,j]=pt$lag[which(pt$acf==pt$acf[which.max(abs(pt$acf))], arr.ind=T)]
        dt$ensolag=lag_fn(x=dt$enso, k=op_ensolag[i,j])
        find_row=which(enso$year==dt[is.na(dt$ensolag),]$year[1] & enso$month==dt[is.na(dt$ensolag),]$month[1])
        #fill lag NA data if the ENSO record is available
        if (!is.na(dt[is.na(dt$ensolag),]$year[1])){
        if (dt[is.na(dt$ensolag),]$year[1] <= 2005) {
            dt[is.na(dt$ensolag),"ensolag"]=enso[(find_row-op_ensolag[i,j]):(find_row-1),]$enso} else {
            dt[is.na(dt$ensolag),"ensolag"]=enso[(find_row-op_ensolag[i,j]):(find_row-2*op_ensolag[i,j]-1),]$enso}}
        mbb=t(replicate(1000, mbfun(formula=o3d~ind+ensolag,data=dt,tau=0.5)))
        op_slope[i,j]=coef(rq(o3d~ind+ensolag, data=dt, tau=0.5))[2]*12
        op_sigma[i,j]=t(apply(mbb, 2, sd, na.rm=TRUE))[2]*12 
    }
print(paste(100*j/24, "% complete", sep=""))
}

#=========================================================
#plot results
#=========================================================
png("omi-mls_enso_lag_2004-2021.png", height=5, width=7, units ="in", res =300)
par(mar=c(2.5, 2.5, 1.7, 0.5), mgp=c(2.4, 0.8, 0), las=1)
image.plot(lon,lat, op_ensolag,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, 
    col=c(rgb(40/255,27/255,13/255), brewer.pal(11,"BrBG"), rgb(0,30/255,0)), breaks=seq(-6.5,6.5,by=1))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Optimized lag (with the peak correlation)",cex=1.3, line=0.5, at=0)
dev.off()


png("omi-mls_enso_cor_2004-2021.png", height=8, width=7, units ="in", res =300)
par(mar=c(2, 2, 0.5, 0.5), mgp=c(2.4, 0.8, 0), las=1, mfrow=c(2,1),  oma=c(0.5,0.5,1.5,0.5))
image.plot(lon,lat, op0_enso,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=colorschemes$GreentoMagenta.16, xaxt="n",
    breaks=seq(-0.8,0.8,length.out=17), lab.breaks=round(seq(-0.8, 0.8, length.out=17),2))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Zero-lag ENSO correlation",cex=1.3, line=0.5, at=0)

image.plot(lon,lat, op_enso,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=colorschemes$GreentoMagenta.16, 
    breaks=seq(-0.8,0.8,length.out=17), lab.breaks=round(seq(-0.8, 0.8, length.out=17),2))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Peak ENSO correlation",cex=1.3, line=0.5, at=0)
dev.off()


png("omi-mls_trends_50th_2004-2021.png", height=8, width=7, units ="in", res =300)
par(mar=c(2, 2, 0.5, 0.5), mgp=c(2.4, 0.8, 0), las=1, mfrow=c(3,1),  oma=c(0.5,0.5,1.5,0.5))
pv1=matrix(0,nrow=72,ncol=24)
pv2=matrix(0,nrow=72,ncol=24)
for (i in 1:72){
    for (j in 1:24){
        pv1[i,j]=ifelse(abs(ops_slope[i,j])>=2*ops_sigma[i,j] & abs(ops_slope[i,j])<3*ops_sigma[i,j],1,NA)
        pv2[i,j]=ifelse(abs(ops_slope[i,j])>=3*ops_sigma[i,j],1,NA)
    }
}
pp1=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv1)))
pp2=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv2)))

image.plot(lon,lat,10*ops_slope,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=difcol, xaxt="n",
    breaks=seq(-5,5,length.out=21), lab.breaks=round(seq(-5, 5, length.out=21),2))
plot(newworld, add=TRUE, lwd=1.5)
points(pp1[,1:2],pch=3,col='azure4')
points(pp2[,1:2],pch=8,col='azure4')
mtext("50th trends (wo/ ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)

pv1=matrix(0,nrow=72,ncol=24)
pv2=matrix(0,nrow=72,ncol=24)
for (i in 1:72){
    for (j in 1:24){
        pv1[i,j]=ifelse(abs(op0_slope[i,j])>=2*op0_sigma[i,j] & abs(op0_slope[i,j])<3*op0_sigma[i,j],1,NA)
        pv2[i,j]=ifelse(abs(op0_slope[i,j])>=3*op0_sigma[i,j],1,NA)
    }
}
pp1=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv1)))
pp2=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv2)))

image.plot(lon,lat,10*op0_slope,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=difcol, xaxt="n",
    breaks=seq(-5,5,length.out=21), lab.breaks=round(seq(-5, 5, length.out=21),2))
plot(newworld, add=TRUE, lwd=1.5)
points(pp1[,1:2],pch=3,col='azure4')
points(pp2[,1:2],pch=8,col='azure4')
mtext("50th trends (zero-lag ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)

pv1=matrix(0,nrow=72,ncol=24)
pv2=matrix(0,nrow=72,ncol=24)
for (i in 1:72){
    for (j in 1:24){
        pv1[i,j]=ifelse(abs(op_slope[i,j])>=2*op_sigma[i,j] & abs(op_slope[i,j])<3*op_sigma[i,j],1,NA)
        pv2[i,j]=ifelse(abs(op_slope[i,j])>=3*op_sigma[i,j],1,NA)
    }
}
pp1=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv1)))
pp2=na.omit(data.frame(lon=rep(lon, 24), lat=rep(lat, each=72), sig=as.vector(pv2)))

image.plot(lon,lat,10*op_slope,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=difcol, 
    breaks=seq(-5,5,length.out=21), lab.breaks=round(seq(-5, 5, length.out=21),2))
plot(newworld, add=TRUE, lwd=1.5)
points(pp1[,1:2],pch=3,col='azure4')
points(pp2[,1:2],pch=8,col='azure4')
mtext("50th trends (peak ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)
dev.off()


png("omi-mls_snr_50th_2004-2021.png", height=8, width=7, units ="in", res =300)
par(mar=c(2, 2, 0.5, 0.5), mgp=c(2.4, 0.8, 0), las=1, mfrow=c(3,1),  oma=c(0.5,0.5,1.5,0.5))
pt=ops_slope/ops_sigma
pt[pt>=4.5]=4.5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=colorschemes$BluetoDarkOrange.18, xaxt="n",
    breaks=seq(-4.5,4.5,length.out=19), lab.breaks=round(seq(-4.5, 4.5, length.out=19),2))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend SNRs (wo/ ENSO correlation)",cex=1., line=0.5, at=0)

pt=op0_slope/op0_sigma
pt[pt>=4.5]=4.5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=colorschemes$BluetoDarkOrange.18, xaxt="n",
    breaks=seq(-4.5,4.5,length.out=19), lab.breaks=round(seq(-4.5, 4.5, length.out=19),2))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend SNRs (zero-lag ENSO correlation)",cex=1., line=0.5, at=0)

pt=op_slope/op_sigma
pt[pt>=4.5]=4.5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=colorschemes$BluetoDarkOrange.18, 
    breaks=seq(-4.5,4.5,length.out=19), lab.breaks=round(seq(-4.5, 4.5, length.out=19),2))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend SNRs (peak ENSO correlation)",cex=1., line=0.5, at=0)
dev.off()

  
png("omi-mls_sigma_50th_2004-2021.png", height=8, width=7, units ="in", res =300)
par(mar=c(2, 2, 0.5, 0.5), mgp=c(2.4, 0.8, 0), las=1, mfrow=c(3,1),  oma=c(0.5,0.5,1.5,0.5))
pt=ops_sigma
pt[pt>=0.1]=0.1
image.plot(lon,lat,pt*10,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=c(0,brewer.pal(9,"YlGn")), xaxt="n",
    breaks=seq(0,1,length.out=11), lab.breaks=c(seq(0, 0.9, length.out=10),'>1'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend uncertainties (wo/ ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)

pt=op0_sigma
pt[pt>=0.1]=0.1
image.plot(lon,lat,pt*10,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=c(0,brewer.pal(9,"YlGn")), xaxt="n",
    breaks=seq(0,1,length.out=11), lab.breaks=c(seq(0, 0.9, length.out=10),'>1'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend uncertainties (zero-lag ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)

pt=op_sigma
pt[pt>=0.1]=0.1
image.plot(lon,lat,pt*10,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, col=c(0,brewer.pal(9,"YlGn")), 
    breaks=seq(0,1,length.out=11), lab.breaks=c(seq(0, 0.9, length.out=10),'>1'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("50th trend uncertainties (peak ENSO correlation) [DU/decade]",cex=1., line=0.5, at=0)
dev.off()


png("omi-mls_rchange_ 50th_2004-2021.png", height=7, width=10, units ="in", res =300)
par(mar=c(2, 2, 0.5, 0.5), mgp=c(2.4, 0.8, 0), las=1, mfrow=c(2,2),  oma=c(0.5,0.5,1.5,0.5))
pt=(op0_slope-ops_slope)/ops_slope
pt[pt>=5]=5
pt[pt<=-5]=-5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, xaxt="n",
    col=difcol, breaks=seq(-5,5, length.out=21),  lab.breaks=c('<-5',seq(-4.5, 4.5, length.out=19),'>5'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Trend change [%] (zero-lag v no ENSO correlation)",cex=1., line=0.5, at=0)

pt=(op0_sigma-ops_sigma)/ops_sigma
pt[pt>=5]=5
pt[pt<=-5]=-5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, xaxt="n",
    col=difcol, breaks=seq(-1,1, length.out=21),  lab.breaks=c('<-1',seq(-0.9, 0.9, length.out=19),'>1'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Uncertainty change [%] (zero-lag v no ENSO correlation)",cex=1., line=0.5, at=0)

pt=(op_slope-op0_slope)/op0_slope
pt[pt>=5]=5
pt[pt<=-5]=-5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, 
    col=difcol, breaks=seq(-5,5, length.out=21),  lab.breaks=c('<-5',seq(-4.5, 4.5, length.out=19),'>5'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Trend change [%] (peak v zero-lag ENSO correlation)",cex=1., line=0.5, at=0)

pt=(op_sigma-op0_sigma)/op0_sigma
pt[pt>=5]=5
pt[pt<=-5]=-5
image.plot(lon,lat,pt,
    xlim=c(-180,180),ylim=c(-60,60), xlab="",ylab="", cex.lab=1, horizontal=F, 
    col=difcol, breaks=seq(-1,1, length.out=21),  lab.breaks=c('<-1',seq(-0.9, 0.9, length.out=19),'>1'))
plot(newworld, add=TRUE, lwd=1.5)
mtext("Uncertainty change [%] (peak v zero-lag ENSO correlation)",cex=1., line=0.5, at=0)
dev.off()

pt=(op0_slope-ops_slope)/ops_slope
pt[pt>=5]=5
pt[pt<=-5]=-5
freq=hist(pt,freq=T,breaks=seq(-5,5, length.out=21))$counts
sum(freq[10:11])/(72*24)
sum(freq[9:12])/(72*24)
#pt=(op_slope-op0_slope)/op0_slope
#coord=which(pt == min(pt), arr.ind=TRUE)
#dt=sat[coord[1],coord[2],]



