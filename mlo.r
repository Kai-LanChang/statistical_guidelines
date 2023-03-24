
########################################################################
#read data
########################################################################
hd=read.csv("/Users/kai-lanchang/Dropbox/data/All_Surface_Ozone_1957-2022_DRAFT.csv")
ymd= data.frame(do.call('rbind', strsplit(as.character(hd[,1]),'/',fixed=TRUE)))
ymdt=data.frame(do.call('rbind', strsplit(as.character(ymd[,3]),' ',fixed=TRUE)))
hd=data.frame(Year=ymdt[,1], Month=ymd[,1], Day=ymd[,2], Time=ymdt[,2], hd[,-1])
hd$Year= as.numeric(hd$Year)
hd$Month= as.numeric(hd$Month)
hd$Day= as.numeric(hd$Day)

dd=aggregate(x=hd, by=list(Day=hd$Day, Month=hd$Month,Year=hd$Year), FUN=function(x) {mean(x,na.rm=T)})
dd=dd[,-c(1,2,3)]
for (i in 1:nrow(dd)){
td=hd[hd$Year==dd[i,1] & hd$Month==dd[i,2] & hd$Day==dd[i,3],]
	for (k in 5:20){
	dd[i,k]=ifelse(sum(!is.na(td[,k])) > 11, mean(td[,k], na.rm = T), NA_real_)
}}

mlo=hd[,c("Year","Month","Day","Time","MLO")]
mlo1=mlo[mlo$Time=="8:00",]
mlo2=mlo[mlo$Time=="9:00",]
mlo3=mlo[mlo$Time=="10:00",]
mlo4=mlo[mlo$Time=="11:00",]
mlo5=mlo[mlo$Time=="12:00",]
mlo6=mlo[mlo$Time=="13:00",]
mlo7=mlo[mlo$Time=="14:00",]
mlo8=mlo[mlo$Time=="15:00",]
mlo=rbind(mlo1,mlo2,mlo3,mlo4,mlo5,mlo6,mlo7,mlo8)
for (i in 1:nrow(dd)){
td=mlo[mlo$Year==dd[i,1] & mlo$Month==dd[i,2] & mlo$Day==dd[i,3],]
	dd[i,14]=ifelse(sum(!is.na(td[,5])) > 3, mean(td[,5], na.rm = T), NA_real_)
}

md=aggregate(x=dd, by=list(Month=dd$Month,Year=dd$Year), FUN=function(x) {mean(x,na.rm=T)})
md=md[,-c(1,2)]
for (i in 1:nrow(md)){
td=dd[dd$Year==md[i,1] & dd$Month==md[i,2],]
	for (k in 5:20){
	md[i,k]=ifelse(sum(!is.na(td[,k])) > 14, mean(td[,k], na.rm = T), NA_real_)
}}

ny=66
date=cbind(Year=rep(seq(1957,2022), each=12), Month=rep(seq(1,12),ny), ind=seq(1,12*ny))
md=merge(date, md, by=c("Year","Month"), all=T)
md=md[,c(1:3,13)]
md=na.omit(md)
colnames(md)=c("year","month","x","y")
md=md[md$year>=1974 & md$year<2022,]
write.csv(md, "mlo.csv", row.names=F)

########################################################################
#read mlo data
########################################################################
dir="/Users/kai-lanchang/Dropbox/2.CurrentProjects/2023-SotC_BAMS/"
library(quantreg)

mlo=read.csv("mlo.csv")
head(mlo)

base=mlo[mlo$year>=2000 & mlo$year<=2020,]
seasonality=predict(lm(y~sin(2*pi*month/12)+cos(2*pi*month/12)+sin(2*pi*month/6)+cos(2*pi*month/6), data=base), newdata=data.frame(month=1:12))
mlo=merge(mlo, data.frame(month=1:12, yd=seasonality), by="month")
mlo$yd=mlo$y-mlo$yd
mlo=mlo[order(mlo$x),]

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
fit=coef(rq(yd~x, data=mlo, tau=0.5))*12 #intercept and slope
op = t(replicate(1000, mbfun(formula=yd~x,data=mlo,tau=0.5)))
fit_se=t(apply(op, 2, sd, na.rm=TRUE))*12 #MBB standard error for intercept and slope
fit_pv=2*pt(q=abs(fit/fit_se), df=nrow(mlo)-2, lower.tail=FALSE) #MBB p value for intercept and slope

op=t(replicate(1000, mbfun(formula=yd~x,data=mlo,tau=seq(0.1,0.9,by=0.1))[2,]))
