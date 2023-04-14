########################################################################
#moving block bootstrap for quantile regression
########################################################################
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
 samp_data=data[unlist(blocks[bn]), ]  
 mod=rq(formula, data=samp_data, tau=tau)
 coef(mod)
}

set.seed(2013)
fit=coef(rq(yd~x, data=mlo, tau=0.5))*12 #intercept and slope
op=t(replicate(1000, mbfun(formula=yd~x,data=mlo,tau=0.5)))
fit_se=apply(op, 2, sd, na.rm=TRUE)*12 #MBB standard error for intercept and slope
fit_pv=2*pt(q=abs(fit/fit_se), df=nrow(mlo)-2, lower.tail=FALSE) #MBB p value for intercept and slop
#multiple quantiles 
#op=t(replicate(1000, mbfun(formula=yd~x,data=mlo,tau=seq(0.1,0.9,by=0.1))[2,]))
