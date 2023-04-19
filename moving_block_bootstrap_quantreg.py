#This code is maintained by Niklas Selke
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.formula.api as smf

#read data
mlo = pd.read_csv("mlo.csv")
mlo.head()

#derive the seasonal cycle
base = mlo[(mlo["year"] >= 2000) & (mlo["year"] <= 2020)]
seasonality = smf.ols("y~np.sin(2*np.pi*month/12)+np.cos(2*np.pi*month/12)+np.sin(2*np.pi*month/6)+np.cos(2*np.pi*month/6)", base).fit(method="qr").predict(pd.DataFrame({"month": range(1, 13)}))

#calculate anomalies
mlo = pd.merge(mlo, pd.DataFrame({"month": range(1, 13), "yd": seasonality}), on="month")
mlo["yd"] = mlo["y"] - mlo["yd"]
mlo = mlo.sort_values("x")


##moving block bootstrap function
def mbfun(formula, data, tau):
    n = len(data)
    b = np.ceil(n**0.25).astype(int)
    nblocks = np.ceil(n/b).astype(int)
    blocks = np.array([list(range(i, i+b)) for i in range(n-b+1)])
    bn = np.random.choice(np.arange(len(blocks)), nblocks, replace=True)
    samp_data = data.iloc[blocks[bn].flatten()]
    mod = smf.quantreg(formula, data=samp_data).fit(q=tau)
    return mod.params


np.random.seed(2013)
#median intercept and slope
fit = smf.quantreg("yd~x", data=mlo).fit(q=0.5).params*12

#MBB standard error for intercept and slope
op = [mbfun(formula="yd~x", data=mlo, tau=0.5) for _ in range(1000)]
fit_se = np.nanstd(op, axis=0)*12

#p value for intercept and slope
fit_pv = 2*scipy.stats.t.sf(x=abs(fit/fit_se), df=len(mlo)-2)
