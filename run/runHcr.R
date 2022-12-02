library(FLasher)
library(FLBRP)
library(FLAssess)
library(mpb)
library(ggplotFL)
library(plyr)
library(magrittr)

#library(rdrop2)
#drop_auth()

setwd("/home/laurie/Desktop/projects/pelagics/R")

source('../R/hcrICES.R')

load("../data/inputs/ices.RData")

stks=names(ices)

stk="whb"

load(paste("../data/om/",stk,".RData",sep=""))
load(paste("../data/om/",stk,"Fwd.RData",sep=""))
load(paste("../data/om/",stk,"Err.RData",sep=""))

om =get(paste(stk,"Fwd",sep=""))
eql=get(paste(stk,"R",  sep=""))
err=get(paste(stk,"Err",sep=""))
refpts=FLPar(attributes(ices[[stks[grep(stk,stks)]]])$benchmark)

nits=dim(stock.n(om))[6]

#"whb.27.1-91214" "mac.27.nea"     "her.27.3a47d"  

mp=om
m_ply(data.frame(year=2021:2001), function(year) {
  res=FLBRP(window(mp,end=year),nyear=3)
  
  stock.wt(mp)[   ,ac(year)]=stock.wt(    res)
  catch.wt(mp)[   ,ac(year)]=catch.wt(    res)
  landings.wt(mp)[,ac(year)]=landings.wt( res)
  discards.wt(mp)[,ac(year)]=discards.wt( res)
  
  harvest(    mp)[,ac(year)]=catch.sel(   res)
  catch.n(    mp)[,ac(year)]=catch.sel(   res)
  landings.n( mp)[,ac(year)]=landings.sel(res)
  discards.n( mp)[,ac(year)]=discards.sel(res)
  
  m(          mp)[,ac(year)]=m(           res)
  mat(        mp)[,ac(year)]=mat(         res)
  
  mp<<-mp})
mp=propagate(iter(mp,seq(nits)),nits)

om=propagate(iter(om,seq(nits)),nits)

fit.new[[1]][[2]][names(fit.new[[1]][[1]])=="logR"][-43]
fit.new[[1]][[1]][names(fit.new[[1]][[1]])=="logR"][-43]

par0=array(c(ftar=refpts["Fmsy"],btrig=refpts["Btrigger"],fmin=refpts["Fmsy"],Blim=refpts["Blim"]),
           dim=c(4),
           dimnames=list(c("ftar","btrig","fmin","blim")))

par1=array(c(ftar=refpts["Fmsy"],btrig=refpts["Btrigger"],fmin=FLPar(0.05),blim=refpts["Blim"]),
           dim=c(4),
           dimnames=list(c("ftar","btrig","fmin","blim")))

start   =2001
end     =2021
interval=1

err =rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2021)),var(err[-1,1])^0.5)

srDev=rec(om)%=%1

#Reference scenarios with steepness=1
sims=list("ICES"=list(om,NULL))
sims[["fmsy"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2021))%=%0.32,
                           sr  =rec(om)),
                     NULL)
sims[["0.6"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2021))%=%0.6,
                          sr =rec(om)),
                    NULL)
sims[["sim0.0"]]=hcrICES(om,eql,rec(om)%=%1,
                      par0,
                      start=start,end=end,interval=interval)

# HCRs 
sims[["sim0"]]=hcrICES(om,eql,srDev,
                      par0,
                      start,end,interval,
                      err=err,
                      bndTac=c(0,Inf))
sims[["sim1"]]=hcrICES(om,eql,srDev,
                      par1,
                      start,end,interval,
                      err=err,
                      bndTac=c(0,Inf))

# HCRs with bounds 
sims[["sim0_bnd"]]=hcrICES(om,eql,srDev,
                      par0,
                      start,end,interval,  
                      err=err,
                      bndTac=c(0.8,1.25))
sims[["sim1_bnd"]]=hcrICES(om,eql,srDev,
                      par1,
                      start,end,interval,
                      err=err,
                      bndTac=c(0.80,1.25))

plot(subset(sims[["sim2_bnd"]][[2]],hcrYrs==2003)[,"ssb"],
     ssb(sims[["sim2_bnd"]][[1]])[,"2002"])

x=subset(sims[["sim2_bnd"]][[2]],hcrYrs==2003)[,"ssb"]
y=ssb(sims[["sim2_bnd"]][[1]])[,"2002"]


#######################################################
## Scenarios with random recruitment ##################
#######################################################
srDevAr=rlnoise(nits,iter(rec(om),1)%=%0,0.6,0.8)
sims[["acf1"]]=hcrICES(om,eql,srDevAr,
                      par1,
                      start,end,interval,
                      err=err,
                      bndTac=c(0,Inf))

### Beverton and Holt with steepness=1
### Assess Error = 0.5
### 2013-2021
## Scenarios with random recruitment
srMn=fmle(as.FLSR(om1,model="geomean"),
         control=list(silent=TRUE))

eqlMn=FLBRP(om1,sr=srMn) 

sims[["Box13-1"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par1,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2021)),0.5),
                          bndTac=c(0,Inf))
sims[["Box13-1-bnd"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par1,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2021)),0.5),
                          bndTac=c(0.8,1.25))

sims=llply(sims, function(x) x)

plot(FLStocks(llply(sims,function(x) x[[1]])))

save(sims,par0,par1,file=paste("../data/results/",stk,"Sims.RData",sep=""),compress="xz")


### Beverton and Holt with steepness estimated

#Reference scenarios
sim2=list("ICES"=list(om,NULL))
sim2[["fmsy"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2021))%=%0.32,
                        sr=eql,residuals=propagate(exp(residuals(sr)),1000)),
                    NULL)
sim2[["0.6"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2021))%=%0.6,
                       sr=eql,residuals=propagate(exp(residuals(sr)),1000)),
                   NULL)
sim2[["sim0.0"]]=hcrICES(om,eql,srDev,
                         par0,
                         start=start,end=end,interval=interval)

# HCRs 
sim2[["sim0"]]=hcrICES(om,eql,srDev,
                       par0,
                       start,end,interval,
                       err=err,
                       bndTac=c(0,Inf))
sim2[["sim1"]]=hcrICES(om,eql,srDev,
                       par1,
                       start,end,interval,
                       err=err,
                       bndTac=c(0,Inf))

# HCRs with bounds 
sim2[["sim0_bnd"]]=hcrICES(om,eql,srDev,
                           par0,
                           start,end,interval,  
                           err=err,
                           bndTac=c(0.8,1.25))
sim2[["sim1_bnd"]]=hcrICES(om,eql,srDev,
                           par1,
                           start,end,interval,
                           err=err,
                           bndTac=c(0.80,1.25))

## Scenarios with rsandom recruitment
srDev=rlnoise(nits,iter(rec(om),1)%=%0,0.6,0.8)
sim2[["acf1"]]=hcrICES(om,eql,srDevAr,
                       par1,
                       start,end,interval,
                       err=err,
                       bndTac=c(0,Inf))

sim2=llply(sim2, function(x) x)

plot(FLStocks(llply(sim2,function(x) x[[1]])))

save(sim2,par0,par1,file=paste("../data/results/",stk,"Sim.RData",sep=""),compress="xz")
