library(FLasher)
library(FLBRP)
library(FLAssess)
library(mpb)
library(ggplotFL)
library(plyr)
library(magrittr)

source('/home/laurie/Desktop/projects/hindcast/R/hcrICES.R')

dirDat="/home/laurie/Desktop/projects/hindcast/data/om"
load(file.path(dirDat,"mac.RData"))
load(file.path(dirDat,"herFwd.RData"))
om =herFwd
eql=herR

nits=dim(stock.n(om))[6]

refpts=FLPar(c(
  blim    =1500000,
  bpa     =2250000,
  flim    =0.88,
  fpa     =0.53,
  msytrig =2250000,
  fmsy    =0.32,
  bmgtlow =1500000,
  bmgt    =2250000,
  fmgstlow=0.05,
  fmgt    =0.32),units="NA")

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
om=propagate(iter(om,seq(nits)),nits)
mp=propagate(iter(mp,seq(nits)),nits)

fit.new[[1]][[2]][names(fit.new[[1]][[1]])=="logR"][-43]
fit.new[[1]][[1]][names(fit.new[[1]][[1]])=="logR"][-43]



par0=array(c(ftar=refpts["fmsy"],btrig=refpts["msytrig"],fmin=refpts["fmsy"],blim=refpts["blim"]),
           dim=c(4),
           dimnames=list(c("ftar","btrig","fmin","blim")))

par1=array(c(ftar=refpts["fmsy"],btrig=refpts["msytrig"],fmin=FLPar(0.05),blim=refpts["blim"]),
           dim=c(4),
           dimnames=list(c("ftar","btrig","fmin","blim")))

par2=array(c(0.20,2250000,
             0.05,1500000,
             0.32,5200000,
             0.20,4000000),
          dim=c(4,2),
          dimnames=list(c("ftar","btrig","fmin","blim"),c("lower","upper")))[c(1,3,2,4),1:2]

start   =2001
end     =2018
interval=1

err =rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2018)),0.3)

srDev=rec(om)%=%1

#Reference scenarios with steepness=1
sims=list("ICES"=list(om,NULL))
sims[["fmsy"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2018))%=%0.32,
                           sr  =rec(om)),
                     NULL)
sims[["0.6"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2018))%=%0.6,
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
sims[["sim2"]]=hcrICES(om,eql,srDev,
                      par2,
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
sims[["sim2_bnd"]]=hcrICES(om,eql,srDev,
                      par2,
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
sims[["acf2"]]=hcrICES(om,eql,srDevAr,
                      par2,
                      start,end,interval,
                      err=err,
                      bndTac=c(0,Inf))

### Beverton and Holt with steepness=1
### Assess Error = 0.5
### 2013-2018
## Scenarios with random recruitment
srMn=fmle(as.FLSR(om1,model="geomean"),
         control=list(silent=TRUE))

eqlMn=FLBRP(om1,sr=srMn) 

sims[["Box13-1"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par1,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2018)),0.5),
                          bndTac=c(0,Inf))
sims[["Box13-2"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par2,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2018)),0.5),
                          bndTac=c(0,Inf))
sims[["Box13-1-bnd"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par1,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2018)),0.5),
                          bndTac=c(0.8,1.25))
sims[["Box13-2-bnd"]]=hcrICES(om,eqlMn,srDev=exp(residuals(srMn)),
                          par2,
                          start=2013,end,interval,
                          err=rlnorm(nits,FLQuant(0,dimnames=list(year=2001:2018)),0.5),
                          bndTac=c(0.8,1.25))

sims=llply(sims, function(x) x)

plot(FLStocks(llply(sims,function(x) x[[1]])))

save(sims,par0,par1,par2,file="/home/laurence/Desktop/Dropbox/PELAC/bluewhiting/data/sims.RData",compress="xz")



### Beverton and Holt with steepness estimated

#Reference scenarios
sim2=list("ICES"=list(om,NULL))
sim2[["fmsy"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2018))%=%0.32,
                        sr=eql,residuals=propagate(exp(residuals(sr)),1000)),
                    NULL)
sim2[["0.6"]]=list(fwd(om,fbar=FLQuant(1,dimnames=list(year=2001:2018))%=%0.6,
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
sim2[["sim2"]]=hcrICES(om,eql,srDev,
                       par2,
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
sim2[["sim2_bnd"]]=hcrICES(om,eql,srDev,
                           par2,
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
sim2[["acf2"]]=hcrICES(om,eql,srDevAr,
                       par2,
                       start,end,interval,
                       err=err,
                       bndTac=c(0,Inf))

sim2=llply(sim2, function(x) x)

plot(FLStocks(llply(sim2,function(x) x[[1]])))

save(sim2,par0,par1,par2,file="/home/laurence/Desktop/Dropbox/PELAC/bluewhiting/data/sim2.RData",compress="xz")
