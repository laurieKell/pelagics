require(FLasher)
require(FLBRP)
require(FLAssess)
require(mpb)
require(ggplotFL)
require(plyr)
require(magrittr)

#require(rdrop2)
#drop_auth()

dropboxdir="~/Dropbox/pelagics"

#setwd("/home/laurie/Desktop/projects/pelagics/R")

load(file.path(dropboxdir,paste("data/runs/",stk,".RData",sep="")))

source("../R/hcrICESV2.R")

stkid=c("whb.27.1-91214","mac.27.nea","her.27.3a47d")

#stk     ="whb"

nits    =100
start   =2001
end     =2021
interval=1

load(paste(file.path(dropboxdir,paste("data/om/",stk,"OM.RData",sep=""))))

ices=get(paste(stk,"Ctc", sep=""))
fmsy=get(paste(stk,"Fmsy",sep=""))
err =get(paste(stk,"Err", sep=""))
ftar=get(paste(stk,"Ftar",sep=""))
refs=get(paste(stk,"Refs",sep=""))
eql =FLBRP(ices)

par=as(model.frame(mcf(FLQuants(ftar=ftar,
                                btrig=refs[["msybtrigger"]][,ac(start:end)],
                                fmin=ftar%=%0.05,
                                blim=refs[["blim"]][,ac(start:end)])),drop=T)[,-1],"FLPar")
par["btrig"][is.na(par["btrig"])]=min(par["btrig"],na.rm=T)
par["blim"][ is.na(par["blim"])] =min(par["blim"],na.rm=T)

err     =rlnorm(nits,FLQuant(0,dimnames=list(year=start:end)),var(err[-1,1])^0.5)

implErr=FLQuants(NULL)
implErr[["10%"]]=propagate(FLQuant(0.10,dimnames=list(year=(start-1):end)),nits)
implErr[["20%"]]=propagate(FLQuant(0.20,dimnames=list(year=(start-1):end)),nits)
implErr[["30%"]]=propagate(FLQuant(0.30,dimnames=list(year=(start-1):end)),nits)

implErr[["10% final"]]=propagate(FLQuant(c(rep(0,12),rep(0.10,10)),dimnames=list(year=(start-1):end)),nits)
implErr[["20% final"]]=propagate(FLQuant(c(rep(0,12),rep(0.20,10)),dimnames=list(year=(start-1):end)),nits)
implErr[["30% final"]]=propagate(FLQuant(c(rep(0,12),rep(0.30,10)),dimnames=list(year=(start-1):end)),nits)

implErr[["random final"]]=implErr[[4]]
implErr[["random final"]][,ac(2012:2021)][] = sample(c(0.10,0.20,0.30),nits*10,T)
# implErr[["random final"]][,ac(2012:2021)][] = sample(runif(nits, 0.1, 0.3),nits*10,T)

if (!FALSE){
#Scenario
sims=list("ICES"=list(ices,NULL))
sims[["Fmsy"]]=list(fmsy,NULL)

## HCRs ########################################################################
sims[["HCR1"]]     =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0,Inf))
sims[["HCR2"]]     =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0,Inf),
                            bndWhen ="blim")
sims[["HCR1 bnd"]] =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25))
sims[["HCR1 +10%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0,Inf),
                            implErr=implErr[["10%"]])
sims[["HCR1 +20%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0,Inf),
                            implErr=implErr[["20%"]])
sims[["HCR1 +30%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0,Inf),
                            implErr=implErr[["30%"]])
sims[["HCR1 +10% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            implErr=implErr[["10%"]])
sims[["HCR1 +20% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            implErr=implErr[["20%"]])
sims[["HCR1 +30% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            implErr=implErr[["30%"]])

sims[["HCR2 bnd"]] =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim")
sims[["HCR2 +10%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["10%"]])
sims[["HCR2 +20%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["20%"]])
sims[["HCR2 +30%"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["30%"]])
sims[["HCR2 +10% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["10%"]])
sims[["HCR2 +20% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["20%"]])
sims[["HCR2 +30% bnd"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                            bndTac=c(0.8,1.25),
                            bndWhen ="blim",
                            implErr=implErr[["30%"]])
sims[["HCR1 bnd cap"]] =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                bndTac=c(0.8,1.25),
                                bndCap=0.4)
sims[["HCR2 bnd cap"]] =hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                bndTac=c(0.8,1.25),
                                bndWhen ="blim",
                                bndCap=0.4)

sims[["HCR1 +10% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["10%"]],
                                    bndCap=0.4)
sims[["HCR1 +20% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["20%"]],
                                    bndCap=0.4)
sims[["HCR1 +30% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["30%"]],
                                    bndCap=0.4)

sims[["HCR2 +10% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["10%"]],
                                    bndWhen ="blim",
                                    bndCap=0.4)
sims[["HCR2 +20% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["20%"]],
                                    bndWhen ="blim",
                                    bndCap=0.4)
sims[["HCR2 +30% bnd cap"]]=hcrICES(fmsy,eql,rec(fmsy),par,start,end,interval,lag=lag,err=err,
                                    bndTac=c(0.8,1.25),
                                    implErr=implErr[["30%"]],
                                    bndWhen ="blim",
                                    bndCap=0.4)
}
#plot(FLStocks(llply(sims,function(x) x[[1]])))

save(sims,refs,file=paste(file.path(dropboxdir,paste("data/runs/",stk,".RData",sep=""))))