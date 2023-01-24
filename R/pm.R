pmYear<-function(x){
  FLQuants(catch=FLQuantPoint(catch(x)),
           ssb  =FLQuantPoint(ssb(  x)),
           rec  =FLQuantPoint(rec(  x)),
           fbar =FLQuantPoint(fbar(  x)),
           cCum =FLQuantPoint(apply(catch(x),c(6), cumsum)),
           aav  =FLQuantPoint(as.FLQuant(adply(catch(x), c(1,3:6),  function(y) 
             data.frame(year=names(y)[-length(y)],
                        data=(y[-1]-y[-length(y)])/y[-length(y)])))),
           aav2 =FLQuantPoint(as.FLQuant(adply(catch(x), c(1,3:6),  function(y) 
             data.frame(year=names(y)[-length(y)],
                        data=(y[-1]-y[-length(y)])/y[-1])))))}

pmIter<-function(x){
  # Median total catch over whole time period, by iter
  p1=FLPar(catch.50=aaply(catch(x),6,median))
  
  # Median inter-annual variability over whole time period, by iter 
  p2=FLPar(iav.50=aaply(catch(x), 6,  function(y) var((y[-1]-y[-length(y)])/y[-length(y)])))
  
  # The median Inter-Annual Variability per iteration
  p3=FLPar(aav=array(aaply(catch(x), 6,  function(y) median(var((y[-1]-y[-length(y)])/y[-length(y)]))),c(1,dim(x)[6])))
  
  rbind(p1,p2,p3)}

stab<-function(x,u=0.2500,l=-0.2000){
  # The number of years when the stability mechanism was applied
  aav  =as.FLQuant(adply(x, c(1,3:6),  function(y) 
    data.frame(year=names(y)[-length(y)],
               data=(y[-1]-y[-length(y)])/y[-length(y)])))

  FLQuants("lower"=apply(as.FLQuant((round(aav,4)==l)),2,sum),
           "upper"=apply(as.FLQuant((round(aav,4)==u)),2,sum))}

stab<-function(x){
  # The number of years when the stability mechanism was applied
  aav  =as.FLQuant(adply(x, c(1,3:6),  function(y) 
    data.frame(year=names(y)[-length(y)],
               data=(y[-1]-y[-length(y)])/y[-length(y)])))
  aav}
  
pm<-function(x){
  
  list(
    # Median total catch over whole time period
    catch.50=aaply(catch(x),6,median),
    
    # Median inter-annual variability over whole time period 
    iav.50=aaply(catch(x), 6,  function(y) var((y[-1]-y[-length(y)])/y[-length(y)])),
    
    # Median stock size by year (and variability)
    ssb.50 =aaply(ssb(x),2,median),
    ssb.var=aaply(ssb(x),2,var),
    
    # Median recruitment by year (and variability)
    rec.50 =aaply(rec(x),2,median),
    rec.var=aaply(rec(x),2,var),
    
    # Median catch by year (and variability)
    catch.50 =aaply(catch(x),2,median),
    catch.var=aaply(catch(x),2,var),
    
    # The number of years when the stability mechanism was applied
    stab=NULL,
    
    # The median Inter-Annual Variability per iteration
    iav.50.iter=median(aaply(catch(x), 6,  function(y) var((y[-1]-y[-length(y)])/y[-length(y)])))
  )}