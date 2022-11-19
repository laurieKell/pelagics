
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