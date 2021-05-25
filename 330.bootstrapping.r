rm(list=ls())
## COMMON LIBRARIES AND FUNCTIONS
source("100.common-variables.r")
source("101.common-functions.r")

source("300.variables.r")
source("301.functions.r")

## SCRIPT SPECIFIC LIBRARIES
##
##
library("doSNOW")
library("foreach") ## loaded by doSNOW

## SCRIPT SPECIFIC FUNCTIONS

## SCRIPT CODE
##
##
if( 1 ) {
    Print.Disclaimer( )

    ##
    ## Local variables
    ##
    
    for( lset in BOOT.SET ) {
        cat( "=====", lset, "=====\n" )

        PATHS.LIST <- Create.Folders( Tag=lset )
        
        HOLDER <- Load.Subset.Wrapper( Tag=lset, LSubset=TRUE, LModel=TRUE, LFit=TRUE )

        CLUSTER <- makeCluster(BOOT.OPT$Number.Cluster)
        registerDoSNOW(CLUSTER)

        PROGRESS <- txtProgressBar(max = BOOT.OPT$Number.Replicates, style = 3)
        OPTS <- list(progress=function(n) setTxtProgressBar(PROGRESS, n))

        FOREACH.OBJ <- foreach(n=1:BOOT.OPT$Number.Replicates, .options.snow=OPTS, .packages=c("gamlss")) ## .combine=rbind,

        BOOT.EXTRACT <- FOREACH.OBJ %dopar% Boot.Function(n=n,Base.Seed=BOOT.OPT$Seed,Holder=HOLDER)

        saveRDS(object=BOOT.EXTRACT,file=file.path(PATHS.LIST$BOOT.EXTRACT,sprintf("s%05i+n%05i.rds",BOOT.OPT$Seed,BOOT.OPT$Number.Replicates)))
        
        close(PROGRESS)
        stopCluster(CLUSTER)

        

    }


    
}


print( warnings() ) ## print any warnings from the code
