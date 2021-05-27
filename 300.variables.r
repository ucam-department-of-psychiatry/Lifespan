##
## Set up variables across 3xx-scripts
##

##
## Libraries (needed by functions in this file)
##
##

library("gamlss") ## need gamlss.family objects in all scripts

##
## Global random seed
##
set.seed( seed = 12345 ) ## Ensure reproducible code blocks
warning("Have set a fixed random seed, all runs will be identical unless this is changed (may be good or bad, depending on what you want)")


##
## Bootstrapping options
##
BOOT.OPT <- list(Seed=666,
                 Number.Cluster=6, ## higher if run on a cluster
                 Number.Replicates=500 ## should be 1000+ for production code
                 )

##
## whether to store certain output or not
##
STORE.OPT <- list(Fit.Summary=FALSE,    ## Not strictly necessary
                  Boot.Summary=FALSE,   ## Good for investigating failed bootstraps, but does increase file sizes
                  Boot.Index=FALSE,     ## Necessary for investigating failed bootstraps, but does increase file sizes
                  Expanded.Index=FALSE, ## Necessary for investigating expanded variability, but does increase file sizes
                  Fit.Full=FALSE,       ## These are very big files, recommend to leave as FALSE, only store the EXTRACT object
                  Fit.Overwrite=c("stop","overwrite","skip")[2]
                  )
## When fitting a model-object, we check if the matching fit-extract already exists. If it does:
## stop:= stop() with a fatal error
## overwrite:= recalcualte the fit and overwrite the existing fit-extract
## skip:= skip this model (give a warning)
##

##
## gamlss control variables
##
GAMLSS.CTRL <- gamlss.control(n.cyc=200) ## long tail for some iterations


##
## Must define FITTING.SET, NOVEL.SET, etc.
##
if( 0 ) {
    ##
    ## Automatically discover all subset directories (using the '-' separator)
    ## NOTE: This may be prone to error if hyphens appear in the names of directories/subsets/tag/etc
    ##

    FITTING.SET <- dir( path=file.path(RDS.DIR), pattern=".*-.*" )

    BOOT.SET <- FITTING.SET
    
    NOVEL.SET <- FITTING.SET

    DERIVED.SET <- FITTING.SET

    COPYMODEL.SET <- list() ## used by 321-script to copy a 'best' model to other subsets for fitting
    
    ## DO NOT USE THE ABOVE IF THE SIMULATION SCENARIOS ARE PRESENT, you will run multiple unnecessary calculations
    stop("Manually specify the xxx.SET objects to determine which scripts run on which data/subset folders!")
    
} else if ( 1 ) {
    ##
    ## 2021-04-01 settings
    ##

    ## using the testing subsets (for speed, change below to full)
    FITTING.SET <- c( "omega-Wild__.n0000", "omega-Wand__.n0000" )
    FITTING.SET <- c( "20210205full-GMVTransformed", "20210205full-sGMVTransformed", "20210205full-WMVTransformed", "20210205full-VentriclesTransformed" )
    
    NOVEL.SET <- FITTING.SET

    BOOT.SET <- FITTING.SET

    DERIVED.SET <- FITTING.SET
    
}
