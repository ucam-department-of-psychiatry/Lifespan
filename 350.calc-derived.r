rm(list=ls())
## COMMON LIBRARIES AND FUNCTIONS
source("100.common-variables.r")
source("101.common-functions.r")

source("300.variables.r")
source("301.functions.r")

## SCRIPT SPECIFIC LIBRARIES

## SCRIPT SPECIFIC FUNCTIONS

## SCRIPT CODE
##
##
if( 1 ){
    Print.Disclaimer( )

    for( lset in DERIVED.SET ) {

        cat( "=====", lset, "=====\n" )

        PATHS.LIST <- Create.Folders( Tag=lset )

        PRIMARY <- Load.Subset.Wrapper( Tag=lset, LSubset=TRUE, LModel=TRUE, LFit=TRUE, LBoot=TRUE, LData=TRUE )

        ##
        ## Using FIT.EXTRACT (fitted to SUBSET), make predictions on DATA
        ##
        ## NOTE: Currently we do not include bootstrap variability in DATA.PRED, seems too much information (although it is important of course)
        PRIMARY$DATA.PRED <- Apply.Param(NEWData=PRIMARY$DATA,
                                         Reference.Holder=PRIMARY,
                                         FITParam=PRIMARY$FIT.EXTRACT$param,
                                         Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                         Prefix="",
                                         Add.Moments=FALSE, ## does not make sense for observations
                                         Add.Normalise=TRUE,
                                         Add.Derivative=FALSE,  ## does not make sense for observations
                                         MissingToZero=TRUE, NAToZero=TRUE,
                                         verbose=FALSE )

        ## Since SUBSET.PRED is entirely contained within DATA.PRED we do not explicitly calculate it again
        ## However, we do not have a convenint back mapping from SUBSET to DATA, perhaps we should have saved a logical oer outcome inclusion vector.
        

        ##
        ## Using predictions on DATA, summarise individuals with longitudinal observations (not included in SUBSET, and hence not within the fitting process)
        ##
        PRIMARY$LONG.SUMMARY <- Make.Longitudinal( Holder=PRIMARY )




        
        ##
        ## Calculate Population Curves
        ##
        PRIMARY$POP.CURVE.LIST <- list()
        attr(PRIMARY$POP.CURVE.LIST,"xlim") <- range( PRIMARY$SUBSET[,PRIMARY$MODEL$covariates$X], na.rm=TRUE )
        PRIMARY$POP.CURVE.LIST[[PRIMARY$MODEL$covariates$X]] <- seq(attr(PRIMARY$POP.CURVE.LIST,"xlim")[1],
                                                                    attr(PRIMARY$POP.CURVE.LIST,"xlim")[2],
                                                                    length.out= 2^12 )
        if( !is.null(PRIMARY$MODEL$covariates$BY) ) {
            for( lBY in PRIMARY$MODEL$covariates$BY ) {
                PRIMARY$POP.CURVE.LIST[[ lBY ]] <- factor( x=1:nlevels(PRIMARY$SUBSET[[ lBY ]]), labels=levels(PRIMARY$SUBSET[[ lBY ]]) )
            }
        }
        ## NOTE: we generate population curves over Cov$X stratified by Cov$BY, all other model covariates will be zeroed [makes sense for contr.sum factors only]
        PRIMARY$POP.CURVE.RAW <- do.call( what=expand.grid, args=PRIMARY$POP.CURVE.LIST )

        ##
        ## Using FIT.EXTRACT and BOOT.EXTRACT, obtain population curves (across the range of DATA [which may extrapolate outside the range of SUBSET])
        ##
        PRIMARY$POP.CURVE.PRED <- Apply.FitAndBoot(NEWDATA=PRIMARY$POP.CURVE.RAW, FIT=PRIMARY$FIT.EXTRACT, BOOT=PRIMARY$BOOT.EXTRACT,
                                                   Boot.wre=TRUE, VERBOSE=TRUE,
                                                   Prefix="",
                                                   Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                                   Add.Moments=TRUE,
                                                   Add.Normalise=FALSE, ## no outcome, so does not make sense
                                                   Add.Derivative=TRUE,
                                                   MissingToZero=TRUE, NAToZero=TRUE
                                                   )                                               


        ##
        ## Calculate Study Curves
        ##
        PRIMARY$STUDY.CURVES <- list()
        for( lSTUDY in levels(PRIMARY$SUBSET[,PRIMARY$MODEL$covariates$RANEF[1]]) ) {
            cat(lset,":",PRIMARY$MODEL$covariates$RANEF[1],"curve for",lSTUDY,"\n")
            lWHICH <- which( PRIMARY$SUBSET[ , PRIMARY$MODEL$covariates$RANEF ] == lSTUDY )
            if( length(lWHICH) < 2 ) {
                cat("Study",lSTUDY,"does not have enough observations within the SUBSET to define a study-specific curve.\n")
                next
            }
            PRIMARY$STUDY.CURVES[[lSTUDY]] <- list()
            attr(PRIMARY$STUDY.CURVES[[lSTUDY]],"xlim") <- range( PRIMARY$SUBSET[ lWHICH, PRIMARY$MODEL$covariates$X ], na.rm=TRUE )
            PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST <- list()
            PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[PRIMARY$MODEL$covariates$X]] <- seq(attr(PRIMARY$STUDY.CURVES[[lSTUDY]],"xlim")[1],
                                                                                     attr(PRIMARY$STUDY.CURVES[[lSTUDY]],"xlim")[2],
                                                                                     length.out= 2^9 )
            PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[PRIMARY$MODEL$covariates$RANEF]] <- lSTUDY

            if( !is.null(PRIMARY$MODEL$covariates$BY) ) {
                for( lBY in PRIMARY$MODEL$covariates$BY ) {
                    PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[ lBY ]] <- factor( x=1:nlevels(PRIMARY$SUBSET[[ lBY ]]), labels=levels(PRIMARY$SUBSET[[ lBY ]]) )
                }
            }
            if( !is.null(PRIMARY$MODEL$covariates$OTHER) ) {
                for( lOTHER in PRIMARY$MODEL$covariates$OTHER ) {

                    if( any(class( PRIMARY$SUBSET[[ lOTHER ]] )=="factor") ) {
                        lTAB <- table(PRIMARY$SUBSET[ lWHICH, lOTHER ])
                        PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[ lOTHER ]] <- factor( names(which.max(lTAB))[1], levels=levels(PRIMARY$SUBSET[ , lOTHER ]) )
                        cat("Setting",lOTHER,"to",PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[ lOTHER ]],"\n")
                    } else if( any(class( PRIMARY$SUBSET[[ lOTHER ]] )=="numeric") ) {
                        PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[ lOTHER ]] <- median( PRIMARY$SUBSET[ lWHICH, lOTHER ], na.rm=TRUE )
                        cat("Setting",lOTHER,"to",PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST[[ lOTHER ]],"\n")
                    } else {
                        cat("Do nothing with",lOTHER,"covariate\n")
                    }
                }
            }
            PRIMARY$STUDY.CURVES[[lSTUDY]]$RAW <- do.call( what=expand.grid, args=PRIMARY$STUDY.CURVES[[lSTUDY]]$LIST )

            ##
            ## Using FIT.EXTRACT and BOOT.EXTRACT, obtain study curves (across the range of study within DATA)
            ##
            PRIMARY$STUDY.CURVES[[lSTUDY]]$PRED <- Apply.FitAndBoot(NEWDATA=PRIMARY$STUDY.CURVES[[lSTUDY]]$RAW,
                                                                    FIT=PRIMARY$FIT.EXTRACT, BOOT=PRIMARY$BOOT.EXTRACT,
                                                                    Boot.wre=TRUE,
                                                                    Prefix="",
                                                                    Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                                                    Add.Moments=TRUE,
                                                                    Add.Normalise=FALSE, ## no outcome, so does not make sense
                                                                    Add.Derivative=TRUE,
                                                                    MissingToZero=TRUE, NAToZero=TRUE,
                                                                    VERBOSE=FALSE )                                               

        }
        

        
        ##
        ## Save HOLDER+DERIVED
        ##
        saveRDS(object=PRIMARY, file=file.path(PATHS.LIST$PATH,"DERIVED.rds"))
        
    }
}

print( warnings() ) ## print any warnings from the code
