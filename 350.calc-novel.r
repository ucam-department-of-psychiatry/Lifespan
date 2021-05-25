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
   
    for( lset in NOVEL.SET ) {
        cat( "=====", lset, "=====\n" )

        PATHS.LIST <- Create.Folders( lset )

        PRIMARY <- Load.Subset.Wrapper( Tag=lset, LSubset=TRUE, LModel=TRUE, LFit=TRUE, LBoot=TRUE )

        ##
        ## Load the novel datasets
        ##

        FILES <- list.files( PATHS.LIST$NOVEL, pattern="*.rds" )

        for( lFILE in FILES ) {
            cat( lset, lFILE, "\n" )
            
            NOVEL <- Load.Subset.Wrapper( Tag=lset, LModel=TRUE, LFit=TRUE )

            NOVEL$NOVEL.TAG <- sub("\\.rds","",lFILE)

            NOVEL$DATA <- readRDS( file.path(PATHS.LIST$NOVEL,sprintf("%s.rds",NOVEL$NOVEL.TAG)) )

            NOVEL$DATA.PRED <- Apply.Param(NEWData=NOVEL$DATA,
                                           FITParam=PRIMARY$FIT.EXTRACT$param,
                                           Reference.Holder=PRIMARY,
                                           Pred.Set=NULL, Prefix="", Add.Moments=FALSE, Add.Normalise=FALSE, Add.Derivative=FALSE, MissingToZero=TRUE,
                                           verbose=FALSE )

            NOVEL$SUBSET <- NOVEL$DATA.PRED[ attr(NOVEL$DATA.PRED,"logical.selectors")$REFIT.VALID, ]
            ##
            ##

            EXPANDED <- Calc.Expanded(NewData=NOVEL$SUBSET,
                                      Cur.Param=PRIMARY$FIT.EXTRACT$param,
                                      Missing=attr(NOVEL$DATA.PRED,"missing.levels") )
            cat( sprintf("Refitting using %i novel observations\n",NROW(NOVEL$SUBSET)) )            

            EXPANDED.PATH <- file.path( PATHS.LIST$NOVEL, NOVEL$NOVEL.TAG )

            if( !dir.exists(EXPANDED.PATH) ) {
                dir.create(EXPANDED.PATH)
            }

            saveRDS(object=list(param=EXPANDED,summary=NULL),
                    file=file.path(EXPANDED.PATH,"FIT.EXPANDED.rds"))


            NOVEL$DATA.PRED <- Apply.Param(NEWData=NOVEL$DATA,
                                         Reference.Holder=PRIMARY,
                                         FITParam=EXPANDED,
                                         Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                         Prefix="",
                                         Add.Moments=FALSE, ## does not make sense for observations
                                         Add.Normalise=TRUE,
                                         Add.Derivative=FALSE,  ## does not make sense for observations
                                         MissingToZero=TRUE, NAToZero=TRUE,
                                         verbose=FALSE )
            
            saveRDS(object=NOVEL$DATA.PRED,
                    file=file.path(EXPANDED.PATH,"DATA.PRED.rds"))
            



            
            ##
            ## Do the same for all the bootstrap replicates
            ##
            BOOT.EXPANDED <- list()
            for( IDX in 1:length(PRIMARY$BOOT.EXTRACT) ) {
                
                ##if( (IDX %% 25) == 0 ) { cat("     ",IDX,"\n") }
                
                if( is.list(PRIMARY$BOOT.EXTRACT[[IDX]]$param) ) {
                    SEED <- PRIMARY$BOOT.EXTRACT[[IDX]]$base + PRIMARY$BOOT.EXTRACT[[IDX]]$offset+0
                    set.seed(seed=SEED)

                    if(is.null(PRIMARY$MODEL$stratify)) {
                        INDEX <- sample(1:NROW(NOVEL$SUBSET),NROW(NOVEL$SUBSET),TRUE)
                    } else {
                        INDEX <- unsplit(lapply(split(1:NROW(NOVEL$SUBSET),NOVEL$SUBSET[PRIMARY$MODEL$stratify]),
                                                function(X){sample(x=X,size=length(X),replace=TRUE)}),
                                         NOVEL$SUBSET[PRIMARY$MODEL$stratify])
                    }

                    EXPANDED <- Calc.Expanded(NewData=NOVEL$SUBSET[INDEX,],
                                              Cur.Param=PRIMARY$BOOT.EXTRACT[[IDX]]$param,
                                              Missing=attr(NOVEL$SUBSET,"missing.levels") )

                    BOOT.EXPANDED[[IDX]] <- list(param=EXPANDED,summary=NULL,seed=SEED)
                    if( STORE.OPT$Expanded.Index ) {
                        BOOT.EXPANDED[[IDX]]$index <- INDEX
                    } else {
                        BOOT.EXPANDED[[IDX]]$index <- NULL
                    }
                }
            }
            saveRDS(object=BOOT.EXPANDED,file=file.path(EXPANDED.PATH,"BOOT.EXPANDED.rds"))
        }
    }
}

print( warnings() ) ## print any warnings from the code

