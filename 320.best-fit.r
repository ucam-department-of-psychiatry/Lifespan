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
if( 1 ) { ## this block fits all models (and associated bfp-versions)
    Print.Disclaimer( )
    
    for( lset in FITTING.SET ) {
        cat( "=====", lset, "=====\n" )

        PATHS.LIST <- Create.Folders( lset )
        
        FILES <- list.files( PATHS.LIST$FIT.EXTRACT, pattern="*.rds")
        
        EXTRACT <- list()
        for( lFILE in FILES ) {
            EXTRACT[[lFILE]] <- readRDS(file.path(PATHS.LIST$FIT.EXTRACT,lFILE))
        }

        Logical.Completed <- sapply(EXTRACT,function(X){ if(is.null(X$diagnostics)){TRUE}else{FALSE} })
        Logical.NotFP <- sapply(EXTRACT,function(X){ if(is.null(X$param)){NA}else{!attr(X$param,"model")$inc.fp} })

        if( 1 ) {
            RUNS <- data.frame(row.names=FILES)
            RUNS$Completed <- sapply(EXTRACT,function(X){ if(is.null(X$diagnostics)){TRUE}else{FALSE} })
            RUNS$BIC <- sapply(EXTRACT,function(X){ if(is.null(X$param)){NA}else{X$param$BIC} })
            RUNS$bfp.version <- grepl("bfpNA",row.names(RUNS))
            tmpA <- grep("bfpNA",row.names(RUNS),value=TRUE)
            RUNS[tmpA,"bfp.version.of"] <- row.names(RUNS)[match(sub("bfpNA","fp",tmpA),row.names(RUNS))]
            RUNS$relative.BIC <- RUNS$BIC - min(RUNS$BIC,na.rm=TRUE)

            table(RUNS$bfp.version)

            print(RUNS[with(RUNS,order(Completed,bfp.version,-1*relative.BIC)),])

            ##cat(unlist(attr(EXTRACT[["baseFO330R110.GGalt.bfpNA.rds"]]$param,"model")[c("mu","sigma","nu")]),sep="\n")

            saveRDS(object=RUNS[with(RUNS,order(Completed,bfp.version,-1*relative.BIC)),],
                    file=file.path(PATHS.LIST$FIT.EXTRACT,"run-summary.rds"))
            
        }
        
        WHICH <- which(Logical.Completed & Logical.NotFP)
        cat("The following models failed to retrn... check?\n")
        cat(FILES[which((!Logical.Completed))],"\n",sep="\n")
        
        BICs <- sapply( EXTRACT[WHICH], function(X){X$param$BIC} )
        BEST.FILE <- names(BICs)[which.min(BICs)]

        sink(file=file.path(PATHS.LIST$FIT.EXTRACT,"BIC-output.txt"), split=TRUE)
        cat( lset, "best (by BIC) is", BEST.FILE, "\n")
        print( BICs - min(BICs) )
        sink(file=NULL)
        
        

        for( LAB in c("FIT.FULL","FIT.EXTRACT","MODEL") ) {

            if( file.exists(file.path(PATHS.LIST$PATH,sprintf("%s.rds",LAB))) ) {
                STR <- sprintf("[%s : %s] Removing existing 'best' model, is this right?",lset,LAB)
                cat(STR,"\n")
                warning(STR)
                file.remove(file.path(PATHS.LIST$PATH,sprintf("%s.rds",LAB)))
            }
            
            if( file.exists(file.path( PATHS.LIST[[ LAB]], BEST.FILE )) ) {
                if( LINK.TYPE=="symlink" ) {
                    file.symlink(from=file.path(LAB,BEST.FILE),
                                 to=file.path(PATHS.LIST$PATH,sprintf("%s.rds",LAB))
                                 )
                } else if ( LINK.TYPE=="copy" ) {
                    file.copy(from=file.path(PATHS.LIST$FIT.EXTRACT,BEST.FILE),
                              to=file.path(PATHS.LIST$PATH,sprintf("%s.rds",LAB))
                              )
                } else {
                    stop("Not implemented hardlinks")                    
                }
            }
        }
    }
}


print( warnings() ) ## print any warnings from the code
