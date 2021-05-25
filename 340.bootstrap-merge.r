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
if( 1 ) {
    Print.Disclaimer( )

    
    for( lset in BOOT.SET ) {
        cat( "=====", lset, "=====\n" )


        PATHS.LIST <- Create.Folders( Tag=lset )

        if( file.exists( file.path( PATHS.LIST$PATH, "BOOT.EXTRACT.rds" ) ) ) {
            cat(lset,"Deleting existing merged file.\n")
            file.remove( file.path( PATHS.LIST$PATH, "BOOT.EXTRACT.rds" ) )
        }

        FILES <- list.files(path=PATHS.LIST$BOOT.EXTRACT, pattern="*.rds",full.names=TRUE )

        BOOTS <- lapply(FILES,readRDS)

        MERGED <- Reduce( c, BOOTS )

        print( table(sapply( MERGED, function(X){ !is.null(X$param) } )) )
        
        saveRDS(object=MERGED,file=file.path(PATHS.LIST$PATH,"BOOT.EXTRACT.rds"))
    }
}


print( warnings() ) ## print any warnings from the code
