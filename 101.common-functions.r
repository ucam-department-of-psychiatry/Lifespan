##
## Load modified gamlss functions
##
## NOTE: We have written alternative GG() family - to avoid computation issues
##       We have written alternative bfp() function - to avoid NA issue
##
source("102.gamlss-recode.r")

##
## Disclaimer and version
##
##
Print.Disclaimer <- function( ) {

    cat("
##### Disclaimer
##
## This code is provided as is. Users are responsible for checking and
## validating the output.
##
## This code make multiple assumptions about the format of the
## dataset, the type of model being fitted, and the ultimate important
## information to be extracted from the fitting process. 
##
##
##### The code is associated with the publication:
##
## Bethlehem, Seidlitz, White et al. (2022). Nature. https://www.nature.com/articles/s41586-022-04554-y
##
##### Code written by:
##
## Dr Simon R. White ( sw539@cam.ac.uk ), University of Cambridge, UK
## Dr Richard A.I. Bethlehem ( rb643@cam.ac.uk ), University of Cambridge, UK
##
##### Version
##
## v11.r1 (2022-05-13)
##
")
    
    invisible(TRUE)
}

##
## Create folder structure and single source of truth for path names
##
Create.Folders <- function(Tag, Anchor=RDS.DIR) {
    OBJ <- list(Tag=Tag,
                Anchor=Anchor,
                PATH=file.path( Anchor, Tag )
                )

    if( !dir.exists(OBJ$PATH) ) {
        dir.create(OBJ$PATH)
    }
    for( Folder in c("MODEL","FIT.FULL","FIT.EXTRACT","BOOT.EXTRACT","NOVEL","LOO.EXTRACT") ) {
        OBJ[[Folder]] <- file.path(OBJ$PATH,Folder)
        if( !dir.exists( OBJ[[Folder]] ) ) {
            dir.create( OBJ[[Folder]] )
        }
    }
    return( OBJ )
}

##
## Check attributes on cleaned DATA/SUBSET objects
##
Check.Attributes <- function( Data ) {

    if( is.null( attr(Data,"columns") ) ) stop("Missing columns attribute") ## used within setup-scripts

    if( is.null( attr(Data,"tag") ) ) stop("Missing tag attribute") ## used to refer to DATA

    if( is.null( attr(Data,"Transformations") ) ) stop("Missing Transformations attribute") ## used to convert Cov$X/etc into original form


    if( !all(c("INDEX.ID","INDEX.OB","INDEX.TYPE") %in% names(Data)) ) stop("Missing pipeline INDEX.xx columns")
}
