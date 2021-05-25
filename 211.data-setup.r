rm(list=ls())
## COMMON LIBRARIES AND FUNCTIONS
source("100.common-variables.r")
source("101.common-functions.r")

source("200.variables.r")

## SCRIPT SPECIFIC LIBRARIES

## SCRIPT SPECIFIC FUNCTIONS

## SCRIPT CODE
##
##
for( TESTING in c(1) ) {
    Print.Disclaimer( )

    ##
    ## Load and clean real data
    ## (note, most cleaning done by RB in other prior to running this script)
    ##

    
    RAW.DATA <- read.csv(file.path( RAW.DIR,"allnew_050221.csv"),stringsAsFactors = TRUE) # load data

    if( TESTING==1 ) {
        DATA.TAG <- "20210205full"
    } else {
        ## This is only for testing+development, a smaller version of the real dataset used to validate the code
        DATA.TAG <- "20210205part"
        RAW.DATA <- RAW.DATA[ sample(1:NROW(RAW.DATA),floor(NROW(RAW.DATA)/5)) ,]
    }

    ##
    ## TOPLEVEL FILTERING
    ##
    RAW.DATA <- RAW.DATA[ RAW.DATA$study!="iADNI", ]
    cat("Dropping iADNI from DATA\n")
    
    ##
    ## Perform any transformations necessary for fitting
    ## 
    ## NOTE: for fp()/bfp() functions, Cov$X must be POSITIVE (apply a shift) and SMALL (apply a scaling factor)
    ##
    TRANSFORMATIONS <- list()
    TRANSFORMATIONS[[ "X" ]] <- list("OriginalName"="age_days",
                                     "TransformedName"="AgeTransformed",
                                     "toTransformed"=function(Z) { log(Z) }, ## must manually scale X-variable for numerical stability within bfpNA()
                                     "toOriginal"=function(Z) { exp(Z) }
                                     )
    RAW.DATA[,TRANSFORMATIONS[["X"]][["TransformedName"]]] <- TRANSFORMATIONS[["X"]][["toTransformed"]]( RAW.DATA[, TRANSFORMATIONS[["X"]][["OriginalName"]] ] )

    TRANSFORMATIONS[[ "Y" ]] <- list("OriginalName"="OUTCOME", ## We will update these in a moment
                                     "TransformedName"="OUTCOMETransformed",
                                     "toTransformed"=function(Z) { Z/10000 }, ## must manually scale X-variable for numerical stability within bfpNA()
                                     "toOriginal"=function(Z) { Z*10000 }
                                     )

    for( OUTCOME in c("GMV","WMV","sGMV","Ventricles") ) {
        RAW.DATA[,sprintf("%sTransformed",OUTCOME)] <- TRANSFORMATIONS[["Y"]][["toTransformed"]]( RAW.DATA[, OUTCOME ] )
    }
        
    ##
    ## Clean columns and factors
    ## NOTE: Use Additional-element to retain columns not used in the modelling (model columns are Outcomes+Covariates)
    ##
    COLUMNS <- list(Outcomes=c("GMVTransformed","WMVTransformed","sGMVTransformed","VentriclesTransformed"),
                    Covariates=c("AgeTransformed","study","fs_version","sex","dx","participant"),
                    Additional=c("session","run","country","site","age_days")
                    )
    if( 1 ) {
        Sufficiently.Complete.Columns <- names( which( sapply( RAW.DATA, function(X){ sum(!is.na(X))/length(X) } ) > 0.75 ) )
        cat("The following columns are 75% complete within RAW.DATA, but not saved within DATA or SUBSET - check:\n")
        print(Sufficiently.Complete.Columns[ !Sufficiently.Complete.Columns %in% unique(c(unlist(COLUMNS),sub("Transformed","",COLUMNS$Outcomes))) ])

        cat("We note a columns called 'age', which is almost equal to 'age_days', but contains NAs. For those NAs what is the mean of 'age_days'?\n")
        AGG <- aggregate( age_days ~ study, data=droplevels(RAW.DATA[is.na(RAW.DATA$age),]),FUN=function(X){c(mean=mean(X),n=length(X))}, simplify=TRUE )
        AGG$age_years <- AGG$age_days[,"mean"]/365.25
        print(AGG)
        cat("Are we happy with the above list of studies with missing age (but available age_days)?\n")
    }

    
    DATA <- RAW.DATA[,unlist(COLUMNS)] ## only keep selected columns
    rm(RAW.DATA)

    DATA$sex.original <- DATA$sex
    DATA$sex <- factor( as.character(DATA$sex.original), levels=c("Female","Male") )
    warning("Assumptions regarding sex coding")
    
    DATA$dx <- relevel(DATA$dx,"CN")


    ##
    ## Manually drop specific unused levels
    ##
    DATA <- droplevels(DATA) # make sure there are no left over levels


    ##
    ## Add INDEX.TYPE column, a script-defined column to separate individuals used to fit the model (ie healthy controls)
    ##
    DATA$INDEX.TYPE   <- factor(DATA$dx=="CN",levels=c(TRUE,FALSE),labels=c("CN","notCN"))
    
    ##
    ## Create unique identifier (INDEX.ID) for each person across all studies
    ## NOTE: This column will be used in later scripts, so it must exist!
    ##
    DATA$INDEX.ID <- factor( paste( DATA$study, DATA$participant, DATA$sex, sep="|" ) )
    warning("Created a bespoke INDEX.ID which \"should\" uniquely identify each individual")


    ##
    ## Check on missing data
    if( 1 ) {
        cat("Proportion of complete values by DATA columns:\n")
        sapply( DATA, function(X){ sum(!is.na(X))/length(X) } )
    }
    
    ##
    ## Reorder dataset
    ## NOTE: This assumes first scan (by age) corresponds to first scan of interest. This may not be true!
    ##
    DATA <- DATA[order(DATA$INDEX.ID,DATA$AgeTransformed),]
    
    DATA$INDICATOR.OB <- (DATA$run==1)
    COLUMNS$Additional <- append( COLUMNS$Additional, "INDICATOR.OB" )
    warning("Within real datasets using ( run==1 ) to select first run within each session")

    ##
    ## Create identifer for 'first' scan (INDEX.OB) within each person
    ## NOTE: This column will be used to subset to cross-sectional data, so it must exist!
    ##
    DATA$INDEX.OB <- NA
    DATA$INDEX.OB[ which(DATA$INDICATOR.OB) ] <- Reduce(c,lapply(rle(as.numeric(DATA$INDEX.ID[which(DATA$INDICATOR.OB)]))$lengths,function(k){1:k}))
    warning("Created a bespoke INDEX.OB which \"should\" identify repeat observations of individuals (akin to \"run\" variable)")



    ##
    ## Need to add these new columns to the 'to keep' list
    ##
    COLUMNS$Index <- c("INDEX.ID","INDEX.OB","INDEX.TYPE")
    


    ##
    ## Adding Outcome specific exclusion criteria
    ##
    DATA[,"sGMVTransformed.DROP"] <- ifelse(DATA$fs_version=="FSInfant", TRUE, FALSE )

    if( 0 ) {
        FTAB <- ftable(xtabs( ~ study + INDEX.TYPE, data=DATA[with(DATA,(INDEX.OB==1)&(!is.na(WMVTransformed))),] ))
        ORDER <- order(FTAB[,1])
        structure(.Data=FTAB[ORDER,],dim=dim(FTAB),class="ftable",row.vars=list(study=attr(FTAB,"row.vars")$study[ORDER]),col.vars=list(INDEX.TYPE=attr(FTAB,"col.vars")$INDEX.TYPE))
    }

    
    COLUMNS$Drop <- c("sGMVTransformed.DROP")
    


    ##
    ## Attaching some attributes
    ##
    attr(DATA,"columns") <- COLUMNS    
    
    attr(DATA,"tag") <- DATA.TAG
    
    attr(DATA,"Transformations") <- TRANSFORMATIONS
    

    ##
    ## Sanity checking dataset
    ##
    if( 1 ) {
        print( xtabs( ~ addNA(session) + addNA(INDEX.OB), data=DATA ) )
        warning("Must ensure session+run Vs INDEX.OB mis-matches are valid (see commented out code to investigate mis-matches)")

        cat("\n\n")
    }


    ##
    ## Save dataset in RDS format for use in later scripts
    ##
    DATA.PATH <- file.path( RDS.DIR, DATA.TAG )
    if( !dir.exists( DATA.PATH ) ) {
        dir.create( DATA.PATH )
    }

    TOSAVE <- DATA[ , unlist(COLUMNS) ]

    attributes(TOSAVE) <- c( attributes(TOSAVE), attributes(DATA)[c("columns","tag","Transformations")] )
    
    Check.Attributes(TOSAVE)

    saveRDS(object=TOSAVE,
            file=file.path( DATA.PATH, "DATA.rds"))

    

    ##
    ## ==============================
    ## Below this point we generate SUBSETs, MODELs, NOVEL datasets derived from DATA
    ##

    

    for( OUTCOME in COLUMNS$Outcomes ) {


        PATHS.LIST <- Create.Folders( Tag=sprintf("%s-%s", DATA.TAG, OUTCOME ) )
        
        ##
        ## Generate subsets (by outcome[=column] and included/excluded[=rows])
        ## NOTE: excluded implicitly means cross-sectional, ie. only 'first' observation
        ##

        WHICH <- list()
        
        WHICH$BASELINE.CONTROL <- with(DATA, (INDEX.TYPE==levels(INDEX.TYPE)[1]) & (INDEX.OB==1))
        warning("Current SUBSET is based on INDEX.ID and INDEX.OB assumptions, these might not be the correct way to subset the data")

        ## Following is outcome-specific code
        if(!is.null(COLUMNS$Drop)){
            MATCH <- match(x=sprintf("%s.DROP",OUTCOME), table=COLUMNS$Drop )
            if( !is.na(MATCH) ) {
                WHICH$KEEP <- !DATA[,COLUMNS$Drop] ## above we specify in terms of which rows to drop, so we must negate to keep those we want to KEEP
                cat("Outcome specific subsetting:",OUTCOME," (dropping ",sum(!WHICH$KEEP,na.rm=TRUE)," rows)\n")
            } else {
                WHICH$KEEP <- rep(TRUE,NROW(DATA))
            }
        } else {
            WHICH$KEEP <- rep(TRUE,NROW(DATA))
        }

        

        ##
        ## Check for NAs in Outcome, Covariates, Index and Drop columns (not Additional, since they do not impact fitting by definition)
        WHICH$VALID <- Reduce(`&`, lapply( DATA[c(OUTCOME,unlist( attr(DATA,"columns")[c("Covariates", "Index", "Drop")] ) )], function(X){!is.na(X)} ) )

        WHICH.COLUMNS <- c( OUTCOME, unlist( attr(DATA,"columns")[c("Covariates", "Additional", "Index", "Drop")] ) ) ## note explicitly including OUTCOME

        
        
        SUBSET <- droplevels( DATA[ Reduce( `&`, WHICH ), WHICH.COLUMNS ] )

        cat( "Subset", PATHS.LIST$PATH, "has", NROW(SUBSET), "rows.\n")

        attributes(SUBSET) <- c( attributes(SUBSET), attributes(DATA)[c("columns", "tag","Transformations")] )

        attr(SUBSET,"DATA.WHICH.LIST") <- WHICH

        ## Set the per-SUBSET trasnformation names for the Y-variable
        ## NOTE: we could imagine doing the outcome transformations within this for-loop
        ##       but since it is common to all outcomes, might as well do it above
        ##       (hence the need for thse lines below to 'fix' the Y-variable name)
        attr(SUBSET,"Transformations")[["Y"]]$OriginalName <- sub("Transformed","",OUTCOME)
        attr(SUBSET,"Transformations")[["Y"]]$TransformedName <- OUTCOME
        
        saveRDS(object=SUBSET, file=file.path(PATHS.LIST$PATH,"SUBSET.rds"))


        ##
        ## Generate model sets
        ##
        cat("Generating models...\n")

        
        ## NOTE: FAMILY.SET allows us to explore multiple gamlss outcome distributions, later scripts will select the 'best' (via AIC/BIC/etc)
        FAMILY.SET <- c("GGalt")
        FP.SET <- matrix(c(1,1,0,
                           2,1,0,
                           2,0,1,
                           2,1,1,
                           2,2,1,
                           2,2,2,
                           3,1,0,
                           3,1,1,
                           3,2,1,
                           3,1,2,
                           3,2,2,
                           3,3,0,
                           3,3,1,
                           3,3,2),
                         byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))

        RANDOM.SET <- matrix(c(0,0,0,
                               1,0,0,
                               0,1,0,
                               1,1,0
                               ),
                             byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))
        row.names(RANDOM.SET) <- LETTERS[1:NROW(RANDOM.SET)]
        RANDOM.STR <- c(""," + random(study)")
        
        for( lFAM in FAMILY.SET ) { ## loop to search multiple outcome distributions

            for( iFP in 1:NROW(FP.SET) ) {

                for( iRAND in 1:NROW(RANDOM.SET) ) {

                    MODEL.NAME <- paste0("baseFO",paste0(FP.SET[iFP,],collapse=""),"R",paste0(RANDOM.SET[iRAND,],collapse=""))
                    
                    MODEL <- list(covariates=list("Y"=OUTCOME, ## The Outcome
                                                  "X"="AgeTransformed", ## The main X-variable (continuous) for plotting against Y
                                                  "ID"="participant", ## Subject-level ID, will be superceded by INDEX.ID in later scripts
                                                  "BY"="sex", ## factor columns to stratify plots (and implicitly within the model)
                                                  "OTHER"="fs_version", ## other variables (note: in later scripts if missing these will be set to zero)
                                                  ## ie.  Y ~ f( X, BY, OTHER )
                                                  "COND"="dx", ## should be all equal to base case in fitted SUBSET
                                                  "RANEF"="study"),
                                  family=lFAM,
                                  contrasts=list("fs_version"="contr.sum"), ## (*1*) Make intercept interpretation as mean fs_version
                                  stratify=c("study","sex"),
                                  mu   =if(FP.SET[iFP,"mu"]>0){
                                            sprintf("%s ~ 1 + sex + fs_version + fp(AgeTransformed,npoly = %i)%s",
                                                    OUTCOME,
                                                    FP.SET[iFP,"mu"],
                                                    RANDOM.STR[RANDOM.SET[iRAND,"mu"]+1])
                                        } else {
                                            sprintf("%s ~ 1 + sex + fs_version%s",
                                                    OUTCOME,
                                                    RANDOM.STR[RANDOM.SET[iRAND,"mu"]+1])
                                        },
                                  sigma=if(FP.SET[iFP,"sigma"]>0){
                                            sprintf("%s ~ 1 + sex + fp(AgeTransformed,npoly = %i)%s",
                                                    OUTCOME,
                                                    FP.SET[iFP,"sigma"],
                                                    RANDOM.STR[RANDOM.SET[iRAND,"sigma"]+1])
                                        } else {
                                            sprintf("%s ~ 1 + sex%s",
                                                    OUTCOME,
                                                    RANDOM.STR[RANDOM.SET[iRAND,"sigma"]+1])
                                        },
                                  nu   =if(FP.SET[iFP,"nu"]>0){
                                            sprintf("%s ~ 1 + fp(AgeTransformed,npoly = %i)%s",
                                                    OUTCOME,
                                                    FP.SET[iFP,"nu"],
                                                    RANDOM.STR[RANDOM.SET[iRAND,"nu"]+1])
                                        } else {
                                            sprintf("%s ~ 1%s",
                                                    OUTCOME,
                                                    RANDOM.STR[RANDOM.SET[iRAND,"nu"]+1])
                                        },
                                  inc.fp=TRUE)

                    saveRDS(object=MODEL,file=file.path(PATHS.LIST$MODEL,sprintf("%s.%s.fp.rds",MODEL.NAME,lFAM)))
                }
            }
        }


        ##
        ## Generate novel-testing-datasets
        ##
        cat("Generating novel datasets...\n")

        CLONE.SET <- c("ADNI","NSPN")
        for( TO.CLONE in CLONE.SET ) {

            NOVEL.NAME <- sprintf("%s-CLONE",TO.CLONE)
            nWHICH <- which(SUBSET$study==TO.CLONE)
            if( length(nWHICH)>50 ) {
                NOVEL <- SUBSET[ nWHICH, ] ## we're going to use the HCP study as a basis
                NOVEL$study <- factor(NOVEL.NAME)       ## MadeUp should end up the same as HCP... by definition

                if(NROW(NOVEL)>0) {
                    ## Want a sequence of novel datasets moving toward complete clone
                    HALF <- if(NROW(NOVEL)>400){floor(NROW(NOVEL)/2)}else{Inf}
                    SIZE.SEQ <- unique(pmin(NROW(NOVEL),c(5,10,15,20,25,30,35,40,45,50,100,150,200,HALF,Inf)))

                    for( lSIZE in SIZE.SEQ ) {

                        SAMPLE <- sort(sample(1:NROW(NOVEL),lSIZE,replace=FALSE))

                        saveRDS(object=NOVEL[ SAMPLE, ],
                                file=file.path(PATHS.LIST$NOVEL,sprintf("%s.n%04i.rds",NOVEL.NAME,lSIZE)))
                    }
                }
            }
        }
    }
}

print( warnings() ) ## print any warnings from the code
