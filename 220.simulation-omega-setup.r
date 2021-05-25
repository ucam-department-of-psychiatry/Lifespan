rm(list=ls())
## COMMON LIBRARIES AND FUNCTIONS
source("100.common-variables.r")
source("101.common-functions.r")

source("200.variables.r")
source("201.functions.r")

## SCRIPT SPECIFIC LIBRARIES

## SCRIPT SPECIFIC FUNCTIONS

## SCRIPT CODE
##
##
if( 1 ) {
    Print.Disclaimer( )

    ##
    ## Set random seed for consistent simulations (over-ride any global setting)
    ##
    set.seed(seed=222)

    ##
    ## This is the omega simulation scenario
    ##
    DATA.TAG <- "omega"
    

    ##
    ## Define simulated dataset structure
    ##
    DATA.spec <- data.frame(Study=LETTERS[1:22],
                            N=c(rep(500,10),rep(1000,5),rep(500,3),rep(1000,2),200,200),
                            novel=c(rep(FALSE,20),rep(TRUE,2)), ## designate two studies for the novel-script
                            r=c(rep(1,15),rep(5,5),1,5),
                            t0=c(seq(0.0,0.8,length=10),    seq(0,0.8,length=5),    seq(0.2,0.8,length=3),    seq(0.2,0.8,length=2),    0.05,0.05),
                            t1=c(seq(0.0,0.8,length=10)+0.1,seq(0,0.8,length=5)+0.1,seq(0.2,0.8,length=3)+0.1,seq(0.2,0.8,length=2)+0.1,0.15,0.15),
                            row.names=LETTERS[1:22]
                            )
    ##
    ## Generate study random-effects
    ##
    RE.SD.STUDY <- 0.075
    if( 1 ) {
        ## Define a finite set of normal steps, in terms of standard normal
        RE.SET <- seq(-2.5,2.5,by=0.25)

        DATA.spec$MU.ranef.Study <- sample(x=RE.SET * RE.SD.STUDY, ## MUST RESCALE STD-NORMAL TO RE.SD
                                           size=NROW(DATA.spec),replace=TRUE,prob=dnorm(x=RE.SET))

        DATA.spec["U","MU.ranef.Study"] <- -0.25 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD
        DATA.spec["V","MU.ranef.Study"] <-  3.00 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD

    } else {
        ## OR fully random, but less control on simulation effects
        DATA.spec$MU.ranef.Study <- rnorm(n=NROW(DATA.spec),mean=0,sd=RE.SD.STUDY)
    }
    ##
    ## Build subject-level data structure (DATA0)
    ##
    RE.SD.ID <- 0.1
    DATA.parts <- list()
    for( iTYPE in 1:2 ) {
        lCOUNT <- if( iTYPE==1 ) {DATA.spec$N} else {floor(DATA.spec$N/2)} ## implies a 1:2 ratio for CN:non-CN
        DATA.parts[[iTYPE]]          <- data.frame(Study=factor(rep(DATA.spec$Study,times=lCOUNT)))
        DATA.parts[[iTYPE]][,"Grp"]  <- factor(sample(1:2,NROW(DATA.parts[[iTYPE]]),TRUE),1:2,c("F","M"))
        DATA.parts[[iTYPE]][,"Type"] <- factor( rep( iTYPE, NROW(DATA.parts[[iTYPE]]) ), 1:2, c("CN","notCN"))

        DATA.parts[[iTYPE]][,"seq"]  <- Reduce(f=c,sapply( rle(as.numeric(DATA.parts[[iTYPE]]$Study))$lengths, function(X){seq(from=1,to=X)} ))
        DATA.parts[[iTYPE]][,"ID"]   <- factor(sprintf("ID%i%04i",iTYPE,DATA.parts[[iTYPE]][,"seq"]))

        DATA.parts[[iTYPE]][,"t0.rand"] <- runif(n=NROW(DATA.parts[[iTYPE]]))
        DATA.parts[[iTYPE]][,"t0"]      <- DATA.spec$t0[(DATA.parts[[iTYPE]]$Study)]
        DATA.parts[[iTYPE]][,"t1"]      <- DATA.spec$t1[(DATA.parts[[iTYPE]]$Study)]
        DATA.parts[[iTYPE]][,"r"]       <- DATA.spec$r[(DATA.parts[[iTYPE]]$Study)]
        DATA.parts[[iTYPE]][,"time0"]      <- with( DATA.parts[[iTYPE]], ((t1 - t0)*t0.rand) + t0 )

        DATA.parts[[iTYPE]][,"MU.ranef.ID"] <- rnorm(n=NROW(DATA.parts[[iTYPE]]),mean=0,sd=RE.SD.ID)
        DATA.parts[[iTYPE]][,"MU.ranef.Study"] <- DATA.spec$MU.ranef.Study[DATA.parts[[iTYPE]]$Study]


        DATA.parts[[iTYPE]][,"SIZE"] <- ifelse((iTYPE==1) & (DATA.parts[[iTYPE]][,"Study"]%in%c("U","V")),
                                               DATA.parts[[iTYPE]][,"seq"],
                                               0)
    }
    DATA.A <- Reduce(rbind,DATA.parts)

    attr(DATA.spec,"re.sd")  <- list(Study=RE.SD.STUDY, ID=RE.SD.ID)
    
    ##
    ## Build dataset
    ##
    DATA <- DATA.A[rep(1:NROW(DATA.A),times=DATA.A$r),]

    DATA[,"obs"] <- Reduce( f=c, lapply( rle(as.numeric(DATA$ID))$lengths, function(X){1:X} ) ) - 1

    DATA[,"time"] <- DATA[,"time0"] + ( 0.04*DATA[,"obs"])

    DATA[,"Time"] <- ( 80 * DATA[,"time"] )

    TRANSFORMATIONS <- list()
    TRANSFORMATIONS[[ "X" ]] <- list("OriginalName"="Time",
                                     "TransformedName"="TimeTransformed",
                                     "toTransformed"=function(X) { X/10 }, ## must manually scale X-variable for numerical stability within bfpNA()
                                     "toOriginal"=function(X) { 10 * X }
                                     )
    DATA[,TRANSFORMATIONS[["X"]][["TransformedName"]]] <- TRANSFORMATIONS[["X"]][["toTransformed"]]( DATA[, TRANSFORMATIONS[["X"]][["OriginalName"]] ] )

    ##
    ## Generate outcome (including random-effects)
    ##
    if( 1 ) {
        ## Generate Wand
        ##
        Func.mu.fixed.1 <- function( X ) {
            X <- X/8 ## time -> Time -> TimeTransformed transformation
            OUT <- log(-1*(0.4-X)*(0.5-X)+1.8)
            return(OUT)
        }
        Func.mu.fixed.2 <- function( X ) {
            X <- X/8 ## time -> Time -> TimeTransformed transformation
            OUT <- log(-1*(0.35-X)*(0.3-X)+1.55)
            return(OUT)
        }

        Func.mu.fixed.3 <- function( X ) {
            X <- X/8 ## time -> Time -> TimeTransformed transformation
            OUT <- 1 + (0.5*X)
            return(OUT)
        }
        Func.mu.fixed.4 <- function( X ) {
            X <- X/8 ## time -> Time -> TimeTransformed transformation
            OUT <- 0.75 + (0.75*X)
            return(OUT)
        }

        
        if( 0 ) {
            SEQ <- seq(0,1,length.out=256)
            plot( x=SEQ, y=exp(Func.mu.fixed.1(SEQ)), type="l", col="black", ylim=c(0,2) )
            lines( x=SEQ, y=exp(Func.mu.fixed.2(SEQ)), col="red" )

            lines( x=SEQ, y=Func.mu.fixed.3(SEQ), col="black", lty=2 )
            lines( x=SEQ, y=Func.mu.fixed.4(SEQ), col="red", lty=2 )
        }



        attr(DATA.spec,"truth") <- list("Wand"=list(family="GGalt",
                                                    MU=list(TypeBase=Func.mu.fixed.1,TypeOther=Func.mu.fixed.2),
                                                    SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                                                    NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
                                                    ),
                                        "Wild"=list(family="NO",
                                                    MU=list(TypeBase=Func.mu.fixed.3,TypeOther=Func.mu.fixed.4),
                                                    SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                                                    NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
                                                    )
                                        )

        DATA[,"RAND"] <- runif(NROW(DATA),min=0,max=1) ## common random number
        

        
        TRUTH.COLUMNS <- c("RAND")
        
        for( LAB in names(attr(DATA.spec,"truth")) ) {

            FAMILY <- get(attr(DATA.spec,"truth")[[LAB]]$family)

            ARGS.FULL <- list(p=DATA[,"RAND"])
            
            for( lP in names(FAMILY()$parameters) ) {
                
                if( toupper(lP) %in% names(attr(DATA.spec,"truth")[[LAB]]) ) {
                    DATA[,sprintf("%s.%s.fixef",LAB,toupper(lP))] <- ifelse(DATA$Type=="CN",
                                                                            attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeBase(DATA[,"TimeTransformed"]),
                                                                            attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeOther(DATA[,"TimeTransformed"]) )
                    TRUTH.COLUMNS <- append(TRUTH.COLUMNS,sprintf("%s.%s.fixef",LAB,toupper(lP)))

                    ASPECTS <- c(sprintf("%s.%s.fixef",LAB,toupper(lP)), sprintf("%s.ranef.ID",toupper(lP)), sprintf("%s.ranef.Study",toupper(lP)) )

                    DATA[,sprintf("%s.%s",LAB,toupper(lP))] <- rowSums( DATA[ , ASPECTS[ ASPECTS %in% names(DATA) ], drop=FALSE ] )

                    TRUTH.COLUMNS <- append( TRUTH.COLUMNS, sprintf("%s.%s",LAB,toupper(lP)) )
                    
                    ARGS.FULL[[lP]] <- FAMILY()[[sprintf("%s.linkinv",lP)]]( DATA[,sprintf("%s.%s",LAB,toupper(lP))] )

                } else {
                    stop("Must specify true-functional forms of all gamlss-components!")
                }
            }
            DATA[,LAB] <- do.call( what=get(sprintf("q%s",attr(DATA.spec,"truth")[[LAB]]$family)), args=ARGS.FULL )

        }
    }

    ##
    ## Add special columns (INDEX.ID and INDEX.OB) used by later scripts
    ##
    DATA$INDEX.ID <- factor( paste(DATA$Study,DATA$ID,sep="_") )
    DATA$INDEX.OB <- as.integer(DATA$obs+1)
    DATA$INDEX.TYPE <- DATA$Type


    ##
    ## Select only a few columns
    ##
    COLUMNS <- list(Outcomes=names(attr(DATA.spec,"truth")),
                    Covariates=c("Study","Grp","Type","ID","TimeTransformed"),
                    Additional=c("t0.rand", "t0", "t1", "r", "time0","time","Time",
                                 TRUTH.COLUMNS),
                    Index=c("INDEX.ID","INDEX.OB","INDEX.TYPE"), ## Would add 'Type' here, to be consistent with real-dataset, but Type is within Covariates
                    Drop=NULL
                    )


    

    ##
    ## Add data spec as an attribute
    ##
    attr(DATA,"spec") <- DATA.spec

    attr(DATA,"columns") <- COLUMNS

    attr(DATA,"tag") <- DATA.TAG

    attr(DATA,"Transformations") <- TRANSFORMATIONS

    ##
    ## For simulations, we will generate clone and entirely new novel datasets
    ##
    COMMON.SET <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
                    "M", "N", "O", "P", "Q", "R", "S", "T")
    SIZE.SEQ <- c(5,10,15,20,25,30,35,40,45,50,100,150,200)
    attr(DATA,"SCENARIOS") <- list(list(label="__",base=COMMON.SET,added=NULL,sizes=c(0)),
                                   list(label="u_",base=COMMON.SET,added=c("U"),sizes=SIZE.SEQ),
                                   list(label="_v",base=COMMON.SET,added=c("V"),sizes=SIZE.SEQ),
                                   list(label="uv",base=COMMON.SET,added=c("U","V"),sizes=SIZE.SEQ)
                                   )
    

    ##
    ## Save dataset in RDS format for use in later scripts
    ##
    
    DATA.PATH <- file.path( RDS.DIR, DATA.TAG )
    if( !dir.exists( DATA.PATH ) ) {
        dir.create( DATA.PATH )
    }
    TOSAVE <- DATA[ , unlist(attr(DATA,"columns")) ]

    attributes(TOSAVE) <- c( attributes(TOSAVE), attributes(DATA)[c("columns","tag","Transformations","spec","SCENARIOS")] )
    
    Check.Attributes(TOSAVE)
    
    saveRDS(object=TOSAVE, file=file.path( DATA.PATH, "DATA.rds"))
    


    ##
    ## ==============================
    ## Below this point we generate SUBSETs, MODELs, NOVEL datasets derived from DATA
    ##




    
    for( OUTCOME in attr(DATA,"columns")$Outcomes ) {

        ##
        ## Within the simulated datasets, we will also explore the novel-expanded estimation
        ##
        for( iSCENARIO in 1:length(attr(DATA,"SCENARIOS")) ) {

            lSCENARIO <- attr(DATA,"SCENARIOS")[[iSCENARIO]]

            for( iSIZE in 1:length(lSCENARIO$sizes) ) {

                lSIZE <- lSCENARIO$sizes[iSIZE]

                lTAG <- sprintf("[%s / %s / %i]", OUTCOME, lSCENARIO$label, lSIZE)
                
                cat( lTAG, "\n" )
        
                


                PATHS.LIST <- Create.Folders( Tag=sprintf("%s-%s%s.n%04i", DATA.TAG, OUTCOME, lSCENARIO$label, lSIZE ) )
                
                ##
                ## Generate subsets (by outcome[=column] and included/excluded[=rows])
                ## NOTE: excluded implicitly means cross-sectional, ie. only 'first' observation
                ##

                WHICH <- list()
                
                WHICH$STUDIES <- with(DATA, Study %in% c(lSCENARIO$base,lSCENARIO$added) )
                
                WHICH$SAMPLE.SIZE <- with(DATA, SIZE <= lSIZE )
                
                WHICH$BASELINE.CONTROL <- with(DATA, (INDEX.TYPE==levels(INDEX.TYPE)[1]) & (INDEX.OB==1))

                ## Following is outcome-specific code
                if(!is.null( attr(DATA,"columns")$Drop )){
                    MATCH <- match(x=sprintf("%s.DROP",OUTCOME), table=attr(DATA,"columns")$Drop )
                    if( !is.na(MATCH) ) {
                        WHICH$KEEP <- !DATA[, attr(DATA,"columns")$Drop ]
                        ## above we specify in terms of which rows to drop, so we must negate to keep those we want to KEEP
                        cat("Outcome specific subsetting:",OUTCOME," (dropping ",sum(!WHICH$KEEP,na.rm=TRUE)," rows)\n")
                    } else {
                        WHICH$KEEP <- TRUE
                    }
                } else {
                    WHICH$KEEP <- TRUE
                }

                
                ##
                ## Check for NAs in Outcome, Covariates, Index and Drop columns (not Additional, since they do not impact fitting by definition)
                WHICH$VALID <- Reduce(`&`, lapply( DATA[unlist( attr(DATA,"columns")[c(OUTCOME, "Covariates", "Index", "Drop")] )], function(X){!is.na(X)} ) )

                


                WHICH.COLUMNS <- c( OUTCOME, unlist( attr(DATA,"columns")[c("Covariates", "Additional", "Index", "Drop")] ) ) ## note we have omitted "Outcomes" element

                SUBSET <- droplevels( DATA[ Reduce( `&`, WHICH ), WHICH.COLUMNS ] )

                cat( "Subset", PATHS.LIST$Tag, "has", NROW(SUBSET), "rows.\n")

                attributes(SUBSET) <- c( attributes(SUBSET), attributes(DATA)[c("spec", "columns", "tag", "Transformations", "SCENARIOS")] )

                attr(SUBSET,"DATA.WHICH.LIST") <- WHICH
                
                attr(SUBSET,"iSCENARIO") <- iSCENARIO
                attr(SUBSET,"iSIZE") <- iSIZE


                Check.Attributes( SUBSET )
                ##
                saveRDS(object=SUBSET, file=file.path( PATHS.LIST$PATH, "SUBSET.rds"))


                ##
                ## Generate model sets
                ##

                ## Only generate models for the primary SUBSET
                ## For all the derived subsets, as extra scenarios, we will use the "best" model selected
                ## NOTE: There will be some code in the fitting script to link/copy the relevant model, then fit these scenario-subsets
                ##
                if( is.null(lSCENARIO$added) ) {
                

                    if( OUTCOME=="Wand" ) {

                        ## NOTE: FAMILY.SET allows us to explore multiple gamlss outcome distributions, later scripts will select the 'best' (via AIC/BIC/etc)
                        FAMILY.SET <- c("GGalt") ## c("GGalt","BCCG","GIG")
                        FP.SET <- matrix(c(1,0,0,
                                           1,1,0,
                                           1,0,1,
                                           2,0,0,
                                           2,1,0,
                                           2,0,1,
                                           2,1,1
                                           ),
                                         byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))
                        

                        for( lFAM in FAMILY.SET ) { ## loop to search multiple outcome distributions

                            for( iFP in 1:NROW(FP.SET) ) {

                                MODEL.NAME <- paste0("base",paste0(FP.SET[iFP,],collapse=""))

                                MODEL <- list(covariates=list("Y"=OUTCOME,
                                                              "X"="TimeTransformed",
                                                              "ID"="ID",
                                                              "BY"="Grp",
                                                              "OTHER"=NULL,
                                                              "COND"="Type", ## should be all equal to base case in fitted SUBSET
                                                              "RANEF"="Study"),
                                              family=lFAM,
                                              contrasts=list("Grp"="contr.sum"), ## (*1*)
                                              stratify=c("Study","Grp"),
                                              mu   =if(FP.SET[iFP,"mu"]>0){
                                                        sprintf("%s ~ 1 + fp(TimeTransformed,npoly=%i) + Grp + random(Study)",OUTCOME,FP.SET[iFP,"mu"])
                                                    } else {
                                                        sprintf("%s ~ 1 + Grp + random(Study)",OUTCOME)
                                                    },
                                              sigma=if(FP.SET[iFP,"sigma"]>0){
                                                        sprintf("%s ~ 1 + fp(TimeTransformed,npoly=%i) + Grp",OUTCOME,FP.SET[iFP,"sigma"])
                                                    } else {
                                                        sprintf("%s ~ 1 + Grp",OUTCOME)
                                                    },
                                              nu   =if(FP.SET[iFP,"nu"]>0){
                                                        sprintf("%s ~ 1 + fp(TimeTransformed,npoly=%i)",OUTCOME,FP.SET[iFP,"nu"])
                                                    } else {
                                                        sprintf("%s ~ 1",OUTCOME)
                                                    },
                                              inc.fp=TRUE)
                                saveRDS(object=MODEL,file=file.path(PATHS.LIST$MODEL,sprintf("%s.%s.fp.rds",MODEL.NAME,lFAM)))
                            }
                            
                        }
                    } else if ( OUTCOME=="Wild" ) {

                        FAMILY.SET <- c("NO") ## c("GGalt","BCCG","GIG")
                        for( lFAM in FAMILY.SET ) { ## loop to search multiple outcome distributions
                            MODEL.NAME <- "base"
                            MODEL <- list(covariates=list("Y"=OUTCOME,
                                                          "X"="TimeTransformed",
                                                          "ID"="ID",
                                                          "BY"="Grp",
                                                          "OTHER"=NULL,
                                                          "COND"="Type", ## should be all equal to base case in fitted SUBSET
                                                          "RANEF"="Study"),
                                          family=lFAM,
                                          contrasts=list("Grp"="contr.sum"), ## (*1*)
                                          stratify=c("Study","Grp"),
                                          mu   =sprintf("%s ~ 1 + fp(TimeTransformed,npoly=1) + Grp + random(Study)",OUTCOME),
                                          sigma=sprintf("%s ~ 1 + Grp",OUTCOME),
                                          nu   =sprintf("%s ~ 1",OUTCOME),
                                          inc.fp=TRUE)
                            saveRDS(object=MODEL,file=file.path(PATHS.LIST$MODEL,sprintf("%s.%s.fp.rds",MODEL.NAME,lFAM)))
                        }
                    }
                }

                ##
                ## Generate novel-testing-datasets
                ##
                if( is.null(lSCENARIO$added) ) {
                    
                    for( CLONE.ORIGINAL in c("B","Q") ) {
                        NOVEL.NAME <- sprintf("%s-CLONE",CLONE.ORIGINAL)
                        NOVEL <- SUBSET[ SUBSET$Study==CLONE.ORIGINAL, ] ## we're going to use the HCP study as a basis
                        NOVEL$Study <- factor( NOVEL.NAME )       ## MadeUp should end up the same as HCP... by definition

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
                ##
                ## Generate novel-scenario-testing
                ## NOTE: only for primary scenario, comparing : gamlss(SUBSET+NOVEL(1)), gamlss(SUBSET+NOVEL(2)), etc
                ##                                  with      : expand(gamlss(SUBSET),NOVEL(1)), expand(gamlss(SUBSET),NOVEL(2)), etc
                if( is.null(lSCENARIO$added) ) {

                    for( jSCENARIO in 1:length(attr(DATA,"SCENARIOS")) ) {
                        kSCENARIO <- attr(DATA,"SCENARIOS")[[jSCENARIO]]
                        if( is.null(kSCENARIO$added) ) {
                            ## if no added studies, then there is no novel data
                            next
                        }
                        for( jSIZE in 1:length(kSCENARIO$sizes) ) {
                            kSIZE <- kSCENARIO$sizes[jSIZE]

                            WHICH$SAMPLE.SIZE <- DATA$SIZE <= kSIZE
                            WHICH$STUDIES <- DATA$Study %in% kSCENARIO$added

                            
                            NOVEL <- droplevels( DATA[ Reduce(`&`,WHICH) , WHICH.COLUMNS ] )

                            NOVEL.TAG <- sprintf("%s%s.n%04i.rds", OUTCOME, kSCENARIO$label, kSIZE )
                            cat( "Novel", NOVEL.TAG, "has", NROW(NOVEL), "rows.\n")
                            
                            saveRDS( object=NOVEL, file=file.path(PATHS.LIST$NOVEL,NOVEL.TAG) )
                        }
                    }
                }
            }
        }
    }    
}


print( warnings() )

    
