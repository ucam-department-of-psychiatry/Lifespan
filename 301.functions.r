##
## Libraries (needed by functions in this file)
##
##

library("gamlss") ## need gamlss.family objects in all scripts

##
## Analysis Functions
##
##

gamlssWrapper <- function( Model, Data, Ctrl, ... ) { ##
    ## The problem:
    ##   (a) The gamlss() function returns an object with the contrast functions, not their names
    ##   (b) The model.matrix() function requires the name of the contrast functions
    ## In order for Add.Pred() to work, we need Extract.Param() to know the names of the contrast functions
    ## Hence, we need this wrapper function to attach the name of the functions as an attribute we can access
    ## While we're doing this, we might as well attach the Model object directly
    ## See (*1*) as a comment in other functions relating to this point
    if(is.null(Model$contrasts)){
        CONTRASTS.LIST <- NULL
    }else{
        CONTRASTS.LIST <- list()
        for ( LAB in names(Model$contrasts) ) {
            if( is.factor(Data[[LAB]]) ) {
                CONTRASTS.LIST[[LAB]] <- get(Model$contrasts[[LAB]])( levels(Data[[LAB]]) )
            } else {
                stop("Have provided contrasts for a non-factor column, ",LAB)
            }
        }
    }
    
    
    TMP <- gamlss(formula=as.formula(Model$mu),
                  sigma.formula=as.formula(Model$sigma),
                  nu.formula=as.formula(Model$nu),
                  tau.formula=as.formula(Model$tau),
                  family=Model$family,
                  contrasts=CONTRASTS.LIST, ## gamlss() requires the contrasts as a matrix or function
                  data=Data,
                  control=Ctrl, ... )
    attr( TMP, "model" )  <- Model
    return(TMP)
}

Extract.Wrapper <- function( Holder, Store.Full=FALSE, ... ) {

    OBJECT <- list()
    if( STORE.OPT$Fit.Summary ) {
        OBJECT$summary <-  Extract.Summary(Data=Holder$SUBSET, Stratified=Holder$MODEL$stratify )
    } else {
        OBJECT$summary <- NULL
    }

    FIT.FULL <- try(gamlssWrapper(Model=Holder$MODEL,
                                  Data=Holder$SUBSET[,unlist(Holder$MODEL$covariates[c("Y","X","BY","OTHER","RANEF")])],
                                  Ctrl=GAMLSS.CTRL,
                                  ... ),
                    silent=TRUE )

    if( any(class(FIT.FULL)=="try-error") ) {
        OBJECT$FIT.FULL <- NULL
        OBJECT$param <- NULL
        OBJECT$diagnostics <- "try-error"
        OBJECT$iter <- -1
    } else if( ! FIT.FULL$converged ) {
        OBJECT$FIT.FULL <- NULL
        OBJECT$param <- NULL
        OBJECT$diagnostics <- "not-converged"
        OBJECT$iter <- FIT.FULL$iter
    } else {
        if( Store.Full ) {
            OBJECT$FIT.FULL <- FIT.FULL
        }
        OBJECT$param <- Extract.Param( FIT.FULL, Holder$SUBSET )
        OBJECT$diagnostics <- NULL
        OBJECT$iter <- FIT.FULL$iter
    }
    return( OBJECT )
}

Extract.Summary <- function( Data=NULL, Stratified=NULL, Index=NULL ) {
    if( is.null(Data) ) {stop("No Data argument")}
    LIST <- list()
    if( !is.null(Index) ) {
        LIST$reps <- table(table(Index))
    } else {
        LIST$reps <- table(table(1:NROW(Data)))
    }
    if( !is.null(Stratified) ) {
        LIST$xtab <- xtabs(formula=~.,data=Data[Stratified])
    }
    class(LIST) <- c("SummaryObj","list")
    return(LIST)
}

Extract.Param <- function( object, HSubset ) {
    ## This function extracts the necessary information from a GAMLSS model fit object

    ## Extend this to include any information from a fitted model object that may be important
    ## However, it does not include predictions based on "new data", see Add.Pred() function
    
    if( !is.gamlss(object) ) { stop("Object is not a gamlss model") }

    PARAM <- list()
    PMODEL <- attr(object,"model")
    PARAM$family <- object$family[1]
    if( PMODEL$family != PARAM$family ) {stop("Something has gone wrong with the family argument")}
    PARAM$contrasts <- PMODEL$contrasts
    PARAM$stratify <- PMODEL$stratify
    PARAM$ID <- list()

    for( LAB in object$parameters ) {
        cat(LAB,"\n")
        
        PARAM[[LAB]] <- list()
        PARAM[[LAB]]$fixef <- coef( object, what=LAB )
        
        PARAM[[LAB]]$equ <- list()
        lEQU <- (rev(as.character(object[[sprintf("%s.formula",LAB)]]))[1])
        PARAM[[LAB]]$equ$fixef <- sub("random\\(.*\\)","1",lEQU)
        if( grepl("random",lEQU) ) {
            PARAM[[LAB]]$equ$ranef <- sub(".*random\\((.*)\\).*","\\1",lEQU)
        } else {
            PARAM[[LAB]]$equ$ranef <- NA
        }
        
        PARAM[[LAB]]$ID.sigma <- NULL ## since we are not fitting a within-subject repeated measure model we do not get a (1|Subject) random-effect variance
        CLASS <- sapply(getSmo(object,what=LAB,which=0),class)
        
        if( any(CLASS=="random") ) {
            lWHICH <- which(CLASS=="random")[1] ## CURRENT VERSION ONLY PERMITS ONE RANDOM-EFFECT PER EQUATION [ENHANCEMENT PLAN]
            if( is.null(getSmo(object,what=LAB,which=lWHICH)$sigb) ) {stop("Issue in extracting random-effect details from getSmo()")}
            PARAM[[LAB]]$sigma <- getSmo(object,what=LAB,which=lWHICH)$sigb
            PARAM[[LAB]]$ranef <- getSmo(object,what=LAB,which=lWHICH)$coef
            PARAM[[LAB]]$ranef.TYPE <- rep("gamlss",length(PARAM[[LAB]]$ranef)) ## added to later identify whether model-fitted or novel-mle (point) estimate
        }
        if( sum(CLASS=="random")>1 ) {
            PARAM[[LAB]]$warn <- TRUE
        } else {
            PARAM[[LAB]]$warn <- FALSE
        }
        if( any(CLASS=="lm") ) {
            if( length(which(CLASS=="lm")) > 1 ) {stop("Code only allows one lm-smoother per equation")}
            lWHICH <- which(CLASS=="lm")[1] ## CURRENT VERSION ONLY PERMITS ONE LM-SMOOTHER PER EQUATION [ENHANCEMENT PLAN]
            PARAM[[LAB]]$fp.smoother <- getSmo(object,what=LAB,which=lWHICH)[c("coefficients","power")]
        }
    }
    PARAM$AIC <- AIC( object )
    PARAM$BIC <- BIC( object )

    PARAM$N <- object$N
    PARAM$df.fit <- object$df.fit
    PARAM$df.residual <- object$df.residual

    if( 1 ) {
        ## Extra information added to ParamObj, required for out-of-sample estimation when we do not have access to SUBSET
        FactorColumns <- names(which(sapply(HSubset[unlist(attr(object,"model")$covariates[c("COND","X","BY","OTHER")])],is.factor)))
        attr(PARAM,"levels") <- lapply(HSubset[FactorColumns],levels)
        attr(PARAM,"transformations") <- attr(HSubset,"Transformations")

        if(!is.null(attr(object,"model")$contrasts)){
            FactorContrasts <- list()
            for ( LAB in names(attr(object,"model")$contrasts) ) {
                if( is.factor(HSubset[[LAB]]) ) {
                    FactorContrasts[[LAB]] <- get(attr(object,"model")$contrasts[[LAB]])( levels(HSubset[[LAB]]) )
                } else {
                    stop("Have provided contrasts for a non-factor column, ",LAB)
                }
            }
            attr(PARAM,"contrasts") <- FactorContrasts
        } else {
            attr(PARAM,"contrasts") <- NULL
        }
    }

    
    attr(PARAM,"model") <- attr(object,"model")
    class(PARAM) <- c("ParamObj","list")
    return(PARAM)
}


##
## Fitting Functions
##
##
Make.bfpNA.model.from.extract <- function( Param ) {
    LOCAL <- attr(Param,"model")
    if( (!is.null(LOCAL$inc.fp)) && (LOCAL$inc.fp==FALSE) ) {stop("Input model does not include fp(), so cannot convert to a bfp/bfpNA() model")}

    for( LAB in names(get(Param$family)()$parameters) ) {
        if( !is.null(Param[[LAB]]$fp.smoother) ) {
            if( !grepl("fp",LOCAL[[LAB]]) ) {stop("Fit.Extract includes a lm-smoother, but equation does not include fp()")}
            REPLACE <- sprintf("powers=c(%s)",paste0(Param[[LAB]]$fp.smoother$power,collapse=","))
            
            LOCAL[[LAB]] <- sub("npoly[ ]*=[ ]*[0-9]*",REPLACE,sub("fp","bfpNA",LOCAL[[LAB]]))
        }
    }
    LOCAL$inc.fp <- FALSE
    return( LOCAL )
}



Make.Longitudinal <- function( Holder ) {

    lPRED <- Holder$DATA.PRED

    Who.Long <- as.character(unique(lPRED[which(lPRED[,"INDEX.OB"]>1),"INDEX.ID"]))

    lCOV <- Holder$MODEL$covariates

    LONG <- data.frame(row.names=Who.Long,stringsAsFactors=FALSE)
    LONG[,sprintf("%s.first",lCOV$COND)] <- factor(NA,levels(lPRED[,lCOV$COND]))
    LONG[,sprintf("%s.last",lCOV$COND)] <- factor(NA,levels(lPRED[,lCOV$COND]))
    LONG[,lCOV$RANEF[1]] <- factor(NA,levels(lPRED[,lCOV$RANEF[1]]))
    for( VAR in lCOV$BY ) {
        LONG[,VAR] <- factor(NA,levels(lPRED[,VAR]))
    }
    TEN.PERCENT <- floor(length(Who.Long)/10)
    for( IDX in 1:length(Who.Long) ) {
        if( (IDX %% TEN.PERCENT) == 0 ) {cat(sprintf("%i [%s] / %i\n",IDX,WHO,length(Who.Long)))}
        WHO <- Who.Long[IDX]
        WHICH <- which(lPRED[,"INDEX.ID"]==WHO)

        Y.Values <- lPRED[WHICH,sprintf("%s",lCOV$Y)]
        if( all(is.na(Y.Values)) ) {
            LONG[WHO,sprintf("%s.max",lCOV$Y)]      <- NA
            LONG[WHO,sprintf("%s.rng",lCOV$Y)]      <- NA
            LONG[WHO,sprintf("%s.iqr",lCOV$Y)]      <- NA
        } else {
            LONG[WHO,sprintf("%s.max",lCOV$Y)]      <- max( Y.Values, na.rm=TRUE)
            LONG[WHO,sprintf("%s.rng",lCOV$Y)]      <- diff( range( Y.Values, na.rm=TRUE) )
            LONG[WHO,sprintf("%s.iqr",lCOV$Y)]      <- IQR( Y.Values ,na.rm=TRUE)
        }
        ##
        Q.Values <- lPRED[WHICH,sprintf("%s.q.wre",lCOV$Y)]
        if( all(is.na(Q.Values)) ) {
            LONG[WHO,sprintf("%s.q.max",lCOV$Y)]      <- NA
            LONG[WHO,sprintf("%s.q.rng",lCOV$Y)]      <- NA
            LONG[WHO,sprintf("%s.q.iqr",lCOV$Y)]      <- NA
        } else {
            LONG[WHO,sprintf("%s.q.max",lCOV$Y)]      <- max( Q.Values, na.rm=TRUE)
            LONG[WHO,sprintf("%s.q.rng",lCOV$Y)]      <- diff( range( Q.Values, na.rm=TRUE) )
            LONG[WHO,sprintf("%s.q.iqr",lCOV$Y)]      <- IQR( Q.Values, na.rm=TRUE)
        }
        ## add baseline age
        LONG[WHO,sprintf("%s.first",lCOV$X)]   <- min(lPRED[WHICH,lCOV$X],na.rm=TRUE)
        ## age spanned
        LONG[WHO,sprintf("%s.rng",lCOV$X)]        <- diff( range(lPRED[WHICH,lCOV$X],na.rm=TRUE) )
        ## number of time points
        LONG[WHO,"nobs"]       <- length(WHICH)

        LONG[WHO,sprintf("%s.first",lCOV$COND)] <- lPRED[WHICH[which.min(lPRED[WHICH,lCOV$X])],lCOV$COND]
        LONG[WHO,sprintf("%s.last",lCOV$COND)] <- lPRED[WHICH[which.max(lPRED[WHICH,lCOV$X])],lCOV$COND]
        LONG[WHO,sprintf("%s.nchanges",lCOV$COND)] <- sum(lPRED[WHICH[-1],lCOV$COND] != lPRED[rev(rev(WHICH)[-1]),lCOV$COND])

        ##
        ## Next we add covariate information, to assist with plotting
        LONG[WHO,lCOV$RANEF[1]] <- lPRED[WHICH[1],lCOV$RANEF[1]] ## by definition, INDEX.ID is unique across all individuals across all studies
        for( VAR in lCOV$BY ) {
            LONG[WHO,VAR] <- lPRED[WHICH[1],VAR] ## not checking if individuals change BY-variable over observations, assume these are constant?
        }

    }

    LONG[,"INDEX.ID"] <- factor(Who.Long,levels=levels(lPRED[,"INDEX.ID"]))
    
    return( LONG )
}

Boot.Function <- function( n, Base.Seed, Holder, Create.INIT=TRUE ) {
    source("100.common-variables.r")
    source("101.common-functions.r")

    source("300.variables.r")
    source("301.functions.r")
    ## must source file within BOOT.FUNCTION, otherwise contents not available witihn %dopar%

    set.seed( seed=Base.Seed + n )
    if(is.null(Holder$MODEL$stratify)) {
        INDEX <- sample(1:NROW(Holder$SUBSET),NROW(Holder$SUBSET),replace=TRUE) ## unstratified bootstrap
    } else {
        SPLIT <- split(1:NROW(Holder$SUBSET),Holder$SUBSET[Holder$MODEL$stratify])
        LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
        INDEX <- unsplit(LAPPLY,Holder$SUBSET[Holder$MODEL$stratify])
    }
    if( is.null(INDEX) ) { return(NA) }
    INDEX <- sort(INDEX)

    
    Holder$SUBSET <- Holder$SUBSET[INDEX,] ## replace SUBSET with bootstrapped-SUBSET

    ##
    ## In the following code we are creating an initialisation block to use within the bootstrapping
    ## NOTE: this must be done within this function, since it depends on the bootstrap sample drawn
    if( Create.INIT ) {
        tPRED <- Apply.Param(NEWData=Holder$SUBSET, FITParam=Holder$FIT.EXTRACT$param,
                             Pred.Set=NULL, Prefix="", Add.Moments=FALSE, Add.Normalise=FALSE, Add.Derivative=FALSE, MissingToZero=TRUE,
                             verbose=FALSE )

        INIT <- structure( list(parameters=names(Holder$FAMILY$parameters)), class="gamlss" )
        for( lPAR in INIT$parameters ) {
            INIT[[sprintf("%s.fv",lPAR)]] <- unname(as.numeric(Holder$FAMILY[[sprintf("%s.linkinv",lPAR)]]( tPRED[,sprintf("%s.pop",lPAR)] )))
        }
    } else {
        INIT <- NULL
    }

    EXTRACT <- Extract.Wrapper( Holder, Store.Full=FALSE, start.from=INIT ) ## DO NOT store full fitting object
    
    if( STORE.OPT$Boot.Index ) {
        EXTRACT$index <- INDEX
    } else {
        EXTRACT$index <- NULL
    }
    EXTRACT$base <- Base.Seed
    EXTRACT$offset <- n
    
    return( EXTRACT )
}



##
##
##
Apply.Param <- function(NEWData, FITParam,
                        Reference.Holder=NULL,
                        MissingToZero=TRUE, NAToZero=TRUE, Prefix="",
                        Pred.Name="PRED", Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                        Add.Normalise=FALSE,
                        Add.Moments=TRUE, Add.Derivative=TRUE,
                        Only.New=FALSE,
                        verbose=FALSE ) {
    if( 0 ) {
        
        
        NEWData=NEWDATA
        FITParam=FIT$param

        NEWData=PRIMARY$DATA
        FITParam=PRIMARY$FIT.EXTRACT$param
        Reference.Holder=PRIMARY
        
        Prefix=""
        Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975)
        MissingToZero=TRUE
        NAToZero=TRUE
        verbose=FALSE
        Pred.Name="PRED"

        Add.Moments=FALSE
        Add.Normalise=TRUE ## no outcome, so does not make sense
        Add.Derivative=FALSE
        ## or
        Add.Moments=TRUE
        Add.Normalise=FALSE ## no outcome, so does not make sense
        Add.Derivative=TRUE        


    }

    ##
    ## check arguments
    if( is.null(NEWData) || is.null(FITParam) ) { stop("Must supply valid NEWData and FITParam objects") }

    if( !is.null(Reference.Holder) ) {
        NEWData <- ValidateCleanInput(IN=NEWData,
                                      Reference.Subset=Reference.Holder$SUBSET,
                                      Reference.Model=Reference.Holder$MODEL,
                                      Reference.Param=Reference.Holder$FIT.EXTRACT$param )
        ## above line checks that input conforms to reference
        ##
    }
    
    ##
    ## vvv not sure this is necessary, but a precaution against mixing with model.frame(), model.matrix(), etc.
    Saved.Attributes <- attributes(NEWData)
    attributes(NEWData)[ !names(attributes(NEWData)) %in% c("names","row.names","class") ]  <- NULL

    if( Only.New ) {
        stop("This isn't quite right at the moment, need to fix this to perhaps Only.Modified?")
        IN.NAMES <- names(NEWData)
    }

    In.Model <- attr(FITParam,"model")
    
    Model.Columns <- unlist(In.Model$covariates[c("X","BY","OTHER","RANEF")])
    
    FAMILY <- get(FITParam$family)()
    
    Model.Parameters <- names(FAMILY$parameters)

    Any.Missing.Columns <- Model.Columns[!c(Model.Columns %in% names(NEWData))]
    if( length(Any.Missing.Columns) > 0 ) {
        NEWData[ , Any.Missing.Columns ] <- NA
        if( verbose ) { cat("Missing columns in NEWData from within Model. Set to NA. [", Any.Missing.Columns,"]\n") }
    }

    for( wPARAM in Model.Parameters ) {
        Model.Formula <- as.formula(sprintf(" ~ %s",FITParam[[ wPARAM ]]$equ$fixef))
        Model.Frame <- model.frame(formula=Model.Formula, data=NEWData, na.action=na.pass )
        Model.Matrix <- model.matrix( Model.Formula, Model.Frame, contrasts.arg=FITParam$contrasts, na.action=na.pass )

        if( MissingToZero==TRUE ) {
            lASSIGN <- attr(Model.Matrix,"assign")
            lAllNA <- apply(Model.Matrix,2,function(X){all(is.na(X))})
            if( length(lASSIGN)!=length(lAllNA) ) {stop("FATAL ERROR: lASSIGN different length from lAllNA")}
            if( sum(lAllNA) > 0 ) {
                lLABELS <- paste( attr(terms(Model.Formula),"term.labels")[unique(lASSIGN[lAllNA])], collapse=", " )
                if( verbose ) { 
                    cat(wPARAM,"Some columns of the MODEL MATRIX are all NAs (",lLABELS,"), setting to zero in MODEL MATRIX.\n")
                    cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
                }
                Model.Matrix[ , which(lAllNA) ]  <- 0
            }
        }
        if( NAToZero==TRUE ) {
            local.sum <- sum(apply(Model.Matrix,1,function(X){any(is.na(X))}))
            if( local.sum > 0 ) {
                ARR.IND <- which( is.na(Model.Matrix), arr.ind=TRUE )

                lCOLUMNS <- paste( attr(terms(Model.Formula),"term.labels")[unique(attr(Model.Matrix,"assign")[unique(ARR.IND[,"col"])])], collapse=", " )
                STR <- sprintf("MODEL MATRIX has %i rows (across columns: %s ) with NAs, all set to zero",local.sum,lCOLUMNS)
                if( verbose ) { 
                    cat(wPARAM,STR,"\n")
                    cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
                }
                Model.Matrix[ ARR.IND ]  <- 0
            }
        }
        if( !all( colnames(Model.Matrix) %in% names(FITParam[[ wPARAM ]]$fixef) ) ) {stop("FITParam inconsistent with INPUT")}
        Fit.Fixef <- matrix( FITParam[[ wPARAM ]]$fixef[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta") )

        NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] <- as.vector(Model.Matrix %*% Fit.Fixef)

        if( !is.na( FITParam[[ wPARAM ]]$equ$ranef ) ) {
            ##
            ## In this block we find the matching random-effect (from FITParam) and insert it into NEWData
            ## NOTE: This code block assumes only one random-effect (per gamlss-component) [ENHANCEMENT PLAN: expand to multiple random-effects]
            ## 
            POSITION.IN.MAP <- as.numeric(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ])
            MAPPING <- match( levels(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ]), names(FITParam[[ wPARAM ]]$ranef) )
            NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- FITParam[[ wPARAM ]]$ranef[ MAPPING[POSITION.IN.MAP] ]

            NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] + NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)]
        } else {
            NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- Inf ## to differentiate from missing random-effect levels (NAs), use Inf if no ranef at all
            NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)]
        }
    }


    ARGUMENTS <- list()
    for( lTYPE in c("pop","wre") ) {
        ARGUMENTS[[ lTYPE ]] <- list()
        for( wPARAM in Model.Parameters ) {
            LinkInvFun <- FAMILY[[sprintf(sprintf("%s.linkinv",wPARAM))]]
            ComponentValues <- NEWData[,sprintf("%s%s.%s",Prefix,wPARAM,lTYPE)]
            
            ARGUMENTS[[ lTYPE ]][[ wPARAM ]] <- LinkInvFun( ComponentValues )
        }
    }

    ##
    ## Select rows (within arguments) that are valid, ie. not NAs
    KEEP <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
    SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ KEEP[[IDX]] ]) } ) } ), names(KEEP) )
    WHICH <- lapply( KEEP, which )

    ##
    ## in below we make versions that account for NAs in outcome (OUT.NAME) column
    OUT.NAME <- In.Model$covariates$Y
    if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
        out.KEEP.ARGS <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
        out.KEEP <- lapply( out.KEEP.ARGS, function(X) { (X) & (!is.na(NEWData[,OUT.NAME])) } )
        out.SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ out.KEEP[[IDX]] ]) } ) } ), names(out.KEEP) )
        out.WHICH <- lapply( out.KEEP, which )

        KEEP.C <- Reduce(`&`,out.KEEP)
        SHORT.C <- lapply( ARGUMENTS, function(Y){ lapply(Y,function(X){ X[KEEP.C] } ) } )
        WHICH.C <- which( KEEP.C )
    }

    for( lTYPE in c("pop","wre") ) {
        if( length(WHICH[[lTYPE]])>0 ) {
            if( !is.null(Pred.Set) ) {
                for( lLEVEL in names(Pred.Set) ) {
                    NEWData[WHICH[[lTYPE]],sprintf("%s%s.%s.%s",Prefix,Pred.Name,lLEVEL,lTYPE)] <- do.call(what=get(paste0("q",FAMILY$family[1])),
                                                                                                       args=c(SHORT[[lTYPE]],list(p=Pred.Set[lLEVEL])))
                }
            }
            if( Add.Moments ) {
                NEWData[WHICH[[lTYPE]],sprintf("%s%s.mean.%s",Prefix,Pred.Name,lTYPE)] <- do.call(what=FAMILY$mean,args=SHORT[[lTYPE]])
                NEWData[WHICH[[lTYPE]],sprintf("%s%s.variance.%s",Prefix,Pred.Name,lTYPE)] <- do.call(what=FAMILY$variance,args=SHORT[[lTYPE]])
                ##
                ## [ENHANCEMENT PLAN: Add 3rd and 4th order moments (skew and kurtosis)
                ##
            }
        
            if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
                
                NEWData[out.WHICH[[lTYPE]],sprintf("%s%s.q.%s",Prefix,OUT.NAME,lTYPE)] <- do.call(what=get(paste0("p",FAMILY$family[1])),
                                                                                                  args=c(out.SHORT[[lTYPE]],
                                                                                                         list(q=NEWData[out.WHICH[[lTYPE]],OUT.NAME,drop=TRUE])))
            }
        }
    }

    if( Add.Normalise && (OUT.NAME %in% names(NEWData)) ) {
        if( (length(WHICH.C)>0) ) {
            NEWData[WHICH.C,sprintf("%s%s.normalised",Prefix,OUT.NAME)] <- do.call(what=get(paste0("q",FAMILY$family[1])),
                                                                                   args=c(SHORT.C[["pop"]],
                                                                                          list(p=NEWData[WHICH.C,sprintf("%s%s.q.wre",Prefix,OUT.NAME)]  )))
        }
    }

    if( (Add.Derivative) & (!Only.New) ) {
        NEWData[,sprintf("%sorig.pos.D",Prefix)] <- 1:NROW(NEWData) ## store original row order from input NEWData
        
        if( "INDEX.IDX" %in% names(NEWData) ) {
            if( any(is.na(NEWData[,"INDEX.IDX"])) ) {stop("Cannot have NAs within INDEX.ID for Add.Derivative")}
            cat("Unexpected to have INDEX.ID within a call to Add.Derivative?\n")
            cat("NOTE: Code does not currently accommodate Subject-level random-effects, if/when it does, this would make sense for a subject-specific deriviative.")
            if( "INDEX.OB" %in% names(NEWData) ) {
                NEWData <- NEWData[ order( NEWData[,"INDEX.IDX"], NEWData[,"INDEX.OB"] ), ]
                stop("The Add.Derivative code assumes the observations are equally spaced, and that the spacing is 'small' relative to the rate of change (for the gradient to be approximated. Using Add.Derivative like this seems suspicious.")
            } else {
                stop("Must recalculate INDEX.OB, which seems suspicious.")
            }
        } else {
            lNAMES <- unlist(In.Model$covariates[c("RANEF","COND","BY","OTHER","X")])
            kNAMES <- lNAMES[ lNAMES %in% names(NEWData) ]
            if( any(duplicated(NEWData[,kNAMES])) ) { stop("Duplicated rows (based on Model$covariates included in input NEWData). This makes no sense for Add.Derivative.") }
            ORDER <- do.call(what=order, args=NEWData[ kNAMES ] ) ## using default behaviour, NAs placed last
            NEWData <- NEWData[ ORDER, ]

            lNAMES <- unlist(In.Model$covariates[c("RANEF","COND","BY","OTHER")]) ## NOTE WE HAVE REMOVED 'X'-covariate from this line
                                                                                  ## want to split by the other covariates only
            kNAMES <- lNAMES[ lNAMES %in% names(NEWData) ]
            lLIST <- as.list(NEWData[ kNAMES ])
            kLIST <- lLIST[ !sapply(lLIST,function(Z){all(is.na(Z))}) ]
            SPLIT <- split( x=1:NROW(NEWData), f=kLIST )
            NEWData[,sprintf("%sIDX.D",Prefix)] <- Reduce(f=c,x=lapply( SPLIT[lengths(SPLIT)>0], function(X){ 1:length(X) } ))
        }

        NEWData[,sprintf("%sdx.D",Prefix)] <- NEWData[c(1:NROW(NEWData)),In.Model$covariates$X] - NEWData[c(NA,1:(NROW(NEWData)-1)),In.Model$covariates$X]
        NEWData[which(NEWData[,sprintf("%sIDX.D",Prefix)]==1),sprintf("%sdx.D",Prefix)] <- NA
        lDX <- NEWData[,sprintf("%sdx.D",Prefix)]
        for( lTYPE in c("pop","wre") ) {
            if( !is.null(Pred.Set) ) {
                for( lLEVEL in names(Pred.Set) ) {
                    lLAB <- sprintf("%s%s.%s.%s",Prefix,Pred.Name,lLEVEL,lTYPE)
                    if( lLAB %in% names(NEWData) ) {
                        lDY <- NEWData[ c(1:NROW(NEWData)), lLAB ] - NEWData[ c(NA,1:(NROW(NEWData)-1)), lLAB ]
                        NEWData[,sprintf("%s.D",lLAB)] <- lDY / lDX
                    }
                }
            }
            if( Add.Moments ) {
                lLAB <- sprintf("%s%s.mean.%s",Prefix,Pred.Name,lTYPE)
                if( lLAB %in% names(NEWData) ) {
                    lDY <- NEWData[ c(1:NROW(NEWData)), lLAB ] - NEWData[ c(NA,1:(NROW(NEWData)-1)), lLAB ]
                    NEWData[,sprintf("%s.D",lLAB)] <- lDY / lDX
                }
                lLAB <- sprintf("%s%s.variance.%s",Prefix,Pred.Name,lTYPE)
                if( lLAB %in% names(NEWData) ) {
                    lDY <- NEWData[ c(1:NROW(NEWData)), lLAB ] - NEWData[ c(NA,1:(NROW(NEWData)-1)), lLAB ]
                    NEWData[,sprintf("%s.D",lLAB)] <- lDY / lDX
                }
            }
            if( (OUT.NAME %in% names(NEWData)) ) {
                STR <- "Add.Derivative is ON, but NEWData includes the outcome column?"
                if( verbose ) { cat(STR,"\n") }
            }
        }
    } else if( (Add.Derivative) & (Only.New) ) {
        stop("Add derivative may re-order NEWData, making Only.New very dangerous" )
    }
    

    if( Only.New ) {
        NEWData <- NEWData[,names(NEWData)[!names(NEWData)%in%IN.NAMES]]
    }
    attributes(NEWData) <- c( attributes(NEWData), Saved.Attributes[!names(Saved.Attributes) %in% c("names","row.names","class")] )
    
    return( NEWData )
}


Apply.FitAndBoot <- function( NEWDATA=NULL, FIT=NULL, BOOT=NULL, Boot.wre=FALSE, QUANTILES=c(0.025,0.975), VERBOSE=TRUE, ... ) {


    ##
    ## We'll be passing all arguments down to two calls to Apply.Param()
    if( is.null(NEWDATA) || is.null(FIT) ) stop("Need NEWDATA, FIT.EXTRACT and BOOT.EXTRACT-list objects")

    CURVE <- Apply.Param(NEWData=NEWDATA, FITParam=FIT$param, verbose=VERBOSE, ... )
    
    if( Boot.wre ) {
        BOOT.COLUMNS <- grep("\\.pop|\\.wre",names(CURVE),value=TRUE)
    } else {
        BOOT.COLUMNS <- grep("\\.pop",names(CURVE),value=TRUE)
    }
    
    if( !is.null( BOOT ) ) {
        BOOT.KEEP <- sapply( BOOT, function(X){!is.null(X$param)})
        BOOT.LIST <- lapply( BOOT[BOOT.KEEP], function(X) {
            as.list(Apply.Param(NEWData=CURVE, FITParam=X$param, verbose=FALSE, ... )[BOOT.COLUMNS] )
        })
        BOOT.CI <- setNames(lapply( BOOT.COLUMNS, function(X) {
            apply( sapply(BOOT.LIST,function(Y){ Y[[X]] } ), 1, quantile, prob=QUANTILES, na.rm=TRUE )
        }),BOOT.COLUMNS)

        for( LAB in BOOT.COLUMNS ) {
            CURVE[,sprintf("%s.lower",LAB)] <- BOOT.CI[[ LAB ]][1,]
            CURVE[,sprintf("%s.upper",LAB)] <- BOOT.CI[[ LAB ]][2,]
        }
    } else {
        warning("BOOT.EXTRACT was null?")
    }
    return( CURVE )
}

Save.Extracted <- function( Extracted, Paths, Tag, Save.Full=FALSE ) {
    if( Save.Full && (!is.null(Extracted$FIT.FULL)) ) {
        saveRDS(object=Extracted$FIT.FULL,file=file.path( Paths$FIT.FULL, Tag ))
    }
    
    saveRDS( object=Extracted[ names(Extracted)[ !names(Extracted)%in%c("FIT.FULL") ] ], file=file.path( Paths$FIT.EXTRACT, Tag ) )
}

Load.Subset.Wrapper <- function(Tag,
                                LSubset=FALSE,
                                LModel=FALSE,
                                LData=FALSE,
                                LFit=FALSE,
                                LFitFull=FALSE,
                                LBoot=FALSE,
                                Anchor=RDS.DIR ) {

    OBJ <- list()

    OBJ$SUBSET.TAG <- Tag
    OBJ$PATH <- file.path( Anchor, Tag )

    if( LSubset ) {
        OBJ$SUBSET <- readRDS(file=file.path(OBJ$PATH,"SUBSET.rds"))
    }
    if( LData ) {
        if( LSubset ) {
            OBJ$DATA.TAG <- attr(OBJ$SUBSET,"tag")
            OBJ$DATA <- readRDS(file=file.path(Anchor,OBJ$DATA.TAG,"DATA.rds"))
        } else {
            stop("Require a SUBSET to be loaded, to access attr(.,'tag').")
        }
    }    
    if( LModel ) {
        OBJ$MODEL <- readRDS(file=file.path(OBJ$PATH,"MODEL.rds"))

        OBJ$FAMILY <- get(OBJ$MODEL$family)()
    }
    if( LFitFull ) {
        OBJ$FIT.FULL <- readRDS(file=file.path(OBJ$PATH,"FIT.FULL.rds"))
    }
    if( LBoot ) {
        OBJ$BOOT.EXTRACT <- readRDS(file=file.path(OBJ$PATH,"BOOT.EXTRACT.rds"))
    }
    if( LFit ) {
        OBJ$FIT.EXTRACT <- readRDS(file=file.path(OBJ$PATH,"FIT.EXTRACT.rds"))
    }
    return(OBJ)
}

Find.Models.To.Fit <- function( Paths.List ) {

    Model.Set.All <- list.files( Paths.List$MODEL )
    Model.Set.notbfp <- grep("bfpNA",Model.Set.All,value=TRUE,invert=TRUE) ## DO NOT FIT bfpNA models, they are derived from the matching fp model

    for( IDX in 1:length(Model.Set.notbfp) ) {
        if( file.exists( file.path(Paths.List$FIT.EXTRACT,Model.Set.notbfp[IDX]) ) ) {
            if( STORE.OPT$Fit.Overwrite=="stop" ) {
                stop(Paths.List$Tag,Model.Set.notbfp[IDX],"Model already has an existing fit-extract output.")
            }
            if( STORE.OPT$Fit.Overwrite=="overwrite" ) {
                STR <- sprintf("%s %s Model already has an existing fit-extract output [OVERWRITE].",Paths.List$Tag, Model.Set.notbfp[IDX])
                cat(STR,"\n")
                warning(STR)
            }
            if( STORE.OPT$Fit.Overwrite=="skip" ) {
                STR <- sprintf("%s %s Model already has an existing fit-extract output [SKIP].",Paths.List$Tag, Model.Set.notbfp[IDX])
                warning(STR)
                cat(STR,"\n")
                Model.Set.notbfp[IDX] <- NA
            }
        }
    }
    return( na.omit( Model.Set.notbfp ) )
}

ValidateCleanInput <- function( IN=NULL, Reference.Subset=NULL, Reference.Model=NULL, Reference.Param=NULL ) {

    if( is.null(IN) || is.null(Reference.Subset) || is.null(Reference.Model) || is.null(Reference.Param) ) {
        stop("Input argument is NULL")
    }

    ##
    ## Check for pipeline columns
    ##
    if( ! all( attr(Reference.Subset,"columns")$Index %in% names(IN) ) ) {
        stop("One of COLUMNS$Index is missing")
    }

    ##
    ## Check for model columns
    ##
    if( ! all(unlist(Reference.Model$covariates) %in% names(IN)) ) {
        stop("IN is missing columns included in the model")
    }
    for( localCOL in unlist(Reference.Model$covariates) ) {
        if( !all(class( Reference.Subset[,localCOL] ) == class( IN[, localCOL ] ) ) ) {
            cat(localCOL, "has different class:", class( Reference.Subset[,localCOL] ), "!=", class( IN[, localCOL ] ), "\n" )
            stop("IN and reference SUBSET have different class(es) for a columns")
        }
    }
    
    ##
    ## Check factor() levels
    ##
    COLUMNS.TO.CHECK <- unlist(Reference.Model$covariates[ !(names(Reference.Model$covariates) %in% c("Y","COND","RANEF","ID")) ])
    ## NOTE: in above we exclude
    ##       [Y] handle that specially below
    ##       [ID column] will clearly differ in factor levels
    ##       [RANEF column] will clearly differ, since these are new studies
    ##       [COND column] by definition, SUBSET is only base-group/healthy-controls/controls/etc, so of NOVEL includes a 'patient' group it will differ
    ## NOTE: we take is as a given that SUBSET conforms to the MODEL

    FACTOR.COLUMNS <- COLUMNS.TO.CHECK[which(sapply( IN[ COLUMNS.TO.CHECK ], is.factor ))]
    
    for( LAB in FACTOR.COLUMNS ) {
        if( ! is.factor(Reference.Subset[,LAB]) ) {
            stop("FATAL: NOVEL dataset has a Model Covariate as a factor, but in the  SUBSET data (which was used to fit the model) it is not a factor")
        }
        if( (nlevels(IN[,LAB])==nlevels(Reference.Subset[,LAB])) && all(levels(IN[,LAB])==levels(Reference.Subset[,LAB])) ) {
            ## do nothing
        } else {
            IN[,sprintf("%s.original",LAB)] <- IN[,LAB] ## save copy of column
            
            if( any( !(levels(IN[,LAB]) %in% levels(Reference.Subset[,LAB]) ) ) ) {
                TEMP <- levels(IN[,LAB])[ !(levels(IN[,LAB]) %in% levels(Reference.Subset[,LAB]) ) ]
                STR <- sprintf("[%s %s] Levels of %s (%s) in INPUT do not exist in Reference SUBSET, they have been set to NA (since they are not in the FIT.EXTRACT model",
                               attr(Reference.Subset,"tag"),Reference.Model$covariates$Y,LAB,paste(TEMP,collapse=", "))
                cat(STR,"\n")
                warning(STR)
            }
            IN[,LAB] <- factor( x=levels(IN[,LAB])[as.numeric(IN[,LAB])], levels=levels(Reference.Subset[,LAB]) )
            ## Do the above to ensure consistency in factor levels (of NOVEL with SUBSET)
        }
    }

    nonFACTOR.COLUMNS <- COLUMNS.TO.CHECK[which(sapply( IN[ COLUMNS.TO.CHECK ], function(X){!is.factor(X)} ))]
    
    LOGICAL <- list()
    LOGICAL[["factor.notNA"]] <- Reduce( f=`&`, lapply(IN[FACTOR.COLUMNS],function(X){!is.na(X)}) )
    LOGICAL[["nonfactor.notNA"]] <- Reduce( f=`&`, lapply(IN[nonFACTOR.COLUMNS],function(X){!is.na(X)}) )
    LOGICAL[["Y.notNA"]] <- (!is.na(IN[,Reference.Model$covariates$Y]))
    LOGICAL[["OB"]] <- (IN[,"INDEX.OB"]==1)
    LOGICAL[["TYPE"]] <- (IN[,"INDEX.TYPE"]==levels(IN[,"INDEX.TYPE"])[1])

    FIND.FITTED <- Find.Fitted.Levels( Reference.Param )
    FIND.MISSING <- Find.Missing.Levels( IN, FIND.FITTED, SEPARATOR="|" )

    LOGICAL[["MISSING.RANEF"]] <- FIND.MISSING$Obs
    LOGICAL[["KNOWN.RANEF"]] <- ( ! FIND.MISSING$Obs )
    
    LOGICAL[["REFIT.VALID"]] <- Reduce( `&`, LOGICAL[ c("Y.notNA","factor.notNA","nonfactor.notNA","OB","MISSING.RANEF","TYPE") ] )
    LOGICAL[["PRED.VALID"]] <- Reduce( `&`, LOGICAL[ c("factor.notNA","nonfactor.notNA","OB","KNOWN.RANEF","TYPE") ] )
    LOGICAL[["NORM.VALID"]] <- Reduce( `&`, LOGICAL[ c("Y.notNA","factor.notNA","nonfactor.notNA","OB","KNOWN.RANEF","TYPE") ] )
    
    attr(IN,"missing.levels") <- FIND.MISSING
    attr(IN,"logical.selectors") <- LOGICAL
    
    return(IN)
}



Find.Fitted.Levels <- function( Param ) {
    GAMLSS.PARAMS <- names(get(Param$family)()$param)
    FITTED.LEVELS <- list()
    for( lparam in GAMLSS.PARAMS ) {
        if( is.null( Param[[ lparam ]]$equ$ranef) || is.na( Param[[ lparam ]]$equ$ranef ) ) {
            ## no random-effect within this gamlss parameter
        } else {
            FITTED.LEVELS[[lparam]] <- list()
            ##
            ## next few lines assume parameter$ranef is a single vector, and parameter$equ$ranef is a one-length vector
            ## since we are constrained to one random-effect per gamlss parameter
            ## [ENHANCEMENT PLAN: expand this]
            ## [WARNING: this expansion of the code capability will NOT be backward compatible...]
            if( length(Param[[ lparam ]]$equ$ranef)>1 ) { stop("Implies multiple random-effects within the same gamlss parameter") }
            if( is.list(Param[[ lparam ]]$ranef) ) { stop("Implies multiple random-effects within the same gamlss parameter") }
            FITTED.LEVELS[[lparam]][[ Param[[ lparam ]]$equ$ranef ]]  <- Param[[ lparam ]]$ranef
        }
    }
    return(FITTED.LEVELS)
}

Find.Missing.Levels <- function( Novel, Fitted.Lvls, SEPARATOR="|" ) {

    MISSING.VECTOR <- numeric(length=0)
    MISSING.LEVELS <- list()
    NOVEL.OBSERVATIONS <- list()
    COUNTER <- 0
    for( lparam in names(Fitted.Lvls) ) {
        MISSING.LEVELS[[lparam]] <- list()
        for( lranef in names(Fitted.Lvls[[lparam]]) ) {
            ##
            ##MISSING.LEVELS[[lparam]][[lranef]] <- levels(INTERIM[,lranef])[ !levels(INTERIM[,lranef]) %in% names(FITTED.LEVELS[[lparam]][[lranef]]) ]
            ## ^^ in line above, must remove levels() aspect, since we have not (and cannot) do droplevels, well we can, but only on RANEF columns
            ## vv replaced with below, using unique() to find actual levels in INTERIM(ie NOVEL) dataset
            in.NOVEL <- as.character(unique(Novel[,lranef]))
            in.FITTED <- names(Fitted.Lvls[[lparam]][[lranef]])
            MISSING.LEVELS[[lparam]][[lranef]] <- in.NOVEL[ !in.NOVEL %in% in.FITTED ]
            
            lLEN <- length(MISSING.LEVELS[[lparam]][[lranef]])
            if( lLEN>0 ) {
                MISSING.VECTOR[ COUNTER+(1:lLEN) ]  <- rep(0,lLEN)
                names(MISSING.VECTOR)[ COUNTER+(1:lLEN) ] <- sprintf("%s%s%s%s%s",lparam,SEPARATOR,lranef,SEPARATOR,MISSING.LEVELS[[lparam]][[lranef]])
                attr(MISSING.LEVELS[[lparam]][[lranef]],"index") <- COUNTER+(1:lLEN)
                COUNTER <- COUNTER + lLEN
            }
            NOVEL.OBSERVATIONS[[sprintf("%s%s%s",lparam,SEPARATOR,lranef)]] <- ( Novel[,lranef] %in% MISSING.LEVELS[[lparam]][[lranef]] )
        }
    }
    MISSING.OBS <- Reduce(f=`|`,x=NOVEL.OBSERVATIONS)
    return( list(Vector=MISSING.VECTOR,Levels=MISSING.LEVELS,Obs=MISSING.OBS) )
}



##
## MLE refitting functions
##
Ranef.MLE.Func <- function( theta, Param, Missing, Novel, Prefix="", Return="optim" ) {
    if( missing(theta)|missing(Param)|missing(Missing)|missing(Novel) ) {stop("Mandatory argmument(s) missing")}
    if(length(theta)!=sum(sapply(Missing,lengths))){stop("Problem: length(theta) != number of Missing levels")}

    LL.ranef <- list()
    for( lIDX in 1:length(Missing) ) {
        LAB <- names(Missing)[lIDX]
        if( length(Missing[[lIDX]])> 1 ) { stop("FATAL ERROR: MLE function currently assumes a single random-effect per gamlss parameter") }
        JDX <- attr(Missing[[lIDX]][[1]],"index") ## NOTE: making assumption of a single random-effect within each gamlss parameter

        LL.ranef[[LAB]] <- dnorm(x=theta[JDX],mean=0,sd=Param[[LAB]]$sigma,log=TRUE) ## this can be a vector of length>1

        lMATCH <- match( Novel[ , names(Missing[[lIDX]])[1] ], Missing[[lIDX]][[1]] )

        if( any(is.na(lMATCH)) ) { stop("This should not happen") }
        Novel[,sprintf("%s%s.wre",Prefix,LAB)]  <-  Novel[,sprintf("%s%s.pop",Prefix,LAB),drop=TRUE] + theta[JDX[lMATCH]]
    }
    
    DIST <- get( Param$family )()
    lARGS <- list()
    CHECK <- DIST$parameters
    for( LAB in names(DIST$parameters) ) {
        lARGS[[LAB]] <- DIST[[sprintf("%s.linkinv",LAB)]]( Novel[,sprintf("%s.wre",LAB)] )
        CHECK[[LAB]] <- DIST[[sprintf("%s.valid",LAB)]]( lARGS[[LAB]] )
    }
    if( !all(unlist(CHECK)) ) {stop("Failed distribution parameter checks")}

    LL.out <- do.call( what=get(paste0("d",Param$family)), args=c(lARGS,list(x=Novel[,attr(Param,"model")$covariates$Y],log=TRUE)))

    if( Return=="optim" ) {
        -1 * ( sum(LL.out) + sum(unlist(LL.ranef)) )
    } else if ( Return=="LL" ) {
        sum(LL.out)
    } else if ( Return=="LL+RE" ) {
        ( sum(LL.out) + sum(unlist(LL.ranef)) )
    }
}

Add.New.Ranefs <- function( new.vector, Fit.Extract, Missing ) {
    if( missing(new.vector)|missing(Fit.Extract)|missing(Missing) ) {stop("Mandatory arguement(s) missing")}
    for( lIDX in 1:length(Missing) ) {
        LAB <- names(Missing)[lIDX]

        JDX <- attr(Missing[[lIDX]][[1]],"index")

        Fit.Extract[[LAB]]$ranef <- append( Fit.Extract[[LAB]]$ranef, setNames( new.vector[JDX], Missing[[lIDX]][[1]] ) )

        Fit.Extract[[LAB]]$ranef.TYPE <- append( Fit.Extract[[LAB]]$ranef.TYPE, rep("newmle",length(JDX)) )
    }
    return( Fit.Extract )
}







Calc.Expanded <- function( NewData, Cur.Param, Missing, Prefix="" ) {

    OPT <- optim(par=Missing$Vector, fn=Ranef.MLE.Func,
                 Param=Cur.Param, Missing=Missing$Levels, Novel=NewData, Prefix=Prefix,
                 method=if(length(Missing$Vector)==1){"Brent"}else{"Nelder-Mead"},
                 lower=if(length(Missing$Vector)==1){-1000}else{-Inf},
                 upper=if(length(Missing$Vector)==1){ 1000}else{ Inf})
    
    if( OPT$convergence > 0 ) {
        return(NULL)
    } else {
        ##
        ## Append new levels to fit object
        ##
        EXPANDED <- Add.New.Ranefs( new.vector=OPT$par, Fit.Extract=Cur.Param, Missing=Missing$Levels )
        return(EXPANDED)
    }
}



Fit.Function <- function(idx=1, List) {
    source("100.common-variables.r")
    source("101.common-functions.r")

    source("300.variables.r")
    source("301.functions.r")

    PATHS.LIST <- Create.Folders( List[[idx]]$subset )

    MODEL.SET <- Find.Models.To.Fit( PATHS.LIST )

    HOLDER <- Load.Subset.Wrapper( Tag=List[[idx]]$subset, LSubset=TRUE ) ## NOTE: LModel=FALSE, since we have not yet determined the 'best' model

    HOLDER$MODEL <- readRDS( file.path( PATHS.LIST$MODEL, List[[idx]]$model ) )
    
    EXTRACT <- Extract.Wrapper( HOLDER, Store.Full=TRUE ) ## store full fitting object to use as initial point of bfpNA() re-fit

    Save.Extracted( EXTRACT, PATHS.LIST, List[[idx]]$model, Save.Full=STORE.OPT$Fit.Full )
    
    if( !is.null(EXTRACT$param) ) {
        if( (HOLDER$MODEL$inc.fp ) || (grepl(".fp.",List[[idx]]$model,fixed=TRUE)) ) {
            lbfp <- sub("\\.fp\\.",".bfpNA.",List[[idx]]$model)
            cat( List[[idx]]$subset, lbfp, "start ...\n" )
            
            HOLDER$MODEL <- Make.bfpNA.model.from.extract( EXTRACT$param )
            
            saveRDS( HOLDER$MODEL, file.path( PATHS.LIST$MODEL, lbfp ) )
            
            EXTRACT.bfp <- Extract.Wrapper( HOLDER, Fit.Full=STORE.OPT$Fit.Full, start.from=EXTRACT$FIT ) ## helpful to start.from, improves convergence speed

            Save.Extracted( EXTRACT.bfp, PATHS.LIST, lbfp, Save.Full=STORE.OPT$Fit.Full )
        }
    }

    cat( List[[idx]]$subset, List[[idx]]$model, "done", "\n" )
}
