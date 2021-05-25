rm(list=ls())
## COMMON LIBRARIES AND FUNCTIONS
source("100.common-variables.r")
source("101.common-functions.r")

source("300.variables.r")
source("301.functions.r")

source("500.plotting-variables.r")
source("501.plotting-functions.r")


## SCRIPT SPECIFIC LIBRARIES

## SCRIPT SPECIFIC FUNCTIONS

## SCRIPT CODE
##
##
if( 1 ){
    Print.Disclaimer( )
    
    for( wTAG in DERIVED.SET ) {

        cat( "=====", wTAG, "=====\n" )
        
        ##
        ## readRDS all the elements within the subset-folder
        ##
        PATHS.LIST <- Create.Folders( wTAG )

        if( !file.exists( file.path( PATHS.LIST$PATH, "DERIVED.rds" ) ) ) {
            cat(sprintf("Warning: %s does not have a DERIVED.rds object, run 350.calc-derived script.",wTAG),"\n")
            next
        }
        PRIMARY <- readRDS( file=file.path( PATHS.LIST$PATH, "DERIVED.rds" ) )
        
        


        ##
        ## Some validation of a simulated dataset (only in a very specific case where we can exactly compute a comparable variance)
        ##
        if( (PRIMARY$MODEL$family=="NO") && (!is.null(attr(PRIMARY$SUBSET,"spec"))) ) {
            ##
            ## what follows only applies under Gaussian/Normal outcome distribution!
            
            ##
            ## USES CURVE OBJECT FROM ABVE BLOCK
            ##
            
            ## investigate variance (within NO-family this is simpler to conceptualise)
            ##
            cat("Within NO()-class we can directly compare the 'effective' true and estimated variance.\n")            
            
            ## First: the truth from the generation script
            unlist( attr(attr(PRIMARY$SUBSET,"spec"),"re.sd") )
            ##
            ## Note: within the simulation we have study and subject(=ID) random-effects, but the model only has study
            ##       hence, the model will fit these as a sort of combined effect

            ## Under normal outcome NO(mu=0+ranef(R)+ranef(S),sigma=T) ==> baseline observations are NO(mu=0,sigma=sqrt(S^2+R^2+T^2))
            TRUE.SIGMA <- unique(PRIMARY$SUBSET[,sprintf("%s.SIGMA.fixef",PRIMARY$MODEL$covariates$Y)])
            if( length(TRUE.SIGMA)==1 ) {
                TRUE.EFFECTIVE.VAR <- sqrt( sum( exp(TRUE.SIGMA)^2, unlist( attr(attr(PRIMARY$SUBSET,"spec"),"re.sd") )^2 ) )
                ## Within NO-fit, we have the same, except only two variance terms: NO-sigma and mu-ranef
                EST.EFFECTIVE.VAR <- sqrt( exp( PRIMARY$FIT.EXTRACT$param[["sigma"]]$fixef["(Intercept)"] )^2 + (PRIMARY$FIT.EXTRACT$param[["mu"]]$sigma)^2 )

                cat("NOTE: This is efective variance, assuming a 'sigma ~ 1' and hence, Effective = NO-sigma + ranef(ID) + ranef(Study)\n")
                cat("      Within Effective there are three parts, but the estimated model only includes two parts (which partition the variance).\n")
                print(c("True Var"=unname(TRUE.EFFECTIVE.VAR),"Est Var"=unname(EST.EFFECTIVE.VAR)))
            } else {
                cat("Model includes varying sigma-component, so we cannot calculate a direct comparator (ie. the effective variance).\n")
            }
        }


        ##
        ## Plot DATA
        ##
        if( 1 ) {
            cat("Plot DATA...")
            
            XLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$X],na.rm=TRUE)
            YLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$Y],na.rm=TRUE)
            COLOURS.STUDY <- rainbow( nlevels( PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF] ) )
            COLOURS.TYPE <- c("black","red")
            
            png(filename=file.path(PATHS.LIST$PATH,"DATA.baseline.png"),width=PNG.OPT$width,height=PNG.OPT$height)
            layout(matrix(1:2,nrow=1))
            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y )
            title(main=wTAG,line=3)
            title(main=sprintf("data (baseline, control only), all data studies%s",""),line=1,font.main=1)        
            points(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                   data=PRIMARY$DATA[ with(PRIMARY$DATA, (INDEX.OB==1) & (INDEX.TYPE==levels(INDEX.TYPE)[1]) ), ],
                   col=COLOURS.STUDY[ get(PRIMARY$MODEL$covariates$RANEF[1]) ], pch=19 )
            LevLIM <- min( 8, nlevels(PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF[1]]) )
            legend("topright",legend=levels(PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF[1]])[1:LevLIM],col=COLOURS.STUDY[1:LevLIM],pch=19,bg="white",cex=0.6)

            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y )
            title(main=wTAG,line=3)
            title(main=sprintf("data (baseline, control Vs other), all data studies%s",""),line=1,font.main=1)
            points(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                   data=PRIMARY$DATA[ with(PRIMARY$DATA, (INDEX.OB==1) ), ], col=COLOURS.TYPE[INDEX.TYPE], pch=19 )
            legend("topright",legend=levels(PRIMARY$DATA[,"INDEX.TYPE"]),col=COLOURS.TYPE[1:nlevels(PRIMARY$DATA[,"INDEX.TYPE"])],pch=19,bg="white",cex=1)
            dev.off()
            

            cat(" done\n")
        }

        ##
        ## Plot longitudinal DATA
        ##
        if( 1 ) {

            XLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$X],na.rm=TRUE)
            YLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$Y],na.rm=TRUE)
            COLOURS.TYPE <- c("black","red")
            
            png(filename=file.path(PATHS.LIST$PATH,sprintf("DATA.longitudinal-sample%s.original.png","")),width=PNG.OPT$width,height=PNG.OPT$height)
            layout(1)
            ##
            ## Below: How to create an axis displaying the original scale (ie. don't change the data, change the axis labels)
            ##
            AXIS.AT <- pretty( PRIMARY$DATA[, PRIMARY$MODEL$covariates$X] )
            if( "X" %in% names(attr(PRIMARY$SUBSET,"Transformations")) ) {
                XLAB <- attr(PRIMARY$SUBSET,"Transformations")[["X"]]$OriginalName
                AXIS.LAB <- attr(PRIMARY$SUBSET,"Transformations")[["X"]]$toOriginal( AXIS.AT )
            } else {
                XLAB <- PRIMARY$MODEL$covariates$X
                AXIS.LAB <- AXIS.AT
            }

            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=XLAB, ylab=PRIMARY$MODEL$covariates$Y, xaxt="n" )
            axis(side=1,at=AXIS.AT,labels=AXIS.LAB)
            title(main=wTAG,line=3)
            title(main=sprintf("Sample of longitudinal data within DATA (with 'long' follow-up)%s",""),line=1,font.main=1)

            Span.max <- max(PRIMARY$LONG[,sprintf("%s.rng",PRIMARY$MODEL$covariates$X)])
            tWHO <- PRIMARY$LONG[ ( PRIMARY$LONG[,sprintf("%s.rng",PRIMARY$MODEL$covariates$X)] > Span.max/5 ), "INDEX.ID" ]
            for( lID in tWHO ) {
                lines(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                      data=na.omit(PRIMARY$DATA[ (PRIMARY$DATA$INDEX.ID==lID), c("INDEX.TYPE",unlist(PRIMARY$MODEL$covariates[c("Y","X")])) ]),
                      col=COLOURS.TYPE[unique(INDEX.TYPE)] )
            }
            legend("topleft",levels(PRIMARY$DATA$INDEX.TYPE),lty=1,col=COLOURS.TYPE)
            dev.off()

        }
        ##
        ## Plot longitudinal DATA
        ##
        if( 1 ) {

            png(filename=file.path(PATHS.LIST$PATH,sprintf("BOXPLOT.longitudinal-q-iqr%s.png","")),width=PNG.OPT$width,height=PNG.OPT$height)
            layout(1)
            PRIMARY$LONG[,"INDEX.TYPE"] <- factor((PRIMARY$LONG[,sprintf("%s.first",PRIMARY$MODEL$covariates$COND)] == levels(PRIMARY$DATA[,"INDEX.TYPE"])[1] ), c(TRUE,FALSE), levels(PRIMARY$DATA[,"INDEX.TYPE"]) )
            oldpar <- par(mar=c(10,5,4,1))
            boxplot(as.formula(sprintf(" %s.q.iqr ~ INDEX.TYPE + %s",PRIMARY$MODEL$covariates$Y,PRIMARY$MODEL$covariates$RANEF[1])),
                    data=droplevels(PRIMARY$LONG), las=2 )
            par(oldpar)
            dev.off()
        }
        ##
        ## Plot 
        ##
        if( 1 ) {
            COLOURS.TYPE <- c("white","red")
            
            png(filename=file.path(PATHS.LIST$PATH,sprintf("BOXPLOT.quantiles%s.png","")),width=PNG.OPT$width,height=PNG.OPT$height)

            boxplot(as.formula(sprintf(" %s.q.wre ~ INDEX.TYPE + %s",PRIMARY$MODEL$covariates$Y,PRIMARY$MODEL$covariates$RANEF[1])),
                    data=droplevels(PRIMARY$DATA.PRED[ PRIMARY$DATA.PRED[,"INDEX.OB"]==1, ]), las=2, col=COLOURS.TYPE )
            legend("topright",levels(PRIMARY$DATA[,"INDEX.TYPE"]),fill=COLOURS.TYPE,bg="white")
            dev.off()
        }
        

        ##
        ## Plot SUBSET with fitted population curves
        ##        
        if( 1 ) {
            cat("Plot SUBSET with POPULATION CURVE...")

            XLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$X],na.rm=TRUE)
            YLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$Y],na.rm=TRUE)
            COLOURS.STUDY <- rainbow( nlevels( PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF] ) )
            COLOURS.TYPE <- c("black","red")

            COLOURS.FIT <- "purple"

            QUANT.SET <- c("l025"=2,"m500"=1,"u975"=2)
            WHICH <- Reduce(f=`&`,x=lapply(2:length(PRIMARY$POP.CURVE.LIST),
                                           function(X){ PRIMARY$POP.CURVE.PRED[,names(PRIMARY$POP.CURVE.LIST)[X]]==PRIMARY$POP.CURVE.LIST[[X]][1]  } ))
            WHICH.LABEL <- paste(sapply( 2:length(PRIMARY$POP.CURVE.LIST),
                                        function(X){ sprintf("%s = %s",names(PRIMARY$POP.CURVE.LIST)[X],PRIMARY$POP.CURVE.LIST[[X]][1])  } ),
                                 collapse=", ")
            
            png(filename=file.path(PATHS.LIST$PATH,"PLOT.fit.png"),width=PNG.OPT$width,height=PNG.OPT$height)
            layout(matrix(1:3,nrow=1))
            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y ) ## add log="x" ?
            title(main=wTAG,line=3)
            title(main=sprintf("fit + subset(baseline, control only), all subset studies [which: %s]",WHICH.LABEL),line=1,font.main=1)            
            points(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                   data=PRIMARY$SUBSET,
                   col=alpha(COLOURS.STUDY,0.25)[ get(PRIMARY$MODEL$covariates$RANEF[1]) ], pch=19 )

            lines(as.formula(sprintf("PRED.m500.pop ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,],
                  lwd=4, lty=1, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.l025.pop.lower ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,],
                  lwd=4, lty=2, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.u975.pop.upper ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,],
                  lwd=4, lty=2, col=COLOURS.FIT )
            legend("topright",c("Predicted median (50th-quantile)","Predicted 95th-quantile","Predicted 5th-quantile"),
                   col=COLOURS.FIT,lty=c(1,2,2),lwd=2,title="Population",bg="white")
            
            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y ) ## add log="x"?
            title(main=wTAG,line=3)
            title(main=sprintf("fit + data (baseline, other only), all data studies [Which: %s]",WHICH.LABEL),line=1,font.main=1)            
            points(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                   data=PRIMARY$DATA[ with(PRIMARY$DATA, (INDEX.OB==1) & (INDEX.TYPE==levels(INDEX.TYPE)[2]) ), ], col=alpha("red",0.25), pch=19 )

            lines(as.formula(sprintf("PRED.m500.pop ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=1, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.l025.pop.lower ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=2, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.u975.pop.upper ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=2, col=COLOURS.FIT )

            legend("topright",c("Predicted median (50th-quantile)","Predicted 95th-quantile","Predicted 5th-quantile"),
                   col=COLOURS.FIT,lty=c(1,2,2),lwd=2,title="Population",bg="white")


            plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y ) ## add , log="x"?
            title(main=wTAG,line=3)
            title(main=sprintf("fit + normalised-data (baseline, other only), all data studies [Which: %s]",WHICH.LABEL),line=1,font.main=1)            
            points(as.formula(sprintf("%s.normalised ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                   data=PRIMARY$DATA.PRED[ with(PRIMARY$DATA.PRED, (INDEX.OB==1) & (INDEX.TYPE==levels(INDEX.TYPE)[2]) ), ], col=alpha("red",0.25), pch=19 )

            lines(as.formula(sprintf("PRED.m500.pop ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=1, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.l025.pop.lower ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=2, col=COLOURS.FIT )
            lines(as.formula(sprintf("PRED.u975.pop.upper ~ %s",PRIMARY$MODEL$covariates$X)),
                  data=PRIMARY$POP.CURVE.PRED[WHICH,], lwd=4, lty=2, col=COLOURS.FIT )

            
            dev.off()

            cat(" done\n")
        }


        ##
        ## Plot moments and derivatives
        ##        
        if( 1 ) {
            cat("Plot moments...")
            
            XLIM <- range(PRIMARY$POP.CURVE.PRED[,PRIMARY$MODEL$covariates$X],na.rm=TRUE)
            YLIM <- list()
            for( lMOMENT in c("mean","variance") ) {
                YLIM[[lMOMENT]] <- list()
                for( lD in c("",".D") ) {
                    LAB <- sprintf("PRED.%s.pop%s",lMOMENT,lD)
                    YLIM[[lMOMENT]][[LAB]] <- range(PRIMARY$POP.CURVE.PRED[,sprintf(paste0(LAB,"%s"),c("",".lower",".upper"))],
                                                    na.rm=TRUE)
                }
            }

            BY.COLOURS <- setNames( c("black","red","blue","purple")[1:length(PRIMARY$POP.CURVE.LIST[[2]])], PRIMARY$POP.CURVE.LIST[[2]] )

            if( length(PRIMARY$POP.CURVE.LIST)>2 ) {
                LOGICAL <- Reduce(f=`&`,x=lapply(3:length(PRIMARY$POP.CURVE.LIST),
                                               function(X){ PRIMARY$POP.CURVE.PRED[,names(PRIMARY$POP.CURVE.LIST)[X]]==PRIMARY$POP.CURVE.LIST[[X]][1]  } ))
                LOGICAL.LABEL <- paste(sapply( 3:length(PRIMARY$POP.CURVE.LIST),
                                              function(X){ sprintf("%s = %s",names(PRIMARY$POP.CURVE.LIST)[X],PRIMARY$POP.CURVE.LIST[[X]][1])  } ),
                                     collapse=", ")
            } else {
                LOGICAL <- rep(TRUE,NROW(PRIMARY$POP.CURVE.PRED))
                LOGICAL.LABEL <- ""
            }
            
            png(filename=file.path(PATHS.LIST$PATH,"PLOT.moments.png"),width=PNG.OPT$width,height=PNG.OPT$height)
            layout(matrix(1:(length(YLIM)*2),nrow=2))
            for( lMOMENT in c("mean","variance") ) {
                for( lD in c("",".D") ) {
                    LAB <- sprintf("PRED.%s.pop%s",lMOMENT,lD)            
                    plot(0, type="n", xlim=XLIM, ylim=YLIM[[lMOMENT]][[LAB]], xlab=PRIMARY$MODEL$covariates$X, ylab=sprintf("%s%s",lMOMENT,lD) ) ## add , log="x"
                    abline(h=0,col="grey")
                    title(main=wTAG,line=3)
                    title(main=sprintf("k-th moment of fit [Which: %s]",LOGICAL.LABEL),line=1,font.main=1)            

                    for( lBY in PRIMARY$POP.CURVE.LIST[[2]] ) {## loop over the second element of the BY variable
                        lWHICH <- LOGICAL & ( PRIMARY$POP.CURVE.PRED[,names(PRIMARY$POP.CURVE.LIST)[2]]==lBY )
                        lines(as.formula(sprintf("PRED.%s.pop%s ~ %s",lMOMENT,lD,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$POP.CURVE.PRED[lWHICH,],
                              lwd=2, lty=1, col=BY.COLOURS[lBY] )
                        lines(as.formula(sprintf("PRED.%s.pop%s.lower ~ %s",lMOMENT,lD,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$POP.CURVE.PRED[lWHICH,],
                              lwd=2, lty=2, col=BY.COLOURS[lBY] )
                        lines(as.formula(sprintf("PRED.%s.pop%s.upper ~ %s",lMOMENT,lD,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$POP.CURVE.PRED[lWHICH,],
                              lwd=2, lty=2, col=BY.COLOURS[lBY] )
                    }
                }
            }
            oldpar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
            legend("center",names(BY.COLOURS),col=BY.COLOURS,lty=1,title=names(PRIMARY$POP.CURVE.LIST)[2])
            par(oldpar)
            dev.off()

            cat(" done\n")
        }


        ##
        ## Plot SUBSET with fitted population curves
        ##        
        if( 1 ) {
            cat("Plot SUBSET with STUDY CURVE...")

            COLOURS.TYPE <- c("black","red")
            COLOURS.FIT <- c("purple","blue")

            for( IDX in 1:length(PRIMARY$STUDY.CURVES) ) {
                LAB <- names(PRIMARY$STUDY.CURVES)[IDX]
                
                WHICH <- Reduce(f=`&`,x=lapply(3:length(PRIMARY$STUDY.CURVES[[LAB]]$LIST),
                                               function(X){ PRIMARY$STUDY.CURVES[[LAB]]$PRED[,names(PRIMARY$STUDY.CURVES[[LAB]]$LIST)[X]]==PRIMARY$STUDY.CURVES[[LAB]]$LIST[[X]][1]  } ))
                WHICH.LABEL <- paste(sapply( 3:length(PRIMARY$STUDY.CURVES[[LAB]]$LIST),
                                            function(X){ sprintf("%s = %s",names(PRIMARY$STUDY.CURVES[[LAB]]$LIST)[X],PRIMARY$STUDY.CURVES[[LAB]]$LIST[[X]][1])  } ),
                                     collapse=", ")
                
                png(filename=file.path(PATHS.LIST$PATH,sprintf("PLOT.study-curve-%03i.png",IDX)),width=PNG.OPT$width,height=PNG.OPT$height)
                layout(1)
                XLIM <- range( PRIMARY$DATA[ which(PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF]==LAB), PRIMARY$MODEL$covariates$X], na.rm=TRUE )
                YLIM <- range( PRIMARY$DATA[ which(PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF]==LAB), PRIMARY$MODEL$covariates$Y], na.rm=TRUE )
                YLIM <- range( PRIMARY$DATA[ , PRIMARY$MODEL$covariates$Y], na.rm=TRUE )
                plot(0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X, ylab=PRIMARY$MODEL$covariates$Y ) ## add , log="x"
                title(main=wTAG,line=3)
                title(main=sprintf("fit + data(baseline, control only), within %s studies [which: %s]",LAB,WHICH.LABEL),line=1,font.main=1)            
                points(as.formula(sprintf("%s ~ %s",PRIMARY$MODEL$covariates$Y, PRIMARY$MODEL$covariates$X)),
                       data=PRIMARY$DATA[ (PRIMARY$DATA[,PRIMARY$MODEL$covariates$RANEF]==LAB) & (PRIMARY$DATA[,"INDEX.OB"]==1),],
                       col=alpha(COLOURS.TYPE,0.25)[ get("INDEX.TYPE") ], pch=19 )

                for( JDX in 1:2 ) {
                    lre <- c("pop","wre")[JDX]
                    lines(as.formula(sprintf("PRED.m500.%s ~ %s",lre,PRIMARY$MODEL$covariates$X)),
                          data=PRIMARY$STUDY.CURVES[[LAB]]$PRED[WHICH,],
                          lwd=4, lty=1, col=COLOURS.FIT[JDX] )
                    if( 0 ) {
                        lines(as.formula(sprintf("PRED.l025.%s.lower ~ %s",lre,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$STUDY.CURVES[[LAB]]$PRED[WHICH,],
                              lwd=4, lty=2, col=COLOURS.FIT[JDX] )
                        lines(as.formula(sprintf("PRED.u975.%s.upper ~ %s",lre,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$STUDY.CURVES[[LAB]]$PRED[WHICH,],
                              lwd=4, lty=2, col=COLOURS.FIT[JDX] )
                    }
                    if( 1 ) {
                        lines(as.formula(sprintf("PRED.m500.%s.lower ~ %s",lre,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$STUDY.CURVES[[LAB]]$PRED[WHICH,],
                              lwd=4, lty=2, col=COLOURS.FIT[JDX] )
                        lines(as.formula(sprintf("PRED.m500.%s.upper ~ %s",lre,PRIMARY$MODEL$covariates$X)),
                              data=PRIMARY$STUDY.CURVES[[LAB]]$PRED[WHICH,],
                              lwd=4, lty=2, col=COLOURS.FIT[JDX] )
                    }                    
                }
                legend("topright",c("Population curve","Study-specific curve"),
                       col=COLOURS.FIT,lty=c(1,1),lwd=2,bg="white")

                dev.off()
            }            
            


            cat(" done\n")
        }
        
        
        ##
        ## Compare fixed-effect curve and random-effect estimates (ONLY FOR SIMULATION SETTING WITH KNOWN TRUTH)
        ##
        if( !is.null( attr(PRIMARY$SUBSET,"spec") ) ) {
            cat("Plot random-effects...")

            for( LAB in names(PRIMARY$FAMILY$parameters) ) {
                ##
                ## Dealing with a simulated subset (with a known truth)
                TRUE.FUNC <- attributes(attr(PRIMARY$SUBSET,"spec"))$truth[[PRIMARY$MODEL$covariates$Y]][[ toupper(LAB) ]]$TypeBase
                PRIMARY$POP.CURVE.PRED[,sprintf("truth.%s.pop",LAB)] <- TRUE.FUNC( PRIMARY$POP.CURVE.PRED[,PRIMARY$MODEL$covariates$X] )
            }


            for( wPARAM in names(PRIMARY$FAMILY$parameters) ) {

                local.columns <- grep(sprintf("%s.pop",wPARAM),names(PRIMARY$POP.CURVE.PRED),fixed=TRUE,value=TRUE)
                YLIM <- range( PRIMARY$POP.CURVE.PRED[,local.columns], na.rm=TRUE )
                XLIM <- range(PRIMARY$DATA[,PRIMARY$MODEL$covariates$X],na.rm=TRUE)
                
                tLAB <- sprintf("%s.ranef.Study",toupper(wPARAM))
                if( tLAB %in% names(attr(PRIMARY$SUBSET,"spec")) ) {
                    RANEF.TRUTH <- attr(PRIMARY$SUBSET,"spec")[,tLAB,drop=FALSE]
                    RANEF.GAMLSS <- PRIMARY$FIT.EXTRACT$param[[wPARAM]]$ranef
                    RLIM <- range( c(RANEF.GAMLSS, RANEF.TRUTH[ names(RANEF.GAMLSS), ]) )

                    ##RANEF.BOOT <- apply(sapply(PRIMARY$BOOT.EXTRACT, function(X){X$param[[wPARAM]]$ranef} ),1,quantile,probs=c(0.1,0.9))
                } else {
                    RANEF.TRUTH <- numeric(0)
                }

                cWHICH <- Reduce(f=`&`,x=lapply(2:length(PRIMARY$POP.CURVE.LIST),
                                                function(X){ PRIMARY$POP.CURVE.PRED[,names(PRIMARY$POP.CURVE.LIST)[X]]==PRIMARY$POP.CURVE.LIST[[X]][1]  } ))
                cWHICH.LABEL <- paste(sapply( 2:length(PRIMARY$POP.CURVE.LIST),
                                             function(X){ sprintf("%s = %s",names(PRIMARY$POP.CURVE.LIST)[X],PRIMARY$POP.CURVE.LIST[[X]][1])  } ),
                                     collapse=", ")
                
                png(filename=file.path(PATHS.LIST$PATH,sprintf("PLOT.%s-component.png",wPARAM)),width=PNG.OPT$width,height=PNG.OPT$height)
                layout(matrix(c(1,2, 4,3),nrow=2,byrow=TRUE),heights=c(3,2))
                plot( 0, type="n", xlim=XLIM, ylim=YLIM, xlab=PRIMARY$MODEL$covariates$X,
                     ylab=PRIMARY$MODEL$covariates$Y, main=sprintf("Comparison of true and fitted [Which: %s] %s-component",cWHICH.LABEL,wPARAM) )
                lines( as.formula(sprintf("%s.pop ~ %s",wPARAM,PRIMARY$MODEL$covariates$X)), data=PRIMARY$POP.CURVE.PRED[cWHICH,], col="black" )
                lines( as.formula(sprintf("%s.pop.upper ~ %s",wPARAM,PRIMARY$MODEL$covariates$X)), data=PRIMARY$POP.CURVE.PRED[cWHICH,], col="black", lty=2 )
                lines( as.formula(sprintf("%s.pop.lower ~ %s",wPARAM,PRIMARY$MODEL$covariates$X)), data=PRIMARY$POP.CURVE.PRED[cWHICH,], col="black", lty=2 )
                lines( as.formula(sprintf("truth.%s.pop ~ %s",wPARAM,PRIMARY$MODEL$covariates$X)), data=PRIMARY$POP.CURVE.PRED[cWHICH,], col="red" )
                legend("topleft",c("best-fit","true"), col=c("black","red"), lty=1, bg="white" )
                
                if( length(RANEF.TRUTH) > 0 ) {
                    plot( 0, type="n", xlab="Truth", ylab="Best", xlim=RLIM, ylim=RLIM, main=sprintf("Compare %s-component study-level ranef",wPARAM) )
                    text(y=RANEF.GAMLSS, x=RANEF.TRUTH[ names(RANEF.GAMLSS), ], labels=names(RANEF.GAMLSS) )
                    abline(a=0,b=1,col="grey")

                    plot( 0, type="n", xlab="Truth", ylab="Best", xlim=RLIM, ylim=RLIM, main=sprintf("Compare %s-component study-level ranef (with bootstrap)",wPARAM) )
                    points(y=RANEF.GAMLSS, x=RANEF.TRUTH[ names(RANEF.GAMLSS), ], pch=19 )
                    abline(a=0,b=1,col="grey")
                    W.B.null <- sapply(PRIMARY$BOOT.EXTRACT,function(X){!is.null(X$param)})
                    E.B <- sapply(PRIMARY$BOOT.EXTRACT[W.B.null],function(X){X$param[[wPARAM]]$ranef})
                    if( !is.null(dim(E.B)) ) {
                        Q.B <- apply(E.B,1,quantile,prob=c(0.05,0.95))
                        arrows(x0=RANEF.TRUTH[ names(RANEF.GAMLSS), ], y0=Q.B[1,], y1=Q.B[2,], code=3, length=0.1, angle=90 )
                    }
                } else {
                    plot.new()
                    plot.new()
                }

                plot( 0, type="n", xlim=XLIM, ylim=RLIM, xlab=PRIMARY$MODEL$covariates$X,
                     ylab=PRIMARY$MODEL$covariates$RANEF[1], main=sprintf("Comparison of fitted %s-component by x-variable",wPARAM) )
                abline(h=0,col="grey")
                abline(v=min(PRIMARY$SUBSET[,PRIMARY$MODEL$covariates$X],na.rm=TRUE),col="grey",lty=2)
                abline(v=max(PRIMARY$SUBSET[,PRIMARY$MODEL$covariates$X],na.rm=TRUE),col="grey",lty=2)
                for( LAB in names(PRIMARY$FIT.EXTRACT$param[[wPARAM]]$ranef) ) {
                    sWHICH <- which( PRIMARY$SUBSET[, PRIMARY$MODEL$covariates$RANEF[1] ] == LAB )
                    X.pt <- min( PRIMARY$SUBSET[ sWHICH, PRIMARY$MODEL$covariates$X ], na.rm=TRUE )
                    points(x=X.pt,y=RANEF.GAMLSS[LAB],pch=19)
                    if( length(RANEF.TRUTH) > 0 ) {
                        points(x=X.pt,y=RANEF.TRUTH[LAB,1],pch=19,col="red")
                    }
                    if( exists("Q.B") ) {
                        arrows(x0=X.pt, y0=Q.B[1,LAB], y1=Q.B[2,LAB], code=3, length=0.1, angle=90 )
                    }
                }
                legend("topleft",c("Minimum","Maximum"),lty=2,col="grey",title="Subset",bg="white")
                legend("topright",c("Estimate","Truth (if known)"),pch=19,col=c("black","red"),bg="white")
                dev.off()


            }
            cat(" done\n")
        }


        ##
        ## Plot CLONE study estimates
        ##
        NOVEL.REFITS <- list.dirs( PATHS.LIST$NOVEL, full.names=FALSE )
        CLONE.LIST <- unique(sub("(.*)-CLONE.*","\\1",grep("-CLONE",NOVEL.REFITS,value=TRUE)))
        if( length(CLONE.LIST)>0 ) {
            cat("Plot CLONES...")
            
            With.Ranef <- sapply(PRIMARY$FIT.EXTRACT$param[names(PRIMARY$FAMILY$parameters)],function(X){!is.null(X$ranef)})
            QUANTILES <- c(0.05,0.95)
            
            for( CLONE.ORIGINAL in CLONE.LIST ) {
                CLONE.NAME <- sprintf("%s-CLONE",CLONE.ORIGINAL)
                for( lPAR in names(With.Ranef)[With.Ranef] ) {
                    CLONE.FILES <- grep(sprintf("%s-CLONE",CLONE.ORIGINAL),NOVEL.REFITS,value=TRUE)

                    CLONE.SUMMARY <- data.frame(row.names=c("TRUTH",CLONE.ORIGINAL,CLONE.FILES))
                    CLONE.SUMMARY[CLONE.ORIGINAL,"EST"] <- PRIMARY$FIT.EXTRACT$param[[lPAR]]$ranef[[CLONE.ORIGINAL]]
                    unlist(sapply( PRIMARY$BOOT.EXTRACT, function(X){X$param[[lPAR]]$ranef[[CLONE.ORIGINAL]]} ))
                    temp.list2vector <- unlist(sapply( PRIMARY$BOOT.EXTRACT, function(X){X$param[[lPAR]]$ranef[[CLONE.ORIGINAL]]} ))
                    CLONE.SUMMARY[CLONE.ORIGINAL,c("LO","UP")] <- quantile(temp.list2vector,probs=QUANTILES )

                    if( !is.null(attr(PRIMARY$SUBSET,"spec")) ) {
                        CLONE.SUMMARY["TRUTH","EST"] <- attr(PRIMARY$SUBSET,"spec")[CLONE.ORIGINAL,sprintf("%s.ranef.Study",toupper(lPAR))]
                    }

                    for( lFILE in CLONE.FILES ) {
                        
                        FIT.LOCAL  <- readRDS(file.path(PATHS.LIST$NOVEL,lFILE,"FIT.EXPANDED.rds"))
                        BOOT.LOCAL <- readRDS(file.path(PATHS.LIST$NOVEL,lFILE,"BOOT.EXPANDED.rds"))
                        CLONE.SUMMARY[lFILE,"SIZE"] <- as.numeric(sub(".*\\.n(.*)","\\1",lFILE))
                        CLONE.SUMMARY[lFILE,"EST"] <- FIT.LOCAL$param[[lPAR]]$ranef[[CLONE.NAME]]
                        temp.list2vector <- unlist(lapply( BOOT.LOCAL, function(X){X$param[[lPAR]]$ranef[[CLONE.NAME]]} ))
                        CLONE.SUMMARY[lFILE,c("LO","UP")] <- quantile(temp.list2vector,probs=QUANTILES)
                    }

                    XLIM <- range(CLONE.SUMMARY$SIZE,na.rm=TRUE)
                    YLIM <- range(CLONE.SUMMARY[,c("EST","LO","UP")],na.rm=TRUE)
                    png(filename=file.path(PATHS.LIST$PATH,sprintf("CLONE.%s-component.%s.png",lPAR,CLONE.NAME)),width=PNG.OPT$width,height=PNG.OPT$height)
                    plot( 0, type="n", xlim=XLIM, ylim=YLIM, ylab=CLONE.NAME, xlab="Cloned data sample size (log)", log="x" ) ## log(n) on the x-axis
                    title(main=sprintf("Random-effect estimate for %s-component of cloned study %s (%s)",lPAR,CLONE.ORIGINAL,CLONE.NAME))
                    abline(h=0,col="grey")
                    abline(h=CLONE.SUMMARY["TRUTH","EST"],col="red")
                    abline(h=CLONE.SUMMARY[CLONE.ORIGINAL,"EST"],col="purple")
                    abline(h=CLONE.SUMMARY[CLONE.ORIGINAL,c("LO","UP")],col="purple",lty=2)
                    points( EST ~ SIZE, data=CLONE.SUMMARY, col="black" )
                    with(CLONE.SUMMARY, arrows( x0=SIZE, y0=LO, y1=UP, code=3, length=0.1, angle=90, col="black" ) )
                    if( !is.null(attr(PRIMARY$SUBSET,"spec")) ) {
                        legend("bottomright",c("Truth","gamlss-fit","gamlss model uncertainty"),lty=c(1,1,2),col=c("red","purple","purple"),bg="white")
                    } else {
                        legend("bottomright",c("gamlss-fit","gamlss model uncertainty"),lty=c(1,2),col=c("purple","purple"),bg="white")
                    }
                    dev.off()
                }
            }
            cat(" done\n")
        }
    }
}


print( warnings() ) ## print any warnings from the code

