rm(list=ls())
source("100.common-variables.r")
source("101.common-functions.r")

source("300.variables.r")
source("301.functions.r")


## 310-script
PATHS <- Create.Folders( "omega-Wand__.n0000" )
HOLDER <- Load.Subset.Wrapper( Tag="omega-Wand__.n0000", LSubset=TRUE )
HOLDER$MODEL <- readRDS( file.path( PATHS$MODEL, "base200.GGalt.fp.rds" ) )
EXTRACT <- Extract.Wrapper( HOLDER, Store.Full=TRUE ) ## store full fitting object to use as initial point of bfpNA() re-fit [expect 12 iterations]
summary(EXTRACT$FIT.FULL) ## standard summary of a GAMLSS fit
getSmo(EXTRACT$FIT.FULL,what="mu")
getSmo(EXTRACT$FIT.FULL,what="mu")$power

HOLDER$MODEL <- Make.bfpNA.model.from.extract( EXTRACT$param )
saveRDS( HOLDER$MODEL, file.path( PATHS$MODEL, "base200.GGalt.bfpNA.rds" ) )
EXTRACT.bfp <- Extract.Wrapper( HOLDER, Fit.Full=FALSE, start.from=EXTRACT$FIT ) ## helpful to start.from, improves convergence speed [expect 5 iterations]
Save.Extracted( EXTRACT.bfp, PATHS, "base200.GGalt.bfpNA.rds", Save.Full=FALSE )

## 320-script
EXTRACT.bfp$param$BIC ## compare BIC on all fitted models
file.copy(from=file.path(PATHS$FIT.EXTRACT,"base200.GGalt.bfpNA.rds"),to=file.path(PATHS$PATH,"MODEL.rds"))
file.copy(from=file.path(PATHS$FIT.EXTRACT,"base200.GGalt.bfpNA.rds"),to=file.path(PATHS$PATH,"FIT.EXTRACT.rds"))
## or
file.symlink(from=file.path("MODEL","base200.GGalt.bfpNA.rds"),to=file.path(PATHS$PATH,"MODEL.rds"))
file.symlink(from=file.path("FIT.EXTRACT","base200.GGalt.bfpNA.rds"),to=file.path(PATHS$PATH,"FIT.EXTRACT.rds"))

## 330-script (and 340-script)
HOLDER <- Load.Subset.Wrapper( Tag="omega-Wand__.n0000", LSubset=TRUE, LModel=TRUE, LFit=TRUE )

BOOT <- list()
BOOT[[1]] <- Boot.Function(n=1,Base.Seed=12345,Holder=HOLDER)
BOOT[[2]] <- Boot.Function(n=2,Base.Seed=12345,Holder=HOLDER)
BOOT[[3]] <- Boot.Function(n=3,Base.Seed=12345,Holder=HOLDER)
for( NUM in 4:100 ) { ## 100s of bootstrap replicates required
    BOOT[[NUM]] <- Boot.Function(n=NUM,Base.Seed=12345,Holder=HOLDER)
}

Reduce(rbind,lapply(BOOT,function(X){X$param$mu$fixef}))
Reduce(rbind,lapply(BOOT,function(X){X$param$sigma$fixef}))

apply( Reduce(rbind,lapply(BOOT,function(X){X$param$mu$fixef})), 2, quantile, probs=c(0.05,0.95), na.rm=TRUE )
apply( Reduce(rbind,lapply(BOOT,function(X){X$param$sigma$fixef})), 2, quantile, probs=c(0.05,0.95), na.rm=TRUE )
apply( Reduce(rbind,lapply(BOOT,function(X){X$param$nu$fixef})), 2, quantile, probs=c(0.05,0.95), na.rm=TRUE )

saveRDS(object=BOOT,file=file.path(PATHS$PATH,"BOOT.EXTRACT.rds"))


## 350-novel-script
PRIMARY <- Load.Subset.Wrapper( Tag="omega-Wand__.n0000", LSubset=TRUE, LModel=TRUE, LFit=TRUE, LBoot=TRUE, LData=TRUE )

dim(PRIMARY$DATA)   ## Note: PRIMARY$DATA and PRIMARY$SUBSET are different,
dim(PRIMARY$SUBSET) ##       the latter contains only observations used for fitting the model
table(PRIMARY$SUBSET$Study) ## Studies U and V were not included in the orginal set

NOVEL <- list()
NOVEL$DATA <- dim(readRDS(file=file.path(PATHS$NOVEL,"Wandu_.n0200.rds")))
## or
NOVEL$DATA <- PRIMARY$DATA[ with(PRIMARY$DATA, which(Study=="U" & INDEX.OB==1 & INDEX.TYPE=="CN") ), ]

NOVEL$DATA.PRED <- Apply.Param(NEWData=NOVEL$DATA,
                               FITParam=PRIMARY$FIT.EXTRACT$param,
                               Reference.Holder=PRIMARY,
                               Pred.Set=NULL, Prefix="", Add.Moments=FALSE, Add.Normalise=FALSE, Add.Derivative=FALSE, MissingToZero=TRUE,
                               verbose=FALSE )
PRIMARY$MODEL ## in our selected model only mu has a random-effect
summary(NOVEL$DATA.PRED) ## see that mu.wre is NA, but sigma.wre and nu.wre are not (as there are no missing random-effects)

attr(NOVEL$DATA.PRED,"missing.levels") ## Apply.Param() returns information on missing random-effects

NOVEL$SUBSET <- NOVEL$DATA.PRED[attr(NOVEL$DATA.PRED,"logical.selectors")$REFIT.VALID,]
EXPANDED <- Calc.Expanded(NewData=NOVEL$SUBSET,
                          Cur.Param=PRIMARY$FIT.EXTRACT$param,
                          Missing=attr(NOVEL$DATA.PRED,"missing.levels") )

tail(data.frame(EXPANDED$mu$ranef,EXPANDED$mu$ranef.TYPE)) ## U-specific random-effects added

EXPANDED.PATH <- file.path( PATHS$NOVEL, "U" )

if( !dir.exists(EXPANDED.PATH) ) {
    dir.create(EXPANDED.PATH)
}

saveRDS(object=list(param=EXPANDED,summary=NULL),
        file=file.path(EXPANDED.PATH,"FIT.EXPANDED.rds"))



## 350-derived-script
PRIMARY <- Load.Subset.Wrapper( Tag="omega-Wand__.n0000", LSubset=TRUE, LModel=TRUE, LFit=TRUE, LBoot=TRUE, LData=TRUE )

PRIMARY$DATA.PRED <- Apply.Param(NEWData=PRIMARY$DATA, Reference.Holder=PRIMARY, FITParam=PRIMARY$FIT.EXTRACT$param,
                                 Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                 Add.Moments=FALSE, Add.Normalise=TRUE, Add.Derivative=FALSE,
                                 MissingToZero=TRUE, NAToZero=TRUE )

PRIMARY$LONG.SUMMARY <- Make.Longitudinal( Holder=PRIMARY )

head(PRIMARY$DATA.PRED)
tail(PRIMARY$DATA.PRED)

range(PRIMARY$DATA[,"TimeTransformed"]) ## whole dataset
range(PRIMARY$DATA[PRIMARY$DATA$Study=="E","TimeTransformed"]) ## only study E

PRIMARY$CURVE <- Apply.Param(NEWData=expand.grid(list(
                                 TimeTransformed=seq(0,9,length.out=2^4),
                                 Grp=c("Female","Male")
                             )),
                             FITParam=PRIMARY$FIT.EXTRACT$param )


PRIMARY$CURVE.E <- Apply.Param(NEWData=expand.grid(list(
                                 TimeTransformed=seq(0,9,length.out=2^8),
                                 Grp=c("Female","Male"),
                                 Study="E"
                             )),
                             FITParam=PRIMARY$FIT.EXTRACT$param )

RANGE <- range(PRIMARY$DATA[PRIMARY$DATA$Study=="E","TimeTransformed"])
plot( PRED.m500.pop ~ TimeTransformed, data=subset(PRIMARY$CURVE,Grp=="Female"), type="l", ylim=c(0,2.5) )
lines( PRED.m500.wre ~ TimeTransformed, data=subset(PRIMARY$CURVE.E,Grp=="Female"&TimeTransformed<RANGE[1]), col="red", lwd=2, lty=2 )
lines( PRED.m500.wre ~ TimeTransformed, data=subset(PRIMARY$CURVE.E,Grp=="Female"&TimeTransformed>RANGE[2]), col="red", lwd=2, lty=2 )
lines( PRED.m500.wre ~ TimeTransformed, data=subset(PRIMARY$CURVE.E,Grp=="Female"&TimeTransformed<RANGE[2]&TimeTransformed>RANGE[1]), col="red", lwd=4, lty=1 )
abline(v=RANGE,col="red",lty=2)
legend("bottomright",c("Population","Study E","Study E (extrapolated)"),lty=c(1,1,2),col=c("black","red","red"),title="50th Centile")

plot( PRED.variance.pop ~ TimeTransformed, data=subset(PRIMARY$CURVE,Grp=="Female"), type="l", ylim=c(0,0.05) )
lines( PRED.variance.pop ~ TimeTransformed, data=subset(PRIMARY$CURVE,Grp=="Male"), col="purple" )
legend("topleft",c("Female","Male"),lty=1,col=c("black","purple"),title="Population variance")


names(PRIMARY$LONG.SUMMARY)
tail(PRIMARY$DATA.PRED[ PRIMARY$DATA.PRED$Study=="V", ])

BP <- boxplot( Wand.q.iqr ~ Study + Type.first, data=droplevels(na.omit(PRIMARY$LONG.SUMMARY[,c("Wand.q.iqr","Study","Type.first")])) )
