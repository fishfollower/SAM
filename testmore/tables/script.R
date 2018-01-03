set.seed(234079)

library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"AR"
nscodConf$keyCorObs[3,]<-c(0,1,1,-1,-1)

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

## forec <- forecast(fit,fval=c(.1,.2,.3,4))

## resid <- residuals(fit)

## resp <- procres(fit)

retro <- retro(fit,year=3,ncores=1)

td <- tempdir()
f <- file.path(td,"SAMtable_%03d.tab")

txt <- c(
    ## Calibration talbe - if this changes something else is wrong
    'matrix(1:9,3,3)',
    ## sam plots
    'ssbtable(fit)',
    'tsbtable(fit)',
    'fbartable(fit)',
    'rectable(fit)',
    'catchtable(fit)',
    'ntable(fit)',
    'faytable(fit)',
    'partable(fit)',
    'modeltable(fit)',
    'ypr(fit)',
    ## samset plots
    'modeltable(retro)'
     )

for(i in 1:length(txt)){
    tt <- capture.output(print(eval(parse(text = txt[i]))))
    con <- file(sprintf(f,i),"w") 
    writeLines(tt,con,sep=ifelse( Sys.info()["sysname"]=="Windows","\n","\r\n"))
    close(con)
}

lf <- list.files(td,pattern="^SAMtable_[[:digit:]]{3}.tab$",full.names=TRUE)

cat(tools::md5sum(lf),"\n",sep="\n", file="res.out")
