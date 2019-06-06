set.seed(234079)

library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"AR"
nscodConf$keyCorObs[3,]<-c(0,1,1,-1,-1)

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

forec <- forecast(fit,fval=c(.1,.2,.3,4))

resid <- residuals(fit)

resp <- procres(fit)

retro <- retro(fit,year=3,ncores=1)

td <- tempdir()
f <- file.path(td,"SAMplot_%03d.pdf")

txt <- c(
    ## Calibration plot - if this changes something else is wrong
    'plot(1:10,col=c("black","red",rgb(0.5,0.9,0.1,0.700000)))',
    ## sam plots
    'plot(fit)',
    'corplot(fit)',
    'fbarplot(fit)',
    'ssbplot(fit)',
    'tsbplot(fit)',
    'recplot(fit)',
    'catchplot(fit)',
    'parplot(fit)',
    'srplot(fit)',
    'fitplot(fit)',
    'dataplot(fit)',
    ## samset plots
    'fbarplot(retro)',
    'ssbplot(retro)',
    'tsbplot(retro)',
    'recplot(retro)',
    'catchplot(retro)',
    'parplot(retro)',
    ## samforecast plots
    'plot(forec)',
    'fbarplot(forec)',
    'ssbplot(forec)',
    'tsbplot(forec)',
    'recplot(forec)',
    'catchplot(forec)',
    ## samres plots
    'plot(resid)',
    'corplot(resid)',
    ## Process residuals
    'plot(resp)'
     )

for(i in 1:length(txt)){
    pdf(sprintf(f,i),
        encoding="ISOLatin1",
        useDingbats=FALSE,
        compress=FALSE,
        onefile=TRUE,
        paper="a4",
        colormodel="cmyk",
        useKerning=FALSE,
        version="1.4")
    eval(parse(text = txt[i]))
    dev.off()
}

lf <- list.files(td,pattern="^SAMplot_[[:digit:]]{3}.pdf$",full.names=TRUE)

## Remove some meta data and change to windows line endings
invisible(sapply(lf, function(ff){
    con <- file(ff,"r")
    a <- readLines(ff)
    close(con)
    b<-a[-2][!grepl("^(/.+Date|/Producer)",a[-2])]
    ## Fix rounding on alpha levels
    indx <- grepl("/(CA|ca) [[:digit:]]\\.[[:digit:]]{1,}",b)
    b[indx] <- unlist(lapply(strsplit(b[indx],"[[:space:]]"),function(x){
        paste(c(x[1],formatC(as.numeric(x[2]),1,format="f"),x[-(1:2)]),collapse=" ")
    }))
    ## Write new pdf
    con <- file(ff,"w") 
    writeLines(b,con,sep=ifelse( Sys.info()["sysname"]=="Windows","\n","\r\n"))
    close(con)
}))

res <- cbind(txt,tools::md5sum(lf))
oo <- options(width=1000)
cat(capture.output(prmatrix(res, rowlab = rep("",nrow(res)), collab = c("",""), quote=FALSE)), sep = "\n", file="res.out")
options(oo)
