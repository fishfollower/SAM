

a <- read.csv("D:/Downloads/NSASfleet.csv",header=T,stringsAsFactors=F)
Afleet <- subset(a,Fleet=="A")
Bfleet <- subset(a,Fleet=="B")
Cfleet <- subset(a,Fleet=="C")
Dfleet <- subset(a,Fleet=="D")
CA <- CB <- CC <- CD <- t(matrix(NA,nrow=length(1997:2016),ncol=c(10),dimnames=list(year=1997:2016,age=0:9)))
CWA <- CWB <- CWC <- CWD <- t(matrix(NA,nrow=length(1997:2016),ncol=c(10),dimnames=list(year=1997:2016,age=0:9)))
CA[] <- as.numeric(Afleet$numbers)*1000
CB[] <- as.numeric(Bfleet$numbers)*1000
CC[] <- as.numeric(Cfleet$numbers)*1000
CD[] <- as.numeric(Dfleet$numbers)*1000

CWA[] <- as.numeric(Afleet$weight)
CWB[] <- as.numeric(Bfleet$weight)
CWC[] <- as.numeric(Cfleet$weight)
CWD[] <- as.numeric(Dfleet$weight)

write.table(t(CA),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cn_A.dat",row.names=F,col.names=F)
write.table(t(CB),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cn_B.dat",row.names=F,col.names=F)
write.table(t(CC),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cn_C.dat",row.names=F,col.names=F)
write.table(t(CD),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cn_D.dat",row.names=F,col.names=F)

write.table(t(CWA),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cw_A.dat",row.names=F,col.names=F)
write.table(t(CWB),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cw_B.dat",row.names=F,col.names=F)
write.table(t(CWC),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cw_C.dat",row.names=F,col.names=F)
write.table(t(CWD),file="D:/Repository/SAM/SAMMULTI/testmore/NSAS/split/cw_D.dat",row.names=F,col.names=F)