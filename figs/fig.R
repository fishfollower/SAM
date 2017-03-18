library(stockassessment)
source("../stockassessment/tests/nsher/script.R", chdir = TRUE)

png("fbar.png", 600, 600)
fbarplot(fit, partial=FALSE)
dev.off()
png("ssb.png", 600, 600)
ssbplot(fit)
dev.off()
png("rec.png", 600, 600)
recplot(fit)
dev.off()
png("catch.png", 600, 600)
catchplot(fit)
dev.off()

res<-residuals(fit)

png("res.png", 600, 600)
plot(res, bubble=.5)
dev.off()

resp<-procres(fit)
png("resp.png", 600, 600)
plot(resp, bubble=.5)
dev.off()

retro <- retro(fit,year=10)

png("retro.png", 600, 1200)
plot(retro)
dev.off()

lo<-leaveout(fit)

png("lo.png", 600, 1200)
plot(lo)
dev.off()
