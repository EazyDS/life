require(demography)
require(MortalitySmooth)
source("~/Dropbox/Research/Rpackages/nicefigs.R")

usa <- hmd.mx("USA","Rob.Hyndman@buseco.monash.edu.au","ezekiel;")
usa1950 <- extract.years(usa,years=1950:2010)


savepdf("us2003")
plot(usa1950, years=2003, series='male',type="p",pch=1)
endpdf()

smus <- smooth.demogdata(usa1950)


savepdf("smus2003")
plot(usa1950, years=2003, series='male',type="p", pch=1, col="gray")
lines(smus, years=2003, series='male')
endpdf()

# CDE method

Ext <- usa$pop$male
Dxt <- usa$rate$male * Ext
fitBIC <- Mort2Dsmooth(x=usa$age, y=usa$year, Z=Dxt, offset=log(Ext))

savepdf("cdesmooth")
par(mfrow=c(1,2))
plot(fitBIC$x,log(usa$rate$male[,"2003"]),xlab="Year",ylab="Log death rate",
     main="USA: male death rates 2003",col="gray")
lines(fitBIC$x,log(fitBIC$fitted.values[,"2003"]/Ext[,"2003"]))
lines(smus,year=2003,series='male',lty=2)
legend("topleft",lty=1:2,legend=c("CDE smoothing", "HU smoothing"))

plot(fitBIC$y,log(Dxt["65",]/Ext["65",]),xlab="Year",ylab="Log death rate",
     main="USA: male death rates age 65",col="gray")
lines(fitBIC$y,log(fitBIC$fitted.values["65",]/Ext["65",]))
lines(smus$year,log(smus$rate$male["65",]),lty=2)
legend("bottomleft",lty=1:2,legend=c("CDE smoothing", "HU smoothing"))
dev.off()

## LEE_CARTER

lc.female <- lca(usa, "female")
forecast.lc.female <- forecast(lc.female, h=20)

plot(usa, series="female")  
plot(lc.female)
plot(forecast.lc.female,"c") 
plot(forecast.lc.female)      

lcnone.female <- lca(usa,"female",adjust="none")

plot(lcnone.female$kt)
lines(lc.female$kt,col="blue")

lm.male <- lca(usa,"male",adjust="e0",years=1950:max(usa$year))
forecast.lm.male <- forecast(lm.male,h=20,jumpchoice = "actual")

bms.total <- bms(usa,"total",minperiod = 20, breakmethod = "bms")
forecast.bms.total <- forecast(bms.total,h=40)

bms.total$kt
plot(bms.total$kt)

bms.total <- lca(usa,"total",adjust="dxt", minperiod = 20, breakmethod = "bms")
forecast.bms.total <- forecast(bms.total,h=40)

usa.sm <- smooth.demogdata(usa)
fdm.male <- fdm(usa.sm,"male",order=3)
forecast.fdm.male <- forecast.fdm(fdm.male,h=30)
plot(forecast.fdm.male)
plot(forecast.fdm.male,"c")

fdm.male <- fdm(usa.sm,"male",method="rapca")
forecast.fdm.male <- forecast.fdm(fdm.male,h=30)
plot(forecast.fdm.male)

fdm.male <- fdm(usa.sm,"male", method="classical", weight=TRUE, beta = 0.1)
forecast.fdm.male <- forecast.fdm(fdm.male,h=20)
plot(forecast.fdm.male)


usa1960 <- extract.years(extract.ages(usa, 0:100, combine.upper=TRUE), 1960:2010)
smus1960 <- smooth.demogdata(usa1960)
usa.pr <- coherentfdm(smus1960,weight=TRUE,beta=0.05)
usa.pr.f <- forecast(usa.pr,h=25,ic="bic")

# compare with independent forecasts
usa.ind.fem <- fdm(smus1960,"female",weight=TRUE,beta=0.5)
usa.ind.mal <- fdm(smus1960,"male",weight=TRUE,beta=0.5)
#usa.ind <- fdm.ind(smus1960,weight=TRUE,beta=0.5)
usa.ind.fem.f <- forecast(usa.ind.fem,h=25,ic="bic")
usa.ind.mal.f <- forecast(usa.ind.mal,h=25,ic="bic")

usa.ind.f <- list(usa.ind.fem.f,usa.ind.mal.f)  ##is this useful??

par(mfrow=c(1,2))
plot(sex.ratio(smus1960),ylab="Sex ratio of rates: M/F",main="a. Coherent forecasts",ylim=c(0.7,3.5),lty=2,las=1,font.lab=2)
lines(sex.ratio(usa.pr.f))
plot(sex.ratio(smus1960),ylab="Sex ratio of rates: M/F",main="b. Independent forecasts",ylim=c(0.7,3.5),lty=2,las=1,font.lab=2)
lines(sex.ratio(usa.ind.f))   ########### fails


# Life expectancy forecasts
usa.pr <- coherentfdm(smus1960,weight=TRUE,beta=0.05)
usa.pr.f <- forecast(usa.pr,h=25)
e0.fcast.m <- e0(usa.pr.f, PI=TRUE,series="male")
e0.fcast.f <- e0(usa.pr.f, PI=TRUE,series="female")
plot(e0.fcast.m,ylim=c(67,85),col="blue",fcol="blue")
par(new=TRUE)
plot(e0.fcast.f,ylim=c(67,85),col="red",fcol="red")


