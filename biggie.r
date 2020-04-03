# mark potosnak
# 2017-11-21
# 2018-06-20
# 2018-08-22
# 2019-07-31 modifying for new GitHub to make cleaner
# works with GitHub files
# before running need to: 
#  1) get hourly data: see RawDataReduction directory
#  2) get epa data: see EPA directory
#  3) generate calibration data: see RawDataReduction directory
#  4) CO and Solar contain necessary data, but no scripts

# EPA data pull should work with my password
# requires calibration file from SPEC, which is on GitHub

# this version checks to see if the hrf agrees with my calculations

# you'll need to make sure you have in a reasonable working directory
if(exists("My.Setwd")) {
   # My.Setwd("GitHub/AoT-AQ")
}

# look at hrf instead of epa (but still called epa in plots)
hrf <- TRUE

# compare co to Max Berkelhammer's instrument
max.co <- FALSE

# apply by month co zero 
CO.by.month <- FALSE

# dates, will include data from last day
# times & dates here are GMT
# this is start of good collocation data
# START <- "2018-09-01"
# START <- "2018-12-01"
# START <- "2018-03-01"
START <- "2018-03-01"
START <- "2020-03-01"

# should be about current date
# END   <- "2019-03-31"
# END   <- "2018-12-31"
# END   <- "2019-06-30"
# END   <- "2018-05-10"
END   <- "2020-03-31"

# History of nodes at ComEd, from photos
# Jun 8 072 and 01C
# Apr 12 072 and 070
# Feb 8 072 and 070
# 072 has been installed since the beginning, except down for a bit
# 01C installed Jun 9
# vsn <- c("072", "01C", "070", "057")
# vsn <- c("072", "01C")
# vsn <- c("072", "028", "04C", "076", "09D", "092", "067")
# 15 nodes with most data from Adam's run
# vsn <- c("072", "08B", "02D", "08F", "07C",
vsn <- c("072", "08B", "08F", "07C",
         "08C", "086", "086", "088", "004",
         "09C", "089", "02A", "073", "07F")
vsn <- c("072")

# source all functions, note vsn needs to be set first
source("func.r")

# read all the hourly data into one object
all.dat <- Combine(START, END, vsn, x.time)

# get EPA data
epa <- EPA()
# convert ozone from ppm to ppb
epa[,"epa.o3.concentration"] <- epa[,"epa.o3.concentration"]*1e3

# apply the calibration for the chemsense board
all.dat <- Calibration(all.dat)

# what species to look at for main sensor that has been running the entire time
# vsn.plot <- "04C"
# vsn.plot <- "09D"
# vsn.plot <- "072"
for(vsn.plot in vsn) {

# write out a file for David
FullHourly <- cbind(all.dat[[vsn.plot]], epa)
write.csv(FullHourly, file=paste("FullHourly", vsn.plot, ".csv", sep=""))

if(hrf) {
   # all chemsense species
   species.list <- c("o3.concentration", "no2.concentration", 
      'reducing_gases.concentration', 'co.concentration', 'h2s.concentration', 
      'so2.concentration', 'oxidizing_gases.concentration')
} else {
   # only species for which we have epa data
   species.list <- c("o3.concentration", "no2.concentration", 
      'so2.concentration')
}

if(max.co) {
   species.list <- c(species.list, "co.concentration")
   co <- read.csv("CO/FinalCOData.csv")
   co.time <- ISOdate(co$Year, co$Month, co$Day, co$Hour, co$Minute, 
      tz="Etc/GMT+6")
   DT.format <- "%Y-%m-%d %H:00"
   co.fac <- as.factor(strftime(co.time, format=DT.format, tz="GMT"))
   co.hr <- tapply(co$CO..ppb., co.fac, mean, na.rm=TRUE)
   co.time <- strptime(levels(co.fac), format=DT.format, tz="GMT")
   co.fin <- co.hr[match(strftime(x.time,   format=std.str, tz="GMT"), 
             strftime(co.time, format=std.str, tz="GMT"))]

   epa <- cbind(epa, 
      epa.co.concentration=co.fin)
   # rownames(epa) <- strftime(x.time, format=std.str, tz="GMT")
}
for(species in species.list) {

   # extract chemical species concentration to plot
   y.aot <- all.dat[[vsn.plot]][,species]
   if(hrf) {
      # says epa, but really hrf
      y.epa <- all.dat[[vsn.plot]][,paste(species, "hrf", sep=".")]*1e3
   } else {
      # normal, really epa data
      y.epa <- epa[,paste('epa', species, sep=".")]
   }

   # something crazy with data from Memorial Day heat wave--kill data
   look <-  x.time >= strptime("2018-05-20 00:00", format=std.str, tz="GMT") & 
            x.time <= strptime("2018-05-31 23:00", format=std.str, tz="GMT")
   # get rid of values that are just way too high
   # look <- look | y.aot > 1e4
   y.aot[look] <- NA

   My.Plot(paste(species, "AoT and EPA"))
   plot(x.time, y.aot, type='l', col=2,
        ylim=c(0, max(y.aot, na.rm=TRUE)), xlab="Date", ylab=labs[[species]], 
        main=paste(labs[[species]], " ", 
                  strftime(min(x.time), "%Y %b %d"), " - ", 
                  strftime(max(x.time), "%Y %b %d"), " ", sep=""))
   lines(x.time, y.epa, col=4)
   legend('topleft', c(vsn.plot, "EPA"), lty=1, col=c(2,4))
   dev.off()

   look <- !is.na(y.aot) & !is.na(y.epa)

   My.Plot(paste(species, "AoT vs EPA"))
   plot(y.epa, y.aot, xlab="EPA", ylab="AoT",
      main=paste(labs[[species]], " ", 
         strftime(min(x.time[look]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look]), "%Y %b %d"), " ", sep=""))
   fit <- lm(y.aot ~ y.epa)
   rsq <- round(summary(fit)$r.sq, 2)
   text(par('usr')[1], par('usr')[4], bquote(R^2 == .(rsq)), adj=c(-0.5,1.5))
   abline(fit)
   abline(0,1,col=4)
   legend('bottomright', c("Linear fit", "1:1"), lty=1, col=c(1,4), bty='n')

   dev.off()

   My.Plot(paste(species, "resid"))
   plot(x.time[look], resid(fit), col=2, xlab="Date", ylab=labs[[species]],
      type='l', 
      main=paste("Residual ", labs[[species]], " ", 
         strftime(min(x.time[look]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look]), "%Y %b %d"), " ", sep=""))
   if(TRUE) {
      solar <- read.csv("Solar/2018_DGD.csv")
      solar.date <- ISOdate(solar$Year, solar$Month, solar$Day)
      lines(solar.date, solar$A*2)
      day.fact <- factor(strftime(x.time[look], format="%j"), levels=c(1:365))
      day.resid <- tapply(resid(fit), day.fact, mean, na.rm=TRUE)
   }
   dev.off()

   My.Plot(paste(species, "vs temp"))
   plot(all.dat[[vsn.plot]][,"chem.temp"], y.aot, 
      xlab="Board Temperature (deg C)", 
      ylab=labs[[species]],
      main=paste("Concentration ", labs[[species]], " ", 
         strftime(min(x.time), "%Y %b %d"), " - ", 
         strftime(max(x.time), "%Y %b %d"), " ", sep=""))
   dev.off()
   
   My.Plot(paste(species, "resid vs temp"))
   plot(all.dat[[vsn.plot]][look,"chem.temp"], resid(fit), 
      xlab="Board Temperature (deg C)", 
      ylab=labs[[species]],
      main=paste("Residual ", labs[[species]], " ", 
         strftime(min(x.time[look]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look]), "%Y %b %d"), " ", sep=""))
   dev.off()
   
   doy <- as.numeric(x.time - ISOdatetime(2018, 1, 1, 0, 0, 0, tz="GMT") + 1)
   hod <- as.factor(floor((doy[look]%%1)*24))

   My.Plot(paste(species, "resid vs hour of day"))
   plot(hod, resid(fit), xlab="Hour of day [GMT]",
      ylab=labs[[species]],
      main=paste("Residual ", labs[[species]], " ", 
         strftime(min(x.time[look]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look]), "%Y %b %d"), " ", sep=""))
   dev.off()

   My.Plot(paste(species, "AoT and EPA night"))
   look3 <- doy%%1 < 6/24 | doy%%1 >= 22/24
   plot(x.time[look3], y.aot[look3], type='p', col=2,
      ylim=c(0, max(y.aot[look3], na.rm=TRUE)), xlab="Date", 
      ylab=labs[[species]], 
      main=paste("Night ", labs[[species]], " ", 
         strftime(min(x.time[look3]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look3]), "%Y %b %d"), " ", sep=""))
   points(x.time[look3], y.epa[look3], col=4)
   legend('topleft', c(vsn.plot, "EPA"), lty=1, col=c(2,4))
   dev.off()

   # look at absolute vapor pressure
   svp <- 0.61365*exp(17.502*all.dat[[vsn.plot]][,'chem.temp']/(240.97 +
      all.dat[[vsn.plot]][,'chem.temp']))
   if(is.na(match("sht25.humidity.hrf", colnames(all.dat[[vsn.plot]])))) {
      vp <- rep(NA, nrow(all.dat[[vsn.plot]]))
   } else {
      vp <- svp*all.dat[[vsn.plot]][,"sht25.humidity.hrf"]
   }

   # only print if there is enough vp data available 
   if(sum(!is.na(vp[look])) > 5) {

      My.Plot(paste(species, "resid vs abs humidity"))
      plot(vp[look], resid(fit), 
         xlab="Absolute humidity (kPa)", 
         ylab=labs[[species]],
         main=paste("Residual ", labs[[species]], " ", 
            strftime(min(x.time[look]), "%Y %b %d"), " - ", 
            strftime(max(x.time[look]), "%Y %b %d"), " ", sep=""))
      dev.off()
   }

   My.Plot(paste(species, "hist"))
   hist(resid(fit), breaks=20, main="Residuals",  
      xlab=labs[[species]])
   dev.off()

   look.may <- x.time[look] > ISOdate(2018,11,5) & x.time[look] < 
      ISOdate(2018,11,10)
   look.may <- x.time[look] > ISOdate(2018,5,1) & x.time[look] < 
      ISOdate(2018,5,7)
   if(sum(look.may) > 5) {
      My.Plot(paste(species, "resid May"))
      par(mar=c(4,4,2,8))
      plot(x.time[look][look.may], resid(fit)[look.may], 
         xlab="Date", ylab=labs[[species]], type='l',
         main=paste("Residual ", labs[[species]], " ",
            strftime(min(x.time[look][look.may]), "%Y %b %d"), " - ", 
            strftime(max(x.time[look][look.may]), "%Y %b %d"), " ", sep=""))
      par(new=TRUE)
      plot(x.time[look][look.may], 
         all.dat[[vsn.plot]][look,'chem.temp'][look.may],
         axes=FALSE, type='l', col=4, xlab="", ylab="")
      axis(4)
      mtext(expression(Temperature~(degree~C)), 4, line=2.5, col=4)
      if(sum(!is.na(vp[look])) > 5) {
         par(new=TRUE)
         plot(x.time[look][look.may], vp[look][look.may],
            axes=FALSE, type='l', col=3, xlab="", ylab="")
         axis(4, line=4)
         mtext("Vapor pressure (kPa)", 4, line=6.5, col=3)
      }
      # see if PM 2.5 from EPA correlates (no, it doesn't)
      # par(new=TRUE)
      # plot(x.time[look][look.may], epa[look,"epa.pm2.5"][look.may],
         # axes=FALSE, type='l', col=5, xlab="", ylab="")
      dev.off()
   }

}

if(sum(look.may) > 5) {
   # look at a spike, o3 and no2, raw and hrf, epa
   My.Plot("may o3 and no2")
   par(mfrow=c(2,1),mar=c(0,4,4,4))
   plot(x.time[look][look.may], 
      all.dat[[vsn.plot]][look,'o3.concentration'][look.may], 
      xlab="Date", ylab="Ozone", type='l',
      main=paste("Spike ", " ",
         strftime(min(x.time[look][look.may]), "%Y %b %d"), " - ", 
         strftime(max(x.time[look][look.may]), "%Y %b %d"), " ", sep=""), 
      axes=FALSE)
   axis(2)
   box()
   lines(x.time[look][look.may], epa[look,'epa.o3.concentration'][look.may], 
      col=2)
   par(new=TRUE)
   plot(x.time[look][look.may], 
      all.dat[[vsn.plot]][look,'o3.concentration.raw'][look.may],
      axes=FALSE, type='l', col=4, xlab="", ylab="")
   axis(4)
   mtext(expression(Raw~current~(nA)), 4, line=2.5, col=4)
   par(mar=c(4,4,0,4))
   plot(x.time[look][look.may], 
      all.dat[[vsn.plot]][look,'no2.concentration'][look.may], 
      xlab="Date", ylab="NO2", type='l')
   lines(x.time[look][look.may], epa[look,'epa.no2.concentration'][look.may], 
      col=2)
   legend('topleft', c("AoT", "EPA", "current"), lty=1, col=c(1, 2, 4), bty='n')
   par(new=TRUE)
   plot(x.time[look][look.may], 
      all.dat[[vsn.plot]][look,'no2.concentration.raw'][look.may],
      axes=FALSE, type='l', col=4, xlab="", ylab="")
   axis(4)
   mtext(expression(Raw~current~(nA)), 4, line=2.5, col=4)
   dev.off()

   # look at all raw currents
   My.Plot("may raw currents")
   gases <- c('no2', 'o3', 'h2s', 'co', 'so2')
   matplot(x.time[look][look.may], all.dat[[vsn.plot]][look,paste(gases, 
      'concentration.raw', sep=".")][look.may,], type='l', lty=1)
   legend('topright', gases, lty=1, col=c(1:5))
dev.off()
}
}
