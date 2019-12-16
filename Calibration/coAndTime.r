# mark potosnak
# 2019-08-01

# calculates calibration values based on observed data

My.Setwd("GitHub/AoT-AQ/Calibration")

# make sure this starts with the default calibration
source("CalFileOffline.r")

# bin size for mins
bin.size <- 2.5

# as of Aug 1 2019, this is the available data
START <- "2018-03-01"
END   <- "2019-06-30"

# try first wil standard node
vsn.plot <- c("04C")
vsn.plot <- c("09D")
vsn.plot <- c("072")

# functions used by main script: only need one that combines data
vsn <- vsn.plot
setwd("..")
source("func.r")

# read in all raw data
all.dat <- Combine(START, END, vsn.plot, x.time)
# don't really use this, except gets board temp
all.dat <- Calibration(all.dat, temp.only = TRUE)
# need this for writing cal data
mac.file <- "RawDataReduction/mac.address.rdata"
load(mac.file)
setwd("Calibration")

# get original cal data
CO.FILE <- "SPECcalibration.csv"
cal.info <- read.csv(CO.FILE, row.names=1)

# only do fit over the range where this is sufficient data
look.temp <- all.dat[[vsn.plot]]$chem.temp > -10 & 
   all.dat[[vsn.plot]]$chem.temp < 35

# make a factor, every 5 degrees
temp.fac <- 
   as.factor(floor(all.dat[[vsn.plot]]$chem.temp[look.temp]/bin.size)*bin.size)

# subtract out the background CO, from Max's measurements (a bit arbitrary)
co.background <- 150 # (from Max's data)
# get cal data
node.mac <- toupper(mac.address[[node.id]])
if(length(node.mac) > 1) {
   warning(paste("Multiple MAC addresses; selecting last one"))
   warning(paste(node.mac, ""))
   node.mac <- node.mac[length(node.mac)]
}
co.span <- as.numeric(cal.info[node.mac, c("Sens.CMO")])

co.noback <- all.dat[[1]]$co.concentration.raw - co.background * co.span

# get %5 value, close to min
min.current <- tapply(co.noback[look.temp], temp.fac, quantile, 0.05)

# diagnostic plot
My.Plot("COcalibration")
plot(all.dat[[1]]$chem.temp, co.noback,
   xlab=My.Labs("Board temperature", "degree*C"),
   ylab=My.Labs("CO current", "nA"),
   main="Node at ComEd site")
x <- as.numeric(levels(temp.fac)) + bin.size/2
lines(x, min.current, lwd=3, col=3)

# do an exponential fit
fit <- lm(log(min.current) ~ x)
x.new <- c(-30:40)
# put coef's into form for fit at 40 deg C
m <- coef(fit)[2]
Izero40 <- exp(m*40+coef(fit)[1])
lines(x.new, Izero40*exp((x.new - 40)*m), lwd=2, col=2)

# read in original calibration file
old.cal <- as.numeric(cal.info[node.mac, c("Izero.CMO", "n.CMO")])
print("old cal")
print(paste("m =", round(1/old.cal[2], 4), "Izero40 =", round(old.cal[1]*1e4)))
lines(x.new, old.cal[1]*exp((x.new - 40)/old.cal[2])*1e3, lwd=2, col=4)

# replace cal info
cal.info[node.mac, c("Izero.CMO", "n.CMO")] <- as.character(c(Izero40/1e3, 1/m))
write.csv(cal.info, CO.FILE)
legend("topleft", c("Factory calibration", "5th percentile", "Exp fit"),
   lty=1, col=c(4, 3, 2))
dev.off()

# first, try to do by quarter
q.fac <- as.factor(paste(as.POSIXlt(x.time)$year + 1900, quarters(x.time)))
coefs <- matrix(NA, ncol=2, nrow=nlevels(q.fac), dimnames=list(levels(q.fac),
   c("A", "m")))
for(foo in levels(q.fac)) {
   look.foo <- foo == q.fac
   min.current <- tapply(co.noback[look.temp][look.foo], temp.fac[look.foo], 
      quantile, 0.05)
   if(sum(!is.na(min.current)) > 2) {
      plot(all.dat[[1]]$chem.temp[look.foo], co.noback[look.foo],
         xlab=My.Labs("Board temperature", "degree*C"),
         ylab=My.Labs("CO current", "nA"),
      main=paste("Quarter = ", foo))
      lines(x, min.current, lwd=3, col=3)
      fit <- lm(log(min.current) ~ x)
      m <- coef(fit)[2]
      Izero40 <- exp(m*40+coef(fit)[1])
      lines(x.new, Izero40*exp((x.new - 40)*m), lwd=2, col=2)
      lines(x.new, old.cal[1]*exp((x.new - 40)/old.cal[2])*1e3, lwd=2, col=4)
      legend("topleft", c("Factory calibration", "5th percentile", "Exp fit"),
         lty=1, col=c(4, 3, 2))
      coefs[foo,] <- c(Izero40, m)
      # locator(1)
   } else {
      coefs[foo,] <- c(NA, NA)
   }
}
# now, trim this to where data actually exists
q.look <- !is.na(coefs[,"A"])
coefs <- coefs[q.look,]
q.names <- levels(q.fac)[q.look]
i <- 1
My.Plot("Changing fits")
plot(NA, NA, xlim=c(0, 35), ylim=c(0, 30e3), 
   xlab=My.Labs("Board temperature", "degree*C"),
   ylab=My.Labs("CO current", "nA"))
lines(x.new, old.cal[1]*exp((x.new - 40)/old.cal[2])*1e3, lwd=2, col=i)
for(foo in q.names) {
   Izero40 <- coefs[foo, 'A']
   m <- coefs[foo, 'm']
   i <- i + 1
   lines(x.new, Izero40*exp((x.new - 40)*m), lwd=2, col=i)
}
legend('topleft', c("Factory", q.names), lty=1, col=c(1:i))
dev.off()
   
My.Plot("Coefs")
q.x <- c(1:length(q.names))
fit <- lm(log(coefs[,1]) ~ q.x)
par(mfrow=c(2,1), mar=c(4,4,2,2), mgp=c(2.5, 1, 0))
plot(q.x, coefs[,2], xlab="", axes=FALSE,
   ylab=My.Labs("m", "degree*C^{-1}"))
# simply do mean
m.all <- mean(coefs[,2])
abline(h=m.all)
axis(1, q.names, at=q.x)
box()
axis(2)
plot(q.x, coefs[,1], axes=FALSE, xlab="", ylab=My.Labs("Zero current", "nA"))
y <- exp(predict(fit))
lines(q.x, match(q.names, names(y), y))
axis(1, q.names, at=q.x)
box()
axis(2)
dev.off()

# now, use average exponential for entire fit, but get new Izero40 value
# do by month, not quarter, since just doing average
m.fac <- as.factor(strftime(x.time, format="%Y %m"))
Izero <- c(1:nlevels(m.fac))
m.time <- c(1:nlevels(m.fac))
names(Izero) <- levels(m.fac)
names(m.time) <- levels(m.fac)
# for compatiability with original data, use Izero40
for(foo in levels(m.fac)) {
   look.foo <- foo == m.fac
   min.current <- tapply(co.noback[look.temp][look.foo], temp.fac[look.foo], 
      quantile, 0.05)
   if(sum(!is.na(min.current)) > 2) {
      # day sinc Jan 1 2018, starting with 1 = Jan 1 2018
      m.time[foo] <- mean(x.time[look.temp][look.foo] - 
         ISOdate(2018, 1, 1, 0, tz="GMT") + 1, na.rm=TRUE)
      plot(all.dat[[1]]$chem.temp[look.foo], co.noback[look.foo],
         xlab=My.Labs("Board temperature", "degree*C"),
         ylab=My.Labs("CO current", "nA"),
      main=paste("Month = ", foo))
      lines(x, min.current, lwd=3, col=3)
      # min.current = Izero40*exp((temp - 40)*m)
      # Izero40 = min.current/exp((temp - 40)*m), simply use mean
      Izero40 <- mean(min.current/exp(m.all*(x-40)), na.rm=TRUE)
      lines(x.new, Izero40*exp((x.new - 40)*m.all), lwd=2, col=2)
      lines(x.new, old.cal[1]*exp((x.new - 40)/old.cal[2])*1e3, lwd=2, col=4)
      legend("topleft", c("Factory calibration", "5th percentile", "Exp fit"),
         lty=1, col=c(4, 3, 2))
      Izero[foo] <- Izero40
      # locator(1)
   } else {
      m.time[foo] <- NA
      Izero[foo] <- NA
   }
}
My.Plot("Izero estimate")
par(mgp=c(2.5, 1, 0))
plot(m.time, Izero, xlab="Days starting Jan 1 2018", 
   ylab=My.Labs("Zero current at 40*degree*C", "nA"))
abline(v=366, lty=2)
text(366, par('usr')[4], "2019", adj=c(-0.2,1.1))
text(366, par('usr')[4], "2018", adj=c(1.1,1.1))
# do linear fit for 2018, excluding Dec
look <- m.time < 350
fit <- lm(Izero ~ m.time, sub=look)
abline(fit)
# just do average for 2019
look <- m.time > 370
Izero.2019 <- mean(Izero[look], na.rm=TRUE)
# get the intersection
day.intersect <- (Izero.2019 - coef(fit)[1])/coef(fit)[2]
# plot
segments(day.intersect, Izero.2019, par('usr')[2], Izero.2019)
dev.off()

# replace again cal info
# note that will only be valid for 2019, so need to do adjustments for 2018
cal.info[node.mac, c("Izero.CMO", "n.CMO")] <- 
   as.character(c(Izero.2019/1e3, 1/m.all))
print("new cal")
print(paste("m =", round(m.all, 4), "Izero40 =", round(Izero.2019)))
write.csv(cal.info, CO.FILE)

# save details of the fit
Izero.estimate <- list(
   slope.2018 = coef(fit)[2]/1e3,
   intercept.2018 = coef(fit)[1]/1e3,
   day.intersect = day.intersect,
   mean.2019 = Izero.2019/1e3)

save(Izero.estimate, file="Izero.rdata")
