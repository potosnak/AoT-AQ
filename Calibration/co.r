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
vsn <- c("072")

# functions used by main script: only need one that combines data
source("../func.r")

# read in all raw data
setwd("..")
all.dat <- Combine(START, END, "072", x.time)
# don't really use this, except gets board temp
all.dat <- Calibration(all.dat)
# need this for writing cal data
mac.file <- "RawDataReduction/mac.address.rdata"
load(mac.file)
setwd("Calibration")

# get original cal data
CO.FILE <- "SPECcalibration.csv"
cal.info <- read.csv(CO.FILE, row.names=1)

# only do fit over the range where this is sufficient data
look <- all.dat[[vsn]]$chem.temp > -10 & all.dat[[vsn]]$chem.temp < 35

# make a factor, every 5 degrees
temp.fac <- as.factor(floor(all.dat[[vsn]]$chem.temp[look]/bin.size)*bin.size)

# subtract out the background CO, from Max's measurements (a bit arbitrary)
co.background <- 150 # (from Max's data)
co.span <- as.numeric(cal.info[toupper(mac.address[[node.id]]), c("Sens.CMO")])
co.noback <- all.dat[[1]]$co.concentration.raw - co.background * co.span

# get %5 value, close to min
min.current <- tapply(co.noback[look], temp.fac, quantile, 0.05)

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
print(paste("m =", round(m, 4), "Izero40 =", round(Izero40)))

# read in original calibration file
old.cal <- as.numeric(cal.info[toupper(mac.address[[node.id]]),
   c("Izero.CMO", "n.CMO")])
lines(x.new, old.cal[1]*exp((x.new - 40)/old.cal[2])*1e3, lwd=2, col=4)

# replace cal info
cal.info[toupper(mac.address[[node.id]]), c("Izero.CMO", "n.CMO")] <- 
   as.character(c(Izero40/1e3, 1/m))
write.csv(cal.info, CO.FILE)
legend("topleft", c("Factory calibration", "5th percentile", "Exp fit"),
   lty=1, col=c(4, 3, 2))
dev.off()
