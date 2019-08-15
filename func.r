# mark potosnak
# 2018-06-21
# 2019-07-31 redoing with new GitHub structure

# functions for analyzing AoT data
# Combine: very simple, just reads hourly data and restricts it
# Calibrate: applies calibration data
# epa: gets EPA data, but most work done in the EPA directory now

# make sure some directories exist
if(!file.exists("Plots")) {
   dir.create("Plots")
}

# file locations
# mac.address file name
mac.file <- "RawDataReduction/mac.address.rdata"
# epa data, once retrieved
epa.file <- "EPA/epa.rdata"
# calibration file, needs to be created
cal.file <- "Calibration/SPECcalibration.csv"

# use this output format consistently
std.str <- "%Y-%m-%d %H:%M"

# following for plotting
x.time <- seq(
   strptime(paste(START, "00:00"), format=std.str, tz="GMT"),
   strptime(paste(END, "23:00"), format=std.str, tz="GMT"),
            by="hour")

labs <- list()
labs[["Alphasense.pm2.5"]] <- "PM2.5 [ug/m3]"
labs[["o3.concentration"]] <- "Ozone [ppb]"
labs[["no2.concentration"]] <- "NO2 [ppb]"
labs[["so2.concentration"]] <- "SO2 [ppb]"
labs[["co.concentration"]] <- "CO [ppb]"

# helper function that gets the date from the rownames
# timestamp is at beginning of hour, consistent with EPA
GetDate <- function(i) { strptime(rownames(all.dat[[i]]), format=std.str, 
                                  tz="GMT")}

# sets the plotting type
My.Plot <- function(title) {
   # pdf(paste("Plots/", title, ".pdf", sep=""))
   png(paste("Plots/", title, ".png", sep=""), res=150, width=750, height=750)
}

# need node_id for each vsn, that's how files are stored
node.info <- read.csv(url(paste(
   "https://raw.githubusercontent.com/waggle-sensor/beehive-server/",
   "master/publishing-tools/projects/AoT_Chicago.complete/nodes.csv", sep="")),
   as.is=TRUE)
# add some by hand, since not in file
node.info <- rbind(node.info, c("001e0610ba13", "AoT Chicago", "01C", 
   "7801 S Lawndale Ave Chicago IL", "41.751142", "-87.71299", 
   "AoT Chicago (S) [CA] {ComEd}"))
node.id <- node.info[match(vsn, node.info[,'vsn']), "node_id"]
names(node.id) <- vsn

# Combine
# on a per node basis, combines daily data into a list
# this works starting from big tarball file

# needs:
# working directory for reading 
# START, END
# vsn 
# note, puts everything into all.dat, which must exist
# note: also needs x.time

Combine <- function(START, END, vsn, x.time) {
   all.dat <- list()
   for (i in vsn) {
      # name of extracted file, from hourly.r script
      file.name <- paste("RawDataReduction/hourly", node.id[i], ".csv", sep="")
      foo <- read.csv(file.name, row.names=1)
      x <- substr(rownames(foo), 1, 10)
      foo <- foo[x >= START & x <= END,]
      # put onto common time stamp
      all.dat[[i]] <- foo[match(strftime(x.time, format=std.str, tz='GMT'), 
         rownames(foo)), ]
      rownames(all.dat[[i]]) <- strftime(x.time, format=std.str, tz='GMT')
   }
   all.dat
}

# calibration

# applies calibration data to hourly data
# needs csv calibration info file
# uses linear interpolation of zero vs temperature
# will do for everything hat is in current all.dat

# needs
# all.dat
# calibration file

Calibration <- function(all.dat) {

   # need mac address information for doing calibrations
   if(is.na(file.info(mac.file)[1])) {
      error("No MAC address available--need for calibraiton")
   } else {
      load(mac.file)
   }


   # using CalFileOffline.r, do inversion ahead of time
   cal.info <- read.csv(cal.file, row.names=1)

   # want in same order as cal data
   chemsense.labs <- paste(c("reducing_gases", "oxidizing_gases", 
      "so2", "h2s", "o3", "no2", "co"), "concentration", sep=".")

   # cal mask is in file, except need to add rows for IAQ and IRR
   cal.labs <- c("IRR", "IAQ", "SO2", "H2S", "OZO", "NO2", "CMO")
   names(chemsense.labs) <- cal.labs

   for(board in names(all.dat)) {
      # get line of cal data
      cal <- as.numeric(cal.info[toupper(mac.address[[node.id[board]]]),])
      names(cal) <- colnames(cal.info)
      # get matrix of cross sensitivities
      # note, need to construct rows for IRR and IAQ
      cal.mat <- matrix(cal[outer(cal.labs, cal.labs, paste, sep='.')], ncol=7)
      row.names(cal.mat) <- cal.labs
      # need to get average temperature of ADC board
      chem.temp <- apply(all.dat[[board]][,
         paste("at", c(0:3), ".temperature.hrf", sep="")], 1, mean)
      all.dat[[board]][,"chem.temp"] <- chem.temp
      # loop through to apply calibration to all species
      # subtrack zero current and apply span factor
      # use Izero40 (zero current measured at 40 deg C) and
      # m (n) value measured for each individual sensor
      # see very big spreadsheet
      for(i in cal.labs) {
         chem.zero <- cal[paste("Izero", i, sep=".")]*1e3* 
            exp((chem.temp - 40)/cal[paste("n", i, sep=".")])
         # substract zero from raw signal and put into new column
         all.dat[[board]][,chemsense.labs[i]] <- 
            all.dat[[board]][,paste(chemsense.labs[i], "raw", sep=".")] - 
            chem.zero
         # apply span
         all.dat[[board]][,chemsense.labs[i]] <- 
            all.dat[[board]][,chemsense.labs[i]]/cal[paste("Sens", i, sep=".")]
      }
      # apply cross sensitivities matrix
      Cross <- function(x){cal.mat%*%x}
      all.dat[[board]][,chemsense.labs] <- 
         t(apply(all.dat[[board]][,chemsense.labs], 1, Cross))
      
      # implement 'advanced' filter
      Filter <- function(x) {
         a <- x[1]; b <- x[2]
         if(is.na(a) | is.na(b)) {return(c(NA, NA))}
         if(a < 0 & b < 0) {return(c(0,0))}
         if(a > 0 & b < 0) {return(c(a + b, 0))}
         if(a < 0 & b > 0) {return(c(0, a + b))}
         return(c(a,b))
      }
      all.dat[[board]][,c("o3.concentration", "no2.concentration")] <- 
         t(apply(all.dat[[board]][,c("o3.concentration", "no2.concentration")],
         1, Filter))
   }
   all.dat
}

# epa
# it's variable, but usually the data posts with a one-month lag
# now, most of work done in a separate script, in a separate directory
# needs epa.file
EPA <- function() {

   # need std.str set and working directory
   # check if a previous file has been saved
   if(file.exists(epa.file)) {load(epa.file)} else {error("no epa file")}

   # save file with specified dates, but then put onto regular time stamp
   epa[match(strftime(x.time,   format=std.str, tz="GMT"), rownames(epa)),]

}
