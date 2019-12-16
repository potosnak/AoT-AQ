# mark potosnak
# 2018-08-14
# 2018-08-24 modifying to remove first 10 min of data after a 5 min 
#            gap, per David Peaslee's advice

# length of gap to search for in seconds
GAP.LENGTH <- 5*60
# length of readings to dump in seconds
DISCARD.LENGTH <- 10*60

# Switching to new tarball method of getting data
# first, need to extract chemsense node data from the tarball
# with Bash on Unbuntu, using script slicer.sh
# of course, that takes a while!

# next, this makes an hourly-averaged file, saved as a text file

# runs independently, in its own directory
if(exists("My.Setwd")) {
   My.Setwd("GitHub/AoT-AQ/RawDataReduction")
}

# mac.address file name
mac.file <- "mac.address.rdata"

# always save mac addresses, since key for getting calibration data
# load old mac address info, will get over written below if necessary
# just in case, save alphasense ID too
if(is.na(file.info(mac.file)[1])) {
   mac.address <- list()
   alpha.id <- list()
} else {
   load(mac.file)
}

# date format for use in R, for aggregating to hour
# trick of using "00" always for minute gets aggregate to hour
DT.format <- "%Y-%m-%d %H:00"
# date format used by waggle
DT.waggle <- "%Y/%m/%d %T"
# currently, hard coded to one node
# just as a check, make sure the node id is correct
# 072, long-term node at ComEd
# node.id <- "001e06113107"
# now, do for multiple nodes

# for(node.id in c("001e06109f62", "001e0611441e")) {
for(node.id in c("001e0610b9e7", "001e0610ef27", "001e06113acb", 
      "001e06113d22", "001e0611441e")) {
   # file name of downloaded data, created by extract.sh script
   DATA.FILE <- paste("node", node.id, ".csv", sep="")

   # column used for mean value by waggle
   wag.cols <- gsub("_", ".", scan(DATA.FILE, nlines=1, what="", sep=","))

   # output file name
   HOURLY.FILE <- paste("hourly", node.id, ".csv", sep="")

   # bring in all the data
   foo <- matrix(scan(DATA.FILE, sep=",", what="", skip=1), byrow=TRUE, 
      ncol=length(wag.cols))
   colnames(foo) <- wag.cols

   # do a couple of sanity checks
   x <- unique(foo[,"subsystem"])
   if(length(x) !=1 || "chemsense" != x) {
      warning("extract.csv does not only contain chemsense info")
   }

   x <- unique(foo[,"node.id"])
   if(length(x) != 1 || node.id != x) {
      warning("extract.csv does not only contain the expected node ID")
   }

   # note, better than converting to POSIXlt, 
   # since following steps are much faster
   foo.time <- as.POSIXct(foo[,"timestamp"], format=DT.waggle, tz="GMT")
   foo.diff <- c(GAP.LENGTH + 1, diff(foo.time))
   look <- foo.diff > GAP.LENGTH
   gaps <- c(1:length(foo.time))[look]

   # not very efficient, but get rid of bad data by looping through gaps
   bad <- rep(FALSE, length(foo.time))
   for(i in gaps) {
      gap.time <- foo.time[i]
      bad <- bad | (foo.time >= gap.time & foo.time <= 
         (DISCARD.LENGTH + gap.time))
   }

   foo <- foo[!bad,]
   foo.time <- foo.time[!bad]

   # each sensor reading is a row, so collect unique sensors into 
   # columns
   sensors <- paste(foo[,"sensor"], foo[,"parameter"], sep=".")
   sensor.labs <- sort(unique(sensors))

   # mac address is a string and used as key for calibration data
   # assume this doesn't change with time, so just save last value
   sensor.labs <- sensor.labs[!grepl("chemsense.id", sensor.labs)]
   mac.address[[node.id]] <- unique(foo[sensors==
     "chemsense.id", "value.hrf"])

   # this will be hourly, since format has "00" for minutes
   time.fac <- as.factor(strftime(foo.time, format=DT.format, tz="GMT"))

   # this will hold hourly data
   sensor.labs.2 <- c(paste(sensor.labs, "raw", sep="."), paste(sensor.labs, 
      "hrf", sep="."))
   foo2 <- matrix(nrow=length(levels(time.fac)), 
      ncol=length(sensor.labs.2), dimnames = list(levels(time.fac), 
      sensor.labs.2))
   for(k in sensor.labs) {
      foo2[,paste(k, "raw", sep=".")] <- tapply(as.numeric(foo[sensors == k,
         "value.raw"]), time.fac[sensors == k], mean)
      foo2[,paste(k, "hrf", sep=".")] <- tapply(as.numeric(foo[sensors == k,
         "value.hrf"]), time.fac[sensors == k], mean)
   }
   # store for later analysis
   write.csv(foo2, file=HOURLY.FILE)

}

# save mac address and alpha id info 
save(mac.address, alpha.id, file=mac.file)

