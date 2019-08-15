# mark potosnak
# 2019-06-20 new script that works with new EPA API (json)
# 2019-07-31 running in new GitHub directory


# runs independently, in its own directory
if(exists("My.Setwd")) {
   My.Setwd("GitHub/AoT-AQ/EPA")
}

# need to update this!
START <- "2018-01-01"
END   <- "2019-06-30"

# file name for output
epa.file <- "epa.rdata"

# use this output format consistently
std.str <- "%Y-%m-%d %H:%M"

# new epa output format is jsonlite
library(jsonlite)

# these are the species measured at the ComEd site
param <- list()
param[["epa.pm2.5"]] <- "88501"
param[["epa.o3.concentration"]] <- "44201"
param[["epa.no2.concentration"]] <- "42602"
param[["epa.so2.concentration"]] <- "42401"

# store in a regular matrix, even if missing data
epa.time <- seq(
   strptime(paste(START, "00:00"), format=std.str, tz="GMT"),
   strptime(paste(END,   "23:00"), format=std.str, tz="GMT"),
            by="hour")

epa <- matrix(nrow=length(epa.time), ncol=length(param))
rownames(epa) <- strftime(epa.time, format=std.str, tz="GMT")
colnames(epa) <- names(param)

# get data via EPA API (version 2)
# https://aqs.epa.gov/aqsweb/documents/data_api.html
username <- "mpotosna@depaul.edu"
passwd <- "taupemallard26"

for(i in names(param)) {
   url.txt <- paste(
      "https://aqs.epa.gov/data/api/sampleData/bySite?", 
      "email=", username, "&key=", passwd, 
      "&param=", 
      param[[i]],
      "&bdate=", gsub("-", "", START),
      "&edate=", gsub("-", "", END),
      "&site=0076",
      "&state=17&county=031", sep="")

   raw <- fromJSON(paste(readLines(url.txt, warn=FALSE)))
   if(raw$Header$status != "Success") {warning("EPA file corrupted")}
   raw <- flatten(raw$Data)  
   raw <- raw[order(paste(raw$date_gmt, raw$time_gmt)),]
   raw$time <- strptime(paste(raw$date_gmt, raw$time_gmt), 
      format="%Y-%m-%d %H:%M", tz="GMT")
   epa[,i] <- raw[match(
      strftime(epa.time, format=std.str, tz="GMT"), 
      strftime(raw$time, format=std.str, tz="GMT")), 
                  'sample_measurement']
}


save(epa, file=epa.file)
