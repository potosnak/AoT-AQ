# mark potosnak
# 2018-10-01

# reads in calibration information provided by SPEC
# extracts and creates mask matrix
# writes out head line
# for each BAD (line in input file)
# extracts appropriate cal constants (Sensitity, Izero, 1/m)
# extracts and creates cross-sensitivity matrix
# applies mask matrix to cross-sensivity matrix
# inverts cross-sensivity matrix
# writes out cal constants & cross-sensivity matrix

# 2018-11-01
# modified to use new calibration file from SPEC, uploaded Oct 2018

# both these files will be in root directory
# name of file from SPEC, except need to convert from Excel to CSV
INPUT.FILE <- "PBAYcumulativeparameters.csv"
OUTPUT.FILE <- "SPECcalibration.csv"

# you'll need to make sure you run in a reasonable working directory
# note, code assumes runs in a subdirectory
if(exists("My.Setwd")) {
   My.Setwd("GitHub/AoT-AQ/Calibration")
}

# read in sensor calibration data
cal.file <- INPUT.FILE
# by default, this will make line 3 the column names
cal.info <- read.csv(cal.file, skip=2) 
cal.header <- read.csv(cal.file, header=FALSE, as.is=TRUE, nrows=2)

# cal mask is in file, except need to add rows for IAQ and IRR
cal.mask <- matrix(c(1, rep(0, 6), 0, 1, rep(0, 5), 
   as.numeric(cal.header[1,47:81])), ncol=7, byrow=TRUE)

# make sure that calibration file is structure as expected
cal.labs <- c("IRR", "IAQ", "SO2", "H2S", "OZO", "NO2", "CMO")
if(!identical(cal.labs, colnames(cal.info)[5:11])) {
   stop("Calibration file has different order of gases")
}

# create matrix that will be output
# columns
# 0 BAD (will be rowname in R)
# 1-7 Sensitivity (Sens.GAS)
# 8-14 Izero (Izero.GAS)
# 15-21 n (n.GAS)
# 22-70 cross-sensitivity matrix, inverted (GAS.GAS)
# note, since rowname, BAD isn't a column and no name
out.labs <- c(
   paste("Sens", cal.labs, sep='.'),
   paste("Izero", cal.labs, sep='.'),
   paste("n", cal.labs, sep='.'),
   outer(cal.labs, cal.labs, paste, sep='.'))

# create the output matrix, with all NA values
output.mat <- matrix(nrow=nrow(cal.info), ncol=length(out.labs), 
   dimnames=list(cal.info[,'B.A.D.'], out.labs))

# do by variable, just so it's clear
# in original SPEC file, had "barcode sensivities, so that's why odd structure
output.mat[,paste("Sens", cal.labs, sep='.')] <- 
   as.matrix(cal.info[,paste(cal.labs, "", sep="")])
output.mat[,paste("Izero", cal.labs, sep='.')] <- 
   as.matrix(cal.info[,paste(cal.labs, "3", sep=".")])
output.mat[,paste("n", cal.labs, sep='.')] <- 
   as.matrix(cal.info[,paste(cal.labs, "5", sep=".")])


# go through the calibration file, board by board, line by line
for(board in cal.info[,'B.A.D.']) {
   # get matrix of cross sensitivities
   # note, need to construct rows for IRR and IAQ, but mask will overwrite most
   # to zero
   cal.mat <- matrix(as.numeric(c(rep(1, 14), 
      cal.info[cal.info[,'B.A.D.'] == board, 47:81])), byrow=TRUE, ncol=7)
   # apply mask and invert
   cal.mat <- t(solve(cal.mask*cal.mat))
   # put into calibration matrix for output
   output.mat[board, outer(cal.labs, cal.labs, paste, sep='.')] <- c(cal.mat)
}

# write the file
write.csv(output.mat, file=OUTPUT.FILE)
