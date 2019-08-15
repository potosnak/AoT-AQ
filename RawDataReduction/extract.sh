#!/bin/bash

# 2018-08-14
# 2019-07-30 modifying to work new monthly files
# mark potosnak

# after the big tarballs are downloaded and (partially) unzipped
# this extracts only the chemsense data from a particular node
# it preserves the header, and writes out a flat file

# node to extract (only one at a time for now)
NODE="001e06113107"
# for now, just extract data for air quality: chemsense
BOARD="chemsense"

# change into directory
cd Downloaded

# only get BOARD info from NODE
# note, could use "awk '/yellow/,0' textfile.txt" to start at a certain date
DIRS=`ls | grep '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}'`
echo $DIRS
for i in $DIRS; do
   if [ ! -f $i/node$NODE.csv ]; then
      echo "Working on $i"
      zcat $i/data.csv.gz | grep $NODE | grep $BOARD > $i/node$NODE.csv
   fi
done

# get the header, hopefully this is always the same!
zcat $i/data.csv.gz | head -1 > node$NODE.csv

# put all the date slices together
cat */node$NODE.csv >> node$NODE.csv

# move final file up one directory, so it gets saved
mv node$NODE.csv ..
