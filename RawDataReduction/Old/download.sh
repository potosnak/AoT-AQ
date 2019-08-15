#!/bin/bash

# 2018-08-15
# mark potosnak

# NOTE! this is obsolete, but use as reference for downloading by hand

# this script downloads the latest AoT data, extracts the compressed flat file
# follows commands suggested on github and deletes extra files to save space
# everything happens in the local directory, where the script is called

# make sure directory for downloading files exists
# note, this doesn't get synced by github, since it's big
# the dash p won't complain if directory already exists
mkdir -p BigTar

cd BigTar

# remove any previously downloaded file, in case not previously deleted
rm AoT_Chicago.complete.latest.tar

# remove any previously downloaded directories, again to save space
rm -rf `ls | grep AoT_Chicago.complete.'[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}$'`

# get the data
wget http://www.mcs.anl.gov/research/projects/waggle/downloads/datasets/AoT_Chicago.complete.latest.tar

# extract the tar file first
tar -xvf AoT_Chicago.complete.latest.tar

# remove the downloaded file to save space
rm AoT_Chicago.complete.latest.tar

