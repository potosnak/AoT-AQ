Run these scripts to download and extract chemsense data

1) download monthly tar files 
   a) https://aot-file-browser.plenar.io/data-sets/chicago-complete
   b) put into Downloaded directory
   c) extract with tar -xvf 
   d) remove original downloaded file to save space
2) run extract.sh, which gets out chemsense data from a specified node
3) run hourlyGAP from R, which removes first 10 min of data
   after a gap of 5 min or longer. This removes data before the 
   instrument has warmed up, which can be noisy/high. Both scripts
   produce a file of the same name, which is used by biggie.r
