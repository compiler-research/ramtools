#
# Script to convert a SAM to a RAM file.
#
#
# Best is to cut and paste the desired line from the script.
#
# Author: Fons Rademakers, 26/8/17.
#
nohup root.exe -b "samtoram.C+(\"/data3/rdm/6148.sam\",\"../data/6148.root\")" >& 6148.log &

#nohup root.exe -b "samtoram.C+(\"../data/6148s10.sam\",\"../data/6148s10.root\")" >& 6148s10.log &

#nohup root.exe -b "samtoram.C+(\"../data/6148s10.sam\",\"../data/6148s10-lz4.root\",true,true,true,ROOT::kLZ4)" >& 6148s10-lz4.log &
