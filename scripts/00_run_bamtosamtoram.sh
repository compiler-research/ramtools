#
# Script to convert a BAM to a SAM to a RAM file.
#
#
# Best is to cut and paste the desired line from the script.
#
# Author: Fons Rademakers, 26/8/17.
#
fileid=$1

bamdir=/eos/genome/local/14007a/original_BAM
samdir=/data3/rdm
ramdir=../data

bam=$bamdir/$fileid.bam
sam=$samdir/$fileid.sam
ram=$ramdir/$fileid.root
log1=$fileid-1.log
log2=$fileid-2.log

cd /data/rdm/root
source  bin/thisroot.sh
cd ../ramtools

date
../samtools-1.7/samtools view -h -o $sam $bam
date
root.exe -b "samtoram.C+(\"$sam\",\"$ram\")"
date

#nohup ../samtools-1.7/samtools view -h -o $sam $bam >& $log1 &
#nohup root.exe -b "samtoram.C+(\"$sam\",\"$ram\")" >& $log2 &

#nohup root.exe -b "samtoram.C+(\"../data/6148s10.sam\",\"../data/6148s10.root\")" >& 6148s10.log &

#nohup root.exe -b "samtoram.C+(\"../data/6148s10.sam\",\"../data/6148s10-lz4.root\",true,true,true,ROOT::kLZ4)" >& 6148s10-lz4.log &
