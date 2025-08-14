set -x
F=$1
A=$2


sudo -v

# Keep-alive: update existing sudo time stamp if set, otherwise do nothing.
while true; do sudo -n true; sleep 60; kill -0 "$$" || exit; done 2>/dev/null &

python tools_perf.py view bam ${F}.bam                       ${F}_views.csv --out ${F}_out/out_${A}_sam

python tools_perf.py view ram ${F}_lzma_split.root           ${F}_views.csv --out ${F}_out/out_${A}_ram_lzma_split
python tools_perf.py view ram ${F}_lzma_split.root           ${F}_views.csv --out ${F}_out/out_${A}_ram_lzma_split_index  --macro ramview_index.C

python tools_perf.py view ram ${F}_lzma_nosplit.root         ${F}_views.csv --out ${F}_out/out_${A}_ram_lzma_nosplit
python tools_perf.py view ram ${F}_lzma_nosplit.root         ${F}_views.csv --out ${F}_out/out_${A}_ram_lzma_nosplit_index  --macro ramview_index.C


python tools_perf.py view ram ${F}_zlib_split.root           ${F}_views.csv --out ${F}_out/out_${A}_ram_zlib_split
python tools_perf.py view ram ${F}_zlib_split.root           ${F}_views.csv --out ${F}_out/out_${A}_ram_zlib_split_index  --macro ramview_index.C

python tools_perf.py view ram ${F}_zlib_nosplit.root         ${F}_views.csv --out ${F}_out/out_${A}_ram_zlib_nosplit
python tools_perf.py view ram ${F}_zlib_nosplit.root         ${F}_views.csv --out ${F}_out/out_${A}_ram_zlib_nosplit_index  --macro ramview_index.C



python tools_perf.py view ram ${F}_lzma_split.root           ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_lzma_split_cache
python tools_perf.py view ram ${F}_lzma_split.root           ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_lzma_split_index_cache  --macro ramview_index.C

python tools_perf.py view ram ${F}_lzma_nosplit.root         ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_lzma_nosplit_cache
python tools_perf.py view ram ${F}_lzma_nosplit.root         ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_lzma_nosplit_index_cache  --macro ramview_index.C


python tools_perf.py view ram ${F}_zlib_split.root           ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_zlib_split_cache
python tools_perf.py view ram ${F}_zlib_split.root           ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_zlib_split_index_cache  --macro ramview_index.C

python tools_perf.py view ram ${F}_zlib_nosplit.root         ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_zlib_nosplit_cache
python tools_perf.py view ram ${F}_zlib_nosplit.root         ${F}_views.csv --cache --out ${F}_out/out_${A}_ram_zlib_nosplit_index_cache  --macro ramview_index.C

# cgexec -g memory:limitmem2 python tools_perf.py view bam ${F}.bam                       ${F}_views.csv --out ${F}_out_mem_${A}limit/out_sam

sudo -K