set -x

F=$1

python tools_perf.py convert bam                   --out ${F}_samtoram ${F}.sam
python tools_perf.py convert ram    -a ROOT::kLZMA --out ${F}_samtoram ${F}.sam ${F}_lzma_split.root
python tools_perf.py convert ram -T -a ROOT::kLZMA --out ${F}_samtoram ${F}.sam ${F}_lzma_nosplit.root
python tools_perf.py convert ram    -a ROOT::kZLIB --out ${F}_samtoram ${F}.sam ${F}_zlib_nosplit.root
python tools_perf.py convert ram -T -a ROOT::kZLIB --out ${F}_samtoram ${F}.sam ${F}_zlib_split.root


python tools_perf.py index bam    --out ${F}_index ${F}.bam
python tools_perf.py index ram -S --out ${F}_index ${F}_*root