trap '
trap "
kill -9 0
wait
" INT TERM KILL
mail -s "lyf workstation jobs done" vkzohj@seu.edu.cn <<< "prune on 230k enhanced 512x512x256 blocks done. There are $(find out -name *swc | wc -l) swc."
cd ../stage5_consensus_weak
bash run.sh
' EXIT

# interrupt
trap '
trap " " INT TERM EXIT KILL
kill -9 0
wait
sleep 10s
mail -s "lyf workstation jobs failed" vkzohj@seu.edu.cn <<< "prune on 1891 enhanced 512x512x256 blocks failed."
find out -name *swc -cmin -10 -delete
' INT TERM KILL


source /PBshare/SEU-ALLEN/Users/zuohan/pylib/venv/bin/activate
python seg_prune.py