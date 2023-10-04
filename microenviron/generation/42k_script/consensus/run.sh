trap '
trap "
kill -9 0
wait
" INT TERM KILL
mail -s "lyf workstation jobs done" vkzohj@seu.edu.cn, yufeng_liu@seu.edu.cn <<< "consensus on 230k enhanced 512x512x256 blocks done. 
Result located at $(pwd)/out. There are $(find out -name *swc | wc -l) swc."
' EXIT

# interrupt
trap '
trap " " INT TERM EXIT KILL
kill -9 0
wait
sleep 10s
mail -s "lyf workstation jobs failed" vkzohj@seu.edu.cn <<< "consensus on 1891 enhanced 512x512x256 blocks failed."
find out -name *swc -cmin -10 -delete
' INT TERM KILL

source /PBshare/SEU-ALLEN/Users/zuohan/pylib/venv/bin/activate
python run.py