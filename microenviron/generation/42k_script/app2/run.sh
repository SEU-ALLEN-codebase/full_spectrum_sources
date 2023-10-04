VAA3D=/PBshare/SEU-ALLEN/Users/zuohan/vaa3d/start_vaa3d.sh
timeLim=300s

[ -f all.txt ] || find ../stage1/out -name *pbd > all.txt

find out -name *swc | awk '{print "../stage1/"$0}' | sed 's/.swc//' > done.txt
sort all.txt done.txt done.txt | uniq -u > unfinished.txt

cat unfinished.txt | xargs -n 1 -P 50 -I {} bash -c "
    img={}
    out=\${img#../stage1/}.swc
    mkdir -p \`dirname \$out\`
    timeout $timeLim xvfb-run -a $VAA3D -x vn2 -f app2 -i \$img -o \$out -p decoy.marker 0 AUTO 0 1 0 0 5 1 || rm \$out
    "