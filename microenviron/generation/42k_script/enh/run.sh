img_dir=/home/vkzohj/public/projects/230k_valid/crop/out
mkdir -p out
find $img_dir -name *v3dpbd > all.txt

find out -name *pbd | sed 's/out//' | awk "{print \"$img_dir/\"\$0}" > done.txt
sort all.txt done.txt done.txt | uniq -u > unfinished.txt

VAA3D=/PBshare/SEU-ALLEN/Users/zuohan/vaa3d/start_vaa3d.sh

cat unfinished.txt | xargs -n 1 -P 20 -I {} bash -c "
    img={}
    out=out/\${img#$img_dir/}
    mkdir -p \`dirname \$out\`
    xvfb-run -d $VAA3D -x imPreProcess -f im_enhancement -i \$img -o \$out > /dev/null
    "