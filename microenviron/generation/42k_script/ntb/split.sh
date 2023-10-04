[ -f all.txt ] || find ~/data/enhanced_soma_images/out -name *pbd > all.txt
find ~/data/enhanced_soma_images/out -name *_neutube.swc | sed 's/_neutube.swc//' > done.txt
sort all.txt done.txt done.txt | uniq -u | shuf > unfinished.txt
wc -l unfinished.txt
mkdir -p list
split -n 3 unfinished.txt list/