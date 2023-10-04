ls list/* | while read list; do
    # qsub -l walltime=48:00:00 -l nodes=1:ppn=100 -N ntb -j oe -o /dev/null -d . -V run.sh
    bash run.sh $list
done
find /home/vkzohj/data/enhanced_soma_images/out -name *_neutube.swc -exec cp -p --parents -r -u {} ~/public/projects/230k/ntb \;
mv home/vkzohj/data/enhanced_soma_images/out .
rm -r home