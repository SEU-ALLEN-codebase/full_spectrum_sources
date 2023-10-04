list=$1
cat $list | xargs -n 1 -P 100 -I {} bash -c '
    timeout 600s xvfb-run -d /home/vkzohj/public/vaa3d_for_neutube/start_vaa3d.sh -x neuTube \
        -f neutube_trace -i {} -p 1 > /dev/null
'