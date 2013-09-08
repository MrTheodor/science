#rdir="Exp8F_01/clusters/exp8f0_1_0/c1_240_240_240_0/models/"
#out="Exp8F_01/clusters/clusters_exp8f0_1_0_c1_240_240_240_0"
echo "input_dir" "output_file" "start_range:stop_range"
rdir="$1"
out="$2"
start_range="`echo $3 | cut -f1 -d:`"
stop_range="`echo $3 | cut -f2 -d:`"
cmd="print ' '.join(map(str, range($start_range, $stop_range+1000, 1000)))"
range="`python -c \"$cmd\"`"
echo $rdir $out $start_range $stop_range
echo $cmd
for f in $range; do
    python tools/calculate_clusters.py -f "${rdir}/jmol_${f}.jmol" -step -box 240,240,240 >> $out; 
done
