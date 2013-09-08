ROOT_DIR=$1

for d in $ROOT_DIR/*; do
	ln -vs "`pwd`/tools/SkAVG" "`pwd`/$d/"
	ln -vs "`pwd`/tools/config_2.txt" "`pwd`/$d/config.txt"
	ln -vs "`pwd`/tools/calc_rdf_avg" "`pwd`/$d/calc_rdf"
done
