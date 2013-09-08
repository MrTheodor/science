#!/bin/bash

if [ "X$1" == "X" ] || [ "X$2" == "X" ] || [ "X$3" == "X" ]; then
	echo "Usage:"
	echo -e "\t./`basename $0` root_path experiment_name out_dir path_to_get"
	exit 1
fi

EXP_PATH="$1"
EXPERIMENT_NAME="$2"
OUT="$3"
PATH_TO_GET="$4"

find $EXP_PATH -name $EXPERIMENT_NAME


find $EXP_PATH -name $EXPERIMENT_NAME | while read a; do
	PREFIX="`echo $a | cut -f2 -d'/'`"
	echo "PREFIX: $PREFIX"
	cp -vv $a/$PATH_TO_GET "$EXPERIMENT_NAME/${PREFIX}_`basename ${PATH_TO_GET}`"
done
