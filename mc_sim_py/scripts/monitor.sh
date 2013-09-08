#!/bin/bash
NAME=$2
ROOT=$1
tail -n1 $ROOT*/${NAME}*/calc/calc_*.csv
for f in $ROOT*/${NAME}*/*.kch; do
    echo $f;
    python -c "import kyotocabinet as kb; db = kb.DB(); db.open('${f}', db.OREADER); print '\n'.join(db.match_prefix(''));" | tail -n 1
done
