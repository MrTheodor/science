#!/bin/bash

DIR="$1"
OUT="$2"

for f in $DIR/*.jmol; do
    jmol -n -w "JPG:$2/$f.jpg" $f
done
