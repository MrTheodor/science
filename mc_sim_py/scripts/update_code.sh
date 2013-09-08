#!/bin/bash
CURRENT_PATH="`pwd`"
CODE="$HOME/code/intraglobular/src"
cd $CODE
git pull
find -name "*.so" | while read a; do ./compile_file.sh "$a"; done
cd $CURRENT_PATH
