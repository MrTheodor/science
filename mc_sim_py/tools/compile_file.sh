##
## compile file

INPUT_PY="$1"
INPUT_C="`echo $1 | cut -f1 -d'.'`.c"
OUTPUT_SO="`echo $1 | cut -f1 -d'.'`.so"

echo "INPUT_PY=$INPUT_PY"
echo "INPUT_C=$INPUT_C"
echo "OUTPUT_SO=$OUTPUT_SO"

rm $INPUT_C
rm $OUTPUT_SO

cython $INPUT_PY -o $INPUT_C
gcc -shared -pthread -fPIC -fwrapv -O3 -Wall -fno-strict-aliasing -I${PYTHON_INCLUDE} -o $OUTPUT_SO $INPUT_C
