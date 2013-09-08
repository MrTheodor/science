#!/bin/bash

FILE=$1
COL=$2

TEMPLATE="
set term x11;
set datafile separator \";\";
plot '$FILE' using 1:$COL with lp
"

echo $TEMPLATE | gnuplot --persist
