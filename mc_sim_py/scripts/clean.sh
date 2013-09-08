#!/bin/bash

EXP_NAME=$1

if [ "X$EXP_NAME" != "X" ]; then
	rm *$EXP_NAME
fi

find code -name *.pyc -delete
