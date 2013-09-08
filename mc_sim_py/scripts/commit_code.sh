#!/bin/bash

COMMENT="$1"

git commit -am "$COMMENT"
git push
