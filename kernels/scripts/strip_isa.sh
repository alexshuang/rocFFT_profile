#!/bin/sh

SRC=$1

sed -i 's;\(.*\S\)\s*// 000.*;\1;g' $SRC
