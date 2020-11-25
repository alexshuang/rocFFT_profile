#!/bin/sh

SRC=$1

sed -i 's;real_type_t<T>;float ;g' $SRC
sed -i 's;\<T\>;float2;g' $SRC
sed -i 's;<float2, .*>;;g' $SRC
sed -i 's;lib_make_vector2<\(.*\)>;\1;g' $SRC



