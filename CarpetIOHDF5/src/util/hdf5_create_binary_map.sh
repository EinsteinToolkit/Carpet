#!/bin/bash

# usage: hdf5_create_map.sh a.file_0.h5 a.file_1.h5 b.file_0.h5
# writes output to a.map b.map

rm -f ${*//.file_*.h5/.map}
while [ $# -gt 0 ] ; do
  h5ls $1 | gawk '{gsub("\\\\","");sub(" +(Dataset|Group).*$","");slurp[++n]=$0}END{for(s in slurp) print slurp[s]}' | ${0%.sh} /dev/stdin $1 ${1%.file_*.h5}.map
  shift
done
