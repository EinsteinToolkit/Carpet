#! /bin/bash

awk '{ if (minval[$5]+0==0 || $9<minval[$5]) { minval[$5]=$9; line[$5]=$0; } }
END { for (p in line) print line[p]; }' |
sort -n -k 5
