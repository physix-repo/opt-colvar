#!/bin/bash
start1=$(head -2 RESTART |tail -1 | awk '{print $1}')
end1=$(tail -1 RESTART | awk '{print $1}')
sed "s/cvmin/${start1}/g"  input.template > tmp
sed "s/cvmax/${end1}/g"  tmp > input
