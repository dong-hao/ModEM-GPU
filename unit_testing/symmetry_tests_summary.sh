#!/bin/bash
#
output=$1
#
#   Generate quick summary output
#
echo "#### SYMMETRY TESTS"
#
line=$(grep '\[d^T L e\]' $output | grep '1 dot product'); echo $line
v1=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:8}
line=$(grep -B 1 '\[e^T L^T d\]' $output ); echo $line
v2=$(echo $line | grep '1 dot product' | sed 's/.*(\(.*\),.*/\1/'); v2=${v2:0:8}
if [[ "$v1" == "$v2" ]] && ! [[ "$v1" =~ "0.000000"* ]] && ! [[ -z "$v1" ]]; then
	echo "#### L: PASSED" 
else
	echo "#### L: FAILED" 
fi
#
line=$(grep '\[b^T S^{-1} b\]' $output | grep '1 dot product'); echo $line
v1=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:12}
line=$(grep '\[b^T (S^{-1})^T b\]' $output | grep '1 dot product'); echo $line
v2=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v2=${v2:1:12}
if [[ "$v1" == "$v2" ]] && ! [[ "$v1" =~ "0.000000"* ]] && ! [[ -z "$v1" ]]; then
	echo "#### S: PASSED" 
else
	echo "#### S: FAILED" 
fi
#
line=$(grep '\[e^T P m\]' $output | grep '1 dot product'); echo $line
v1=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:12}
line=$(grep '\[m^T P^T e\]' $output | grep '1 dot product'); echo $line
v2=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v2=${v2:1:12}
if [[ "$v1" == "$v2" ]] && ! [[ "$v1" =~ "0.000000"* ]] && ! [[ -z "$v1" ]]; then
	echo "#### P: PASSED" 
else
	echo "#### P: FAILED" 
fi
#
line=$(grep '\[d^T Q m\]' $output | grep '1 dot product'); echo $line
v1=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:12}
line=$(grep '\[m^T Q^T d\]' $output | grep '1 dot product'); echo $line
v2=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v2=${v2:1:12}
if [[ "$v1" == "$v2" ]] && ! [[ -z "$v1" ]]; then
	echo "#### Q: PASSED" 
else
	echo "#### Q: FAILED" 
fi
#
line=$(grep '\[d^T J m\]' $output); echo $line
v1=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:8}
line=$(grep '\[m^T J^T d\]' $output); echo $line
v2=$(echo $line | sed 's/.*:\(.*\)\[.*/\1/'); v2=${v2:1:8}
if [[ "$v1" == "$v2" ]] && ! [[ "$v1" =~ "0.000000"* ]] && ! [[ -z "$v1" ]]; then
	echo "#### J: PASSED" 
else
	echo "#### J: FAILED" 
fi
#