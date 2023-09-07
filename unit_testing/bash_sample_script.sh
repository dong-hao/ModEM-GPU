#!/bin/bash
#
v1=$(grep '\[e^T P m\]' level0.out | grep '1 dot product' | sed 's/.*:\(.*\)\[.*/\1/'); v1=${v1:1:16}
echo $v1
v2=$(grep '\[m^T P^T e\]' level0.out | grep '1 dot product' | sed 's/.*:\(.*\)\[.*/\1/'); v2=${v2:1:16}
echo $v2
if [[ "$v1" == "$v2" ]]; then
	echo 'string comparison PASSED'
else
	echo 'string comparison FAILED'
fi
if ! [[ "$v1" =~ *"0.000000"* ]]; then
	echo 'regexp comparison PASSED'
else
	echo 'regexp comparison FAILED'
fi
if [[ "$v1" == "$v2" ]] && ! [[ "$v1" =~ *"0.000000"* ]] && ! [[ -z "$v1" ]]; then
	echo 'logical PASSED'
else
	echo 'logical FAILED'
fi