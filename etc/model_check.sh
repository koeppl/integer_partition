#!/bin/zsh

function die {
echo ""
echo "$1"
echo -n $'\a'
exit 1
}


[[ -x ./naive/integer_partition_demo_naive ]] || die "No executable ./naive/integer_partition_demo_naive found"
[[ -x ./demo/integer_partition_demo ]] || die "No executable ./demo/integer_partition_demo found"

wheel="_/\\|"
wheeli=1
it=1

while [ 1 ]; do 
	dim=$((RANDOM % 9 +1))
	maxz=0
	bounds=()
	for i in $(seq 1 $dim); do
		val=$((RANDOM % 200 +1))
		bounds+=$val
		((maxz=maxz+val))
	done
	z=$((RANDOM % maxz))
	# echo $z $bounds 

timeA=$(date +%s%N | cut -b1-13)
outA=$(./naive/integer_partition_demo_naive $z $bounds)
timeB=$(date +%s%N | cut -b1-13)
outB=$(./demo/integer_partition_demo $z $bounds)
timeC=$(date +%s%N | cut -b1-13)
echo -e -n "\r                         \r$t "
if [[ $? -ne 0 ]]; then
	it=0
echo -n $'\a'
	echo -n ':-( '
fi
echo -n "$wheel[$wheeli] $it $(expr $timeB - $timeA)ms $(expr $timeC - $timeB)ms $outA $outB" 
[[ $outA -ne $outB ]] && die "Differ at $z $bounds"
 ((wheeli = (wheeli % 4) + 1))
 ((++it))
 done

