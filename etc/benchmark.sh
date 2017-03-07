#!/bin/zsh

function die {
echo ""
echo "$1"
echo -n $'\a'
exit 1
}

[[ -r datasets.csv ]] || die "No file datasets.csv found!"
[[ -x ./naive/integer_partition_demo_naive ]] || die "No executable ./naive/integer_partition_demo_naive found"
[[ -x ./demo/integer_partition_demo ]] || die "No executable ./demo/integer_partition_demo found"


function printZeit {
echo -n "$1 ;"
echo $zeit | \
perl -M'List::Util qw(sum max min)' -MPOSIX -0777 -a -ne 'printf "%-7s : %d; "x4, "Min", min(@F), "Max", max(@F), "Average", sum(@F)/@F,  "Median", sum( (sort {$a<=>$b} @F)[ int( $#F/2 ), ceil( $#F/2 ) ] )/2;'
echo ''
}

while read line; do
	grep -q '^#' <<< "$line" && continue
	cols=$(grep -o ';' <<< "$line" | wc -l)
	[[ cols -lt 3 ]] && continue
	typ=$(cut -d';' -f1 <<< "$line")
	no=$(cut -d';' -f2 <<< "$line")
	z=$(cut -d';' -f$(expr $cols + 1) <<< "$line")
	bounds=$(cut -d';' -f3-$cols <<< "$line" | sed 's@;@ @g')
	echo no = $no : typ = $typ : z = $z : bounds = $bounds
	outexact=0
	if [[ $typ != 'F' ]]; then
		zeit=()
		for it in $(seq 1 10); do
			timeA=$(date +%s%3N)
			outexact=$(timeout 60m ./naive/integer_partition_demo_naive $z $(echo $bounds))
			[[ $? -eq 124 ]] && outexact=0
			timeB=$(date +%s%3N)
			zeit+=$(expr $timeB - $timeA)
			[[ "$outexact" = "0" ]] && break
		done
		printZeit "Naive"
#		echo $outs
	fi
#	zeit=()
#	outs=()
#	for it in $(seq 1 10); do
#		timeA=$(date +%s%3N)
#		outs+=$(./demo/integer_partition_demo $z $(echo $bounds))
#		timeB=$(date +%s%3N)
#		zeit+=$(expr $timeB - $timeA)
#	done
#	printZeit "FaulS"
##	echo $outs
#	outs=()
for proc in $(seq 1 $(nproc)); do
		zeit=()
		for it in $(seq 1 10); do
			timeA=$(date +%s%3N)
			out=$(./demo/integer_partition_demo -threads $proc $z $(echo $bounds))
			timeB=$(date +%s%3N)
			zeit+=$(expr $timeB - $timeA)
			if [[ $outexact -gt 0 ]] && [[ $out != $outexact ]]; then
				echo "exact $outexact differs from computed value $out"
				echo "Call was ./demo/integer_partition_demo -threads $proc $z $(echo $bounds)"
			fi
		done
		printZeit "Faul$proc"
	done
	echo $outs
	
done < datasets.csv
