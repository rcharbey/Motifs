#!/bin/sh

data="../../../Graphs/Data/Edgelists"
results_file="../../Results/enum_all_kmotifs.csv"
list_edgelists=`ls $data`

echo "ego;3;4;5;6;7" > $results_file

for edgelist in $list_edgelists
do
	ego=`echo $edgelist | cut -d '.' -f 1`
	echo $ego
	result="$ego"
	for k in 3 4 5 6 7 8
	do
		temp=`./kmotif $k $data/$edgelist`
		result="$result;$temp"
		
	done
	echo $result >> $results_file 
done
