#!/usr/bin/env bash 

metrics=( hellinger unweighted_unifrac unifrac bray_curtis chisq euclidean weighted_unifrac ) 

for metric in ${metrics[@]}; do 
	echo "Submitting $metric"
	beta_diversity.py -i ~/Git/data/hmp_biom.json \
		-o ~/Git/data/hmp.pcoa \
		-t ~/Git/data/hmp_rep_set.tre \
		-m $metric & 
done
