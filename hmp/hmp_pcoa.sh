#!/usr/bin/env bash 

# Compute the princpal coordinates for the HMP data set. Run beta_diversity.py 
# on the HMP data. The output of the analysis is the King data files that can 
# be used to visualize the results from PCoA. 

metrics=( hellinger unweighted_unifrac unifrac bray_curtis chisq euclidean \
	weighted_unifrac ) 

for metric in ${metrics[@]}; do 
	echo "Submitting $metric"
	beta_diversity.py -i ~/Git/data/hmp_biom.json \
		-o ~/Git/data/hmp.pcoa \
		-t ~/Git/data/hmp_rep_set.tre \
		-m $metric & 
done
wait
principal_coordinates.py -i ~/Git/data/hmp.pcoa/ -o ~/Git/data/hmp.beta/

pcoa_files=`find ~/Git/data/hmp.beta/*.txt`
for fname in ${pcoa_files[@]}; do 
	make_3d_plots.py -i $fname -m ~/Git/data/hmp.map -o ~/Git/data/hmp.plots/ & 
	make_2d_plots.py -i $fname -m ~/Git/data/hmp.map -o ~/Git/data/hmp.plots/ & 
done
wait 
