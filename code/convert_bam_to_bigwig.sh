#!/bin/bash
#
## Stop if a command fails
#
 

set -euo pipefail


#### Paths for editing

BAM_DIR="/mnt/c/Users/renee/Other_projects_data/DXR_data/final_data"

OUT_DIR="/mnt/c/Users/renee/Other_projects_data/DXR_data/final_data/bigwig_ind_files"

THREADS=8

mkdir -p "$OUT_DIR"

## Loop through individuals and treatments


for IND in ind1 ind2 ind3 ind4 ind5; do
	for TREATMENT in DOX VEH; do

		SEARCH_DIR="${BAM_DIR}/${IND}/${TREATMENT}"


		if [[ ! -d "$SEARCH_DIR" ]]; then
			echo "Directory not found, skipping: $SEARCH_DIR"
			continue
		fi

		echo "Searching: $SEARCH_DIR"

		find "$SEARCH_DIR" -type f -name "*.bam" | grep -i "H3K27ac" | while read -r BAM; do

		echo "Processing:  $BAM"

		BASENAME=$(basename "$BAM" .bam)

		###clean up basename##
		
		BASENAME=${BASENAME%_final}

		BW_OUT="${OUT_DIR}/${BASENAME}.CPM.bw"

		if [[ -f "$BW_OUT" ]]; then 
			echo " Skipping, already exists: $BW_OUT"
			continue
		fi

		if [[ ! -f "${BAM}.bai" ]]; then
			echo " BAM index not found. creating index...."
			samtools index "$BAM"
		fi

		bamCoverage -b "$BAM" -o "$BW_OUT" --binSize 10 --normalizeUsing CPM --centerReads --numberOfProcessors "$THREADS"

		echo " Finished: $BW_OUT"

	done
done
done

echo "ALL H3K27ac BAMs converted."
