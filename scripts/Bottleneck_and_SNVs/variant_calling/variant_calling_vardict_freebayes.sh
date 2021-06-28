#!/bin/bash
# conda activate ngs
cd ~/work/2020-09-01_COVID_NGS_pipeline/tmp2
cp ../NGS_data_input/reference.fasta ./
~/softwares/samtools-1.11/bin/samtools faidx reference.fasta
~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index reference.fasta
~/softwares/bwa/bwa index reference.fasta
java -jar ~/softwares/picard.jar CreateSequenceDictionary -R reference.fasta -O reference.dict
~/softwares/bioawk/bioawk -c fastx '{print $name"\t1\t"length($seq)"\t"$name}' reference.fasta > reference.bed

while true
do			
	# work for samples in order
	for bam in `ls | grep 'sorted.bam.bai'`
	do
		sample=$(echo $bam | cut -d"_" -f 1)
		bamfile=$sample"_sorted.bam"
		# make sure file1 and file2 exist and are not opening
		if [[ -f "$bamfile" ]] && ! [[ `find $bamfile -mmin -1` ]]; then
		
			# begin varcall
			## freebayes
			~/softwares/freebayes/scripts/freebayes-parallel <(~/softwares/freebayes/scripts/fasta_generate_regions.py reference.fasta.fai 650) 46 -f reference.fasta -F 0.01 -C 1 -p 1 -q 30 -K --min-coverage 5 $bamfile > ../results/"freebayes_"$sample".vcf"
					
			## VarDict
			bedtools makewindows -n 46 -b reference.bed > ref_slide_window.bed
			cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 46 -k "~/softwares/VarDictJava/build/install/VarDict/bin/VarDict -G reference.fasta -f 0.01 -N samplename -th 1 -b $bamfile -R {} " | ~/softwares/VarDictJava/build/install/VarDict/bin/teststrandbias.R | ~/softwares/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N samplename -E -f 0.01 > ../results/"vardict_"$sample".vcf"

			## lofreq
			# /home/hggu/softwares/lofreq/src/lofreq/lofreq call-parallel --pp-threads 46 -f reference.fasta -o ../results/"lofreq_"$sample".vcf" $bamfile

			# mv data to achive
			mv $bamfile ../tmp2_lofreq/$bamfile
			mv $bam ../tmp2_lofreq/$bam
			mv $bam".bai" ../tmp2_lofreq/

			# update individual reports
			Rscript ../scripts/mutations_idv.R $sample

			#cp snvs comparing
			cp /home/hggu/work/2020-09-01_COVID_NGS_pipeline/scripts/snvs_compare.Rmarkdown /srv/shiny-server/index.rmd

		else
			if [[ -f "$bamfile" ]]; then 
				echo $bamfile exist
				if [[ `find $bamfile -mmin -1` ]]; then 
					echo but $bamfile is still busy
				fi
			fi
			
		fi
	done
	
	if `! ls | grep "bam"`; then
		echo No new data
	fi

	sleep 6

done
