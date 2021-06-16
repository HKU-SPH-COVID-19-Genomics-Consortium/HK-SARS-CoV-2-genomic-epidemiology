
# Parameters and scripts for calling single nucleotide polymorphisms (SNPs).

---


```
for

 bam in `ls | grep 'sorted.bam.bai'`
do
	sample=$(echo $bam | cut -d"_" -f 1)
	bamfile=$sample"_sorted.bam"

	# begin varcall
	## freebayes
	~/softwares/freebayes/scripts/freebayes-parallel <(~/softwares/freebayes/scripts/fasta_generate_regions.py reference.fasta.fai 650) 46 -f reference.fasta -F 0.01 -C 1 -p 1 -q 30 -K --min-coverage 5 $bamfile > ../results/"freebayes_"$sample".vcf"
			
	## VarDict
	bedtools makewindows -n 46 -b reference.bed > ref_slide_window.bed
	cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 46 -k "~/softwares/VarDictJava/build/install/VarDict/bin/VarDict -G reference.fasta -f 0.01 -N samplename -th 1 -b $bamfile -R {} " | ~/softwares/VarDictJava/build/install/VarDict/bin/teststrandbias.R | ~/softwares/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N samplename -E -f 0.01 > ../results/"vardict_"$sample".vcf"

	## lofreq
	/home/hggu/softwares/lofreq/src/lofreq/lofreq call-parallel --pp-threads 46 -f reference.fasta -o ../results/"lofreq_"$sample".vcf" $bamfile

	# mv data to achive
	mv $bamfile ../tmp2_lofreq/$bamfile
	mv $bam ../tmp2_lofreq/$bam
	mv $bam".bai" ../tmp2_lofreq/

done
```