bamfile=$1
reference=$2
REF_SNPs_file=$3
nameout=$(basename $bamfile | cut -f1 -d ".")

while read line; do
	contig=$(echo $line | cut -f2 -d " ")
	pos=$(echo $line | cut -f3 -d " ")
	res=$(samtools mpileup --reference $reference -uD -r $contig\:$pos\-$pos $bamfile | bcftools call -c --ploidy 1 | bcftools view -H | cut -f 4-5)
	echo $res | awk '{if(NF == 2){print $0}else{print "- ."}}' | awk '{if($2 == "."){print $1}else{print $2}}' | awk '{if(length($1) == 1){print $1}}'
	done < $REF_SNPs_file 2> /dev/null | tr "\n" " " | sed 's/ //g' | sed 's/$/\n/g' | sed -e "s/^/\>$nameout\t/g" | tr "\t" "\n"
