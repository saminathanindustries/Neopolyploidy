# Bash script to run Salmon

# Paths
nviSalmonIdx='/path/to/salmonIndex'
nviMap='/path/to/MapGenesTranscripts.txt'
fastqFiles='/path/to/WGDL'

# Sample list (when changing, change the above FASTQ files PATH too)
samples='C21 C22 C23 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C51 C52 C53 C54 C55B C56 C57 C58 C59 C60 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C61 C62 C63 C64 C66 C67 C68 C69 C71 C72 C73 C74 C76 C77 C78 C79 C81 C82 C83 C84 C86 C87 C88 C89 C91 C92 C93 C94 C96 C97 C99'

# Read the sample list and run salmon on each sample
for sample in $samples; do
	echo "Currently analyzing the sample: $sample"
	fullName1=(${fastqFiles}/${sample}"_R1.fq.gz")
	#echo "$fullName1"
    fullName2=(${fastqFiles}/${sample}"_R2.fq.gz")
	#echo "$fullName2"
    echo salmon quant -i $nviSalmonIdx --libType A -1 ${fullName1} -2 ${fullName2} -o $sample -p 16 --geneMap $nviMap --dumpEq
	salmon quant -i $nviSalmonIdx --libType A -1 ${fullName1} -2 ${fullName2} -o $sample -p 16 --geneMap $nviMap --dumpEq
done

