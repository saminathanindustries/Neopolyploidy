# Neopolyploidy
# Genome files (Nasonia vitripennis; Nvit_psr_1.1; TaxonomyID 7425) from two sources:
    # NCBI - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009193385.2/

# Lists of samples 
samplesTRA='A1A A2A A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15A A16 A18 A19 A20 A21 A22 A23 A24 A25 A26 A27 A28 A31 A32 A33 A34 A35 A36 A37 A38 A39 A40 A41 A42 A43 A44 A45 A46 A47 A48 A49 A50 A51 A52 A53 A54 A55 A56 A57 A58 A59 A60 A61 A62 A63 A64 A65 A66 A67 A68 A69 A70 A71 A72 A73 A74 A75 A76 A77 A78 A79 A80 A81 A82 A83 A84 A85 A81A A86 A87 A88 A89 A90'
samplesWGDL='C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C21 C22 C23 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C51 C52 C53 C54 C55B C56 C57 C58 C59 C60 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50'
samplesWOM='D1 D1A D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 D13 D14 D15 D16 D17 D18 D19 D20 D21 D22 D23 D24 D25 D26 D27 D28 D29 D30 D31 D32 D33 D34 D35 D36 D37 D38 D39 D40 D41 D42 D43 D44 D45 D46 D47 D48 D49 D50 D51 D56 D61 D63 D64 D66 D67 D68 D69 D81 D82B D83B D84 D85 D81B D84A D84B D85A D85B D86 D87 D88 D89 D90A D86B D87B D88B D90B D90C D91 D92 D93 D94 D95 D96 D97 D98 D99 D100 D101 D102 D103 D104 D105 D106 D107 D108 D109 D110 D110B D111 D112 D113 D114 D115 D116 D117 D118 D119 D120'
samplesWPL='B1A B2A B3 B4A B5D B6A B7 B8 B9 B10 B11 B12A B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26A B27 B28A B30 B31 B32 B33 B34A B35A B36A B37A B38 B39A B40A'

# Merge raw data FASTQ samples with multiple files from multiple lanes (for each set)
for sample in $samplesWPL; do
    cd $sample
    if [[ $(ls -1 *.fq.gz | wc -l) != 2 ]]; then
        cat $sample*"_1.fq.gz" > ../$sample"_R1.fq.gz"
        cat $sample*"_2.fq.gz" > ../$sample"_R2.fq.gz"
    else
        mv $sample*"_1.fq.gz" ../$sample"_R1.fq.gz"
        mv $sample*"_2.fq.gz" ../$sample"_R2.fq.gz"
    fi
    cd ..
done

# FastQC (v0.12.1) for raw data
    $ allSets='TRA WGDL WOM WPL'
    for eachSet in $allSets; do
        cd $eachSet
        fastqc -t 18 *.fq.gz
        mv *fastqc* $eachSet
        cd ..
    done

# MultiQC (v1.21) for the raw data
    for eachSet in $allSets; do
        multiqc $eachSet -o FastqcRaw --title $eachSet
    done

# FastQC summary: Skip first 15 bases in the mapping; There are atleast 20M reads in all samples, some even have 200M, will downsample later if needed. No adapter contamination and no other trimming needed. 

# Salmon (v1.10.0) index the genome (NCBI - GCF_009193385.2; Nvit_psr_1.1)
    $ samtools faidx GCF_009193385.2_Nvit_psr_1.1_genomic.fna
    $ grep "^>" GCF_009193385.2_Nvit_psr_1.1_genomic.fna | cut -d " " -f1 | sed -e 's/>//g' > decoys.txt
    $ cat rna.fna GCF_009193385.2_Nvit_psr_1.1_genomic.fna | gzip > gentrome.fa.gz
    $ salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i salmonIndex

# Extracting identifiers from genome GTF file to make gene-transcript map
    $ awk 'BEGIN{FS=OFS="\t"}{split($9,a,";"); print a[2]"\t"a[1]}' genomic.gtf | grep -v "db_xref" | sed 's/gene_id "//' | sed 's/transcript_id "//' | sed 's/"//g' | sed 's/ //g' | sort | uniq | sed 1d > MapGenesTranscripts.txt

# Updated lists after removing unnecessary ones (TRA:77; WGDL:71; WOM:80; WPL:39; WGDL:70)
samplesTRA='A1A A2A A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15A A16 A18 A19 A20 A21 A22 A23 A24 A25 A26 A27 A28 A31 A32 A33 A34 A35 A36 A37 A38 A39 A40 A51 A52 A53 A54 A55 A56 A57 A58 A59 A60 A61 A62 A63 A64 A65 A66 A67 A68 A69 A70 A71 A72 A73 A74 A75 A76 A77 A78 A79 A80 A81 A82 A83 A84 A85 A86 A87 A88 A89 A90'
samplesWOM='D1 D1A D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 D13 D14 D15 D16 D17 D18 D19 D20 D21 D22 D23 D24 D25 D26 D27 D28 D29 D30 D31 D32 D33 D34 D35 D36 D37 D38 D39 D40 D81 D82B D83B D84 D85 D86 D87 D88 D89 D90A D91 D92 D93 D94 D95 D96 D97 D98 D99 D100 D101 D102 D103 D104 D105 D106 D107 D108 D109 D110 D111 D112 D113 D114 D115 D116 D117 D118 D119 D120'
samplesWPL='B1A B2A B3 B4A B5D B6A B7 B8 B9 B10 B11 B12A B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26A B27 B28A B30 B31 B32 B33 B34A B35A B36A B37A B38 B39A B40A'
samplesWGDL='C21 C22 C23 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C51 C52 C53 C54 C55B C56 C57 C58 C59 C60 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C61 C62 C63 C64 C66 C67 C68 C69 C71 C72 C73 C74 C76 C77 C78 C79 C81 C82 C83 C84 C86 C87 C88 C89 C91 C92 C93 C94 C96 C97 C98 C99'

# Run Salmon on each dataset
    $ bash runSalmon.sh
    
# MultiQC for Salmon results
    $ multiqc TRA/ -o Salmon --title TRA
    $ multiqc WOM/ -o Salmon --title WOM
    $ multiqc WPL/ -o Salmon --title WPL
    $ multiqc WGDL/ -o Salmon --title WGDL

# Gene Ontology (GO) data files (accessed on 12-May-2024):
    # NCBI - https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    $ zcat gene2go.gz | awk '$1 == 7425' > gene2go_7425.tsv
    # Added header manually to "gene2go_7425.tsv" from the original file 
    # Extract the GeneID conversion table from the GTF file
    $ awk 'BEGIN{FS=OFS="\t"}{split($9,a,";"); print a[2]"\t"a[1]}' genomic.gtf | grep "db_xref" | sed 's/gene_id "//' | sed 's/ db_xref "GeneID://' | sed 's/"//g' | sed 's/ //g' > GeneID_GeneNames.tsv
    # Added 'GeneNames' from the GTF file to the "gene2go_7425.tsv" in Excel and saved as "gene2go_7425_updated.tsv"




