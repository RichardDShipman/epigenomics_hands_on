# Epigenomics Tasks

Richard Shipman MAR2025

# Download and setup

```sh
git clone https://github.com/bborsari/epigenomics_uvic
cd epigenomics_uvic
```

# Run the docker container

```sh
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

# 4_EN‐TEx_ATAC‐seq_data_downstream_analyses

# Tasks

1. Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

# Move to the ATAC-seq folder and create folder structure for analysis.

```sh
cd ATAC-seq
mkdir 
mkdir bigBed_files peaks_analysis bed_files annotations
mkdir analysis
mkdir analysis/peaks.analysis
```

2. Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

# Download metadata data files from ENCODE.

```sh
 ../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
```

# Meta Data file structure

```sh
# Top of metadata file, head
head -1 metadata.tsv
# Show info
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

# Use metadata file to gather correct files for analysis
# bigBed narrow
# pseudoreplicated peaks
# assembly GRCh38

```sh
grep -F ATAC-seq metadata.tsv | 
grep -F "bigBed_narrowPeak" |
grep -F "pseudoreplicated_peaks" |
grep -F "GRCh38" |
awk 'BEGIN{FS=OFS="\t"}{print $1}' |
sort -k2,2 -k1,1r > analyses/bigBed.peaks.ids.txt  
```

# Download the bigBed files from the text from above

Download the two bigBed files.

```sh
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

# md5sum data integretity 

Find column in metadata with md5sum hash.

# Which field of the metadata table contains the MD5 hash? 

```sh
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'  
```

md5sum	46

# Cut md5sum hash 46 row from bigBed files

```sh
../bin/selectRows.sh <(cut -f1 analyses/bigBed.*.ids.txt) metadata.tsv | cut -f1,46 > data/bigBed.files/md5sum.txt

# Check the integrity of the downloaded files

```sh
for file_type in bigBed; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

# Check if columns equal eachother

```sh
awk '$2 != $3 {print "true"}' data/bigBed.files/md5sum.txt
```

ENCFF762IFP	f6a97407b6ba4697108e74451fb3eaf4	f6a97407b6ba4697108e74451fb3eaf4

ENCFF287UHP	46f2ae76779da5be7de09b63d5c2ceb9	46f2ae76779da5be7de09b63d5c2ceb9

true
true

# Convert bigBed files to BED files

```sh
cut -f1 analyses/bigBed.peaks.ids.txt | while read filename; do
    bigBedToBed data/bigBed.files/"$filename".bigBed bed_files/"$filename".bed
done
```

# List bedfiles

```sh
ls bed_files/
```

ENCFF287UHP.bed  ENCFF762IFP.bed

3. For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). Hint: have a look at what we did here and here.

# Download annotation files into the annotation directory

Change directory, download and upzip the genome annoation file.

```sh
cd annotations
wget -P annotations "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"
unzip gencode.v24.primary_assembly.annotation.gtf
```

# Define promotor regions around genes

2000 nt plus or minus from promotor regions - transcription start sites TTS

```sh
awk '$3=="gene"' gencode.v24.primary_assembly.annotation.gtf |
grep -F "protein_coding" |
awk 'BEGIN {FS="\t"; OFS="\t"} {
    if($7=="+") { start = $4 - 2000; end = $4 + 2000 }
    else if($7=="-") { start = $5 - 2000; end = $5 + 2000 }
    else { start = $4 - 2000; end = $4 + 2000 }
    print $1, (start<0?0:start), end, $10, 0, $7 
}' | sed 's/\"//g' | awk 'BEGIN {FS="\t"; OFS="\t"} $1!="chrM" {print $0}' 
    > promoter_regions.bed
```

# Set up gene body coordinate file

```sh
grep -F "protein_coding" | 
awk 'BEGIN {FS="\t"; OFS="\t"} { print $1, $4-1, $5, $10, 0, $7 }' |
sed 's/\"//g' | 
awk 'BEGIN {FS="\t"; OFS="\t"} $1!="chrM" {print $0}' > gene_body_coordinates.bed
```

# Change directory to ATAC-seq

```sh
cd ..
```

# Run intersect analysis

```sh
cut -f-2 analyses/bigBed.peaks.ids.txt |
while read filename tissue; do 
  bedtools intersect -a annotations/gene_body_coordinates.bed -b annotations/gencode.v24.protein.coding.non.redundant.TSS.bed -u |
  cut -f7 |
  sort -u > peaks.analysis/genes.with.peaks."$tissue".H3K4me3.txt
done
```

1) the number of peaks that intersect promoter regions:

# change directory to output loaction

```sh
cd analyses/peak.analysis
```

# Check word count

```sh
wc -l peaks.overlap.promoters.stomach.ATAC.txt 

```

44749 

There are a reported 44,749 peaks identified in specific promoter regions associated with stomach tissue.

```sh
wc -l peaks.overlap.promoters.sigmoid_colon.ATAC.txt 
```

47871

There are a reported 47,871 peaks in specific promoter regions associated with sigmoid colon tissue.

2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions):

# Move up two directory

```sh
cd ../..
```
# Number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions):

```sh
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed -v > analyses/peaks.analysis/peaks.out.gene."$tissue".H3K4me3.txt
done
```

34537

# 5_Distal_regulatory_activity

## Study distal regulatory activity

From section 4., you should have obtained a set of ATAC-seq peaks in stomach and sigmoid_colon that lie outside gene coordinates.

We will use these peaks as a starting point to build a catalogue of distal regulatory regions.

## Tasks

Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.

# Make regulatory_elements directory

```sh
mkdir regulatory_elements
```

Task 2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

```sh

```

# How many are they?

Task 3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

```sh

```

Task 4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

```sh
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}'
```

```sh

```

Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

```sh
python ../bin/get.distance.py -h
```

This script takes as input two distinct arguments: 1) --input corresponds to the file gene.starts.tsv (i.e. the file you generated in Task #4); 2) --start corresponds to the 5' coordinate of a regulatory element. Complete the python script so that for a given coordinate --start the script returns the closest gene, the start of the gene and the distance of the regulatory element.

To make sure your script is working fine, run the following command:

```sh
python ../bin/get.distance.py --input gene.starts.tsv --start 980000
```

You should be getting this result:

ENSG00000187642.9	982093 2093

```sh

```


Task 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

```sh
cat regulatory.elements.starts.tsv | while read element start; do 
   python ../bin/get.distance.py ... # to be completed by you; 
done > regulatoryElements.genes.distances.tsv
```

```sh

```

Task 7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.


```sh

```