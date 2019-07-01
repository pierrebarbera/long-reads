#!/bin/bash

# Read curation pipeline, as described in Jamy et al, 2019

# usage: curateReads.sh <sampleName> <path to SILVA SSU or PR2 database>

# Requires that the following softwares are installed:
## mothur
## barrnap 
## vsearch

# Requires the following custom scripts
## polarize_reads.pl
## extract_rrna.pl

# parameters
sample=$1
database=$2


## Convert fastq file to fasta (in Jamy et al, 2019, all fastq files were concatenated and the resulting pooled.fastq file was analysed altogether)
echo "Running mothur's fastq.info() ..."
mothur "#fastq.info(fastq=$sample.fastq, pacbio=T)"


## Filter sequences: discard sequences with >6homopolymers, length<2500 and >6000, and trim based on quality.
echo "Running mothur's trim.seqs() ..."
mothur "#trim.seqs(fasta=$sample.fasta, processors=20, minlength=2500, maxlength=6000, maxhomop=6, qfile=$sample.qual, qwindowsize=50, qwindowaverage=30)"

rm $sample.trim.qual
rm $sample.scrap.qual


## Detect sequences with 18S and 28S and discard sequences with 
### A note on barrnap: this detects 18S and 28S based on HMM built from the SILVA database. However, since it 
### uses NHMMER which currently only supports local alignment, barrnap sometimes gets the ends 'wrong' by a few bases. 
### It may also miss weird rRNAs. RNAmmer uses HMMER in global alignment (which would solve the above problems) - however 
### it does not work with partial eukaryotic sequences

### detect 18S and 28S
echo "Running barrnap ..."
barrnap --threads 50 --reject 0.4 --kingdom euk $sample.trim.fasta > barrnap.SSU.LSU.$sample.gff

### Generate sorted list of single 18S and 28S to compare. Note the sequences common to both and extract them into a new fasta file.

cat barrnap.SSU.LSU.$sample.gff | grep "18S" | cut -f 1 | sort | uniq -u > barrnap.SSU.$sample.list
cat barrnap.SSU.LSU.$sample.gff | grep "28S" | cut -f 1 | sort | uniq -u > barrnap.LSU.$sample.list

comm -12 barrnap.SSU.$sample.list barrnap.LSU.$sample.list > barrnap.SSU.LSU.$sample.comm.list

echo "Extracting seqs with 18S and 28S ..."
mothur "#get.seqs(accnos=barrnap.SSU.LSU.$sample.comm.list, fasta=$sample.trim.fasta)"
### This generates a file labelled $sample.trim.pick.fasta


## Pre-cluster sequences to remove remaining errors and to make number of sequences more manageable
### Polarize reads so they are in same direction. 
echo "polarizing reads ..."
perl polarize_reads.pl barrnap.SSU.LSU.$sample.gff $sample.trim.pick.fasta $sample.trim.pick.pol.fasta

### Cluster at 99% identity
mkdir preclusters
echo "Preclustering ..."
vsearch --cluster_fast $sample.trim.pick.pol.fasta --id 0.99 --clusters preclusters/$sample.c-

### Separate preclusters > 2 sequences
cd preclusters/
mkdir large_preclusters
for i in $sample.c-*; do number=$(grep -c ">" $i); if (($number > 2)); then mv $i large_preclusters/; fi; done 
cd large_preclusters

### Align large preclusters with mafft
echo "Generating consensus for large preclusters ..."
for i in $sample.c-*; do mafft --reorder --thread 50 $i > $i.aligned.fasta; done

### Generate majority rule consensus sequences
for i in *.aligned.fasta; do mothur "#consensus.seqs(fasta=$i, cutoff=51)"; done

### Change header name so that it is informative
for i in *.aligned.cons.fasta; do x=$(echo $i | cut -d '.' -f2); sed -i -E "s/(>)(.*)/\1${x}_\2/" $i; done

### mv consensus sequences fasta files back to the parent folder
for i in *.aligned.cons.fasta; do mv $i ./..; done
cd ..
cat * > $sample.trim.pick.pol.preclusters.fasta

### degap sequences
mothur "#degap.seqs(fasta=$sample.trim.pick.pol.preclusters.fasta)"
cp $sample.trim.pick.pol.preclusters.ng.fasta ./..


## Chimera detection
### Construct abundance table for chimera detection
echo "De novo chimera detection ..."
cd large_preclusters
for i in $sample.c-*.aligned.fasta; do x=$(echo $i | cut -d '.' -f2); count=$(grep -c ">" $i); echo -e ""$x"_conseq\t"$count"" >> $sample.trim.pick.pol.preclusters.ng.count_table; done

mv $sample.trim.pick.pol.preclusters.ng.count_table ./..
cd ..

grep ">m" $sample.trim.pick.pol.preclusters.ng.fasta | while read line; do header=$(echo $line | sed -E 's/>(.*)/\1/'); echo -e "$header\t1" >> $sample.trim.pick.pol.preclusters.ng.count_table; done

mv $sample.trim.pick.pol.preclusters.ng.count_table ./..
cd ..

echo 'Representative_Sequence\ttotal' | cat - $sample.trim.pick.pol.preclusters.ng.count_table > temp && mv temp $sample.trim.pick.pol.preclusters.ng.count_table

### chimera detection
### Note: this step often takes the longest. See Martijn et al. (2019) for justification of values picked for abskew and chunks 
mothur "#chimera.uchime(fasta=$sample.trim.pick.euk.pol.preclusters.ng.fasta, count=$sample.trim.pick.euk.pol.preclusters.ng.count_table, reference=self, chunks=40, abskew=1, chimealns=T)"

mothur "#remove.seqs(fasta=$sample.trim.pick.pol.preclusters.ng.fasta, count=$sample.trim.pick.pol.preclusters.ng.count_table, accnos=$sample.trim.pick.pol.preclusters.ng.denovo.uchime.accnos)"


## Remove prokaryotic sequences
### based on blast results against the SILVA SSU or PR2 database
echo "Removing prok seqs ..."
blastn -query $sample.trim.pick.pol.preclusters.ng.pick.fasta -db $database \
		-evalue 1e-10 -num_threads 50 \
        -out $sample.trim.pick.pol.preclusters.ng.pick_vs_SSUdatabase.blastn \
        -outfmt '6 std salltitles'
        
cat $sample.trim.pick.pol.preclusters.ng.pick_vs_SSUdatabase.blastn | sort -k1,1 -k12,12nr | \
	awk '!seen[$1]++' > tophit.$sample.trim.pick.pol.preclusters.ng.pick_vs_SSUdatabase.blastn
    
cat tophit.$sample.trim.pick.pol.preclusters.ng.pick_vs_SSUdatabase.blastn | 
    grep -v "Eukaryota" | cut -f 1 > $sample.trim.pick.pol.preclusters.ng.pick.prok.list

     
mothur "#remove.seqs(fasta=$sample.trim.pick.pol.preclusters.ng.pick.fasta, \
        count=$sample.trim.pick.pol.preclusters.ng.pick.count_table, \
        accnos=$sample.trim.pick.pol.preclusters.ng.pick.prok.list)"
        
        
## Extract 18S and 28S from each sequence using a perl script (need to re-run barrnap).
echo "Extracting 18S and 28S ..."

barrnap --threads 50 --reject 0.4 \
        --kingdom euk $sample.trim.pick.pol.preclusters.ng.pick.pick.fasta \
        > barrnap.SSU.LSU.$sample.trim.pick.pol.preclusters.ng.pick.pick.gff
        
perl extract_rrna.pl \
     barrnap.SSU.LSU.$sample.trim.pick.pol.preclusters.ng.pick.pick.gff \
     tophit.$sample.trim.pick.pol.preclusters.ng.pick_vs_SSUdatabase.blastn \
     $sample.trim.pick.pol.preclusters.ng.pick.pick.fasta \
     18S.$sample.fasta 28S.$sample.fasta
     
     
## Clustering into OTUs
echo "Aligning 18S and 28S ..."
mafft --reorder --thread 25 18S.$sample.fasta > 18S.$sample.aligned.fasta
mafft --reorder --thread 25 28S.$sample.fasta > 28S.$sample.aligned.fasta


mkdir otus
cd otus

echo "Clustering into OTUs using average neighbour method ..."
mothur "#dist.seqs(fasta=18S.$sample.aligned.fasta, 
         countends=F, cutoff=0.20, processors=18)"
         
mothur "#cluster(column=18S.$sample.aligned.dist, 
         count=../$sample.trim.pick.pol.preclusters.ng.pick.pick.count_table, 
         method=average)"
    
         
echo "Picking OTU representative ..."
### See number of OTUs when clustering at different similarity levels
cat 18S.$sample.aligned.an.list | cut -f 1,2

### Clustering at 97% siilarity (label=0.03 option)
mothur "#get.oturep(column=18S.$sample.aligned.dist, \
  count=../$sample.trim.pick.pol.preclusters.ng.pick.pick.count_table, \
  list=18S.$sample.aligned.an.list, label=0.03, \
  fasta=../18S.$sample.fasta, method=distance)"


### Get 28S OTU representatives based on 18S OTUs
grep ">" 18S.$sample.rep.fasta | tr -d '>' > 28S.$sample.rep.list
mothur "#get.seqs(accnos=28S.$sample.rep.list, fasta=../28S.$sample.fasta)"


## Files containing OTU representatives of the 18S and 28S are:
## 18S.$sample.rep.fasta, and
## 28S.$sample.pick.fasta
## repectively

