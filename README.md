# recommendSVandHAsoftware
This project mainly realizes the automatic recommendation method of Structural variation detection tools and haplotype assembly tools for three generations of sequencing data.

With the rapid development of the third generation sequencing technology, haplotype assembly tools emerge in endlessly. However, when selecting haplotype assembly tools, the impact of upstream Structural variation detection tools on the assembly results needs to be considered, so Structural variation detection tools and haplotype assembly tools are combined for selection. This tool is based on a meta learning framework and constructs a metadatabase containing third-generation sequencing data and its optimal tool combination. When recommending new sequencing data, the closest sample in the metadatabase is found by comparing meta features, and the optimal combination of structural variation detection tools and haplotype assembly tools for the data is determined based on the optimal tool combination of the closest samples.

## configuration environment
| package | version |
|---------|---------|
| linux   | CentOS Linux release 7.7.1908 |
| python  | 3.6 |
| pandas  | 1.1.5 |
| numpy   | 1.19.5 |
| scikit-learn   | 0.19.2 |
| samtools  | 1.6 |
| minimap2  | 2.24 |

##Recommendation
Installation-free mode, copy the code directly and run the following command.

Recommendation running example:

python newSamplerecommendation.py -f sample.fq

Then, the recommendation results will be printed in the console.

## Experiment
Experiment running example

Picky:

 theLAST/bin/lastdb -v -P 1 hg19.lastdb hg19.fa samtools dict -H hg19.fa > hg19.seq.dict
 
 picky/Picky-0.2.a/src/picky.pl script --fastq LongRead.fastq --thread 1 > run.sh
 
 Let “theLAST/bin/lastal”, “picky/Picky-0.2.a/src/picky.pl”, “/refGenome/hg19.lastdb” and “/refGenome/hg19.fa” in “export LASTAL=”，“export PICKY=”，“export LASTALDB=” and “export LASTALDBFASTA=”
 
Nanosv:

 bedtools bamtobed -i sor.bam > sor.bed
 
 /svim_env/bin/NanoSV -t 1 -s /theSAMBAMBA/bin/sambamba -b sor.bed -o sor.vcf sor.bam
 
Pbsv:

 /thePBMM2/bin/pbmm2 align /refGenome/hg19.fa movie1.Q20.fastq ref.movie1.bam --sort --sample sample1 --rg '@RG\tID:movie1' 
 
 or /thePBMM2/bin/pbmm2 align /refGenome/hg19.fa movie1.Q20.fastq ref.movie1.bam --sort --preset CCS --sample sample1 --rg '@RG\tID:movie1'
 
 /thePBSV/bin/pbsv discover ref.movie1.bam ref.sample1.svsig.gz /thePBSV/bin/pbsv call /refGenome/hg19.fa ref.sample1.svsig.gz ref.var.vcf 
 
Flye:

 /theFlye/bin/flye --pacbio-raw LongRead.fastq --out-dir  /flye --threads 4
 
wtdbg2:

 /thewtdbg2/wtdbg2 -x rs -g size -i LongRead.fastq -t 16 -fo dbg
 
 /thewtdbg2/wtpoa-cns -t 16 -i dbg.ctg.lay.gz -fo dbg.raw.fa
 
 minimap2 -t16 -ax map-pb -r2k dbg.raw.fa LongRead.fastq | samtools sort -@4 >dbg.bam
 
 samtools view -F0x900 dbg.bam | /wtdbg2/wtpoa-cns -t 16 -d dbg.raw.fa -i - -fo dbg.cns.fa
 
Canu:

 /thecanu/canu -p project1 -d dir genomeSize=size -pacbio-raw LongRead.fastq
 
Mecat2:

 /MECAT2/Linux-amd64/bin/mecat.pl config config_file.txt
 
 /MECAT2/Linux-amd64/bin/mecat.pl correct config_file.txt
 
 /MECAT2/Linux-amd64/bin/mecat.pl trim config_file.txt
 
 /MECAT2/Linux-amd64/bin/mecat.pl assemble config_file.txt
 
Smartdenovo

 awk 'NR%4==1||NR%4==2' dir LongRead.fastq | sed 's/^@/>/g' > reads.fa
 
 smartdenovo.pl -c 1 reads.fa > wtasm.mak
 
 make -f wtasm.mak
 
