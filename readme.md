<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/csmFInder/blob/master/figures/csmFinder.gif"/></div>

# csmFinder

# Introduction

`csmFinder` is an R package for identifying putative cell-subset specific DNA methylation (pCSM) loci from single-cell or bulk methylomes. For single cell methylomes, a beta mixture model is involved to group the single cells into hyper- and hypo-methylated subsets. For bulk methylomes, a nonparametric Bayesian clustering algorithm is used for grouping the sequence reads into hyper- and hypo-methylated subset. Both of them identify the genomic loci with significant methylation difference bwtween two subsets as pCSM loci. 

# Installation
`csmFinder` needs the following tools to be installed and available in the `PATH` environment:
1.  [R](https://www.r-project.org/)(>=3.4.4)
2.  [python2](https://www.python.org/downloads/) (>=2.7.10), with Numpy and Pandas module, to process the bismark extractor results
3.  [bedtools2](https://github.com/arq5x/bedtools2) (>=2.28.0), to merge the overlapped pCSM segments into pCSM loci

In R console,
```R
install.packages("devtools")
library("devtools")
devtools::install_github("PhanstielLab/bedtoolsr") # to merge the overlopped pCSM segments
devtools::install_github("Gavin-Yinld/csmFinder")
```
# How to Use

## Step 1. Generate 4-CpG segments from bismark extractor results
`csmFinder` takes methylation state per base information of each sequence read in CpG context specificlly as input. Such input file may be obtained from `bismark` pipeline. A typical input file should be in ".gz" compressed format and looks like this:

```
ST-E00523:376:HL7JTCCXY:4:1102:23054:14494_1:N:0:ATGAGCAT + chr1  3023890 Z
ST-E00523:376:HL7JTCCXY:4:1102:23054:14494_1:N:0:ATGAGCAT + chr1  3023859 Z
ST-E00523:376:HL7JTCCXY:4:1107:31730:36592_1:N:0:ATGAGCAT + chr1  3023890 Z
ST-E00523:376:HL7JTCCXY:4:1107:31730:36592_1:N:0:ATGAGCAT + chr1  3023859 Z
ST-E00523:376:HL7JTCCXY:4:1204:29305:15232_1:N:0:ATGAGCAT + chr1  3022537 Z
```
To generate the 4-CpG segment coordinate, another table seperated text file contains the position of C in CpG context
in forward strand of target genome is needed.

```chr1	3021025
chr1	3021077
chr1	3021345
chr1	3021400
chr1	3021720
...
```
For bulk methylome, the 4-CpG segments could be extracted as follow:

```perl
library('csmFinder')

#get the demo datasets
bismark_result=paste(system.file(package="csmFinder"),"extdata/bulk_CpG_extract_file/demo.dataset.gz",sep='/')
CpG_ref=paste(system.file(package="csmFinder"),"extdata/CpG_plus.reference",sep='/')

#generate the 4-CpG segment
segment <- bismark2segment(files=bismark_result,CpG_file=CpG_ref)

#see what the segment looks like
segment[2:5,] 
                                    V1                    V2
2 chr1:3026187_3026278_3026310_3026413 0001:1;0101:2;1101:1;
3 chr1:3026278_3026310_3026413_3026420 0011:3;1011:4;1111:1;
4 chr1:3026864_3026883_3026895_3026926        0111:2;1111:1;
5 chr1:3026883_3026895_3026926_3026929               1111:3;
#The first column denotes the genomic coordinate of the segment, the second column denotes the methylation pattern of the segment. 
#All patterns observed on the segment are listed followed; each pattern information contains the pattern and the read depth supporting the pattern; all pattern information are separated by semicolon.
#For example, "1111:3" means this segment covered by 3 totlly methylated reads in this genomic loci.
```
For single-cell methylomes, `file_type="single-cell"` argument is needed, and the 4-CpG segments could be extracted as follow:

```R
#get the demo datasets
dir <- paste(system.file(package="csmFinder"),"extdata/single_cell_CpG_extract_file",sep='/')
file_list <- paste(dir,list.files(dir),sep='/')

#generate the 4-CpG segment
segment2 <- bismark2segment(files=file_list,file_type="single-cell",CpG_file=CpG_ref)
```
The data format of single-cell methylome analysis is keep pace with our previous study, please see [beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)

## Step 2. Find the candidate pCSM segment
The segment meet the following 2 conditions are considered as candidate pCSM segment:
1. The read depth (for bulk methylome) or number of cells covered the segment (for single-cell methylome) greater than threshold (default: 10)
2. Covered by totally methylated read and totally unmethylated read (totlly methylated cell and totlly unmethylated cell for single-cell analysis ) in this segment at the same time. 
In R console,
```perl
#for bulk methylome
candidate <- find_candidate(segment,depth=10)
#see what the segment looks like
candidate[1:5,]
                                     V1                     V2
8  chr1:3026929_3026936_3026969_3027017         0000:8;1111:5;
16 chr1:3027091_3027102_3027106_3027144 0000:15;1110:3;1111:6;
28 chr1:3031539_3031569_3031573_3031581  0000:1;0111:1;1111:9;
29 chr1:3031569_3031573_3031581_3031586 0000:4;1110:1;1111:11;
30 chr1:3031573_3031581_3031586_3031623         0000:2;1111:9;

#for single-cell methylome
candidate2 <- find_candidate(segment2,data_type="single-cell",depth=10)
```
## Step 3. Identify pCSM segment 
For bulk methylomes, a nonparametric Bayesian clustering algorithm is used for grouping the sequence reads into hyper- and hypo-methylated subset and determining the genomic loci with significant difference bwtween two subsets as so called pCSM loci. For single cell methylomes, a beta mixture model is involved to group the single cells into hyper- and hypo-methylated subsets.
```perl
#for bulk methylome
pcsm_segment <- csmFinder(candidate,data_type='regular')
pcsm_segment
                                     V1                     V2         d
8  chr1:3026929_3026936_3026969_3027017         0000:8;1111:5; 1.0000000
28 chr1:3031539_3031569_3031573_3031581  0000:1;0111:1;1111:9; 0.9750000
29 chr1:3031569_3031573_3031581_3031586 0000:4;1110:1;1111:11; 0.9791667
30 chr1:3031573_3031581_3031586_3031623         0000:2;1111:9; 1.0000000
           pval
8  0.000000e+00
28 2.482661e-05
29 1.766860e-07
30 0.000000e+00
#"d" means the methylation difference between hypo- and hyper-methylated reads.

#for single-cell methylome
pcsm_segment2 <- csmFinder(candidate2,data_type='single-cell')
```
For the illustration of the output of single-cell analysis, please see [beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)
## Step 4. Merge the overlapped pCSM segment to pCSM loci
```R
#for bulk methylome
pcsm_loci <- merge_segment(pcsm_segment,extension=0)
pcsm_loci
    V1      V2      V3
1 chr1 3027075 3027102
2 chr1 3031383 3031539
3 chr1 3031581 3031632
4 chr1 3035642 3035833
5 chr1 3035927 3036246
6 chr1 3037739 3037802

#for single-cell methylome
pcsm_loci2 <- merge_segment(pcsm_segment2,data_type="single-cell",extension=0)
```

