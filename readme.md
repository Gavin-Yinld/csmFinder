<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/csmFInder/blob/master/figures/csmFinder.gif"/></div>

# csmFinder

# Introduction

`csmFinder` is an R package for identifying putative cell-subset specific DNA methylation (pCSM) loci from single-cell or bulk methylomes. For single cell methylomes, a beta mixture model is involved to group the single cells into hyper- and hypo-methylated subsets. For bulk methylomes, a nonparametric Bayesian clustering algorithm is used for grouping the sequence reads into hyper- and hypo-methylated subsets. Both of them identify the genomic loci with significant methylation difference bwtween two subsets as pCSM loci. 

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
...
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
# The first column denotes the genomic coordinate of the segment.
# The second column denotes the methylation pattern of the segment. 
# All pattern information are separated by semicolon.
# Each pattern information contains two parts splited by colon, the pattern and the read depth supporting the pattern. 
# The pattern part denotes the methylation states of all CpG in the segment with 0 representing for unmethylation and 1 for methylation.
# For example, "1111:3" means this segment covered by 3 totlly methylated reads in this genomic loci.
```
For single-cell methylomes, `file_type="single-cell"` argument is needed, and the 4-CpG segments could be extracted as follow:

```R
#get the demo datasets
scDataDir <- paste(system.file(package="csmFinder"),"extdata/single_cell_CpG_extract_file",sep='/')
file_list <- paste(scDataDir,list.files(scDataDir),sep='/')

#generate the 4-CpG segment
scSegment <- bismark2segment(files=file_list,file_type="single-cell",CpG_file=CpG_ref)
```
The data format of single-cell methylome analysis is consistent with our previous study, please see [beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)

## Step 2. Find the candidate pCSM segments
The segments satisfy the following 2 criterions are considered as candidate pCSM segments:
1. The number of reads (for bulk methylome) or number of cells (for single-cell methylome) covering the segment is greater than threshold (default: 10).
2. The candidate pCSM segments should be covered by both totally methylated read and totally unmethylated read (for bulk methylome), or totally methylated cell and totally unmethylated cell (for single-cell methylome). 

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
scCandidate <- find_candidate(scSegment,data_type="single-cell",depth=10)
```
## Step 3. Identify pCSM segments 
For bulk methylomes, a nonparametric Bayesian clustering algorithm is used for grouping the sequence reads into hyper- and hypo-methylated subset and determining the genomic loci with significant difference bwtween two subsets as pCSM segments. For single cell methylomes, a beta mixture model is involved to group the single cells into hyper- and hypo-methylated subsets.
```perl
#for bulk methylome
pcsm_segment <- csmFinder(candidate,data_type='regular',thread=1)
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
scPcsm_segment <- csmFinder(scCandidate,data_type='single-cell',thread=1)
```
For the illustration of the output of single-cell analysis, please see [beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)
## Step 4. Merge the adjacent pCSM segments to pCSM loci
For both bulk methylome and single-cell methylome, the immediately overlapped pCSM segments are merged into pCSM region when `extension=0`. `extension` represents the bases need to be extended, for example, `extension=100` means that extend the pCSM segments 100bp in both of their upstream and downstream, and the immediately overlapped pCSM segments are merged after entension.
```R
#for bulk methylome
pcsm_loci <- merge_segment(pcsm_segment,extension=0)
pcsm_loci

    V1      V2      V3
1 chr1 3026929 3027017
2 chr1 3031539 3031623


#for single-cell methylome
scPcsm_loci <- merge_segment(scPcsm_segment,data_type="single-cell",extension=0)
```

## Further analysis. Co-methylation analysis and engen-pCSM loci extraction
csmFinder is the first step in our pipeline to perfom virtual methylome dissection, downstream analysis for grouping the pCSM loci, extracting eigen-pCSM loci and dissecting methylomes can be found in [coMethy](https://github.com/Gavin-Yinld/coMethy)


