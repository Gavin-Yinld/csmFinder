<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/csmFInder/blob/master/figures/csmFinder.gif"/></div>

# csmFinder

# Introduction

`csmFinder` is an R package for identifying putative cell-subset specific DNA methylation (pCSM) loci from single-cell or bulk methylomes. For single cell methylomes, a beta mixture model is used to group the single cells into hyper- and hypo-methylated subsets. For bulk methylomes, a nonparametric Bayesian clustering algorithm is used to group the pooled sequence reads into hyper- and hypo-methylated subsets. Both methods identify the genomic loci with significant methylation difference between the two subsets as pCSM loci.

# Installation
`csmFinder` needs the following tools to be installed and available in the `PATH` environment:
1.  [R](https://www.r-project.org/)(>=3.4.4)
2.  [python2](https://www.python.org/downloads/) (>=2.7.10), with Numpy and Pandas modules installed, to process the bismark extractor results
3.  [bedtools2](https://github.com/arq5x/bedtools2) (>=2.28.0), to merge the overlapped pCSM segments

In R console,
```R
install.packages("devtools")
library("devtools")
devtools::install_github("PhanstielLab/bedtoolsr") # to merge the overlopped pCSM segments
devtools::install_github("Gavin-Yinld/csmFinder")
```
# How to Use

## Step 1. Generate 4-CpG segments from bismark extractor results
`csmFinder` takes methylation data from bisulfite sequence reads on the genome as input. Such an input file should be in ".gz" compressed format, and may be obtained from the bismark pipeline. An uncompressed example is listed below, with each row representing a bisulfite sequence read:

```
ST-E00523:376:HL7JTCCXY:4:1102:23054:14494_1:N:0:ATGAGCAT + chr1  3023890 Z
ST-E00523:376:HL7JTCCXY:4:1102:23054:14494_1:N:0:ATGAGCAT + chr1  3023859 Z
ST-E00523:376:HL7JTCCXY:4:1107:31730:36592_1:N:0:ATGAGCAT + chr1  3023890 Z
ST-E00523:376:HL7JTCCXY:4:1107:31730:36592_1:N:0:ATGAGCAT + chr1  3023859 Z
ST-E00523:376:HL7JTCCXY:4:1204:29305:15232_1:N:0:ATGAGCAT + chr1  3022537 Z
...
```
To generate the 4-CpG segment coordinates, a table separated text file containing the position of all the cytosines in CpG context in the forward strand of genome is needed.

```chr1	3021025
chr1	3021077
chr1	3021345
chr1	3021400
chr1	3021720
...
```
For bulk methylomes, the 4-CpG segments could be extracted as follows:

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
# All pattern information is separated by semicolon.
# Each pattern information contains two parts split by colon, the pattern and the read depth supporting the pattern. 
# The pattern part denotes the methylation states of all CpG in the segment with 0 representing for unmethylation and 1 for methylation.
# For example, "1111:3" means this segment covered by 3 totally methylated reads in this genomic segment.
```
For single-cell methylomes, the argument `file_type="single-cell"`  is needed for function bismark2segment, and the 4-CpG segments could be extracted as follows:

```R
#get the demo datasets
scDataDir <- paste(system.file(package="csmFinder"),"extdata/single_cell_CpG_extract_file",sep='/')
file_list <- paste(scDataDir,list.files(scDataDir),sep='/')

#generate the 4-CpG segment
scSegment <- bismark2segment(files=file_list,file_type="single-cell",CpG_file=CpG_ref)
```
The data format used in the single-cell methylome analysis is consistent with our previous study ([beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)).

## Step 2. Find candidate pCSM segments
The segments satisfying the following 2 criteria are considered as candidate pCSM segments:
1. To ensure the statistical power, the number of reads (for bulk methylome) or number of cells (for single-cell methylome) covering the segment is greater than the threshold (default: depth=10). Increasing depth cutoff leads to more reliable results, but also less loci retained.

2. The candidate pCSM segments should be covered by both totally methylated reads and totally unmethylated reads (for bulk methylome), or totally methylated cells and totally unmethylated cells (for single-cell methylome).

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
For bulk methylomes, a nonparametric Bayesian clustering algorithm is used to group the sequence reads into hyper- and hypo-methylated subsets and pCSM segments are determined by testing difference between the two subsets. For single cell methylomes, a beta mixture model is used to group the single cells into hyper- and hypo-methylated subsets. 

By default, the loci with methylation difference between two subsets over 0.3 (distance=0.3) and the significance of the methylation difference between two subsets less than 0.05 (pval=0.05) are determined as pCSM loci. Increase the `distance` cutoff or reduce `pval` cutoff leads to increased methylation difference between two subsets. One can speed up the run by increasing the number of computational threads.
```perl
#for bulk methylome
pcsm_segment <- csmFinder(candidate,data_type='regular',distance=0.3,pval=0.05,thread=1)
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
scPcsm_segment <- csmFinder(scCandidate,data_type='single-cell',distance=0.3,pval=0.05,thread=1)
```
For the illustration of the output of the single-cell methylome analysis, please see [beta mixture model](https://github.com/Evan-Evans/Beta-Mixture-Model)
## Step 4. Merge adjacent pCSM segments to pCSM loci
For both bulk and single-cell methylomes, the overlapped pCSM segments are merged into pCSM loci.
```R
#for bulk methylome
pcsm_loci <- merge_segment(pcsm_segment)
pcsm_loci

    V1      V2      V3
1 chr1 3026929 3027017
2 chr1 3031539 3031623


#for single-cell methylome
scPcsm_loci <- merge_segment(scPcsm_segment,data_type="single-cell")
```

## Further analysis to perfom virtual methylome dissection
csmFinder is the first step in our pipeline to perfom virtual methylome dissection, downstream analysis to group the pCSM loci, extract eigen-pCSM loci and dissect methylomes can be found in the [coMethy](https://github.com/Gavin-Yinld/coMethy) package.


