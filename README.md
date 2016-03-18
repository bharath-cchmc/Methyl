---
title: "Methylation Site Identification"
author: "Bharath Manicka Vasagam, CCHMC, The Barski Laboratory"
date: "March 8, 2016"
output: html_document
---

##Objective
  The objective of this code is to to comapre the experimental methylation site output to that of a control of *Mus musculus* and identify the number of 
  1. Number of CPG islands in the experimental data.
  2. Average of the % methylation in the experimental data.

##Packages Used
1. data.table
2. GenomicRanges
3. rtracklayer
4. plyr

## Introduction

The methylation site identification code helps to compare the experimental promoter methylation sites output with that of the control (with all methylation site details) of *Mus musculus*.

## Conversion of files to BED files using LIFTOVER:
The experimental output file was in TXT format. It was converted to a standard BED format as `ChrN:Promoter Start - Promoter End` using R.

Then this csv file was uploaded to liftover <https://genome.ucsc.edu/cgi-bin/hgLiftOver> 

The BED output file from liftover was then data cleaned, compared with the control to find the related overlaps from the **GenomicRanges** package.


```{}
groupA<- data.table(exp_data)

#setting the table key
setkey(groupA, chr, promoter_start, promoter_end)
groupB<- data.table(bgwg)

#setting the table key
setkey(groupB, chr, start, end)

#overlapped files
over <- foverlaps(groupA, groupB, nomatch = 0)
```


## Result:
  Thus after identifying overlaps, number of CPG islands and % average methylation in the experimental data was calculated. To ensure data integrity, manual calculation of CPG islands was performed for few sanple promoter regions.
