# SCEM <img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/SCEM.JPG" align = "right" width="180"/>
Single cell event matrix (SCEM) is an algorithm for event matrix generation from single cell whole genome copy number sequencing data. The event matrix can then be used to generate phylogenetic trees of the evolutionary relationship between individual cells.

<a href="https://zenodo.org/badge/latestdoi/297145258"><img src="https://zenodo.org/badge/297145258.svg" alt="DOI"></a>

When using the algorithm, please cite: 
bioRxiv 2022.04.01.486670

doi: https://doi.org/10.1101/2022.04.01.486670

## Setting up SCEM

SCEM can be installed by running the following lines of code

```
library(devtools)
devtools::install_github('NatalieKAndersson/SCEM')
library("SCEM")
```

The algorithm has the following dependencies. If they are not installed you can install them by using the command

```
install.packages(c("readxl","xlsxjars","xlsx","writexl","phangorn","stats","ade4","ape","tree","tidyverse","colorspace","ggplot2","ggtree","phytools","bbmle","picante","dendextend","ggforce","stringr"))
```
## Example data set
Let's start by analyzing an example data set! The file "Example.xlsx" contains an example of what a scWGS-dataset might look like. Through single cell whole genome sequencing, the copy number of each 1 Mbp (mega base pair) segment of each chromosome in each individual cell is approximated. This results in a matrix with thousands of rows indicating the copy number for each chromosome segment in each analyzed cell. In S. Figure 1 an example of a small proportion of what such a matrix might look like is visualized.

- First column: The bin number.
- Second and third column: The start and end position for that bin.
- Fourth column: The chromosome.
- Fifth column and beyond: The cells. Each column corresponds to a single cell or group of cells, indicated by row two. For example C1_1 contains 20 cells with identical genomic profiles and C1_7 is a single cell having a unique genomic profile.
- Each matrix element is the copy number in that particular chromosomal segment.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/Exampledata.JPG" width="800"/>

## Algorithm
For a detailed explanation of how the algorithm works, see the document "Extended Methods.pdf".

As a first step the algorhtm considers each column by itself. It then identifies all genetic alterations, along with their copy numbers and start and end positions in that column.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/Searchalgorithm.png" width="400"/>

As a second step it compares the genetic alterations found across the columns to identify cases where there are consecutive chromosomal alterations.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/Comparebetween.png" width="400"/>


## Event matrix
The output of the algorithm is an event matrix illustrating which genetic alterations are found in each unique cell group. Each row is a genetic alteration and each column a unique cell group. The matrix elements is either a "1" if the alteration is present or "0" if it is not.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/EM.png" width="200"/>

## Phylogeny
Based on the event matrix, phylogenetic trees can be reconstructed.

Maximum parsimony method
<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/MP_tree.png" width="200"/>

Maximum likelihood method
<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/ML_tree.png" width="200"/>

