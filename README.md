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

Let's start by loading the data into R!

```

#Load your data.
data <- load_data(filename="Example.xlsx",sheetname ="Example") #Extracting all of the CNV:s.
dm <- load_dist(filename="Example.xlsx",sheetname ="Example") #Extracting each bin's start and end position.
names(dm)[1]<-"Positions"
nm <- load_name(filename="Example.xlsx",sheetname ="Example") #Extracting the column names (cell names or name of a group of cells with identical profiles).
nrcells <- load_nrcells(filename="Example.xlsx",sheetname ="Example") #Extracting the number of cells represented by each column.

#Let the algorithm know if any of the cells are from a particular time-point is otherwise more closely related to each other for some reason.
#Example:
samples <- list(1:7,8:16) #Group 1 and 2.

co <- 1 #1 = 100 %. Choosing the cutoff for an event to be considered a stem event and thus be allocated to all cells in the EM.
```

## Algorithm
**See the document "Extended Methods.pdf" and "singlecell_flowchart.pdf" for a detailed explanation of the entire algorithm with examples.**

To make it simple:

1. As a first step the algorhtm considers each column by itself. It then identifies all genetic alterations, along with their copy numbers and start and end positions in that column.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/Searchalgorithm.png" width="400"/>

2. As a second step it compares the genetic alterations found across the columns to identify cases where there are consecutive chromosomal alterations.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/Comparebetween.png" width="400"/>

To run the algorithm, simply write
```
EM <- SCEM(data,dm,nm,nrcells,samples,co) #Generates an event matrix based on the scWGS-dataset.
```

## Event matrix
The output of the algorithm is an event matrix illustrating which genetic alterations are found in each unique cell group. Each row is a genetic alteration and each column a unique cell group. The matrix elements is either a "1" if the alteration is present or "0" if it is not.

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/EM_example.JPG" width="800"/>

Save the event matrix in an excel spreadsheet by writing.
```
write_xlsx(as.data.frame(EM), "EM_final.xlsx")
```

## Phylogeny
Based on the event matrix, phylogenetic trees can be reconstructed.

**Maximum parsimony method**

```
mp_tree <- function(EM_phy){
  MP_tree_pratchet <- pratchet(EM_phy, start = NULL, method = "fitch", maxit = 2000, k = 10, #Funktionen anv?nder sig av Fithchs algoritm. Pratchet = p ratchett (1999).
                               trace = 1, all = FALSE, rearrangements = "TBR",
                               perturbation = "ratchet")
  MP_tree_pratchet <- root(MP_tree_pratchet, outgroup = root, resolve.root = TRUE)
  treeRatchet <- acctran(MP_tree_pratchet, EM_phy)
  plot(treeRatchet)
  return(treeRatchet)}

RMS_mptree <- mp_tree(EM_phy)
RMSmp <- ggplot(RMS_mptree) + geom_tree() + geom_tiplab(size = 4)
RMSmp <- RMSmp +  theme_tree()+xlim(c(0, 10))
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
print(RMSmp+mytheme)
s <- 10
w <- 10
ggsave(RMSmp,file="MP_tree.png",width = w,height = s)
```

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/MP_tree.png" width="500"/>

**Maximum likelihood method**
```
ml_tree <- function(Eventmatrix) {
  dm_h <- dist.hamming(Eventmatrix) #Tar fram en avst?ndsmatris som ska ligga till grund f?r tr?dgenereringen.
  starting_tree <- NJ(dm_h)
  starting_tree <- root(starting_tree, outgroup = root,resolve.root = TRUE)
  Lf <- pml(starting_tree, Eventmatrix) #Obtaining an object of class pml
  Lf_JC <- optim.pml(Lf, model = "ER")#Jukes Cantor model.
  #bs <- bootstrap.pml(Lf_JC, bs=100, optNni=TRUE)
  #treeBS <- plotBS(Lf_JC$tree,bs, type = "phylogram")
  #s <- 15
  #ggsave(treeBS,file="bootstrap_ml.png",width = s,height = s)
  return(Lf_JC)
}

RMS_mltree <- ml_tree(EM_phy)
RMS_mltree <- RMS_mltree$tree
RMSml <- ggplot(RMS_mltree) + geom_tree() + geom_tiplab(size = 4)
RMSml <- RMSml +  theme_tree()+xlim(c(0, 0.4))
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
print(RMSml+mytheme)
s <- 10
w <- 10
ggsave(RMSml,file="ML_tree.png",width = w,height = s)
```

<img src="https://github.com/NatalieKAndersson/SCEM/blob/main/Figures/ML_tree.png" width="500"/>

