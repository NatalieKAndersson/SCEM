########################################
#Single Cell Event Matrix - SCEM       #
#Building an event matrix from CNA data#
########################################
setwd("~/") #Set your working directory (where the scWGS file is).
#Incorporating dependencies
library("readxl")
library("xlsxjars")
library("xlsx")
library("writexl")
library("phangorn")
library("stats")
library("ade4")
library("ape")
library("tree")
library("tidyverse")
library("colorspace")
library("ggplot2")
library("ggtree")
library("phytools")
library("bbmle")
library("picante")
library("dendextend")
library("ggforce")
library("stringr")

######################################
#Loading and preparing the data files#
######################################

load_data <- function(filename, sheetname) {  #Skapar en funktion som fixar till datan.
  data <- as.data.frame(read_xlsx(filename, sheetname)) #H?r l?ser jag in xlsx-filen och sparar den i variabeln data
  View(data)
  subdata <- data[ c(3:nrow(data)), c(5:ncol(data)) ] #H?r plockar jag ut den del av filen som jag faktiskt vill anv?nda i mina ber?kningar sen.
  return(subdata)
}

load_dist <- function(filename, sheetname) {  #Skapar en funktion som fixar till datan.
  data <- as.data.frame(read_xlsx(filename, sheetname)) #H?r l?ser jag in xlsx-filen och sparar den i variabeln data
  subdata <- data[ c(3:nrow(data)), c(1:4) ] #H?r plockar jag ut den del av filen som jag faktiskt vill anv?nda i mina ber?kningar sen.
  return(subdata)
}

load_name <- function(filename, sheetname) {  #Skapar en funktion som fixar till datan.
  data <- as.data.frame(read_xlsx(filename, sheetname)) #H?r l?ser jag in xlsx-filen och sparar den i variabeln data
  subdata <- data[2, c(5:ncol(data)) ] #H?r plockar jag ut den del av filen som jag faktiskt vill anv?nda i mina ber?kningar sen.
  return(subdata)
}

load_nrcells <- function(filename, sheetname) {  #Skapar en funktion som fixar till datan.
  data <- as.data.frame(read_xlsx(filename, sheetname)) #H?r l?ser jag in xlsx-filen och sparar den i variabeln data
  subdata <- data[1, c(5:ncol(data)) ] #H?r plockar jag ut den del av filen som jag faktiskt vill anv?nda i mina ber?kningar sen.
  return(subdata)
}

phydatdata <- function(excelfil){
  patient.matrix <- as.matrix(t(excelfil))
  patient.phydat <- phyDat(patient.matrix,
                           type="USER", #G?r s? att jag sj?lv kan v?lja vilka niv?er vi har i datan.
                           levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,53,24,112,101,"-","?","n"), #Anger niv?erna fr?n 0 till det maximala v?rdet p? CNV.
                           ambiguity='2') #Hanterar de v?rden som inte st?mmer in p? ovanst?ende.
  return(patient.phydat) #returnerar datan i phyDat-formatet
}

phydatevent <- function(excelfil,cellnames,normal){
  patient.matrix <- as.matrix(t(excelfil))
  rownames(patient.matrix) <- cellnames
  if(normal == "Normal"){
    M <- matrix(0, 1, ncol(patient.matrix))
    patient.matrix <- rbind(patient.matrix,M)
    rownames(patient.matrix)[nrow(patient.matrix)] <- "Normal"}
  patient.phydat <- phyDat(patient.matrix,
                           type="USER", #G?r s? att jag sj?lv kan v?lja vilka niv?er vi har i datan.
                           levels=c(0,1), #Anger niv?erna fr?n 0 till det maximala v?rdet p? CNV.
                           ambiguity='0') #Hanterar de v?rden som inte st?mmer in p? ovanst?ende.
  return(patient.phydat) #returnerar datan i phyDat-formatet
}


#Load your data.
data <- load_data(filename="File.xlsx",sheetname ="Sheet1") #Extracting all of the CNV:s.
dm <- load_dist(filename="File.xlsx",sheetname ="Sheet1") #Extracting each bin's start and end position.
names(dm)[1]<-"Positions"
nm <- load_name(filename="File.xlsx",sheetname ="Sheet1") #Extracting the column names (cell names or name of a group of cells with identical profiles).
nrcells <- load_nrcells(filename="File.xlsx",sheetname ="Sheet1") #Extracting the number of cells represented by each column.

#Let the algorithm know if any of the cells are from a particular time-point is otherwise more closely related to each other for some reason.
#Example:
samples <- list(1:7,8:16, 17:22) #Group 1, 2 and 3.

co <- 1 #1 = 100 %. Choosing the cutoff for an event to be considered a stem event and thus be allocated to all cells in the EM.

######################
#Fusing equal columns#
######################
#
# cells <- data[3:nrow(data),]
#
# cell_clusters <- matrix(0,100,ncol(cells))
#
# i <- 1
# c <- 1
# for(i in 1:(ncol(cells)-1)){
#
#   cell1 <- cells[,i]
#   cell1_n <- nm[i]
#
#   j <- i+1
#   for(j in (i+1):ncol(cells)){
#
#     cell2 <- cells[,j]
#     cell2_n <- nm[j]
#
#     celldiff <- as.numeric(cell1)-as.numeric(cell2)
#
#     if(all(celldiff==0) && cell1_n!=cell2_n){
#       #They are equal.
#       print("equal")
#       print(cell1_n)
#       if(cell1_n %in% cell_clusters == FALSE){
#         k <- 1
#         print("Not in cell_clusters")
#
#         for(k in 1:(as.numeric(nrow(cell_clusters))-1)){
#           if(cell_clusters[k,c]==0 && cell1_n %in% cell_clusters == FALSE){
#             print("Adding")
#             print(c)
#             print(k)
#             cell_clusters[k,c] <- as.character(cell1_n)
#             cell_clusters[k+1,c] <- as.character(cell2_n)
#             c <- c+1
#             break
#           }
#           k <- k+1
#         }
#
#       }else{
#
#         if(cell2_n %in% cell_clusters == FALSE){
#
#           l <- 1
#
#           for(l in 1:(as.numeric(ncol(cell_clusters))-1)){
#             cell1_pos <- match(cell1_n,cell_clusters[,l])
#             print("oj")
#
#             if(is.na(cell1_pos)==FALSE){
#
#               m <- 1
#               for(m in 1:length(cell_clusters[,l])){
#
#                 if(cell_clusters[m,l]==0 && cell2_n %in% cell_clusters == FALSE){
#                   cell_clusters[m,l] <- as.character(cell2_n)
#                 }
#                 m <- m+1
#               }
#
#             }
#             k <- k+1
#           }
#
#
#         }
#       }
#
#
#     }
#
#    j <- j+1
#   }
#
#   i <- i+1
# }
#
# cell_clusters_all <- cell_clusters[,cell_clusters[1,]!=0]
#
# j <- 1
# lengths <- matrix(0,1,ncol(cell_clusters_all))
# for(j in 1:ncol(cell_clusters_all)){
#   true <- cell_clusters_all[cell_clusters_all[,j]!="0",j]
#   lengths[1,j] <- length(true)
#   j <- j+1
# }
#
# cell_clusters_all <- rbind(lengths,cell_clusters_all)
# cell_clusters_order <- cell_clusters_all[,order(as.numeric(cell_clusters_all[1,]),decreasing=TRUE)]
# names_clusters <- t(as.matrix(paste(c("A (","B (","C (","D (","E (","F (","G (","H (", "I ("),cell_clusters_order[1,],c(")",")",")",")",")",")",")",")",")"),sep="")))
# cell_clusters_order <- rbind(names_clusters[1:ncol(cell_clusters_order)],cell_clusters_order)
#
# write.xlsx(cell_clusters_order,"test.xlsx")
#
#
# i <- 1
# data_n <- rbind(nm,data)
# for(i in 1:ncol(cell_clusters_order)){
#
#   j <-1
#   for(j in 3:nrow(cell_clusters_order)){
#
#     if(cell_clusters_order[j,i]!="0"){
#
#       col <- match(cell_clusters_order[j,i],data_n[1,])
#       data_n[1,col] <- cell_clusters_order[1,i]
#
#     }
#     j <- j+1
#   }
#  i <- i+1
# }
#
# data_n <- t(as.matrix(data_n))
#
# new_data <- data_n[,!duplicated(c(data_n[1,]))]
#
# write.xlsx(new_data,"test.xlsx")
######################################################
#Finding all events belonging to each individual cell#
######################################################
#Creating an empty matrix that will contain the events for each cell.
j = 1
size <- 2000
d <- matrix(0,size,5*ncol(data))

for(j in 1:ncol(data)){
  e = 1
  i = 1
  for(i in 1:nrow(data)){ #Loopar igenom alla rader f?r varje cell.
    if(as.numeric(data[i,j]) != 2){ #Om v?rdet skiljer sig fr?n 2 har vi en h?ndelse.
      
      if(i==1){ #Hanterar h?ndelser som f?s redan i den f?rsta bin:en.
        d[e,1+5*(j-1)] <- dm[1,1]
        d[e,2+5*(j-1)] <- dm[1,2]
      }
      
      if(i > 1){
        if(i < 2721){
          if(data[i,j] != data[i-1,j]){ #Skiljer sig v?rdet i bin:en fr?n det v?rdet i tidigare bin? Spara ner startpositionen och kr isf.
            d[e,1+5*(j-1)] <- dm[i,1] #H?ndelsens position.
            d[e,2+5*(j-1)] <- dm[i,2]} #Startpositionen.
          
          else if(dm[i,4] != dm[i-1,4]){ #Ny kromosom? Spara isf ner slutpositionen fr?n bin:en innan samt nya startpositioner.
            d[e,3+5*(j-1)] <- dm[i-1,3] #Slutpositionen.
            e <- e+1
            d[e,1+5*(j-1)] <- dm[i,1] #Sparar ner positionen och startpositionen f?r n?sta h?ndelse.
            d[e,2+5*(j-1)] <- dm[i,2]}
          
          if(data[i,j] != data[i+1,j]){ #Skiljer sig v?rdet i bin:en fr?n n?stkommande v?rde? Spara d? ned slutpositionen.
            d[e,3+5*(j-1)] <- dm[i,3]
            d[e,4+5*(j-1)] <- dm[i,1]
            d[e,5+5*(j-1)] <- data[i,j] #CNA,v?rdet
            e <- e+1}
          
          else if(dm[i,4] != dm[i+1,4]){ #?r n?stkommande bin en ny kromosom? Spara ner slutpositionen.
            d[e,3+5*(j-1)] <- dm[i,3]
            d[e,4+5*(j-1)] <- dm[i,1]
            d[e,5+5*(j-1)] <- data[i,j]}
        }
        if(i == 2721){ #Hanterar den sista positionen
          d[e,3+5*(j-1)] <- dm[2721,3]
          d[e,4+5*(j-1)] <- dm[2721,1]
          d[e,5+5*(j-1)] <- data[2721,j]}}
      
    }
    i <- i+1}
  j <- j+1
}

eventpos <- matrix(0,size,3*ncol(data)) #Samlar alla h?ndelsepositionerna f?r varje cell i en matris.
print("hej")
for(j in 1:ncol(data)){
  print(j)
  eventpos[,1+3*(j-1)] <- d[,1+5*(j-1)]
  eventpos[,2+3*(j-1)] <- d[,4+5*(j-1)]
  eventpos[,3+3*(j-1)] <- d[,5+5*(j-1)]
  j <- j+1
}

eventstot <- matrix(0,nrow(eventpos),ncol(data)) #S?tter samman de tre positionerna en enda bin.
for(j in 1:ncol(data)){
  eventstot[,j] <- paste(eventpos[,1+3*(j-1)],eventpos[,2+3*(j-1)],eventpos[,3+3*(j-1)])
  j <- j+1
}

#######################################
#Treating events overlaying each other#
#######################################
# eventstot_sparad <- eventstot
# eventpos_sparad <- eventpos
j = 1
i = 1
event = 35
u <- matrix(0,1,3)
v <- matrix(0,1,3)
x <- matrix(0,1,3)
y <- matrix(0,1,3)
z <- matrix(0,1,3)
w <- matrix(0,1,3)
ev1 <- matrix(0,1,3)
ev2 <- matrix(0,1,3)
samples <- lapply(samples, as.character)
for(j in 1:ncol(data)){ #Hanterar h?ndelser som befinner sig ovanp? andra h?ndelser.
  i = 1
  event = 35
  
  if(j == ncol(data)){
    print("Sista cellen")
    print(data[1:5,j])
  }
  
  cluster <- which(sapply(samples, FUN=function(X) j %in% X)) #Finding the list in which this cell belongs.
  cluster_cells <- unlist(samples[cluster])
  first <- cluster_cells[1] #Start position
  last <- cluster_cells[length(cluster_cells)] #Stop position
  tot_cells <- sum(as.numeric(nrcells[1,first:last]))
  
  for(i in 1:(nrow(eventpos)-4)){
    #Hanterar en händelse med en 2:händelse ovanpå. Ex: 3333223333.
    if(eventpos[i,3+3*(j-1)] == eventpos[i+1,3+3*(j-1)]){
      
      w[1,1] <- eventpos[i,1+3*(j-1)] #Första händelsens startposition.
      w[1,2] <- eventpos[i+1,2+3*(j-1)] #Andra händelsens slutposition.
      w[1,3] <- eventpos[i+1,3+3*(j-1)] #CNA
      
      f <- paste(w[1,1],w[1,2],w[1,3])
      
      chr1 <- dm[match(eventpos[i,2+3*(j-1)],dm$Positions),4]
      chr2 <- dm[match(eventpos[i+1,2+3*(j-1)],dm$Positions),4]
      
      if(f != "0 0 0"){
        if(f != "X_2575 0 0"){
          if(is.na(chr1) == FALSE){
            if(is.na(chr2) == FALSE){
              if(chr1 == chr2){
                if(eventpos[i,3+3*(j-1)] == eventpos[i+1,3+3*(j-1)]){
                  pos1 <- as.numeric(match(eventpos[i,2+3*(j-1)],dm$Positions)) #Slutpositionen f?r den f?rsta h?ndelsen.
                  pos2 <- as.numeric(match(eventpos[i+1,1+3*(j-1)],dm$Positions)) #Startpositionen f?r den andra h?ndelsen.
                  
                  start <- dm[pos1+1,1] #Hittar start- och slutpositionerna f?r 2:h?ndelsen.
                  stop <- dm[pos2-1,1]
                  
                  start2 <- pos1+1 #Hittar start- och slutpositionerna för 2:h?ndelsen.
                  stop2 <- pos2-1
                  pos1_2 <- as.numeric(match(eventpos[i,1+3*(j-1)],dm$Positions)) #Startpositionen f?r den f?rsta h?ndelsen.
                  pos2_2 <- as.numeric(match(eventpos[i+1,2+3*(j-1)],dm$Positions)) #Slutpositionen f?r den andra h?ndelsen.
                  
                  if((stop2-start2)/(pos2_2-pos1_2) < 0.5){ #2-händelsen måste isf täcka mindre än 50 % av den underliggande händelsen. Annars är det inte troligt att de hör ihop.
                    
                    
                    ev1[1,1] <- eventpos[i,1+3*(j-1)] #Första händelsens startposition.
                    ev1[1,2] <- eventpos[i,2+3*(j-1)] #Första händelsens slutposition.
                    ev1[1,3] <- eventpos[i,3+3*(j-1)] #CNA
                    ev1_fused <- paste(ev1[1,1],ev1[1,2],ev1[1,3])
                    
                    ev2[1,1] <- eventpos[i+1,1+3*(j-1)] #Andra händelsens startposition.
                    ev2[1,2] <- eventpos[i+1,2+3*(j-1)] #Andra händelsens slutposition.
                    ev2[1,3] <- eventpos[i+1,3+3*(j-1)] #CNA
                    ev2_fused <- paste(ev2[1,1],ev2[1,2],ev2[1,3])
                    
                    print(ev1_fused)
                    print(ev2_fused)
                    #Nr of cells having the first event.
                    nr_places1 <- 0
                    pl1 <- as.numeric(first)
                    
                    for(pl1 in as.numeric(first):as.numeric(last)){ #Looping through the cells in this cluster
                      
                      if(ev1_fused %in% eventstot[,pl1]){ #If the event exist in a cell in the sample.
                        nr_places1 <- as.numeric(nr_places1)+ as.numeric(nrcells[1,pl1]) #We count this.
                      }
                      
                      pl1 <- pl1+1}
                    
                    #Nr of cells having the second event.
                    nr_places2 <- 0
                    pl2 <- as.numeric(first)
                    
                    for(pl2 in as.numeric(first):as.numeric(last)){ #Looping through the cells in this cluster
                      
                      if(ev1_fused %in% eventstot[,pl2]){ #If the event exist in a cell in the sample.
                        nr_places2 <- as.numeric(nr_places2)+ as.numeric(nrcells[1,pl2]) #We count this.
                      }
                      
                      pl2 <- pl2+1}
                    
                    #print("Här är nr")
                    #print(nr_places1)
                    #print(nr_places2)
                    if(nr_places1 && nr_places2 <= 3){ #If either of the two parts of the "new event" exist in at least three other cells we do not want to create this new event.
                      if(f %in% eventstot){ #Added this 200614 because I had a lot of segments with 11112221111222111. Then it grouped the first part. I rather wanted three separate in this case.
                        #print("Nu är vi i 2:loopen")
                        eventpos[event,1+3*(j-1)] <- w[1,1]
                        eventpos[event,2+3*(j-1)] <- w[1,2]
                        eventpos[event,3+3*(j-1)] <- w[1,3]
                        
                        eventpos[i,1+3*(j-1)] <- 0
                        eventpos[i,2+3*(j-1)] <- 0
                        eventpos[i,3+3*(j-1)] <- 0
                        eventpos[i+1,1+3*(j-1)] <- 0
                        eventpos[i+1,2+3*(j-1)] <- 0
                        eventpos[i+1,3+3*(j-1)] <- 0
                        event <- event+1
                        
                        eventpos[event,1+3*(j-1)] <- start #L?gger till 2:h?ndelsen.
                        eventpos[event,2+3*(j-1)] <- stop
                        eventpos[event,3+3*(j-1)] <- 2
                        event <- event+1}}
                  }}}}}
        }}}
    
    #Kombinerar två på varandra följande event. Andra eventets CNA.
    x <- matrix(0,1,3)
    x[1,1] <- eventpos[i,1+3*(j-1)]
    x[1,2] <- eventpos[i+1,2+3*(j-1)]
    x[1,3] <- eventpos[i+1,3+3*(j-1)]
    
    a <- paste(x[1,1],x[1,2],x[1,3])
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+1,2+3*(j-1)],dm$Positions),4]
    
    #Tittar på hur många celler som förändringen finns i. 200521
    # if(a != "0 0 0"){
    # pl <- 1
    # nr_places <- 0
    # for(pl in 1:ncol(eventstot)){
    #
    #   if(a %in% eventstot[,pl]){
    #   nr_places <- as.numeric(nr_places)+as.numeric(nrcells[1,pl])
    #   }
    #
    # pl <- pl+1
    # }} #sum(as.numeric(nrcells[1,])) This was used below to compare the event. This is the total number of cells.
    
    nr_places <- 0
    if(a != "0 0 0"){
      pl <- as.numeric(first)
      #print("a")
      #print(first)
      #print(last)
      for(pl in as.numeric(first):as.numeric(last)){ #Looping through the cells in this cluster
        
        if(a %in% eventstot[,pl]){ #If the event exist in a cell in the sample.
          nr_places <- as.numeric(nr_places)+ as.numeric(nrcells[1,pl]) #We count this.
        }
        
        pl <- pl+1
      }}
    #if(a != "0 0 0"){
    # if(names(data)[j] == "O_D"){
    # print(names(data)[j])
    # print(a)
    # print(nr_places)}}
    
    
    if(a != "0 0 0"){
      if(a != "X_2575 0 0"){
        if(a %in% eventstot){
          pos1 <- match(x[1,1],dm$Positions)
          pos2 <- match(x[1,2],dm$Positions)
          eventvector <- as.matrix(data[pos1:pos2,j])
          common <- sum(eventvector == x[1,3]) #Beräknar antalet element med det värdet.
          if(common/nrow(eventvector) >= 0.5 && as.numeric(nr_places)/as.numeric(tot_cells) > 1  && chr1==chr2){ #Om den andra händelsen utgör majoriteten av elementen samt om den nya händelsen finns i mer än 1 annan cell (200906).
            eventpos[event,1+3*(j-1)] <- x[1,1] #Lägger in den underliggande händelsen.
            eventpos[event,2+3*(j-1)] <- x[1,2]
            eventpos[event,3+3*(j-1)] <- x[1,3]
            
            eventpos[i+1,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande h?ndelsen som vi ser.
            eventpos[i+1,2+3*(j-1)] <- 0
            eventpos[i+1,3+3*(j-1)] <- 0
            event <- event+1}else if(as.numeric(nr_places)/as.numeric(tot_cells) > 0.5){ #Den andra händelsen är en stamhändelse i detta prov.
              # print("Vi kom in i a")
              # print(names(data)[j])
              # print(eventpos[1:3,1+3*(j-1)])
              # print(eventpos[1:3,2+3*(j-1)])
              # print(eventpos[1:3,3+3*(j-1)])
              # print(a)
              # print(event)
              eventpos[event,1+3*(j-1)] <- x[1,1] #Lägger in den underliggande h?ndelsen.
              eventpos[event,2+3*(j-1)] <- x[1,2]
              eventpos[event,3+3*(j-1)] <- x[1,3]
              
              # print(x)
              # print(i)
              # print(eventpos[i+1,1+3*(j-1)])
              # print(eventpos[i+1,2+3*(j-1)])
              # print(eventpos[i+1,3+3*(j-1)])
              
              eventpos[i+1,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande händelsen som vi ser.
              eventpos[i+1,2+3*(j-1)] <- 0
              eventpos[i+1,3+3*(j-1)] <- 0
              print(eventpos[1:3,1+3*(j-1)])
              print(eventpos[1:3,2+3*(j-1)])
              print(eventpos[1:3,3+3*(j-1)])
              event <- event+1
            }
          
          #Tystade detta 200613.
          # q <- 1
          # g <- 0
          # for(q in 1:ncol(eventstot)){
          #   if(b %in% eventstot[,q]){
          #     g <- g+1
          #   }
          # }
          #
          # if(g/ncol(eventstot) >= co){ #Om den delen ?r en stamh?ndelse.
          #   eventpos[event,1+3*(j-1)] <- y[1,1] #L?gger in den underliggande h?ndelsen.
          #   eventpos[event,2+3*(j-1)] <- y[1,2]
          #   eventpos[event,3+3*(j-1)] <- y[1,3]
          #
          #   eventpos[i,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande h?ndelsen som vi ser.
          #   eventpos[i,2+3*(j-1)] <- 0
          #   eventpos[i,3+3*(j-1)] <- 0
          #   event <- event+1}
          
        }}}
    
    y[1,1] <- eventpos[i,1+3*(j-1)] #Kombinerar två på varandra följande event. Första eventets CNA.
    y[1,2] <- eventpos[i+1,2+3*(j-1)]
    y[1,3] <- eventpos[i,3+3*(j-1)]
    b <- paste(y[1,1],y[1,2],y[1,3])
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+1,2+3*(j-1)],dm$Positions),4]
    
    #Tittar p? hur m?nga celler som f?r?ndringen finns i. 200521
    # if(b != "0 0 0"){
    # pl <- 1
    # nr_places <- 0
    # for(pl in 1:ncol(eventstot)){
    #
    #   if(b %in% eventstot[,pl]){
    #     nr_places <- as.numeric(nr_places)+as.numeric(nrcells[1,pl])
    #   }
    #
    #   pl <- pl+1
    # }}
    
    nr_places <- 0
    if(b != "0 0 0"){
      pl <- as.numeric(first)
      #print("b")
      #print(first)
      #print(last)
      for(pl in as.numeric(first):as.numeric(last)){ #Looping through the cells in this cluster
        
        if(b %in% eventstot[,pl]){ #If the event exist in a cell in the sample.
          nr_places <- as.numeric(nr_places)+as.numeric(nrcells[1,pl]) #We count this.
        }
        
        pl <- pl+1
      }}
    
    # if(b != "0 0 0"){
    #   if(names(data)[j] == "O_D"){
    #     print(names(data)[j])
    #     print(b)
    #     print(nr_places)}}
    
    
    if(b != "0 0 0"){
      if(b != "X_2575 0 0"){
        if(b %in% eventstot){
          pos1 <- match(y[1,1],dm$Positions)
          pos2 <- match(y[1,2],dm$Positions)
          eventvector <- as.matrix(data[pos1:pos2,j])
          common <- sum(eventvector == y[1,3])
          if(common/nrow(eventvector) >= 0.5 && as.numeric(nr_places)/as.numeric(tot_cells) > 1 && chr1==chr2){
            
            eventpos[event,1+3*(j-1)] <- y[1,1] #L?gger in den underliggande h?ndelsen.
            eventpos[event,2+3*(j-1)] <- y[1,2]
            eventpos[event,3+3*(j-1)] <- y[1,3]
            
            eventpos[i,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande h?ndelsen som vi ser.
            eventpos[i,2+3*(j-1)] <- 0
            eventpos[i,3+3*(j-1)] <- 0
            event <- event+1}else if(as.numeric(nr_places)/as.numeric(tot_cells) > 0.5){
              # print("Vi kom in i b")
              # print(names(data)[j])
              # print(b)
              # print(event)
              eventpos[event,1+3*(j-1)] <- y[1,1] #L?gger in den underliggande h?ndelsen.
              eventpos[event,2+3*(j-1)] <- y[1,2]
              eventpos[event,3+3*(j-1)] <- y[1,3]
              
              # print(eventpos[i,1+3*(j-1)])
              # print(eventpos[i,2+3*(j-1)])
              # print(eventpos[i,3+3*(j-1)])
              eventpos[i,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande h?ndelsen som vi ser.
              eventpos[i,2+3*(j-1)] <- 0
              eventpos[i,3+3*(j-1)] <- 0
              event <- event+1
            }
          
          #Tystade 200613
          # q <- 1
          # g <- 0
          # for(q in 1:ncol(eventstot)){
          #   if(b %in% eventstot[,q]){
          #     g <- g+1
          #   }
          # }
          #
          # if(g/ncol(eventstot) >= co){ #Om den delen ?r en stamh?ndelse.
          #   eventpos[event,1+3*(j-1)] <- y[1,1] #L?gger in den underliggande h?ndelsen.
          #   eventpos[event,2+3*(j-1)] <- y[1,2]
          #   eventpos[event,3+3*(j-1)] <- y[1,3]
          #
          #   eventpos[i,1+3*(j-1)] <- 0 #Tar bort den lilla del av den underliggande h?ndelsen som vi ser.
          #   eventpos[i,2+3*(j-1)] <- 0
          #   eventpos[i,3+3*(j-1)] <- 0
          #   event <- event+1}
        }}}
    
    
    u[1,1] <- eventpos[i,1+3*(j-1)] #Hanterar 2 events på ett större event.
    u[1,2] <- eventpos[i+4,2+3*(j-1)]
    u[1,3] <- eventpos[i+4,3+3*(j-1)]
    
    CNA1_1 <- eventpos[i,3+3*(j-1)]
    CNA1_2 <- eventpos[i+2,3+3*(j-1)]
    CNA1_4 <- eventpos[i+4,3+3*(j-1)]
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+4,2+3*(j-1)],dm$Positions),4]
    
    if(CNA1_1 == CNA1_2){
      if(CNA1_1 == CNA1_4){
        Chr1_1 <- as.numeric(dm[match(eventpos[i,2+3*(j-1)],dm$Positions),4])
        Chr1_2 <- as.numeric(dm[match(eventpos[i+2,2+3*(j-1)],dm$Positions),4])
        Chr1_4 <- as.numeric(dm[match(eventpos[i+4,2+3*(j-1)],dm$Positions),4])
        
        t <- paste(u[1,1],u[1,2],u[1,3])
        if(t != "0 0 0"){
          if(t != "X_2575 0 0"){
            if(t %in% eventstot){
              if(Chr1_1 == Chr1_2){
                if(Chr1_1 == Chr1_4){
                  eventpos[event,1+3*(j-1)] <- u[1,1]
                  eventpos[event,2+3*(j-1)] <- u[1,2]
                  eventpos[event,3+3*(j-1)] <- u[1,3]
                  
                  eventpos[i,1+3*(j-1)] <- 0
                  eventpos[i,2+3*(j-1)] <- 0
                  eventpos[i,3+3*(j-1)] <- 0
                  eventpos[i+2,1+3*(j-1)] <- 0
                  eventpos[i+2,2+3*(j-1)] <- 0
                  eventpos[i+2,3+3*(j-1)] <- 0
                  eventpos[i+4,1+3*(j-1)] <- 0
                  eventpos[i+4,2+3*(j-1)] <- 0
                  eventpos[i+4,3+3*(j-1)] <- 0
                  event <- event+1}}}
          }}
      }}
    
    
    z[1,1] <- eventpos[i,1+3*(j-1)] #Kombinerar två events som har ett annat event mellan sig (event ovanpå ett stårre event).
    z[1,2] <- eventpos[i+2,2+3*(j-1)]
    z[1,3] <- eventpos[i+2,3+3*(j-1)]
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+2,2+3*(j-1)],dm$Positions),4]
    
    if(z[1,1]=="2_212"){
      print("Här är denna, mellan.")
      print(z)
    }
    
    c <- paste(z[1,1],z[1,2],z[1,3])
    if(c != "0 0 0"){
      if(c != "X_2575 0 0"){
        
        if(c %in% eventstot && chr1==chr2){
          print("Här är c")
          print(z)
          eventpos[event,1+3*(j-1)] <- z[1,1]
          eventpos[event,2+3*(j-1)] <- z[1,2]
          eventpos[event,3+3*(j-1)] <- z[1,3]
          
          eventpos[i,1+3*(j-1)] <- 0
          eventpos[i,2+3*(j-1)] <- 0
          eventpos[i,3+3*(j-1)] <- 0
          eventpos[i+2,1+3*(j-1)] <- 0
          eventpos[i+2,2+3*(j-1)] <- 0
          eventpos[i+2,3+3*(j-1)] <- 0
          event <- event+1
        }
        else{
          if((eventpos[i+1,2+3*(j-1)] != "0") && (eventpos[i+1,1+3*(j-1)] != "0") && (eventpos[i+2,2+3*(j-1)] != "0") && (eventpos[i,1+3*(j-1)] != "0") && chr1==chr2){
            pos1 <- as.numeric(match(eventpos[i,2+3*(j-1)],dm$Positions)) #Slutpositionen för den första händelsen.
            pos2 <- as.numeric(match(eventpos[i,1+3*(j-1)],dm$Positions))
            pos3 <- as.numeric(match(eventpos[i+2,2+3*(j-1)],dm$Positions)) #Slutpositionen för den första händelsen.
            pos4 <- as.numeric(match(eventpos[i+2,1+3*(j-1)],dm$Positions))
            
            if(is.na(abs(pos4-pos3)) == FALSE){
              if(is.na(abs(pos2-pos1)) == FALSE){
                if(dm[pos1,4] == dm[pos4,4]){ #De måste befinna sig på samma kromosom
                  if(eventpos[i,3+3*(j-1)] == eventpos[i+2,3+3*(j-1)]){ #De måste ha samma CNV
                    if(abs(pos1-pos2) > 5 && abs(pos3-pos4) > 5 && abs(pos4-pos1)/abs(pos3-pos2) < 0.25){ #200906.
                      #if(abs(pos2-pos1)/abs(pos4-pos3) < 0.25){
                      print("Här inne i 0.25")
                      print(z)
                      eventpos[event,1+3*(j-1)] <- z[1,1]
                      eventpos[event,2+3*(j-1)] <- z[1,2]
                      eventpos[event,3+3*(j-1)] <- z[1,3]
                      
                      eventpos[i,1+3*(j-1)] <- 0
                      eventpos[i,2+3*(j-1)] <- 0
                      eventpos[i,3+3*(j-1)] <- 0
                      eventpos[i+2,1+3*(j-1)] <- 0
                      eventpos[i+2,2+3*(j-1)] <- 0
                      eventpos[i+2,3+3*(j-1)] <- 0
                      event <- event+1
                    }}}}}
          }
        }
        
      }}
    
    z[1,1] <- eventpos[i,1+3*(j-1)] #Kombinerar två events som har ett annat event mellan sig (event ovanpå ett större event). Första höndelsens CNV.
    z[1,2] <- eventpos[i+2,2+3*(j-1)]
    z[1,3] <- eventpos[i,3+3*(j-1)]
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+2,2+3*(j-1)],dm$Positions),4]
    
    if(z[1,1]=="2_212"){
      print("Här är denna, första")
      print(z)
    }
    
    pos1 <- as.numeric(match(eventpos[i,2+3*(j-1)],dm$Positions)) #Slutpositionen för den första händelsen. Ändrade från +1 till utan det 200906.
    pos2 <- as.numeric(match(eventpos[i,1+3*(j-1)],dm$Positions)) #Start för den första.
    pos3 <- as.numeric(match(eventpos[i+2,2+3*(j-1)],dm$Positions)) #Slutpositionen för den andra händelsen.
    pos4 <- as.numeric(match(eventpos[i+2,1+3*(j-1)],dm$Positions)) #Ändrade från +1 till +2. Start för den andra.
    
    
    c <- paste(z[1,1],z[1,2],z[1,3])
    if(c != "0 0 0"){
      if(c != "X_2575 0 0"){
        
        if(c %in% eventstot){
          if(eventpos[i,3+3*(j-1)]!=eventpos[i+2,3+3*(j-1)] && abs(pos1-pos2) > abs(pos3-pos4)){
            print("Här är c")
            print(z)
            eventpos[event,1+3*(j-1)] <- z[1,1]
            eventpos[event,2+3*(j-1)] <- z[1,2]
            eventpos[event,3+3*(j-1)] <- z[1,3]
            
            eventpos[i,1+3*(j-1)] <- 0
            eventpos[i,2+3*(j-1)] <- 0
            eventpos[i,3+3*(j-1)] <- 0
            
            event <- event+1
          }
        }
        else{
          if((eventpos[i+1,2+3*(j-1)] != "0") && (eventpos[i+1,1+3*(j-1)] != "0") && (eventpos[i+2,2+3*(j-1)] != "0") && (eventpos[i,1+3*(j-1)] != "0") && chr1==chr2){
            
            print("Kom hit")
            print(z)
            print(pos1)
            print(pos2)
            print(pos3)
            print(pos4)
            print(dm[pos3,4])
            print(dm[pos4,4])
            print(eventpos[i,3+3*(j-1)])
            print(eventpos[i+2,3+3*(j-1)])
            
            if(is.na(abs(pos4-pos3)) == FALSE){
              if(is.na(abs(pos2-pos1)) == FALSE){
                if(dm[pos3,4] == dm[pos4,4]){ #De måste befinna sig på samma kromosom
                  if(eventpos[i,3+3*(j-1)] != eventpos[i+2,3+3*(j-1)]){ #De måste ha olika CNV.
                    if(abs(pos1-pos2) > 5 && abs(pos3-pos4) > 5 && abs(pos4-pos1)/abs(pos3-pos2) < 0.25 && abs(pos1-pos2) > abs(pos3-pos4)){ #200906.
                      #if(abs(pos2-pos1)/abs(pos4-pos3) < 0.25){
                      print("H?r inne i 0.25")
                      print(z)
                      eventpos[event,1+3*(j-1)] <- z[1,1]
                      eventpos[event,2+3*(j-1)] <- z[1,2]
                      eventpos[event,3+3*(j-1)] <- z[1,3]
                      
                      eventpos[i,1+3*(j-1)] <- 0
                      eventpos[i,2+3*(j-1)] <- 0
                      eventpos[i,3+3*(j-1)] <- 0
                      
                      event <- event+1
                    }}
                }}}
          }
        }
        
      }}
    
    z[1,1] <- eventpos[i,1+3*(j-1)] #Kombinerar två events som har ett annat event mellan sig (event ovanpå ett större event). Andra händelsens CNV.
    z[1,2] <- eventpos[i+2,2+3*(j-1)]
    z[1,3] <- eventpos[i+2,3+3*(j-1)]
    
    chr1 <- dm[match(eventpos[i,1+3*(j-1)],dm$Positions),4]
    chr2 <- dm[match(eventpos[i+2,2+3*(j-1)],dm$Positions),4]
    
    if(z[1,1]=="2_212"){
      print("Här är denna, andra.")
      print(z)
    }
    
    pos1 <- as.numeric(match(eventpos[i,2+3*(j-1)],dm$Positions)) #Slutpositionen f?r den f?rsta h?ndelsen. ?ndrade fr?n +1 till utan det 200906.
    pos2 <- as.numeric(match(eventpos[i,1+3*(j-1)],dm$Positions)) #Start f?r den f?rsta.
    pos3 <- as.numeric(match(eventpos[i+2,2+3*(j-1)],dm$Positions)) #Slutpositionen f?r den andra h?ndelsen.
    pos4 <- as.numeric(match(eventpos[i+2,1+3*(j-1)],dm$Positions)) #?ndrade fr?n +1 till +2. Start f?r den andra.
    
    
    c <- paste(z[1,1],z[1,2],z[1,3])
    if(c != "0 0 0"){
      if(c != "X_2575 0 0"){
        
        if(c %in% eventstot){
          if(eventpos[i,3+3*(j-1)]!=eventpos[i+2,3+3*(j-1)] && abs(pos3-pos4) > abs(pos1-pos2)){
            print("Här är c")
            print(z)
            eventpos[event,1+3*(j-1)] <- z[1,1]
            eventpos[event,2+3*(j-1)] <- z[1,2]
            eventpos[event,3+3*(j-1)] <- z[1,3]
            
            eventpos[i+2,1+3*(j-1)] <- 0
            eventpos[i+2,2+3*(j-1)] <- 0
            eventpos[i+2,3+3*(j-1)] <- 0
            event <- event+1
          }
        }
        else{
          if((eventpos[i+1,2+3*(j-1)] != "0") && (eventpos[i+1,1+3*(j-1)] != "0") && (eventpos[i+2,2+3*(j-1)] != "0") && (eventpos[i,1+3*(j-1)] != "0") && chr1==chr2){
            
            print("Kom hit")
            print(z)
            print(pos1)
            print(pos2)
            print(pos3)
            print(pos4)
            print(dm[pos3,4])
            print(dm[pos4,4])
            print(eventpos[i,3+3*(j-1)])
            print(eventpos[i+2,3+3*(j-1)])
            
            if(is.na(abs(pos4-pos3)) == FALSE){
              if(is.na(abs(pos2-pos1)) == FALSE){
                if(dm[pos1,4] == dm[pos4,4]){ #De m?ste befinna sig p? samma kromosom
                  if(eventpos[i,3+3*(j-1)] != eventpos[i+2,3+3*(j-1)]){ #De m?ste ha olika CNV.
                    if(abs(pos1-pos2) > 5 && abs(pos3-pos4) > 5 && abs(pos4-pos1)/abs(pos3-pos2) < 0.25 && abs(pos3-pos4) > abs(pos1-pos2)){ #200906.
                      #if(abs(pos2-pos1)/abs(pos4-pos3) < 0.25){
                      print("H?r inne i 0.25")
                      print(z)
                      eventpos[event,1+3*(j-1)] <- z[1,1]
                      eventpos[event,2+3*(j-1)] <- z[1,2]
                      eventpos[event,3+3*(j-1)] <- z[1,3]
                      
                      eventpos[i+2,1+3*(j-1)] <- 0
                      eventpos[i+2,2+3*(j-1)] <- 0
                      eventpos[i+2,3+3*(j-1)] <- 0
                      event <- event+1
                    }}
                }}}
          }
        }
        
      }}
    
    i <- i+1
  }
  j <- j+1
}

eventstot <- matrix(0,nrow(eventpos),ncol(data)) #Fuserar eventpositionerna till en enda bin.
for(j in 1:ncol(data)){
  eventstot[,j] <- paste(eventpos[,1+3*(j-1)],eventpos[,2+3*(j-1)],eventpos[,3+3*(j-1)])
  j <- j+1
}

events <- as.matrix(unique(c(eventstot))) #Plockar ut alla unika händelser.
events <- as.matrix(events[events != "0 0 0"]) #Tar bort händelser innehållande bara nollor.
EMr <- nrow(events)+1
EMc <- ncol(data)+1
eventmatrix <- matrix(0,EMr,EMc) #Skapar en tom händelsematris.
eventmatrix[2:EMr,1] <- events[,1] #Händelsenamnen placeras i den första kolumnen i händelsematrisen.
eventmatrix[1,2:EMc] <- as.matrix(nm) #Cellnamnen placeras på första raden i händelsematrisen.

#############################################################
#Creating the event matrix from the event data for each cell#
#############################################################

i = 1
j = 1
for(j in 1:ncol(eventstot)){
  for(i in 1:nrow(events)){
    
    if(events[i,1] %in% eventstot[,j] == TRUE){
      eventmatrix[i+1,j+1] <- 1
    }
    else{eventmatrix[i+1,j+1] <- 0}
    
    i <- i+1
  }
  j <- j+1
}

write.table(eventmatrix, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)

################################################################
#Stem events and events consisting of 2:s on top of stem events#
################################################################

i = 2
j = 1
v <- matrix(0,200,3)
erase <- matrix(0,200,3)
stemmatrix <- matrix(0,100,1)
#posvector <- matrix(0,(as.numeric(ncol(eventmatrix))-1),2)
r = 1
s = 1
t = 1
x = 1
for(i in 2:nrow(eventmatrix)){ #Hittar stamhändelserna.
  
  if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) >= co){ #ncol(eventmatrix)-1 stod det innan.
    stemmatrix[s,1] <- eventmatrix[i,1]
    s <- s+1
    eventmatrix[i,2:ncol(eventmatrix)] <- 1
    
    pos1 <- match(word(eventmatrix[i,1],1),dm$Positions)
    pos2 <- match(word(eventmatrix[i,1],2),dm$Positions)
    
    for(j in 1:ncol(eventstot)){ #Loopar igenom varje cell över just det eventet.
      for(k in pos1:pos2){
        if(data[k,j]==2){ #Hanterar händelser med 2:or på stamhhndelsen.
          
          if(v[r,1] == 0){
            v[r,1] <- j
            v[r,2] <- k
            if(k != pos1){
              erase[t,1] <- dm[pos1,1]
              erase[t,2] <- dm[as.numeric(k)-1,1]
              erase[t,3] <- word(eventmatrix[i,1],3)
              t <- t+1}
          }
          else if(v[r,1] != j){
            r <-r+1
            v[r,1] <- j
            v[r,2] <- k
            
            if(k != pos1){
              t <- t+1
              erase[t,1] <- dm[pos1,1]
              erase[t,2] <- dm[as.numeric(k)-1,1]
              erase[t,3] <- word(eventmatrix[i,1],3)
              t <- t+1}
          }
          
          if(k < 2721){
            v[r,3] <- k
            if(data[k+1,j]!=2){r <-r+1
            if(k != pos2){
              erase[t,1] <- dm[as.numeric(k+1),1]
              erase[t,2] <- dm[pos2,1]
              erase[t,3] <- word(eventmatrix[i,1],3)
              t <- t+1
            }
            }}
        }
        k <- k+1
      }
      j <- j+1
    }
  }
  eventstwo <- as.matrix(paste(dm[v[,2],1],dm[v[,3],1],2))
  eventstwo <- as.matrix(eventstwo[eventstwo != "0 0 0"])
  erasetwo <- as.matrix(paste(erase[,1],erase[,2],erase[,3]))
  erasetwo <- as.matrix(erasetwo[erasetwo != "0 0 0"])
  i <- i+1
}

eventmatrixtwo <- matrix(0,nrow(eventstwo),EMc) #Skapar en tom händelsematrisdel.
eventmatrixtwo[,1] <- eventstwo

for(i in 1:nrow(eventstwo)){ #Lägger till tvähändelserna till rätt cell.
  eventmatrixtwo[i, v[i,1]+1] <- 1
}

i = 1
j = 1
for(i in 1:nrow(eventmatrixtwo)){ #Hanterar dubletter av händelser i tvåmatrisen.
  for(j in 1:nrow(eventmatrixtwo)){
    if(i != j){
      if(eventmatrixtwo[i,1] == eventmatrixtwo[j,1]){
        eventmatrixtwo[i,2:ncol(eventmatrixtwo)] <- as.numeric(eventmatrixtwo[i,2:ncol(eventmatrixtwo)])+as.numeric(eventmatrixtwo[j,2:ncol(eventmatrixtwo)])
        eventmatrixtwo[j,] <- 0
      }}
    j <- j+1
  }
  i <- i+1
}

i = 1
j = 2
for(i in 1:nrow(eventmatrixtwo)){ #Normaliserar matrisen.
  for(j in 2:ncol(eventmatrixtwo)){
    if(eventmatrixtwo[i,j] != "0"){
      eventmatrixtwo[i,j] <- as.character(as.numeric(eventmatrixtwo[i,j])/as.numeric(eventmatrixtwo[i,j]))
    }
    j <- j+1
  }
  i <- i+1
}

eventmatrix <- rbind(eventmatrix,eventmatrixtwo) #Fuserar den tidigare EM:n med den nya med tvåhändelserna.

if(length(erasetwo) != 0){
  i = 1
  for(i in 1:nrow(erasetwo)){ #Plockar bort de händelser som ligger omkring 2:händelsen på stammen.
    eventmatrix <- eventmatrix[eventmatrix[,1] != erasetwo[i,1],]
    i <- i+1
  }}

###############################
#Treating isochromosome events#
###############################

i = 2
k = 1
isochr <- matrix(0,4,3)
isochrE <- 0
isochrEM_13 <- matrix(0,10,as.numeric(ncol(eventmatrix)))
isochrEM_31 <- matrix(0,10,as.numeric(ncol(eventmatrix)))

for(i in 2:(as.numeric(nrow(eventmatrix))-2)){
  
  stop <- 0
  chromosome1 <- 0
  chromosome3 <- 0
  
  pos1_1 <- as.numeric(match(word(eventmatrix[i,1],1),dm$Positions)) #Startpositionen för p-händelsen.
  pos1_2 <- as.numeric(match(word(eventmatrix[i,1],2),dm$Positions)) #Slutpositionen för p-händelsen.
  CNA1 <- as.numeric(word(eventmatrix[i,1],3)) #Kopietalet i händelsen.
  chrom1 <- as.numeric(dm[pos1_1,4]) #Kromosomtillhörighet.
  
  if(CNA1 == 1){ #Är kopietalet 1?
    if(pos1_1 != 1){ #Är vi inte på den första binen i kromosomen?
      if(dm[pos1_1,4] != dm[pos1_1-1,4]){ #Är vi i kromosomens början?
        if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) == 1){ #Rör det sig om en stamhändelse?
          j = pos1_1
          for(j in pos1_1:2721){ #Hittar slutpositionen av kromosomen.
            if(dm[j,4] == chrom1){
              stop <- dm[j,1] #Sparar ner slutpositionen
              chromosome1 <- dm[j,4]}
            j <- j+1
          }
          start <- dm[pos1_2+1,1]
          if(paste(start,stop,3) %in% eventmatrix[,1]){
            isochr[k,1] <- start
            isochr[k,2] <- stop
            isochr[k,3] <- 3
            isochrE <- paste(isochr[k,1],isochr[k,2],isochr[k,3])
            isochrEM_13[k,1] <- isochrE
            isochrEM_13[k,2:ncol(isochrEM_13)] <- 1
            k <- k+1}}
      }
    }
  }
  else{ #Om vi är på position 1 av datan.
    j = pos1_1
    if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) == 1){
      for(j in pos1_2:2721){ #Hittar slutpositionen av kromosomen.
        if(dm[j,4] == chrom1){
          stop <- dm[j,1]}
        j <- j+1
      }
      start <- dm[pos1_2+1,1]
      position <- paste(start,stop,3)
      if(position %in% eventmatrix[,1]){
        rowevent <- match(position,eventmatrix[,1])
        if(sum(as.numeric(eventmatrix[rowevent,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) > 0.5){
          isochr[k,1] <- start
          isochr[k,2] <- stop
          isochr[k,3] <- 3
          isochrE <- paste(isochr[k,1],isochr[k,2],isochr[k,3])
          isochrEM_13[k,1] <- isochrE
          isochrEM_13[k,2:ncol(isochrEM_13)] <- 1
          k <- k+1}}
    }
  }
}

if(CNA1 == 3){
  if(dm[pos1_1,4] != 23){
    if(dm[pos1_2,4] != dm[pos1_2+1,4]){ #Vi är i kromosomens slut.
      if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) == 1){
        j = pos1_1
        for(j in pos1_1:1){
          if(dm[j,4] == chrom1){
            start <- dm[j,1]
            chromosome3 <- dm[j,4]}
          j <- j+1}
        stop <- dm[pos1_1-1,1]
        position <- paste(start,stop,1)
        if(position %in% eventmatrix[,1]){
          rowevent <- match(position,eventmatrix[,1])
          if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) > 0.5){
            isochr[k,1] <- start
            isochr[k,2] <- stop
            isochr[k,3] <- 1
            isochrE <- paste(isochr[k,1],isochr[k,2],isochr[k,3])
            isochrEM_31[k,1] <- isochrE
            isochrEM_31[k,2:ncol(isochrEM_31)] <- 1
            k <- k+1}}
      }
    }
  }
  else{ #Om vi är i position 1 av datan.
    j = pos1_1
    if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) == 1){
      for(j in pos1_1:2721){ #Hittar slutpositionen av kromosomen.
        if(dm[j,4] == chrom1){
          stop <- dm[j,1]
          chromosome1 <- dm[j,4]}
        j <- j+1
      }
      start <- dm[pos1_2+1,1]
      position <- paste(start,stop,1)
      if(position %in% eventmatrix[,1]){
        rowevent <- match(position,eventmatrix[,1])
        if(sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) > 0.5){
          isochr[k,1] <- start
          isochr[k,2] <- stop
          isochr[k,3] <- 3
          isochrE <- paste(isochr[k,1],isochr[k,2],isochr[k,3])
          isochrEM_31[k,1] <- isochrE
          isochrEM_31[k,2:ncol(isochrEM_31)] <- 1
          k <- k+1}
      }
    }
  }
  i <- i+1
}

eventmatrix[1,1] <- "Cells"
eventmatrix <- as.matrix(eventmatrix[eventmatrix[,1] != isochrE,]) #Tar bort händelser innehållande bara nollor.
eventmatrix <- rbind(eventmatrix,isochrEM_13)
eventmatrix <- rbind(eventmatrix,isochrEM_31)
eventmatrix <- as.matrix(eventmatrix[eventmatrix[,1] != 1,])
eventmatrix <- as.matrix(eventmatrix[eventmatrix[,1] != 0,])
events <- as.matrix(unique(eventmatrix[,1])) #Plockar ut alla unika händelser.
eventmatrix <- eventmatrix[!duplicated(eventmatrix), ]

eventnumber <- matrix(0,nrow(eventmatrix))
i <- 2
for(i in 2:nrow(eventmatrix)){
  eventnumber[i,1] <- as.numeric(match(word(eventmatrix[i,1],1), dm$Positions))
  i <- i+1
}

eventmatrixx <- cbind(eventnumber,eventmatrix)
eventmatrix <- eventmatrix[order(as.numeric(eventmatrixx[,1])),]

# stemcell <- matrix(0,as.numeric(nrow(eventmatrix)))
# stemcell[1,1] <- "Stemcell"
# i <- 2
# for(i in 2:as.numeric(nrow(eventmatrix))){
#   x <- sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:as.numeric(ncol(nrcells))]))
#   if(x/((sum(as.numeric(nrcells)))) == 1){
#     stemcell[i,1] <- 1
#     print(eventmatrix[i,1])
#   }
#   else{stemcell[i,1] <- 0}
#   i <- i+1
# }
#
# eventmatrix_withX <- cbind(eventmatrix,stemcell)
#eventmatrix <- eventmatrix_withX[eventmatrix_withX[,1] != "X_2575 X_2721 1",]
eventmatrix <- eventmatrix[eventmatrix[,1] != "X_2575 X_2721 1",]

#EM <- eventmatrix[2:as.numeric(nrow(eventmatrix)),2:as.numeric(ncol(eventmatrix))] #Plockar ut h?ndelserna ur matrisen.
eventmatrix_sparad <- eventmatrix
#############################
#Duplications of chromosomes#
#############################
i <- 2
#eventmatrix <- eventmatrix_sparad
for(i in 2:(as.numeric(nrow(eventmatrix))-1)){ #Hittar händelser som stämmer in på kriteriet
  if(sum(as.numeric(eventmatrix[i,2:as.numeric(ncol(eventmatrix))])*as.numeric(nrcells[1,1:ncol(nrcells)]))/(sum(as.numeric(nrcells))) >= 0.5){ #Anger att händelsen måste finnas hos minst hälften av cellerna.
    numbercells <- sum(as.numeric(eventmatrix[i,2:as.numeric(ncol(eventmatrix))])*as.numeric(nrcells[1,1:ncol(nrcells)])) #Antalet celler som har förändringen.
    columns_event <- sum(as.numeric(eventmatrix[i,2:as.numeric(ncol(eventmatrix))]))
    duplication_EM <- matrix(0,6,as.numeric(ncol(eventmatrix))) #Gör en ny liten händelsematris som jag ska lägga in de eventuella nya händelserna i.
    
    pos1 <- match(word(eventmatrix[i,1],1),dm[,1]) #Startpositionen
    pos2 <- match(word(eventmatrix[i,1],2),dm[,1]) #Slutpositionen
    Chr_pos <- dm[pos2,4]
    print("Första")
    print(pos1)
    print(pos2)
    print(Chr_pos)
    if(dm[pos2+1,4] == Chr_pos){ #Händelsen får inte täcka hela kromosomen. Annars kommer vi gå in över nästa kromosom.
      
      s <- 1
      t <- 1
      j <- 2
      pos3 <- pos2+1
      print("pos3")
      print(pos3)
      for(j in 2:(as.numeric(ncol(eventmatrix))-1)){ #Loopar genom cellerna
        
        if(eventmatrix[i,j] == 1){ #Om en cell har händelsen. Allt inom denna är för en enda cell.
          k <- pos2+1
          pos3 <- pos2 + 1 #Startpositionen av den andra händelsen.
          
          # print("pos3")
          # print(pos3)
          
          if(t == 1){
            CNV_x_m <- matrix(0,2,columns_event)
            CNV_x_m[1,t] <- as.numeric(data[pos2,j])}
          
          CNV_x_m[1,t] <- as.numeric(data[pos2,j]) #Kopietalet för den första händelsen.
          t <- t+1
          
          for(k in (pos2+1):2721){ #Lokaliserar slutpositionen av kromosomen.
            if(dm[k,4] != dm[k+1,4]){
              pos4 <- k #Slutpositionen av den andra händelsen.
              #print("Slut")
              #print(pos4)
              if(s == 1){
                submatrix <- matrix(0,(pos4-pos3+2),columns_event)} #Skapar en matris för att analysera den andra händelsen.
              
              if(sum(as.numeric(data[(pos2+1):k,j-1]))/(k-pos2) == data[pos2+1,j-1]){ #Händelsen måste ha en och samma kopietalsförändring.
                
                submatrix[1,s] <- eventmatrix[1,j] #Sparar ner cellen och händelsen i en matris.
                submatrix[2:as.numeric(nrow(submatrix)),s] <- data[pos3:pos4,j-1]
                s <- s+1
                
              }
              break
            }
            else(k <- k+1)
          }
        }
        j <- j+1}
      
      submatrix[1,] <- "0"
      class(submatrix) <- "numeric"
      print("Här")
      
      m <- colSums(submatrix[2:nrow(submatrix),])/(pos4-pos2)
      uy <- unique(m)
      tab <- tabulate(match(m, uy))
      CNV_y <- uy[tab == max(tab)]
      print(submatrix)
      print(pos4)
      print(pos2)
      print(uy)
      print(tab)
      print(max(tab))
      print(CNV_y)
      
      n <- colSums(CNV_x_m[1:2,])
      ux <- unique(n)
      tab <- tabulate(match(n, ux))
      CNV_x <- ux[tab == max(tab)]
      
      l <- 2
      for(l in 2:(as.numeric(ncol(eventmatrix))-1)){ #Loopar genom cellerna
        
        if(as.numeric(eventmatrix[i,l]) == 0){ #Om en cell inte har händelsen
          
          if((sum(as.numeric(data[pos1:pos2,l-1]))/(pos2-(pos1-1))) != CNV_x){ #De får inte ha samma.
            if((sum(as.numeric(data[pos3:pos4,l-1]))/(pos4-pos2)) != CNV_y){
              
              if((sum(as.numeric(data[pos1:pos2,l-1]))/(pos2-(pos1-1))) %% as.numeric(CNV_x) == 0){ #Hittar de celler som har en duplikation av de andra cellernas kopietal.
                if((sum(as.numeric(data[pos3:pos4,l-1]))/(pos4-pos2)) %% as.numeric(CNV_y) == 0){
                  
                  if(paste(dm[pos1,1],dm[pos2,1],CNV_x) %in% eventmatrix[,1]){ #1 original
                    row <- match(paste(dm[pos1,1],dm[pos2,1],CNV_x), eventmatrix[,1])
                    eventmatrix[row,l] <- 1
                  }
                  else{
                    if(as.numeric(CNV_x) != 2){
                      duplication_EM[2,1] <- paste(dm[pos1,1],dm[pos2,1],CNV_x)
                      duplication_EM[2,l] <- 1
                    }
                  }
                  
                  if(paste(dm[pos3,1],dm[pos4,1],CNV_y) %in% eventmatrix[,1]){ #2 original
                    row <- match(paste(dm[pos3,1],dm[pos4,1],CNV_y), eventmatrix[,1])
                    eventmatrix[row,l] <- 1
                  }
                  else{
                    if(as.numeric(CNV_y) != 2){
                      duplication_EM[3,1] <- paste(dm[pos3,1],dm[pos4,1],CNV_y)
                      duplication_EM[3,l] <- 1
                    }
                  }
                  
                  if(paste(dm[pos1,1],dm[pos2,1],(sum(as.numeric(data[pos1:pos2,l-1]))/(pos2-(pos1-1)))) %in% eventmatrix[,1]){ #1 duplikation
                    row <- match(paste(dm[pos1,1],dm[pos2,1],(sum(as.numeric(data[pos1:pos2,l-1]))/(pos2-(pos1-1)))), eventmatrix[,1])
                    eventmatrix[row,l] <- 1
                  }
                  else{
                    duplication_EM[4,1] <- paste(dm[pos1,1],dm[pos2,1],(sum(as.numeric(data[pos1:pos2,l-1]))/(pos2-(pos1-1))))
                    duplication_EM[4,l] <- 1
                  }
                  
                  if(paste(dm[pos3,1],dm[pos4,1],(sum(as.numeric(data[pos3:pos4,l-1]))/(pos4-pos2))) %in% eventmatrix[,1]){ #2 duplikation
                    row <- match(paste(dm[pos3,1],dm[pos4,1],(sum(as.numeric(data[pos3:pos4,l]))/(pos4-pos2))), eventmatrix[,1])
                    eventmatrix[row,l] <- 1
                  }
                  else{
                    duplication_EM[5,1] <- paste(dm[pos3,1],dm[pos4,1],(sum(as.numeric(data[pos3:pos4,l-1]))/(pos4-pos2)))
                    duplication_EM[5,l] <- 1
                  }
                  
                }
              }
            }}}
        
        l <- l+1}
      
      u <- 1
      
      for(u in 1:nrow(duplication_EM)){
        if(duplication_EM[u,1] != 0){
          if(duplication_EM[u,1] %in% eventmatrix[,1] == FALSE){
            print("Duplication matrix")
            print(duplication_EM)
            eventmatrix <- rbind(eventmatrix,duplication_EM[duplication_EM[,1] != "0",])
          }}
        u <- u+1}
    } #Avlutning för en händelse.
    
  }
  i <- i+1}

eventmatrix_duplications <- eventmatrix

#############################
#Treating consecutive events#
#############################
i <- 2
s <- 2
Breakpoints <- matrix(0,nrow(eventmatrix),ncol(eventmatrix))
Breakpoints[1,] <- eventmatrix[1,1:as.numeric(ncol(eventmatrix))]

for(i in 2:nrow(eventmatrix)){
  Breakpointstartrow <- match(word(eventmatrix[i,1],1,2),word(eventmatrix[,1],1,2)) #Väljer ut en händelse och jämför den med alla andra händelser i EM.
  if(Breakpointstartrow != i){ #Om den har jämförts med sig själv ska den inte fortsätta.
    Breakpointstart <- word(eventmatrix[i,1],1)
    Breakpointstop <- word(eventmatrix[i,1],2)
    
    CNV1 <- word(eventmatrix[i,1],3) #Utgångspunkten
    CNV2 <- word(eventmatrix[as.numeric(Breakpointstartrow),1],3) #Den matchade händelsen
    
    if(eventmatrix[Breakpointstartrow,1] %in% Breakpoints[,1] == FALSE){
      Breakpoints[s,] <- eventmatrix[as.numeric(Breakpointstartrow),]
      s <- s+1
    }
    
    if(eventmatrix[i,1] %in% Breakpoints[,1] == FALSE){
      Breakpoints[s,] <- eventmatrix[i,]
      s <- s+1}
    
    # print(c("This is i",i))
    # print(CNV2)
    # print(c(Breakpointstart, Breakpointstop,CNV1, i))
    #
    # Breakpoints[s,] <- eventmatrix[i,]
    # Breakpoints[s+1] <- eventmatrix[as.numeric(Breakpointstart),]
    # s <- s+1
  }
  i <- i+1}

i <- 2
s <- 2
Breakpoints_new <- matrix(0,10,ncol(eventmatrix))
Breakpoints_new[1,] <- eventmatrix[1,1:as.numeric(ncol(eventmatrix))]
Breakpoints_new_cell <- matrix(0,10,ncol(eventmatrix))
eventmatrix_modal <- eventmatrix

for(i in 3:nrow(Breakpoints)){ #Loopar över alla händelser som har delade breakpoints med andra händelser.
  
  if(is.na(word(Breakpoints[i,1],1,2)) == FALSE){
    if(word(Breakpoints[i-1,1],1,2)==word(Breakpoints[i,1],1,2)){
      
      Breakpoints_new[s,] <- Breakpoints[i-1,]
      Breakpoints_new[s+1,] <- Breakpoints[i,]
      s <- s+1
    }
    
    if(word(Breakpoints[i-1,1],1,2) != word(Breakpoints[i,1],1,2)){ #Vi är i slutet av det området.
      j <- 2
      for(j in 2:as.numeric(nrow(Breakpoints_new))){
        Breakpoints_new_cell[j,2:as.numeric(ncol(Breakpoints_new))] <- as.numeric(as.numeric(Breakpoints_new[j,2:as.numeric(ncol(Breakpoints_new))])*as.numeric(nrcells[1,1:ncol(nrcells)]))
        j <- j+1
      }
      #Beräknar modaltalet
      thesum <- rowSums(Breakpoints_new_cell[2:as.numeric(nrow(Breakpoints_new)),2:as.numeric(ncol(Breakpoints_new))]) #Beräknar alla radsummor för den brottpunkten samt tar reda på vilken rad den tillhör.
      
      modal <- which.max(rowSums(Breakpoints_new_cell[2:as.numeric(nrow(Breakpoints_new)),2:as.numeric(ncol(Breakpoints_new))])) #Beräknar alla radsummor för den brottpunkten samt tar reda på vilken rad den tillhör.
      rownumber_modalevent <- match(Breakpoints_new[modal+1,1],eventmatrix[,1]) #På vilken rad ligger denna händelsen?
      
      if(thesum > 10){ #Adding this cutoff so that we do not apply the rule for individual cells.
        k <- 2
        for(k in 2:as.numeric(nrow(Breakpoints_new))){ #Om en cell har en "icke-modal" händelse men inte den modala sätts en 1:a där.
          if(Breakpoints_new[k,1] != "0"){
            rownumber_not_modalevent <- match(Breakpoints_new[k,1],eventmatrix_modal[,1])
            if(rownumber_modalevent != rownumber_not_modalevent){
              
              l <- 2
              for(l in 2:ncol(eventmatrix_modal)){
                if(eventmatrix_modal[rownumber_not_modalevent,l] == 1){
                  
                  if(eventmatrix_modal[rownumber_modalevent,l] == 0){
                    eventmatrix_modal[rownumber_modalevent,l] <- 1
                  }
                }
                l <- l+1}
            }}
          k <- k+1}}
      
      s <- 2 #Återställer breakpointmatrisen.
      Breakpoints_new <- matrix(0,10,ncol(eventmatrix))
      Breakpoints_new[1,] <- eventmatrix[1,1:as.numeric(ncol(eventmatrix))]
      Breakpoints_new_cell <- matrix(0,10,ncol(eventmatrix))
    }}
  i <- i+1
}

eventmatrix_original <- eventmatrix

eventmatrix <- eventmatrix_modal #New EM. Updated with consecutive events.
##################################################
#Ny stamkod - Stamhändelsefördelning inom grupper#
##################################################
#eventmatrix <- eventmatrix_modal
#eventmatrix <- eventmatrix_original
co <- cutoff #Changed from 1 to cutoff 2101030.
new_events <- matrix(0,20,as.numeric(ncol(eventmatrix))+1)

samples <- lapply(samples, as.character)
i <- 1
j <- 1
s <- 1
nr_groups <- length(samples) #The number of groups.
for(i in 2:nrow(eventmatrix)){ #Looping through the events.
  j <- 1
  
  for(j in 1:nr_groups){ #Looping through the groups of cells.
    cluster_cells <- unlist(samples[j])
    first <- as.numeric(cluster_cells[1]) #Start position
    last <- as.numeric(cluster_cells[length(cluster_cells)]) #Stop position
    cells_total <- sum(as.numeric(eventmatrix[i,(first+1):(last+1)])*as.numeric(nrcells[1,as.numeric(first):as.numeric(last)]))
    
    if(cells_total/sum(as.numeric(nrcells[1,first:last])) >= co){
      print("The event.")
      print(eventmatrix[i,1])
      
      #stemcell[i,1] <- 1
      
      k <- first+1
      for(k in (as.numeric(first)+1):(as.numeric(last)+1)){ #Looping though the cells in the group.
        
        if(eventmatrix[i,k] == "0"){ #This cell does not have the stem event.
          print("This cell does not have it")
          print(eventmatrix[1,k])
          
          pos1 <- as.numeric(match(word(eventmatrix[i,1],1),dm$Positions)) #Start position of stem event.
          pos2 <- as.numeric(match(word(eventmatrix[i,1],2),dm$Positions)) #End position of stem event.
          CNV_ny <- data[pos1,k-1] #New event. The CNV in the beginning.
          CNV_ny_end <- data[pos2,k-1] #New event. The CNV in the other end.
          print(CNV_ny)
          print(CNV_ny_end)
          
          eventmatrix[i,k] <- 1 #Adding a 1 to the stem event.
          
          if(CNV_ny == CNV_ny_end){ #The whole stem event has been replaced by another event.
            thenewevent <- paste(word(eventmatrix[i,1],1,2),CNV_ny) #The name of the new event.
            if(thenewevent %in% new_events[,1]){
              print("It already exist in new_events")
              row <- match(thenewevent,new_events[,1])
              new_events[row,k] <- 1 #Adding the new event.
            }else if(thenewevent %in% eventmatrix[,1]){
              print("It already exist in the event matrix")
              row <- match(thenewevent,eventmatrix[,1])
              eventmatrix[row,k] <- 1
              print(eventmatrix[i,1])
              print(eventmatrix[1,k])
            }else{
              print("This is a new event entirely")
              new_events[s,1] <- paste(word(eventmatrix[i,1],1,2),CNV_ny)
              new_events[s,k] <- 1
              s <- s+1}
          }else{ #They are different.
            print("They are different")
            pos_beginning <- "0"
            pos_end <- "0"
            switch <- 0
            CNV_stem <- word(eventmatrix[i,1],3)
            
            if(CNV_ny == "2"){
              pos_beginning <- pos1
              m <- 1
              for(m in pos1:pos2){
                if(switch == 0){
                  if(data[m,k-1] != "2"){
                    pos_end <- m-1
                    switch <- 1
                  }}
                m <- m+1
              }
            }
            
            if(CNV_ny_end == "2"){
              pos_end <- pos2
              m <- 1
              for(m in pos1:pos2){
                if(switch == 0){
                  if(data[m,k-1] == "2"){
                    pos_beginning <- m
                    switch <- 1
                  }}
                m <- m+1
              }
            }
            
            if(CNV_ny == CNV_stem){ #First part has the same CNV as the stem. Remove this event since we already have the stem.
              print("The first part has the same CNV as the stem.")
              switch <- 0
              pos_beginning_remove <- pos1
              m <- 1
              for(m in pos1:pos2){
                if(switch == 0){
                  if(data[m,k-1] != CNV_stem){
                    pos_end_remove <- m-1
                    switch <- 1
                  }}
                m <- m+1
              }
              eventtoremove <- paste(dm[pos_beginning_remove,1],dm[pos_end_remove,1],CNV_stem)
              rowtoremove <- match(eventtoremove,eventmatrix[,1])
              eventmatrix[rowtoremove,] <- 0
            }
            
            if(CNV_ny_end == CNV_stem){ #Last part has the same CNV as the stem. Remove this event since we already have the stem.
              print("The last part has the same CNV as the stem.")
              switch <- 0
              pos_end_remove <- pos2
              m <- 1
              for(m in pos1:pos2){
                if(switch == 0){
                  if(data[m,k-1] == CNV_stem){
                    pos_beginning_remove <- m
                    switch <- 1
                  }}
                m <- m+1
              }
              eventtoremove <- paste(dm[pos_beginning_remove,1],dm[pos_end_remove,1],CNV_stem)
              rowtoremove <- match(eventtoremove,eventmatrix[,1])
              eventmatrix[rowtoremove,] <- 0
              print(eventtoremove)
              print(rowtoremove)
            }
            
            # thenewevent <- paste(dm[pos_beginning,1],dm[pos_end,1],"2") #The name of the new event.
            # print(eventmatrix[1,k])
            # print(thenewevent)
            # print(pos_beginning)
            # print(pos_end)
            # if(pos_beginning != "0" && pos_end != "0"){
            # if(thenewevent %in% new_events[,1]){
            #   print("It already exist in new_events")
            #   row <- match(thenewevent,new_events[,1])
            #   new_events[row,k] <- 1 #Adding the new event.
            #
            # }else if(thenewevent %in% eventmatrix[,1]){
            #   print("It already exist in the event matrix")
            #   row <- match(thenewevent,eventmatrix[,1])
            #   eventmatrix[row,k] <- 1
            #   print(eventmatrix[i,1])
            #   print(eventmatrix[1,k])
            #
            # }else{
            #   print("This is a new event entirely")
            #   new_events[s,1] <- thenewevent
            #   new_events[s,k] <- 1
            #   s <- s+1}}
            
          }
          
        }
        
        k <- k+1
      }
      
    }
    
    j <- j+1
  }
  i <- i+1
}

eventmatrix <- rbind(eventmatrix,new_events[new_events[,1]!="0",1:(ncol(new_events)-1)])

##########################################################
#Treating events with the same cutoff but different CNV:s#
##########################################################
event_order <- matrix(0,nrow(eventmatrix),ncol(eventmatrix)+1)
event_order[1,2:ncol(event_order)] <- eventmatrix[1,]

i <- 2
row <- 2
for(i in 2:nrow(eventmatrix)){
  
  start1 <- word(eventmatrix[i,1],1) #First event.
  end1 <- word(eventmatrix[i,1],2)
  CNV1 <- word(eventmatrix[i,1],3)
  Chr1 <- dm[match(start1,dm[,1]),4]
  
  j <- 2
  for(j in 2:nrow(eventmatrix)){
    start2 <- word(eventmatrix[j,1],1) #Second event.
    end2 <- word(eventmatrix[j,1],2)
    CNV2 <- word(eventmatrix[j,1],3)
    Chr2 <- dm[match(start2,dm[,1]),4]
    
    if(eventmatrix[i,1] != "0" && eventmatrix[j,1] != "0"){
      if(start1 == start2 && end1 == end2 && CNV1 != CNV2){ #Same breakpoints but different CNV.
        
        if(eventmatrix[i,1] %in% event_order[,2] == FALSE){
          event_order[row,2:ncol(event_order)] <- eventmatrix[i,]
          event_order[row,1] <- Chr1
          row <- row+1
        }
        
        if(eventmatrix[j,1] %in% event_order[,2] == FALSE){
          event_order[row,2:ncol(event_order)] <- eventmatrix[j,]
          event_order[row,1] <- Chr2
          row <- row+1
        }
      }
    }
    
    
    j <- j+1
  }
  
  i <- i+1
}

event_order <- event_order[event_order[,2]!="0",]

event_mode <- matrix(0,20,length(unique(word(event_order[2:nrow(event_order),2],1,2)))*2)
events <- unique(word(event_order[2:nrow(event_order),2],1,2))
i <- 1
for(i in 1:(ncol(event_mode)/2)){
  
  event_mode[1,i+(i-1)] <- events[i]
  
  i <- i+1
}

i <- 2
for(i in 2:nrow(event_order)){
  
  matchcolumn <- match(word(event_order[i,2],1,2),event_mode[1,])
  count <- as.numeric(event_mode[1,matchcolumn+1])+1
  event_mode[1,matchcolumn+1] <- count
  
  event_mode[count+1,matchcolumn] <- word(event_order[i,2],3)
  event_mode[count+1,matchcolumn+1] <- sum(as.numeric(event_order[i,3:ncol(event_order)])*as.numeric(nrcells))
  
  word(event_order[i,2],3)
  
  
  
  i <- i+1
}


#Nu ska jag hitta modaltalet för varje. Därefter kommer alla celler som har dessa breakpoints att få denna
#och därefter kommer alla ytterligare förändringar vara en ny händelse.
#Lägg in så att man kan göra en egen matris över vilken ordning som är rimlig.
#eventmatrix_before <- eventmatrix

i <- 1
for(i in 1:(ncol(event_mode)/2)){
  
  row_mode <- which.max(as.numeric(event_mode[2:nrow(event_mode),i+i]))+1
  
  event_stem <- paste(event_mode[1, i+(i-1)],event_mode[row_mode,i+(i-1)])
  event_stem_row <- match(event_stem,eventmatrix[,1])
  
  j <- 2
  for(j in 2:(as.numeric(event_mode[1,i+i])+1)){
    
    if(j != as.numeric(row_mode)){
      
      event_branch <- paste(event_mode[1, i+(i-1)],event_mode[j,i+(i-1)])
      
      event_branch_row <- match(event_branch,eventmatrix[,1])
      
      k <- 2
      for(k in 2:ncol(eventmatrix)){
        
        if(eventmatrix[event_branch_row,k]== "1"){
          
          eventmatrix[event_stem_row,k] <- "1" #It should have the stem event as well.
          
        }
        
        k <- k+1
      }
      
    }
    
    j <- j+1
  }
  
  i <- i+1
}


# i <- 2
# for(i in 2:nrow(event_order)){
#
#   start1 <- word(event_order[i,2],1)
#   end1 <- word(event_order[i,2],2)
#   CNV1 <- word(event_order[i,2],3)
#   rowmatch_start1 <- match(start1,dm[,1])
#   rowmatch_end1 <- match(end1,dm[,1])
#
#   if(rowmatch_start1 == 1){
#     rowmatchnew <- 2
#   }else{
#     rowmatchnew <- rowmatch_start1
#   }
#
#   if(dm[rowmatchnew,4]!=dm[rowmatchnew-1,4]){
#     print("Beginning")
#
#     j <- 2
#     for(j in 2:nrow(event_order)){
#
#       start2 <- word(event_order[j,2],1)
#       end2 <- word(event_order[j,2],2)
#       CNV2 <- word(event_order[j,2],3)
#       rowmatch_start2 <- match(start2,dm[,1])
#       rowmatch_end2 <- match(end,dm[,1])
#
#       if(dm[rowmatch_end2,4]!=dm[rowmatch_end2+1,4] && dm[rowmatch_end1+1,1]==dm[rowmatch_start1,1] && dm[rowmatch_start1,4]==dm[rowmatch_end2,4]){
#        print("Take up the whole chromosome.")
#
#       }
#
#      j <- j+1
#     }
#
#
#   }
#   i <- i+1
# }


############################################
#Lägger till en cell med alla stamhändelser#
############################################
new_events <- matrix(0,40,(as.numeric(ncol(eventmatrix))+1))
co <- 1
stemcell <- matrix(0,as.numeric(nrow(eventmatrix)),1)
stemcell[1,1] <- "Stemcell"

i <- 2
s <- 1
for(i in 2:as.numeric(nrow(eventmatrix))){ #Looping through events.
  print("i")
  print(i)
  x <- sum(as.numeric(eventmatrix[i,2:ncol(eventmatrix)])*as.numeric(nrcells[1,1:as.numeric(ncol(nrcells))]))
  print(x)
  print(sum(as.numeric(nrcells)))
  if(x/((sum(as.numeric(nrcells)))) >= co){ #The event is a stem event.
    stemcell[i,1] <- 1
    for(k in 2:ncol(eventmatrix)){
      print("k")
      print(k)
      if(eventmatrix[i,k] == "0"){ #Finding cells not having the stem event.
        pos1 <- as.numeric(match(word(eventmatrix[i,1],1),dm$Positions)) #Start position of stem event.
        pos2 <- as.numeric(match(word(eventmatrix[i,1],2),dm$Positions)) #End position of stem event.
        CNV_ny <- data[pos1,k-1] #The CNV of the new event.
        print(eventmatrix[i,1])
        print(pos1)
        print(pos2)
        print(CNV_ny)
        thenewevent <- paste(word(eventmatrix[i,1],1,2),CNV_ny) #The name of the new event.
        if(thenewevent %in% eventmatrix[,1] == TRUE){
          therow <- match(thenewevent,eventmatrix[,1])
          eventmatrix[therow,k] <- 1
          eventmatrix[i,k] <- 1
        }else{
          if(thenewevent %in% new_events[,1]){
            row <- match(thenewevent,new_events[,1])
            eventmatrix[i,k] <- 1 #Adding a 1 to the stem event.
            new_events[row,k] <- 1 #Adding the new event.
          }
          else{
            new_events[s,1] <- paste(word(eventmatrix[i,1],1,2),CNV_ny)
            eventmatrix[i,k] <- 1
            new_events[s,k] <- 1
            s <- s+1}
        }
      }
    }
  }
  else{stemcell[i,1] <- 0}
  i <- i+1
}


# i <- 1
# for(i in 2:nrow(eventmatrix)){ #Looping through the events.
# print(i)
#
#     cells_total <- sum(as.numeric(eventmatrix[i,2:(ncol(eventmatrix))])*as.numeric(nrcells[1,]))
#     print(cells_total)
#     print(sum(as.numeric(nrcells[1,])))
#     print(co)
#     if(cells_total/sum(as.numeric(nrcells[1,])) >= co){
#       print(eventmatrix[i,1])
#       stemcell[i,1] <- 1
#     }
#
#   i <- i+1
# }


new_events <- new_events[new_events[,1] != "0",]
#eventmatrix_withX <- cbind(eventmatrix,stemcell)
eventmatrix_final <- cbind(eventmatrix,stemcell)
eventmatrix_final <- rbind(eventmatrix_final,new_events)

#eventmatrix_final_saved <- eventmatrix_final
#eventmatrix_final <- eventmatrix_final_saved

i <- 2
for (i in 2:nrow(eventmatrix_final)) {
  pos1 <- match(word(eventmatrix_final[i,1],1),dm[,1])
  pos2 <- match(word(eventmatrix_final[i,1],2),dm[,1])
  CNA <- word(eventmatrix_final[i,1],3)
  #print(CNA)
  if(pos2-pos1 < 4){ #5-1 = 4. S? 5 bins blir 4.
    if(as.numeric(CNA) <= 5){
      print(CNA)
      eventmatrix_final[i,1] <- "0"
    }
  }
  i <- i+1
}

cellnames <- eventmatrix_final[1,2:as.numeric(ncol(eventmatrix_final))]
eventmatrix_final <- eventmatrix_final[eventmatrix_final[,1] != "0",]
EM <- eventmatrix_final[2:as.numeric(nrow(eventmatrix_final)),2:as.numeric(ncol(eventmatrix_final))] #Plockar ut h?ndelserna ur matrisen.


load_EM <- function(filename, sheetname) {  #Skapar en funktion som fixar till datan.
  data <- as.data.frame(read_xlsx(filename, sheetname)) #H?r l?ser jag in xlsx-filen och sparar den i variabeln data
  subdata <- data[ c(1:nrow(data)), c(2:ncol(data)) ] #H?r plockar jag ut den del av filen som jag faktiskt vill anv?nda i mina ber?kningar sen.
  return(subdata)
}

#################################
#Building the phylogenetic trees#
#################################
root <- "Normal"
normal <- "Normal"
cellnames <- as.vector(colnames(EM))
EM_phy <- phydatevent(EM,cellnames,normal)



#MP tree function
mp_tree <- function(EM_phy){
  MP_tree_pratchet <- pratchet(EM_phy, start = NULL, method = "fitch", maxit = 2000, k = 10, #Funktionen anv?nder sig av Fithchs algoritm. Pratchet = p ratchett (1999).
                               trace = 1, all = FALSE, rearrangements = "TBR",
                               perturbation = "ratchet")
  MP_tree_pratchet <- root(MP_tree_pratchet, outgroup = root, resolve.root = TRUE)
  treeRatchet <- acctran(MP_tree_pratchet, EM_phy)
  plot(treeRatchet)
  #bs <- bootstrap.pml(MP_tree_pratchet, bs=100, optNni=TRUE)
  #treeBS <- plotBS(Lf_JC$tree,bs, type = "phylogram")
  #BStrees <- bootstrap.phyDat(EM_phy, pratchet, bs = 100)
  #treeMP <- plotBS(treeRatchet, BStrees, "phylogram")
  return(treeRatchet)}
mp_tree(EM_phy)

#ML tree function
ml_tree <- function(Eventmatrix) {
  dm_h <- dist.hamming(Eventmatrix) #Tar fram en avst?ndsmatris som ska ligga till grund f?r tr?dgenereringen.
  starting_tree <- NJ(dm_h)
  starting_tree <- root(starting_tree, outgroup = root,resolve.root = TRUE)
  Lf <- pml(starting_tree, Eventmatrix) #Obtaining an object of class pml
  Lf_JC <- optim.pml(Lf, model = "JC")#Jukes Cantor model.
  #bs <- bootstrap.pml(Lf_JC, bs=100, optNni=TRUE)
  #treeBS <- plotBS(Lf_JC$tree,bs, type = "phylogram")
  #s <- 15
  #ggsave(treeBS,file="bootstrap_ml.png",width = s,height = s)
  return(Lf_JC)
}
ml_tree(EM_phy)


#Plotting the trees.
RMS_mptree <- mp_tree(EM_phy)
RMSmp <- ggplot(RMS_mptree) + geom_tree() + geom_tiplab(size = 4)
RMSmp <- RMSmp +  theme_tree()+xlim(c(0, 20))
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
print(RMSmp+mytheme)
s <- 10
w <- 10
ggsave(RMSmp,file="MP_tree.pdf",width = w,height = s)


RMS_mltree <- ml_tree(EM_phy)
RMS_mltree <- RMS_mltree$tree
RMSml <- ggplot(RMS_mltree) + geom_tree() + geom_tiplab(size = 4)
RMSml <- RMSml +  theme_tree()+xlim(c(0, 1))
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
print(RMSml+mytheme)
s <- 10
w <- 10
ggsave(RMSml,file="ML_tree.pdf",width = w,height = s)

#Saving the EM
write_xlsx(as.data.frame(eventmatrix_original), "EM_first.xlsx")
write_xlsx(as.data.frame(eventmatrix_final), "EM_final.xlsx")




#You could also do phylogenies where you color the cell names according to their group belonging.
blue <- c("#6bb5d8","#6fb9e7","#4d8dc6","#2c6f9a","#205575")
red <- c("#ed6d70","#ea5456","#e52421","#a71916")
orange <- c("#e55514","#ed8606","#f6c400")
yellow <- c("#f6c400","#ed8606","#e55514")
#grey <- c("#b9b8b8","#9d9d9c","#706f6f","#3c3c3b") In article.
green <- c("#add3a2","#6abfa4","#497f7a","#2c574a")
brown <- c("#ca9e67","#936037")
purple <- c("#d4bae0","#c9a0dc","#ae87d0","#7851a9","#522d80","#500691","#330066")
grey <- c("#b9b8b8","#9d9d9c","#8a8a8a","#706f6f","#595858","#3c3c3b","#212121")
pink <- c("#f9c0dc","#f29bc2", "#eb619f", "#dd2f89")

#Choose which colors you want. Change this according to which groups you have.
custom_col <- t(as.matrix(c(green[2],green[3],
                            blue[2],blue[3],
                            yellow[1],yellow[2],yellow[3], red[3], red[4], red[5],  #PDX3
                            purple[2],purple[4])))


branchnames <- as.character(nm)
samples <- list(1:54,55:103,104:128,129:150,151:169,170:196) #lp, hp, C1, C4, T4, T5.
branches <- list(lp = branchnames[1:54],hp =branchnames[55:103], #lp and hp.
                 C1 = branchnames[104:128],C4=branchnames[129:150], #Control 1 and 4.
                 T4=branchnames[151:169],T5=branchnames[170:196],
                 Stemcell = c("Stemcell"),Normal=c("Normal"))


RMS_mptree <- mp_tree(EM_phy)
RMS_mptree <- groupOTU(RMS_mptree,branches)
RMSmp <- ggplot(RMS_mptree) + geom_tree() + geom_tiplab(size = 4, aes(color = group))+
  scale_color_manual(values=c(C1=orange[1],C2=orange[2],C3=orange[3],
                              T1_pc= purple[2],T4_pc= purple[3],T5_pc=purple[4],T9_pc=purple[5],
                              T6_reg = blue[1],T7_reg= blue[3],T11_reg= blue[4],
                              T_rel= pink[1],T2_rel= pink[2],T5_rel= pink[3],T17_rel= pink[4],
                              Stemcell = "black",Normal="black"))


RMSmp <- RMSmp +  theme_tree()
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"),legend.pos = "none")
print(RMSmp+mytheme)
s <- 14
ggsave(RMSmp,file="MP_color.pdf",width = 20,height = 24)


RMS_mltree <- ml_tree(EM_phy)
#plot(RMS_mltree)
RMS_mltree <- RMS_mltree$tree
branches <- list(lp = branchnames[1:54],hp =branchnames[55:103], #lp and hp.
                 C1 = branchnames[104:128],C4=branchnames[129:150], #Control 1 and 4.
                 T4=branchnames[151:169],T5=branchnames[170:196],
                 Stemcell = c("Normal"))
branches <- list(C1 = branchnames[1:25],C4=branchnames[26:47],
                 T4=branchnames[48:66],T5=branchnames[67:93],
                 Stemcell = c("Stemcell"))
RMS_mltree <- groupOTU(RMS_mltree,branches)
RMSml <- ggplot(RMS_mltree) + geom_tree() + geom_tiplab(size = 4, aes(color = group))+
  scale_color_manual(values=c(C1=orange[1],C2=orange[2],C3=orange[3],
                              T1_pc= purple[2],T4_pc= purple[3],T5_pc=purple[4],T9_pc=purple[5],
                              T6_reg = blue[1],T7_reg= blue[3],T11_reg= blue[4],
                              T_rel= pink[1],T2_rel= pink[2],T5_rel= pink[3],T17_rel= pink[4],
                              Stemcell = "black",Normal="black"))
RMSml <- RMSml +  theme_tree()
mytheme <- theme(plot.title = element_text(hjust = 0.5, size = (14), color = "black"))
print(RMSml+mytheme)
s <- 14
ggsave(RMSml,file="PDX3_211011_noP_col_ml.pdf",width = 20,height = 24)