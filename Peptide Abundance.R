####This script is to filter peptides following the steps in the paper and conduct differential expression analysis using Limma####
### Read in the data files for HLA-A29 affinity purification 
setwd("/Users/jonaskuiper/Documents/r/march 2020/HLA peptidomics/")

### Setup dependencies ----
library(readr)
library(gplots)
library(limma)
library(qvalue)

# data preprocessing, load all functions following PMID: 25821719
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")
### END: Setup dependencies ----

#create a new function which is the opposite of %in% (will not select terms you specify later in "Quan.info")
`%nin%` = Negate(`%in%`)

# Read in the mass spectrometry datasets for the two biological replicates of the HLA-A29-binding antibody
library(readr)
Replicate_1 <-read_delim("P1146 Elutie set 1 HLA-A29 2400 2402 abc_SILAC Phe10 en Tyr10 Percolator_PeptideGroups.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

Replicate_2<-read_delim("P1146 Elutie set 2 HLA-A29 2421 2423 abc_SILAC Phe10 en Tyr10 Percolator_PeptideGroups.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

#Next, We filter for peptides with Perculator q<0.01 and remove peptides with "quan.info" flags
Replicate_1_filtered<-Replicate_1[Replicate_1$`Percolator q-Value by Search Engine Mascot`<0.01 & Replicate_1$`Quan Info` %nin% c("InconsistentlyLabeled","NoQuanValues","Redundant","IndistinguishableChannels"),]
Replicate_2_filtered<-Replicate_2[Replicate_2$`Percolator q-Value by Search Engine Mascot`<0.01 & Replicate_2$`Quan Info` %nin% c("InconsistentlyLabeled","NoQuanValues","Redundant","IndistinguishableChannels"),]


 
#the number of unique peptides for each replicate
Replicate_1_filtered_unique<-Replicate_1_filtered[!duplicated(Replicate_1_filtered$Sequence),]
Replicate_2_filtered_unique<-Replicate_2_filtered[!duplicated(Replicate_2_filtered$Sequence),]
#venn diagram of the number of unique overlapping peptides
Venn.overlap<-list(Replicate_1= Replicate_1_filtered_unique$Sequence,Replicate_2=Replicate_2_filtered_unique$Sequence)
overlap_replicates<-venn(Venn.overlap) #2315 unique peptides overlap in both experiments

#extract the list of overlapping peptides from both replicates  
overlap<-as.character(attr(x = overlap_replicates, "intersections")$`Replicate_1:Replicate_2`)

A29_1_Rep<-Replicate_1_filtered[Replicate_1_filtered$Sequence %in% overlap,]
A29_2_Rep<-Replicate_2_filtered[Replicate_2_filtered$Sequence %in% overlap,]

#Prepare and normalize Biological replicate 1 as outlined in the r guide (next steps are adapted from PMID: 25821719):
#http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
#We extract the peptide sequence and abundance signal in light and heavy channels
pre_A29_elution_1<-A29_1_Rep[,c("Sequence","Abundances Grouped Light","Abundances Grouped Heavy")]
#we will only compute stats for peptides detected in both heavy and light so remove missing values
A29_elution_1<-pre_A29_elution_1[complete.cases(pre_A29_elution_1),]#1459 peptides without NA
#because we have prefiltered the data and want to keep al indiviual peptides we will use dummy variables for:
A29_elution_1$Quan.Usage<-"Used"
A29_elution_1$Protein.Group.Accessions<-A29_elution_1$Sequence
A29_elution_1$Isolation.Interference<-30
A29_elution_1$Quan.info<-"Unique"
colnames(A29_elution_1)<-c("Sequence","Abundances Grouped Light.x","Abundances Grouped Heavy.x","Quan.Usage","Protein.Group.Accessions","Isolation.Interference","Quan.Info")
cha_1<-c("Abundances Grouped Light.x","Abundances Grouped Heavy.x")

# use the two functions to read and normalize the peptide abundances. 
dat1 <- read.peptides(A29_elution_1, cha_1) 
dat1 <- quantify.proteins(dat1, cha_1)

#prepare the Biological replicate 2
pre_A29_elution_2<-A29_2_Rep[,c("Sequence","Abundances Grouped Light","Abundances Grouped Heavy")]
A29_elution_2<-pre_A29_elution_2[complete.cases(pre_A29_elution_2),] 
A29_elution_2$Quan.Usage<-"Used"
A29_elution_2$Protein.Group.Accessions<-A29_elution_2$Sequence
A29_elution_2$Isolation.Interference<-30
A29_elution_2$Quan.info<-"Unique"

colnames(A29_elution_2)<-c("Sequence","Abundances Grouped Light.y","Abundances Grouped Heavy.y","Quan.Usage","Protein.Group.Accessions","Isolation.Interference","Quan.Info")
cha_2<-c("Abundances Grouped Light.y","Abundances Grouped Heavy.y")

dat2 <- read.peptides(A29_elution_2, cha_2) 
dat2 <- quantify.proteins(dat2, cha_2)

#now we merge both noralized datasets for limma moderated t-statistics
dat3<-merge(dat1, dat2, by="row.names") # 1896 peptides considered for differential expression analysis
row.names(dat3)<-dat3$Row.names

# limma, factorial design; the first 8 columns of ech data set contain the channel measurements
tr <- as.factor(c(2,1,2,1)) #2 is Light(KO) and 1 is Heavy (WT-ERAP2)
ex <- as.factor(c(1,1,2,2)) #1 is experiment 1 and 2 is experiment 2
design <- model.matrix(~ ex+tr) #test groups 1 versus 2 controlling for experiment. 
res.eb.mult_all <- eb.fit.mult(dat3[,c(2,3,6,7)], design) #select ony the columns 2,3,6, and 7
res.eb.mult_all$logP<- -log10(res.eb.mult_all$p.mod) #change the moderated p values to -log10(P) for plotting purpose 
res.eb.mult_all$pep<-row.names(res.eb.mult_all)
res.eb.mult_all$length<-nchar(row.names(res.eb.mult_all))

#load the bindings scores for 1768 8-11 mers (total from limma was 1896 peptides) predicted by HLAthena.tools
HLA_thena_predictions<- read_delim("predictions_20200519.113620.tsv","\t", escape_double = FALSE, trim_ws = TRUE) 

# a total of 1330 peptides are considered to bind to HLA-A29. 
HLA_A29_res.eb.mult <- res.eb.mult_all[res.eb.mult_all$pep %in% HLA_thena_predictions$pep[HLA_thena_predictions$MSi_A2902 > 0.6],]

#write the data to an excelfile
write_table<-HLA_A29_res.eb.mult[,c("pep","length","logFC","p.mod","q.mod")]
write_xlsx(write_table,"Table 1330 HLA-A29 Limma.xlsx")

#to proceed to NMDS extract the 948 9-mers from the dataset and save for use in NMDS.R
HLA_A29_9_mers<-HLA_A29_res.eb.mult$pep[HLA_A29_res.eb.mult$length==9]
save(HLA_A29_9_mers, file="HLA_A29_948_9mers.RData")
