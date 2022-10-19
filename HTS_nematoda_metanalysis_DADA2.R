# Lara Macheriotou #
# Metanalysis of HTS data in DADA2 #

#Load libraries####
library(easypackages)
libraries("ggplot2","RColorBrewer","reshape2","vegan","Hmisc","tidyr","ggpubr","tibble",
          "plyr","knitr","readxl",
          "magrittr","dplyr","stats","ade4","GUniFrac","ape","UpSetR","gridExtra","utils",
          "methods","grDevices","scales","BiodiversityR","Rarity","phangorn","RVAideMemoire",
          "ggdendro","phytools","Biostrings","ggrepel","nlme","methods","mgcv","gridExtra","scales",
          "UpSetR","picante","readxl","car","broom","dada2","OSMscale","xlsx","cowplot",
          "googleway", "ggplot2", "ggrepel", "ggspatial","sf", "rnaturalearth", "rnaturalearthdata",
          "rgeos","maps","mapdata","plotly","ShortRead","phylotools","seqinr","DECIPHER",
          "coil", "iNEXT", "phyloseq", "devtools","ggExtra", "writexl", "FSA", "readr",
          "metagMisc", "ampvis2", "doParallel", "mgcv", "pairwiseAdonis", "lme4", "mvabund", "lsmeans",
          "GUniFrac", "effects", "jtools")
#----------------------------------------------------------------------------------------------------#
#Installations####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("FastqCleaner")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.15")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")

#Memory####
memory.limit(size=1.759219e+20)

#TRUE the maximum amount of memory obtained from the OS is reported
#FALSE the amount currently in use
#NA the memory limit.
memory.size(max=TRUE)
memory.size(max=FALSE)
memory.size(max=NA)

#Linux
install.packages("unix")
library(unix)
rlimit_as(1e19)  #increases to ~19GB

#----------------------------------------------------------------------------------------------------#
#Colour palettes####
#Create a function to interpolate between some colours
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
#Run function asking for n colours
getPalette(20)

#Function to calculate the mean and the standard deviation for each group
#data:data frame
#varname:the name of a column containing the variable to be summariezed
#groupnames:vector of column names to be used as grouping variables
data_summary <- function(data, varname, groupnames){require(plyr)
  summary_func <- function(x, col){c(mean = mean(x[[col]], na.rm=TRUE),
                                     sd = sd(x[[col]], na.rm=TRUE))}
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
#----------------------------------------------------------------------------------------------------#
#Set working directory####
setwd("F:/HTS/Nematoda_metanalysis")
#----------------------------------------------------------------------------------------------------#
#Latitude + longitude####
#Convert all longitude+latitude values to one coordinate system degrees-minutes-seconds (dms)
#Set working directory
setwd("~/1.Post-doc_UGent/2.Nematoda_metanalysis")
#Load dataframe
library(readxl)
lat.long1<-as.data.frame(read_excel("HTS_metanalysis_long_lat.xlsx",sheet = "dms"))

#Convert from dms to decimal
decimal<-degree(lat,long,data = lat.long1,todms = FALSE)

#Combine dataframes
res.1<-cbind(lat.long1,decimal)

#Write to excel sheet
write.xlsx(res.1,file="HTS_metanalysis_long_lat.xlsx",sheetName = "decimal.conv",
           col.names = TRUE, row.names = TRUE,append = TRUE)

#NOTE: Delete first numbered column in excel file after appending!#
#----------------------------------------------------------------------------------------------------#
#Fasta headers####

#Load dataframe with first column for original name, second column for the new name of the sequence
#Easier with excel files as input due to presence of many colons in text file
rename_fasta<-read_excel("rename_fasta.xlsx",sheet = "asv.all",col_names = TRUE)
rename.fasta(infile = "Nematoda_metanalysis_Nematoda_11034.fna",ref_table =rename_fasta,
             outfile = "Nematoda_metanalysis_Nematoda_sq11034.fna")
#----------------------------------------------------------------------------------------------------#
#World map####
#Generate map with sampling stations
#Tutorial: https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html#
#Biogeographic Marine Realms based on species endemicity map ####
#IMPORTANT:Move all files to directory conatining .shp file

#Read shape file
bmre<-read_sf("~/PostDoc_UGent/Nematoda_metanalysis/BMRE_vector_data/MarineRealms.shp")

#Make factors
bmre$Realm<-as.factor(bmre$Realm)

#Load stations coordinates + metadata
stations<-as.data.frame(read_excel("~/PostDoc_UGent/Nematoda_metanalysis/HTS_metanalysis_long_lat.xlsx",
                                   sheet = "decimal.all.final"))
stations<-as.data.frame(read_excel("~/PostDoc_UGent/Nematoda_metanalysis/HTS_metanalysis_long_lat.xlsx",
                                   sheet = "long.lat.33"))

#Datasets as factors
stations$Reference<-as.factor(stations$Reference)
stations$Sample.number<-as.factor(stations$Sample.number)
stations$Realm<-as.factor(stations$Realm)
stations$Type<-as.factor(stations$Type)

#Load world shape file
world.shapefile<-rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))

#Generate base map
world.map<-ggplot() +
  geom_sf(data = world.shapefile, size = .2, fill = "gray80", col = "black") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))
world.map

#Generate 14 colours
mycolors = c(brewer.pal(name="Dark2", n = 4), brewer.pal(name="Paired", n = 12))

#Add bmre data
world.map.f<-world.map+
  geom_sf(data = bmre, aes(), size = .2, alpha=.3)+
  #ggtitle("Sampling stations within Biogeographic Marine Realms")+
  geom_sf_text(data = bmre, aes(label = Realm), colour = "black")+
  coord_sf(expand = FALSE)+
  geom_point(data = stations, aes(x = mean.long, y = mean.lat, 
                                  fill=Reference, shape=Type, colour=Reference), size = 4)+
  #scale_color_manual(values = mycolors)+ scale_shape_manual(values=c(0:14))+
  theme(panel.grid.major = element_line(color = "white", linetype = "dashed", size = 0.5))+
  theme(axis.title = element_blank())+ guides(fill="none")+
  theme(legend.text = element_text(size=8), legend.position = "right", 
        legend.title = element_blank())+
  scale_shape_manual(values = c(21:25))
world.map.f

#Save as pdf
ggsave2("Metanalysis_world_map_bmre.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#Interactive map
world.map.inter<-world.map+
  geom_sf(data = bmre, aes(fill = Realm), size = .2, alpha=.3)+
  ggtitle("Sampling stations within Biogeographic Marine Realms")+
  geom_sf_text(data = bmre, aes(label = Realm), colour = "black")+
  coord_sf(expand = FALSE)+
  geom_point(data = stations, aes(x = long, y = lat, colour=Reference, shape=Reference), size = 2.5)+
  scale_color_manual(values = mycolors)+ scale_shape_manual(values=c(0:14))+
  theme(panel.grid.major = element_line(color = "white", linetype = "dashed", size = 0.5))+
  theme(axis.title = element_blank())+ guides(fill=FALSE)+theme(legend.position = "none")
world.map.inter

#Make interactive map
world.map.inter.f<-ggplotly(world.map.inter)
world.map.inter.f<-ggplotly(world.map.f)
#Save as html: call->export
world.map.inter.f
#----------------------------------------------------------------------------------------------------#
#Rename files####
#Set directory#
setwd("E:/HTS/12.PRJEB23641")

#View files
list.files()

#Load file containing old + new names
#FORMAT: Column Header A1:fastq.file; B1:new name, .txt file
rename.list<-read.table(file="Bik_et_al_2012_MolEcol_rename.txt",header = TRUE)

#Transform columns to vectors
rename.list$fastq.file<-as.vector(rename.list$fastq.file)
rename.list$new.name<-as.vector(rename.list$new.name)

#Rename files
file.copy(from=rename.list$fastq.file,to=rename.list$new.name,overwrite=FALSE)
#----------------------------------------------------------------------------------------------------#
#Merge fastq files####
#Set directory#
setwd("E:/HTS/6.Cowart_et_al_2020_FroMarSci")
combine.list<-read.table(file="Macheriotou_et_al_combine.txt",header = FALSE)

#Make as vector
combine.list$V1<-as.vector(combine.list$V1)

#Read fastq files into 1 R object
fastqs<-ShortRead::readFastq(combine.list$V1)

#Write single fastq output
ShortRead::writeFastq(fastqs,file="Realm13.Sample54.AbyssalPlain.PRJNA623689_R2.fastq.gz",
                      compress = TRUE)
#----------------------------------------------------------------------------------------------------#
#QC Check####
setwd("F:/HTS/Nematoda_metanalysis/4.Haenel_et_al_2017_BioDivDataJour")
path<-"F:/HTS/Nematoda_metanalysis/4.Haenel_et_al_2017_BioDivDataJour"

#View files
list.files(path)

#Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs<-sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs<-sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

#For Roche 454 data (single read)
fnFs<-sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names<-sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visualize the quality profiles of the forward (R1) reads | not aggregated
R1.quality<-plotQualityProfile(fnFs[1:14],aggregate=FALSE)+
  geom_hline(yintercept= 30,colour="blue")+labs(title="R1")+
  geom_vline(xintercept = 200, colour="red")+
  scale_x_continuous(breaks=c(seq(from=0,to=500,by=50)),limits = c(0,300))
R1.quality

#Visualize the quality profiles of the forward (R1) reads | not aggregated Roche 454
R1.quality<-plotQualityProfile(fnFs[1:25],aggregate=FALSE)+
  geom_hline(yintercept= 30,colour="blue")+labs(title="R1")+
  scale_x_continuous(breaks=c(seq(from=0,to=500,by=50)),limits = c(0,500))
R1.quality

#Visualize the quality profiles of the reverse (R2) reads| not aggregated
R2.quality<-plotQualityProfile(fnRs[1:14],aggregate=FALSE)+
  geom_hline(yintercept= 30,colour="blue")+labs(title="R2")+
  scale_x_continuous(breaks=c(seq(from=0,to=500,by=50)),limits = c(0,350))
R2.quality

#Save plots
ggsave(R1.quality,filename ="R1_quality.pdf",units = "in",width = 14,height = 8.5,scale = 1.5)
ggsave(R2.quality,filename ="R2_quality.pdf",units = "in",width = 14,height = 8.5,scale = 1.5)
#----------------------------------------------------------------------------------------------------#
#Mean depth####
#Calculate mean depth by sample
depth<-read.table("Nematoda_metanalysis_samples_depth.txt", header = TRUE)
depth$Sample<-factor(depth$Sample)
mean.depth<-tapply(X=depth$Depth, INDEX = depth$Sample, FUN= mean)
data.frame(mean.depth)

#Mean lat.long
lat.long<-read.table("Nematoda_metanalysis_long_lat.txt", header = TRUE)
lat.long$Realm.sample<-as.factor(lat.long$Realm.sample)
mean.lat<-tapply(X=lat.long$Lat.decimal, INDEX = lat.long$Realm.sample, FUN= mean)
mean.long<-tapply(X=lat.long$Long.decimal, INDEX = lat.long$Realm.sample, FUN= mean)
means.lat.long<-cbind(data.frame(mean.lat, mean.long))
#----------------------------------------------------------------------------------------------------#
#DADA2####
#Change path to the directory containing the fastq files after unzipping
setwd("F:/HTS/Nematoda_metanalysis/3.Macheriotou_et_al/dada/Pam-Moz01")
path<-"F:/HTS/Nematoda_metanalysis/3.Macheriotou_et_al/dada/Pam-Moz01"

#View files
list.files(path)

#Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs<-sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs<-sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

fnFs<-sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs<-sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names<-sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names<-sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#Visualize the quality profiles of the forward reads
R1.quality<-plotQualityProfile(fnFs[1:21],aggregate=FALSE)+geom_vline(xintercept=250, colour="blue")+
  geom_hline(yintercept= 30,colour="blue")+labs(title="R1 after Cutadapt")+
  scale_x_continuous(breaks=c(seq(from=0,to=300,by=50)),limits = c(0,300))
R1.quality

#Visualize the quality profiles of the reverse reads
R2.quality<-plotQualityProfile(fnRs[1:14],aggregate=FALSE)+geom_vline(xintercept=250, colour="blue")+
  geom_hline(yintercept= 30,colour="blue")+labs(title="R2 after Cutadapt")+
  scale_x_continuous(breaks=c(seq(from=0,to=300,by=50)),limits = c(0,300))
R2.quality

#Save plots
ggsave(R1.quality,filename ="Klunder_et_al_2020_FroMarSci.dada2.quality_R1.fastq.pdf",
       units = "in",width = 14,height = 8.5,scale = 1.5)
ggsave(R2.quality,filename ="Klunder_et_al_2020_FroMarSci.dada2.quality_R2.fastq.pdf",
       units = "in",width = 14,height = 8.5,scale = 1.5)

#Place filtered files in filtered/subdirectory
filtFs<-file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs<-file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs)<-sample.names
names(filtRs)<-sample.names

#Quality filter (~fastq_filter in USearch9)
#Truncate forward @250bp, reverse @200bp
filter<-filterAndTrim(fnFs,filtFs,fnRs,filtRs,maxN=0,maxEE=c(2,5),rm.phix=TRUE,truncQ=2,
                      compress=TRUE,matchIDs=FALSE,multithread=FALSE,truncLen=c(250,200))

#Truncate forward @200bp, reverse @200bp
filter<-filterAndTrim(fnFs,filtFs,fnRs,filtRs,maxN=0,maxEE=c(2,5),rm.phix=TRUE,truncQ=2,
                      compress=TRUE,matchIDs=FALSE,multithread=FALSE,truncLen=c(200,200))

#Roche 454: single-read
filter<-filterAndTrim(fnRs,filtRs,maxN=0,maxEE=5,rm.phix=TRUE,truncQ=2,
                      compress=TRUE,matchIDs=FALSE,multithread=FALSE,truncLen=250)

filter<-filterAndTrim(fnFs,filtFs,maxN=0,maxEE=5,rm.phix=TRUE,truncQ=2,
                      compress=TRUE,matchIDs=FALSE,multithread=FALSE,truncLen=250)

#View filter output
filter

#Learn error rates from forward and reverse fastqs by alternating between sample
#inference and error rate estimation until convergence.
errF<-learnErrors(filtFs, multithread=TRUE)
errR<-learnErrors(filtRs, multithread=TRUE)

#Visualize the estimated error rates in forward reads
plotErrors(errF, nominalQ=TRUE)

#Visualize the estimated error rates in reverse reads
plotErrors(errR, nominalQ=TRUE)

#The error rates for each possible transition (A->C, A->G, ...) are shown.
#Points are the observed error rates for each consensus quality score.
#The black line shows the estimated error rates after convergence of the machine-learning algorithm.
#The red line shows the error rates expected under the nominal definition of the Q-score.

#Dereplicate forward and reverse fastq files
derepFs<-derepFastq(filtFs, verbose=TRUE)
derepRs<-derepFastq(filtRs, verbose=TRUE)

#Name the derep-class objects by the sample names
names(derepFs)<-sample.names
names(derepRs)<-sample.names

#Run sample inference on dereplicated forward and reverse reads using calcualted error rates
dadaFs<-dada(derepFs, err=errF, multithread=TRUE)
dadaRs<-dada(derepRs, err=errR, multithread=TRUE)

#Roche 454
dadaFs<-dada(derepFs, err=errF, multithread=TRUE,HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
dadaFs<-dada(derepRs, err=errR, multithread=TRUE,HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

#Inspect dada-class
dadaFs[[1]]
dadaRs[[1]]

#Merge forward and reverse denoised reads
#Merging is performed by aligning the denoised forward reads with the reverse-complement
#of the corresponding denoised reverse reads, and then constructing the merged "contig" sequences
mergers<-mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Inspect the merger data.frame
mergers

#Construct an amplicon sequence variant table (ASV) table
seqtab<-makeSequenceTable(mergers)

#Roche 454
seqtab<-makeSequenceTable(dadaFs)
seqtab<-makeSequenceTable(dadaRs)

#Check dimensions (samples x ASVs)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtab.nochim.all)))
table(nchar(getSequences(brandt2020)))

#Save as R object 
saveRDS(seqtab, "Gabriella_Panto_MSc_seqtab.rds") 

#Merge runs
bik2012<-readRDS("F:/HTS/Nematoda_metanalysis/1.Bik_et_al_2012_MolEcol/dada2/Bik_et_al_2012_MolEcol_seqtab.rds")
brandt2020<-readRDS("F:/HTS/Nematoda_metanalysis/2.Brandt_et_al_2020_MolEcolRes/dada2/Brandt_et_al_2020_MolEcolRes_seqtab.rds")
fonseca2017<-readRDS("F:/HTS/Nematoda_metanalysis/5.Fonseca_et_al_2017_SciRep/dada2/Fonseca_et_al_2017_SciRep_seqtab.rds")
haenel2017<-readRDS("F:/HTS/Nematoda_metanalysis/4.Haenel_et_al_2017_BioDivDataJour/dada2/Haenel_et_al_2017_seqtab.rds")
cowart2020<-readRDS("F:/HTS/Nematoda_metanalysis/6.Cowart_et_al_2020_FroMarSci/dada2/Cowart_et_al_2020_seqtab.rds")
cowart2015<-readRDS("F:/HTS/Nematoda_metanalysis/7.Cowart_et_al_2015_PlosOne/dada2/Cowart_et_al_2015_PlosOne_seqtab.rds")
faria2018<-readRDS("F:/HTS/Nematoda_metanalysis/8.Faria_et_al_2018_MarEnvRes/dada2/Faria_et_al_2018_MarEnvRes_seqtab.rds")
sinniger2016<-readRDS("F:/Nematoda_metanalysis/HTS/9.Sinniger_et_al_2016_FroMarSci/dada2/Sinniger_et_al_2016_FroMarSci_seqtab.rds")
fais2020<-readRDS("F:/HTS/Nematoda_metanalysis/10.Fais_et_al_2020_EstCoaSheSci/dada2/Fais_et_al_2020_EstCoaSheSci_seqtab.rds")
klunder2020<-readRDS("F:/HTS/Nematoda_metanalysis/11.Klunder_et_al_2020_FroMarSci/dada2/Klunder_et_al_2020_FroMarSci_seqtab.rds")
PRJEB23641<-readRDS("F:/HTS/Nematoda_metanalysis/12.PRJEB23641/dada2/PRJEB23641_seqtab.rds")
panto2020<-readRDS("F:/HTS/Nematoda_metanalysis/13.Gabriella_Panto_MSc/Metanalysis/Gabriella_Panto_MSc_seqtab.rds")
macheriotou<-readRDS("F:/HTS/Nematoda_metanalysis/14.Macheriotou_et_al/dada2/Macheriotou_et_al_seqtab.rds")

#Rename matrix rows for datasets with ONLY one sample
rownames(cowart2020)<-("Realm18.Sample23.HydrothermalVent.PRJNA540908")
rownames(cowart2015)<-("Realm18.Sample55.Shelf.NA")
rownames(fais2020)<-("Realm18.Sample57.Intertidal.SRR11266719")
rownames(faria2018)<-("Realm21.Sample24.Intertidal.NA")
rownames(klunder2020)<-("Realm18.Sample64.HydrothermalVent.PRJEB36829")
rownames(panto2020)<-("Realm30.Sample46.Antarctic.NA")
rownames(fonseca2017)<-("Realm30.Sample56.Antarctic")

seqtab.all<-mergeSequenceTables(bik2012, brandt2020, cowart2015, cowart2020, fais2020, 
                                faria2018,fonseca2017,haenel2017, klunder2020, panto2020,
                                PRJEB23641, sinniger2016, macheriotou, orderBy = "nsamples")

#Remove chimeras
seqtab.nochim.all<-removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)

#Check dimensions (samples x non-chimeric ASVs)
dim(seqtab.nochim.all)

#Percent non-chimeric ASVs by abundance
sum(seqtab.nochim.all)/sum(seqtab.all)

#Track reads#
#This step needs to executed within each RData file per dataset in order to combine for all samples
getN<-function(x) sum(getUniques(x))
track<-cbind(filter, sapply(dadaFs, getN), sapply(dadaRs, getN),
             sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track)<-c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track)<-sample.names
track

#Track reads through pipeline for Roche454
getN<-function(x) sum(getUniques(x))
track<-cbind(filter, sapply(dadaRs, getN))
colnames(track)<-c("input", "filtered", "denoisedR")
rownames(track)<-sample.names

#Inspect track
track

#ASV lentgths####
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtab.all)))
table(nchar(getSequences(seqtab.nochim.all)))

#Generate table to feed to ggplot2 for histogram
reads.per.seqlen<-tapply(colSums(seqtab.nochim.all), 
                         factor(nchar(getSequences(seqtab.nochim.all))), sum)

reads.per.seqlen<-tapply(colSums(seqtab), 
                         factor(nchar(getSequences(seqtab))), sum)

reads.per.seqlen<-tapply(colSums(seqtab.all), 
                         factor(nchar(getSequences(seqtab.all))), sum)

#Histogram 
df<-data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)
asv.hist<-ggplot(data=df, aes(x=length, y=count)) + geom_col()+
  scale_x_continuous(expand = c(0,0), breaks = seq(from=100, to=450, by=50), limits=c(100,450))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+theme(plot.margin=unit(c(t=0.5,r=0.5,b=0.5,l=0.5),"cm"))
asv.hist

ggsave(asv.hist,filename ="Nematoda_metanalysis_asv_length_histogram_0-5cm.pdf",
       units = "in",width = 14,height = 8.5,scale = 0.5)

#RDP does not run if ASVs are >50 nucelotides
#Filter seqtab ASVs by length
seqtab.nochim.all.filt<-seqtab.nochim.all[,nchar(colnames(seqtab.nochim.all)) %in% 350:425]
seqtab.nochim.all.filt250<-seqtab.nochim.all[,nchar(colnames(seqtab.nochim.all)) %in% 250:425]

dim(seqtab.nochim.all.filt350)
#----------------------------------------------------------------------------------------------------#
#RDP####
#RDP and Metazoan refDB
#Does not require sequence alignment
#!Check number of taxonomic levels

#18S_Silva_Nematoda_18992.fasta
taxa.all<-assignTaxonomy(seqtab.nochim.all.filt250,
                         "F:/refDBs/18S/converted_18S_Silva_Nematoda_18992.fasta",
                         tryRC = TRUE,multithread=TRUE,minBoot = 80,
                         taxLevels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

#Re-assign taxonomy for Nematoda ASVs
#Transpose
seqtab.nochim.nematoda<-t(seqtab.nochim.nematoda)

#18S_Nematoda_marine_971.fasta
taxa.nema<-assignTaxonomy(seqtab.nochim.nematoda,
                     "F:/refDBs/18S/18S_Nematoda_marine_971.fasta",
                     tryRC = TRUE,multithread=TRUE,minBoot = 80,
                     taxLevels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

#Tax to seqtab.no.chim####
#Check dimensions
dim(seqtab.nochim.all.filt350)
dim(taxa.all)

#To merge these two files, the seqtab.nochim file has to be transponed
seqtab.table<-as.data.frame(t(seqtab.nochim.all.filt350))
nema.tax.seqtab<-as.data.frame(t(seqtab.nochim.nematoda))

taxa.seqtab<-cbind(taxa.all,ASV=rownames(taxa.all))
rownames(taxa.seqtab)<-NULL

nema.tax<-cbind(taxa.nema,ASV=rownames(taxa.nema))
rownames(nema.tax)<-NULL

#Add extra columns containing taxonomy and sequences to the count table
taxa.seqtab<-as.data.frame(taxa.seqtab)
seqtab.table$Phylum<-taxa.seqtab$Phylum
seqtab.table$Class<-taxa.seqtab$Class
seqtab.table$Order<-taxa.seqtab$Order
seqtab.table$Family<-taxa.seqtab$Family
seqtab.table$Genus<-taxa.seqtab$Genus
seqtab.table$Species<-taxa.seqtab$Species
seqtab.table$ASV<-as.vector(colnames(seqtab.nochim.all.filt350))

nema.tax<-as.data.frame(nema.tax)
nema.tax.seqtab$Phylum<-nema.tax$Phylum
nema.tax.seqtab$Class<-nema.tax$Class
nema.tax.seqtab$Order<-nema.tax$Order
nema.tax.seqtab$Family<-nema.tax$Family
nema.tax.seqtab$Genus<-nema.tax$Genus
nema.tax.seqtab$Species<-nema.tax$Species
nema.tax.seqtab$ASV<-as.vector(colnames(seqtab.nochim.nematoda))

#Extract Nematoda only from seqtab.table
nematoda.taxa<-subset(seqtab.table, Phylum=="Nematoda")

#Make Nematoda only seqtab for second round of RDP
#samples as ROWS, ASV as column
seqtab.nochim.nematoda<-nematoda.taxa[,1:64]
#Transpose
seqtab.nochim.nematoda<-t(seqtab.nochim.nematoda)
dim(seqtab.nochim.nematoda)

#Rename seqtab rows####
#Load sample names
stations<-read_excel("Nematoda_metanalysis_samples.xlsx", sheet = "samples")

#replace rownames
rownames(seqtab.nochim.all)<-stations$Samples
rownames(seqtab.nochim.all.filt350)<-stations$Samples
rownames(seqtab.nochim.nematoda)<-stations$Samples

#Write to excel
write_xlsx(seqtab.table, "Nematoda_metanalysis_seqtab_taxa_Silva_Nematoda_18991_sq81216.xlsx")
write_xlsx(as.data.frame(taxa.nema), "Nematoda_metanalysis_RDP_Nematoda_971_sq81216.xlsx")
write_xlsx(nema.tax.seqtab, "Nematoda_metanalysis_nematoda_tax_seqtab_11034.xlsx")

#----------------------------------------------------------------------------------------------------#
#Phyloseq files####
#Transform taxonomic assignments to dataframe
taxa.table<-as.data.frame(taxa)
taxa.table<-taxa.table %>% mutate_if(is.character,as.factor)

taxa.table<-as.data.frame(taxa.250)
taxa.table<-taxa.table %>% mutate_if(is.character,as.factor)

#Filter sequences not assigned to Nematoda
taxa.nematoda<-dplyr::filter(taxa.table,Phylum=="Nematoda")

#Write taxonomic assignments to excel
write_xlsx(taxa.table,"F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_taxa_98854_refDB970.xlsx")
write_xlsx(taxa.table,"F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_taxa_122467_refDB970.xlsx")

#Transpose seqtab.nochim
seqtab.nochim.all.filt350.t<-t(seqtab.nochim.all.filt350)
seqtab.nochim.nematoda.t<-t(seqtab.nochim.nematoda)

#Transform seqtab.nochim.t to dataframe
seqtab.nochim.filt.df<-as.data.frame(seqtab.nochim.filt.t, check.names=FALSE, optional=TRUE)
seqtab.nochim.filt250.df<-as.data.frame(seqtab.nochim.filt250.t, check.names=FALSE, optional=TRUE)

#Export seqtab.nochim as tab-delimited file, get headers from seqtab.nochim.fasta| CHANGE FILE NAME!
write.table(seqtab.nochim.filt.t, "F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_seqtab_98854.txt"
            ,sep="\t", row.names=TRUE, col.names=NA,quote=FALSE)

write.table(seqtab.nochim.filt250.t, 
            "F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_seqtab_122467.txt"
            ,sep="\t", row.names=TRUE, col.names=NA,quote=FALSE)


#Ampvis2
asv.tax.table<-cbind(otu_table(phylo.rare),tax_table(phylo.rare))
write.xlsx2(asv.tax.table, "Nematoda_metanalysis_asv_tax_table_phylo_rare_8548.xlsx")

asv.tax.table<-cbind(otu_table(phylo.rare.pa),tax_table(phylo.rare.pa))
write.xlsx2(asv.tax.table, "Nematoda_metanalysis_asv_tax_table_phylo_rare_pa_8548.xlsx")

#----------------------------------------------------------------------------------------------------#
#Phyloseq####
#Combine to phyloseq object (https://vaulot.github.io/tutorials/Phyloseq_tutorial.html)

#Write phyloseq files
#https://rdrr.io/github/microbiome/microbiome/man/write_phyloseq.html

#Set working directory
setwd("F:/HTS/Nematoda_metanalysis")

#Read in ASV table, taxonomy table, sample table
#!!A1 should read "asv"
#!!Check sample names in seqtab
asv.mat<-read.table("Nematoda_metanalysis_nematoda_seqtab_11034.txt", header = TRUE)
tax.mat<-read_excel("Nematoda_metanalysis_RDP_Nematoda_971_sq11034.xlsx", sheet="nema")
samples.df<-read_excel("Nematoda_metanalysis_samples.xlsx", sheet="samples.63.metadata")

samples.df$Realm.type<-as.factor(samples.df$Realm.type)
samples.df$Type<-as.factor(samples.df$Type)
samples.df$Mean.depth<-as.numeric(samples.df$Mean.depth)

#Phyloseq object requires row.names
#Define row.names from the asv column
asv.mat<-column_to_rownames(asv.mat, var="asv")
tax.mat<-column_to_rownames(tax.mat, var="asv")
samples.df<-column_to_rownames(samples.df, var="Sample")

#Transform to matrix
asv.mat<-as.matrix(asv.mat)
tax.mat<-as.matrix(tax.mat)

#Import reference seqs as DNAString
#!!Fasta headers as "sq1"
asv.seqs<-readDNAStringSet("Nematoda_metanalysis_Nematoda_sq11034.fna")

asv<-otu_table(asv.mat, taxa_are_rows = TRUE)
tax<-tax_table(tax.mat)
samples<-sample_data(samples.df)
asv.seqs<-refseq(asv.seqs)

#Combine to phyloseq
phylo<-phyloseq(asv, tax, samples,asv.seqs)
phylo

#----------------------------------------------------------------------------------------------------#
#Prune samples####
#Zero reads
phylo.prune<-prune_samples(sample_sums(phylo)>1, phylo)
phylo.prune

#Reads per sample####
sums<-data.frame(sample_sums(phylo.prune))
sums<-rownames_to_column(sums)

#Reads per sample histogram
reads.hist<-ggplot(data=sums, aes(x=rowname, y=sample_sums.phylo.prune.)) + geom_col()+
  theme_bw()+theme(plot.margin=unit(c(t=0.5,r=0.5,b=0.5,l=0.5),"cm"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_blank())+ylab("Total reads")
reads.hist

ggsave2("Nematoda_metanalysis_reads_per_sample.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#Remove ASVs that only occur once in the dataset
#Only Keep taxa that were observed at least twice 
phylo.rm1<-filter_taxa(phylo.prune, function(x){sum(x > 0) > 1}, prune = TRUE)
phylo.rm1

#Prune ASVs with total relative abundance across samples of less than 0.000005%
minTotRelAbun = 5e-6
x<-taxa_sums(phylo)
keepTaxa<-(x / sum(x)) > minTotRelAbun
phylo.filt<-prune_taxa(keepTaxa, phylo)
phylo.filt

#----------------------------------------------------------------------------------------------------#
#Rarefaction####
#1.Rarefaction curve
#Rarefy at number of reads where curve begins to plateau
rarecurve(t(otu_table(phylo.prune)), step=50, cex=0.5)

ggsave2("Nematoda_metanalysis_rarefaction_rm1_sq17495.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#2.Rarefaction even depth
#Rarefy @9000 --> Lose 12 samples
#Set.seed(1)
phylo.rare<-rarefy_even_depth(phylo, rngseed=1,sample.size = 7985,
                                             replace=FALSE, trimOTUs = TRUE)

#Phylo.rare
#2486 ASvs were removed because they are no longer present in any sample
phylo.rare

#Prune samples again after removing taxa only observed once
phylo.rare<-prune_samples(sample_sums(phylo.rare)>1, phylo.rare)
phylo.rare

#----------------------------------------------------------------------------------------------------#



#----------------------------------------------------------------------------------------------------#
#Merge by factors####
#Realm number and sample stored as numeric -> convert to character
sample_data(phylo.rare)$Realm<-as.character(sample_data(phylo.rare)$Realm)
sample_data(phylo.rare)$Sample<-as.character(sample_data(phylo.rare)$Sample)
sample_data(phylo.rare)$Type<-as.character(sample_data(phylo.rare)$Type)
sample_data(phylo.rare)$Type<-as.factor(sample_data(phylo.rare)$Type)

#Merge samples
phylo.rare.type<-merge_samples(phylo.rare, "Type")
phylo.rare.gen.type<-merge_samples(phylo.rare.gen, "Type")
phylo.rare.zone<-merge_samples(phylo.rare, "Depth.zone")
phylo.rare.gen.zone<-merge_samples(phylo.rare.gen, "Depth.zone")

#Merged phyloseqs do not contain refseqs
#Append refseqs to merged phyloseq
#Import reference seqs as DNAString
asv.seqs<-readDNAStringSet("Nematoda_metanalysis_rare_sq8548.fna")
asv.seqs<-refseq(asv.seqs)

#Append seqs
phylo.rare.type<-merge_phyloseq(phylo.rare.type, asv.seqs)

#Agglomerate####
#At Genus-level
phylo.rare.type.gen<-tax_glom(phylo.rare.type, taxrank = "Genus")
phylo.rare.gen<-tax_glom(phylo.rare, taxrank = "Genus")
#----------------------------------------------------------------------------------------------------#
#Subset phyloseq####
#Exlcude NA for downstream analysis
phylo.filt.rare.no.na<-subset_taxa(phylo.filt.rare, Genus!="NA")
phylo.rare.no.na<-subset_taxa(phylo.rare, Genus!="NA")


#Step-wise becuase !=c() doesnt work?
phylo.rare.subset<-subset_samples(phylo.rare,Type!="Seamount")
phylo.rare.subset<-subset_samples(phylo.rare.subset,Type!="MudVolcano")
phylo.rare.subset<-subset_samples(phylo.rare.subset,Type!="HydrothermalVent")
phylo.rare.subset<-subset_samples(phylo.rare.subset,Type!="Canyon")

phylo.rare.plain<-subset_samples(phylo.rare,Type=="Plain")
#----------------------------------------------------------------------------------------------------#
#Transform counts####
#Relative abundance
phylo.rare.type.rel<-transform_sample_counts(phylo.rare.type, function(x) x /sum(x))

#Presence-absence
phylo.rare.pa<-phyloseq_standardize_otu_abundance(phylo.rare, method = "pa")

phylo.rare.gen.zone.pa<-phyloseq_standardize_otu_abundance(phylo.rare.gen.zone, method = "pa")

#----------------------------------------------------------------------------------------------------#
#Seqtab to fasta####
#CHANGE FILE NAME!
uniquesToFasta(seqtab.nochim.all.filt250, 
               fout="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_asv_122467.fna", ids=NULL)

uniquesToFasta(seqtab.nochim.nematoda, 
               fout="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_Nematoda_asv_sq21996.fna", ids=NULL)

uniquesToFasta(seqtab.nochim.nematoda, 
               fout="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_Nematoda_asv_sq21996.fna",
               ids=paste0("ASV", seq(length(getSequences(seqtab.nochim))))))

#Save fasta####
writeXStringSet(refseq(phylo.rare),"Nematoda_metanalysis_rare_sq8548.fna")


#Extract depth metadata
mean.depth<-data.frame(sample_data(nematoda.rare.no.na)$Mean.depth)

#OTU table to .txt file
otu.table<-as.data.frame(t(otu_table(nematoda.rm1.rare.type)))
otu.table<-rownames_to_column(otu.table)

write_xlsx(otu.table,"Nematoda_metanalysis_rm1_rare_otu_table.xlsx")

#Taxonomy table as dataframe
taxa<-as.data.frame(tax_table(phylo.rm1.rare))
taxa<-rownames_to_column(taxa)

otu.table<-t((otu_table(nematoda.rm1.rare.type)))
otu.table<-as.data.frame(otu.table)
otu.table<-rownames_to_column(otu.table)

write.xlsx2(otu.table,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_evol_distinct.xlsx", sheetName = "otu.table",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(taxa,
            file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_taxa_rm1_rare_sq13307.xlsx",
            sheetName = "rm1.rare",col.names = TRUE, row.names = TRUE, append = TRUE)
#----------------------------------------------------------------------------------------------------#
#Taxonomy plot####
#Taxonomic plot by TYPE absolute abundance NA included
#GENUS
tax.plot.na<-plot_bar(phylo.rare, fill = "Genus", x="Type")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+
  labs(y="Number of Reads")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand = c(0,0))+ scale_x_discrete(expand = c(0,0))
tax.plot.na


tax.plot.na<-plot_bar(phylo.rare, fill = "Genus", x="Depth.zone")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+
  labs(y="Number of Reads")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand = c(0,0))+ scale_x_discrete(expand = c(0,0))
tax.plot.na

#Save as pdf
ggsave2("Nematoda_metanalysis_tax_rare_sq11480.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#Taxonomic plot by TYPE relative abundance NA
tax.plot.na.rel<-plot_bar(phylo.rare.type.rel, fill = "Genus")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+
  labs(y="Frequency of ASVs")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand = c(0,0))+ scale_x_discrete(expand = c(0,0))
tax.plot.na.rel

#Save as pdf
ggsave2("Nematoda_metanalysis_tax_rare_rel_sq11480.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#FAMILY
tax.plot.na<-plot_bar(phylo.rare.type.rel, fill = "Family", x="Type")+
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+
  labs(y="Relative Abundance")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand = c(0,0))+ scale_x_discrete(expand = c(0,0))
tax.plot.na

#Genera pres-abs per TYPE
tax.plot.gen<-plot_bar(phylo.rare.gen.zone.pa, fill = "Genus")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))+
  labs(y="Number of Genera")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  geom_text(aes(label = stat(y), group = Type), stat = 'summary', fun = sum, vjust = -1)+
  scale_y_continuous(expand = c(0,0),breaks = seq(from=0, to=70, by=10), limits = c(0,70))
tax.plot.gen

#Re-arrange factors by depth
sample_data(phylo.rare.gen.zone.pa)$Depth.zone<-factor(sample_data(phylo.rare.gen.zone.pa)$Depth.zone,
                                               levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))

tax.plot.gen<-plot_bar(phylo.rare.gen.zone.pa)+
  geom_bar(aes(), stat="identity", position="stack")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 0.5, hjust = 0.5))+
  labs(y="Number of Genera")+theme(axis.line=element_line(colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.text = element_text(face = "italic"))+
  geom_text(aes(label = stat(y), group = Type), stat = 'summary', fun = sum, vjust = -1)+
  scale_y_continuous(expand = c(0,0),breaks = seq(from=0, to=70, by=10), limits = c(0,70))
tax.plot.gen

#Save as pdf
ggsave2("Nematoda_metanalysis_rare_type_gen_pa_sq11480.pdf",units = "in",width = 14,height = 8.5,scale = 1)
#----------------------------------------------------------------------------------------------------#
#Ampvis2####

#Load data
myotutable<-read_excel("Nematoda_metanalysis_asv_tax_table_phylo_rare_8548.xlsx", sheet="Sheet1")
myotutable<-read_excel("Nematoda_metanalysis_asv_tax_table_phylo_rare_pa_8548.xlsx")

mymetadata<-read_excel("Nematoda_metanalysis_diversity.xlsx", sheet = "rare.8548")

#Add factors
mymetadata$Depth.zone<-factor(mymetadata$Depth.zone, levels =c("Intertidal", "Shelf", "Bathyal", "Abyssal"))
mymetadata$Type<-recode_factor(mymetadata$Type, 
                                       HydrothermalVent="Hydrothermal vent", MudVolcano="Mud volcano")
mymetadata$Type<-factor(mymetadata$Type, levels = c("Canyon", "Hydrothermal vent", "Mud volcano", 
                                                    "Plain", "Seamount"))
#Correct excel to numeric values!

#Make ampvis2 object
amp<-amp_load(otutable = myotutable, metadata = mymetadata)





#Genus read abundance FACET
amp_heatmap(amp, group_by = "Type", tax_aggregate = "Genus", tax_show = 31,
            showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
            plot_colorscale = "sqrt", facet_by = "Depth_zone",
            plot_legendbreaks = c(seq(from=0, to=60, by=5)))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="% Read Abundance"))+
theme(axis.ticks = element_blank(), axis.text.y = element_text(hjust = 1.1))

#Genus ASV abundance FACET
amp_heatmap(amp, group_by = "Type", tax_aggregate = "Genus", tax_show = 31,
            showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
            plot_colorscale = "sqrt", facet_by = "Depth_zone",
            plot_legendbreaks = c(seq(from=0, to=60, by=5)))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="% ASV Abundance"))+
  theme(axis.ticks = element_blank(), axis.text.y = element_text(hjust = 1.1))

#Genus ASV abundance TYPE
ampvis.type<-amp_heatmap(amp, group_by = "Type", tax_aggregate = "Genus", tax_show = 31,
            showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(seq(from=0, to=60, by=5)))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="Relative read abundance (%)", nrow = 1))+
  theme(axis.ticks = element_blank(), axis.text.y = element_text(hjust = 1.1))
ampvis.type



#Genus ASV abundance DEPTH
ampvis.depth<-amp_heatmap(amp, group_by = "Depth_zone", tax_aggregate = "Genus", tax_show = 31,
            showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(seq(from=0, to=60, by=5)))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="Relative read abundance (%)", nrow = 1))+
  theme(axis.ticks = element_blank(), axis.text.y = element_text(hjust = 1.1))
ampvis.depth


ampvis.depth<-amp_heatmap(amp, group_by = "Depth_zone", tax_aggregate = "Genus", tax_show = 31,
                          showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
                          plot_colorscale = "log10",normalise = FALSE,
                          plot_legendbreaks = c(seq(from=0, to=60, by=5)))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="Relative read abundance (%)", nrow = 1))+
  theme(axis.ticks = element_blank(), axis.text.y = element_text(hjust = 1.1))
ampvis.depth



ggarrange(ampvis.depth, ampvis.type, ncol=2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

#All genera
amp_heatmap(amp, group_by = "Type", tax_aggregate = "Genus",tax_show = "all", 
            showRemainingTaxa = TRUE, tax_empty = "remove", round = 2,
            plot_colorscale = "sqrt",
            plot_legendbreaks = c(seq(from=0, to=60, by=5)), normalise = FALSE)+
  theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))+
  guides(fill=guide_legend(title="% ASV Abundance"))


amp_ordinate(amp, type="pcoa", distmeasure = "bray", sample_colorframe = FALSE, transform = "none",
             sample_color_by = "Depth_zone")


#----------------------------------------------------------------------------------------------------#
#Ggplot extract####
taxplot.data<-ggplot_build(tax.plot.type.no.na.rel)$plot$data[,c(2,3,12)]
rownames(taxplot.data)<-NULL

#Add factors
taxplot.data$Sample<-as.factor(taxplot.data$Sample)
taxplot.data$Genus<-as.factor(taxplot.data$Genus)

#Aggregate by Sample + Genus
taxplot.data<-aggregate(Abundance ~ Sample + Genus, taxplot.data, sum)

#Reshape from long format to wide
taxplot.data.w<-as.data.frame(pivot_wider(taxplot.data, names_from = "Sample", 
                                          values_from = "Abundance"))

#Write to file
write_tsv(taxplot.data.w, file = "Nematoda_metanalysis_tax_type_genus_nematoda.rm1.rare.txt")

#Make rownames of column 1
taxplot.data.w<-column_to_rownames(taxplot.data.w, var = "Genus")

#Remove row names for Upset
rownames(taxplot.data.w)<-NULL
#Use this file for Upset plots

#Transform to otu.table.pa to long format for ggplot
otu.table.pa.l<-pivot_longer(otu.table.pa.df, !Genus,names_to = "Type", values_to = "Count")
#----------------------------------------------------------------------------------------------------#
#Alpha-div####
#The lower and upper hinges correspond to the first and third quartiles.
#The upper whisker extends from the hinge to the largest value no further than 1.5 * IQR from the hinge
#IQR is the inter-quartile range, or distance between the first and third quartiles). 
#The lower whisker extends from the hinge to the smallest value at most 1.5 * IQR of the hinge.
#Data beyond the end of the whiskers are called "outlying" points and are plotted individually.
#!Exclude NA for diversity measures!

#Generate df with diversity indices
alpha.div.df<-estimate_richness(phylo.rare, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.subset, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.gen, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.pa, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.zone, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.plain, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.gen.zone, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))

alpha.div.df<-estimate_richness(phylo.rare.gen.type, measures = c("Observed","Chao1","ACE","Shannon","Simpson"))


#Append to excel
write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.subset",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.gen",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.pa",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.zone",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.plain",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.gen.zone",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(alpha.div.df,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_diversity.xlsx",
            sheetName = "rare.8548.gen.type",col.names = TRUE, row.names = TRUE, append = TRUE)

#Mean + SD####
diversity.alpha<-as.data.frame(read_excel("Nematoda_metanalysis_diversity.xlsx", sheet = "rare.8548"))
diversity.alpha$Depth.zone<-as.factor(diversity.alpha$Depth.zone)

summary.diversity.alpha<-data_summary(diversity.alpha, varname ="Observed", groupnames ="Depth.zone")

mean.diversity.alpha<-as.data.frame(tapply(X=diversity.alpha$Observed, INDEX = diversity.alpha$Depth.zone, FUN=mean)) 

sd.diversity.alpha<-as.data.frame(tapply(X=diversity.alpha$Observed, INDEX = diversity.alpha$Depth.zone, FUN=sd)) 

summary.diversity.alpha<-cbind(mean.diversity.alpha, sd.diversity.alpha)

#----------------------------------------------------------------------------------------------------#
#Alpha-div plots####
#Observed, Simpson, Simpson
#Alpha diversity NA
alpha.div.plot<-plot_richness(phylo.rare.pa, x="Type", 
                         measures=c("Observed","Shannon","Simpson"), color = "Type")+
  geom_boxplot()+ theme(axis.title.x = element_blank(), 
        axis.text.x=element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "none", legend.title = element_blank())+
  theme(strip.text.x =element_text(face="bold",size=12))
alpha.div.plot

alpha.div.plot<-plot_richness(phylo.rare, x="Depth.zone", 
                              measures=c("Observed","Shannon","Simpson"))+
  geom_boxplot()+ theme(axis.title.x = element_blank(), 
                        axis.text.x=element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  theme(legend.position = "none")
alpha.div.plot


#Save as pdf
ggsave2("Nematoda_metanalysis_alpha-div_rare_sq11480.pdf",units = "in",
        width = 14,height = 8.5,scale = 0.95)

alpha.div.subset<-plot_richness(phylo.rare.subset, x="Depth.zone",
                             measures=c("Observed","Shannon","Simpson"),
                             )+ geom_boxplot()+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  theme(legend.position = "none")
alpha.div.subset

#Re-arrange factors by increasing depth
sample_data(phylo.rare)$Depth.zone<-factor(sample_data(phylo.rare)$Depth.zone,
                                           levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))

sample_data(phylo.rare.gen)$Depth.zone<-factor(sample_data(phylo.rare.gen)$Depth.zone,
                                           levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))

sample_data(phylo.rare.pa)$Depth.zone<-factor(sample_data(phylo.rare.pa)$Depth.zone,
                                               levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))


alpha.div<-plot_richness(phylo.rare, x="Depth.zone",
                                measures=c("Observed","Shannon","Simpson"))+ geom_boxplot()+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.4)+
  theme(strip.text.x =element_text(face="bold",size=12))
alpha.div

alpha.div<-plot_richness(phylo.rare, x="Depth.zone",
                                measures=c("Observed","Shannon","Simpson"))+ geom_boxplot()+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  geom_jitter(color="darkgrey", size=3, alpha=0.5, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun=mean, geom="point", shape=18,
                 size=3, color="red")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")
alpha.div

alpha.div<-plot_richness(phylo.rare, x="Depth.zone",
                                measures=c("Observed","Shannon","Simpson"))+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  geom_jitter(color="darkgrey", size=3, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=2), 
                 geom="crossbar", color="black", width=0.5)
alpha.div

#Observed ONLY
alpha.div<-plot_richness(phylo.rare, x="Depth.zone",
                         measures=("Observed"))+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  geom_jitter(color="darkgrey", size=3, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="black", width=0.5)+
  scale_y_continuous(name = "Value")+
  scale_x_discrete(name= "Depth Zone")
alpha.div

alpha.div<-plot_richness(phylo.rare.plain, x="Depth.zone",
                         measures=("Observed"))+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  geom_jitter(color="darkgrey", size=3, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="black", width=0.5)+
  scale_y_continuous(name = "Value")+
  scale_x_discrete(name= "Depth Zone")
alpha.div

alpha.div<-plot_richness(phylo.rare.zone, x="Depth.zone",
                         measures=("Observed"))+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  geom_jitter(color="darkgrey", size=3, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="black", width=0.5)+
  scale_y_continuous(name = "Value")+
  scale_x_discrete(name= "Depth Zone")
alpha.div


alpha.div<-plot_richness(phylo.rare.zone)

alpha.div<-plot_richness(phylo.rare.gen, x="Depth.zone",
                         measures=("Observed"))+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  geom_jitter(color="darkgrey", size=4, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="black", width=0.5)+
  scale_y_continuous(name = "Value")+
  scale_x_discrete(name= "Depth Zone")
alpha.div

#Observed, Shannon, Simpson
alpha.div<-plot_richness(phylo.rare, x="Depth.zone",
                         measures=c("Observed", "Shannon", "Simpson"))+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=0, hjust=0.5, size = 9))+
  theme(legend.position = "none")+
  geom_jitter(color="darkgrey", size=3, position = position_jitter(0))+
  theme(strip.text.x =element_text(face="bold",size=12))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="crossbar", color="black", width=0.5)+
  scale_y_continuous(name = "Value")
alpha.div

#geom_jitter(color="darkgrey", size=3, alpha=0.5)+
#stat_summary(fun=mean, geom="point", shape=22,size=4, color="darkgrey", alpha=1)+

alpha.div<-plot_richness(phylo.rare.pa, x="Depth.zone",
                         measures=c("Observed","Shannon","Simpson"))+ geom_boxplot()+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(angle=0, hjust=0.5, size = 10))+
  axis.text.y=element_text(size = 10)+
  theme(legend.position = "none")+
  theme(strip.text.x =element_text(face="bold",size=12))
alpha.div


alpha.div.subset<-plot_richness(phylo.rare.subset, x="Depth.zone",
                                measures=c("Observed","Shannon","Simpson"))+
  geom_boxplot()+
  theme(axis.title.x = element_blank(), 
        axis.text.x=element_text(hjust=1, size = 9, vjust =0.5))+
  theme(legend.position = "none")
alpha.div.subset

#Alpha diversity GENUS
alpha.div.plot<-plot_richness(phylo.rare.gen, x="Type", 
                              measures=c("Observed","Shannon","Simpson"), color = "Type")+
  geom_boxplot()+ theme(axis.title.x = element_blank(), 
                        axis.text.x=element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  theme(legend.position = "none")
alpha.div.plot

alpha.div.plot<-plot_richness(phylo.rare.gen, x="Depth.zone", 
                              measures=c("Observed","Shannon","Simpson"))+
  geom_boxplot()+ theme(axis.title.x = element_blank(), 
                        axis.text.x=element_text(angle=0, hjust=0.5, size = 9, vjust =1))+
  theme(legend.position = "none")
alpha.div.plot

#Save as pdf
ggsave2("Nematoda_metanalysis_alpha-div_rare_rm1_no_na_sq5115.pdf",units = "in",
        width = 14,height = 8.5,scale = 0.95)
#----------------------------------------------------------------------------------------------------#
#Depth & diversity####
#Load data
div.depth<-read_excel("Nematoda_metanalysis_diversity.xlsx", sheet = "rare.8548")
div.depth<-read_excel("Nematoda_metanalysis_diversity.xlsx", sheet = "rare.8548.subset")
div.depth$Mean.depth<-as.numeric(div.depth$Mean.depth)

#Adding regression line FAILS when shape is defined by Relm.type

#scatterplot ASVs vs Depth
div.depth.plot<-ggplot(div.depth, aes(x=Mean.depth, y=Observed, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  scale_shape_manual(values = c(0:14))+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5500, by=500),limits = c(0,5500))+
  scale_y_continuous(breaks = seq(from=0, to=700, by=100), expand = c(0,0), limits = c(0,700))+
  theme(axis.text.x = element_text(angle=0, hjust=1, size = 9, vjust =0.5))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")+
  theme(legend.position = "right")
div.depth.plot

#Shannon
div.depth.plot<-ggplot(div.depth, aes(x=Mean.depth, y=Shannon, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  scale_shape_manual(values = c(0:14))+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5400, by=200),limits = c(0,5400))+
  scale_y_continuous(breaks = seq(from=2, to=6, by=0.5), expand = c(0,0), limits = c(2,6))+
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  xlab("Mean depth (m)")+ylab("Shannon value")+
  theme(legend.position = "right")
div.depth.plot

#Simpson
div.depth.plot<-ggplot(div.depth, aes(x=Mean.depth, y=Simpson, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  scale_shape_manual(values = c(0:14))+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5400, by=200),limits = c(0,5400))+
  scale_y_continuous(breaks = seq(from=0.7, to=10, by=0.02), expand = c(0,0), limits = c(0.7,1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  xlab("Mean depth (m)")+ylab("Shannon value")+
  theme(legend.position = "right")
div.depth.plot

#Subset to Plains
plain<-subset(div.depth, Type=="Plain")

div.depth.plot.plain<-ggplot(plain, aes(x=Mean.depth, y=Observed, shape=Realm.type))+
  geom_point(size=5)+
  theme_classic()+
  scale_shape_manual(values = c(0:14))+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5500, by=500),limits = c(0,5500))+
  scale_y_continuous(breaks = seq(from=0, to=700, by=100), expand = c(0,0), limits = c(0,700))+
  theme(axis.text.x = element_text(hjust=0.5, size = 9, vjust =0.5))+
  theme(legend.position = c(0.15,0.75))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.plot.plain

#Subset to Offshore&NWAtlantic
onwatlantic<-subset(div.depth, Realm.type=="Offshore&NWAtlantic")

div.depth.plot<-ggplot(onwatlantic, aes(x=Mean.depth, y=Observed, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=2200, by=200),limits = c(0,2200))+
  scale_y_continuous(breaks = seq(from=0, to=700, by=100), expand = c(0,0), limits = c(0,700))+
  theme(axis.text.x = element_text(hjust=0.5, size = 9, vjust =0.5))+
  theme(legend.position = c(0.1,0.8))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.plot

#Subset to Atlantic 
atlantic1<-subset(div.depth, Realm.type==("Offshore&NWAtlantic"))
atlantic2<-subset(div.depth, Realm.type==("OffshoreSAtlantic"))
atlantic3<-subset(div.depth, Realm.type==("NEAtlantic"))
atlantic<-rbind(atlantic1, atlantic2, atlantic3)

div.depth.plot<-ggplot(atlantic, aes(x=Mean.depth, y=Observed, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=2200, by=200),limits = c(0,2200))+
  scale_y_continuous(breaks = seq(from=0, to=700, by=100), expand = c(0,0), limits = c(0,700))+
  theme(axis.text.x = element_text(hjust=0.5, size = 9, vjust =0.5))+
  theme(legend.position = c(0.9,0.25))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.plot

pacific1<-subset(div.depth, Realm.type==("Indo-Pacificseas&IndianOcean"))
pacific2<-subset(div.depth, Realm.type==("Mid-tropicalNPacificOcean"))
pacific3<-subset(div.depth, Realm.type==("Offshoremid-EPacific"))
pacific<-rbind(pacific1, pacific2, pacific3)

div.depth.plot<-ggplot(pacific, aes(x=Mean.depth, y=Observed, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5000, by=500),limits = c(0,5000))+
  scale_y_continuous(breaks = seq(from=0, to=500, by=100), expand = c(0,0), limits = c(0,500))+
  theme(axis.text.x = element_text(hjust=0.5, size = 9, vjust =0.5))+
  theme(legend.position = c(0.2,0.8))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.plot

#Subset to Mediterranean
mediter<-subset(div.depth, Realm.type=="Mediterranean")

div.depth.plot.med<-ggplot(mediter, aes(x=Mean.depth, y=Observed, shape=Realm.type, color=Type))+
  geom_point(size=5)+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=2800, by=200),limits = c(0,2800))+
  scale_y_continuous(breaks = seq(from=0, to=700, by=100), expand = c(0,0), limits = c(0,700))+
  theme(axis.text.x = element_text(hjust=0.5, size = 9, vjust =0.5))+
  theme(legend.position = c(0.1,0.8))+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.plot.med

#Save as pdf
ggsave2(plot=div.depth.plot,"Nematoda_metanalysis_depth_div_rare_sq11480.pdf",units = "in",
        width = 14,height = 8.5,scale = 0.75)

ggarrange(div.depth.plot.med, div.depth.plot.plain, labels = "AUTO", ncol = 1)

#scatterplot ASVs vs Depth + lm regression line
div.depth.lm.plot<-ggplot(div.depth, aes(x=Mean.depth, y=Observed, color=Realm.type))+
  geom_point(size=3)+
  scale_shape_manual(values = c(0:14))+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5400, by=200),limits = c(0,5400))+
  scale_y_continuous(breaks = seq(from=0, to=2000, by=200), expand = c(0,0), limits = c(0,2000))+
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  scale_color_manual(values = getPalette(20))+
  geom_smooth(method="lm",  linetype="dashed",color="black", fill="gray",fullrange=TRUE)+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.lm.plot

#scatterplot ASVs vs Depth + loess local regression line
div.depth.loess.plot<-ggplot(div.depth, aes(x=Mean.depth, y=Observed, color=Realm.type))+
  geom_point(size=3)+
  scale_shape_manual(values = c(0:14))+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_continuous(breaks = seq(from=0, to=5400, by=200),limits = c(0,5400))+
  scale_y_continuous(breaks = seq(from=0, to=1000, by=100), expand = c(0,0), limits = c(0,1000))+
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 9, vjust =0.5))+
  scale_color_manual(values = getPalette(20))+
  geom_smooth(method="loess",span=0.75,linetype="dashed",color="black", fill="gray",fullrange=TRUE)+
  xlab("Mean depth (m)")+ylab("Number of ASvs")
div.depth.loess.plot

#Arrange plots
ggarrange(div.depth.plot, div.depth.lm.plot, div.depth.loess.plot,nrow = 3, 
          common.legend = TRUE, legend = "right")

#Save as pdf
ggsave2("Nematoda_metanalysis_depth_div_rare_rm1_no_na_reg_sq13307.pdf",units = "in",
        width = 14,height = 8.5,scale = 0.95)
#----------------------------------------------------------------------------------------------------#
#Depth regression####

#Check normality of DEPENDENT + INDEPENDENT variables
shapiro.test(div.depth$Observed) #0.7919
shapiro.test(div.depth$Mean.depth) #0.001027

#Non-normal --> Pearson correlation
cor.test(div.depth$Observed, div.depth$Mean.depth, method = "pearson")

#OffshoreNWAtlantic
shapiro.test(onwatlantic$Observed) #0.28
shapiro.test(onwatlantic$Mean.depth) #0.6063

#Normal --> Spearman correlation
cor.test(onwatlantic$Observed, onwatlantic$Mean.depth, method = "spearman")

#Mediterranean
shapiro.test(mediter$Observed) #0.5802
shapiro.test(mediter$Mean.depth) #0.7784

#Normal --> Spearman correlation
cor.test(mediter$Observed, mediter$Mean.depth, method = "spearman")

#Linear regression of diversity and mean depth
#All data
lm.div<-lm(div.depth$Observed~div.depth$Mean.depth)
summary(lm.div)

#OffshoreNWAtlantic
lm.div<-lm(onwatlantic$Observed~onwatlantic$Mean.depth)
lm.div<-lm(plain$Observed~plain$Mean.depth)
summary(lm.div)

#Generalized linear model (GLM) is a flexible generalization of ordinary linear 
#regression that allows for response variables that have error distribution models 
#other than a normal distribution.
#Use POISSON --> used to model COUNT DATA
lm.div<-glm(div.depth$Observed~div.depth$Mean.depth, family = poisson)
summary(lm.div)

#glm equation: ObservedASV = 5.66 + 5.745e-05*Depth

#If you have overdispersion (see if residual deviance is much larger than degrees of freedom),
#you may want to use quasipoisson() instead of poisson().
lm.div<-glm(div.depth$Observed~div.depth$Mean.depth, family = quasipoisson)
summary(lm.div)

#GAM####
#Generalised additive models
#From: http://environmentalcomputing.net/intro-to-gams/

gam.div<-gam(div.depth$Observed~s(div.depth$Mean.depth), family = gamma(link=log))
summary(gam.div)

#Residual plots
par(mfrow = c(2,2))
gam.check(gam.div)
res<-resid(gam.div)
fit<-fitted(gam.div)

#Pick gamma family of function when residuals resemble cone-shape
#Link function: function applied to data value

plot(res~fit)
plot(res~div.depth$Mean.depth)

#AIC
#include the 2 models and obtain the values
aic()


#Model validation with AIC (lower=better) or ANOVA
#----------------------------------------------------------------------------------------------------#
#dbRDA####
#Option 1: Univariate
#Response variable: #ASVs (continuous)
#Explanatory variables: Mean depth (continuous), Realm (discreet), Type (discreet)

#Add factors
div.depth$Realm.type<-as.factor(div.depth$Realm.type)
div.depth$Type<-as.factor(div.depth$Type)
div.depth$Depth.zone<-as.factor(div.depth$Depth.zone)

#dbRDA
depth.dbRDA<-capscale(div.depth$Observed~div.depth$Depth.zone+div.depth$Realm.type+div.depth$Type,
                      distance = "euclidean")

depth.dbRDA<-capscale(div.depth$Observed~div.depth$Mean.depth,
                      distance = "euclidean")

plot(depth.dbRDA)
summary(depth.dbRDA)
anova(depth.dbRDA,by="axis", permu=200)

#Option 2: Multivariate
#Response variable: #ASVs table (continuous)
#Explanatory variables: Mean depth (continuous), Realm (discreet), Type (discreet)
asv.table<-as.data.frame(otu_table(phylo.rare))

depth.dbRDA<-capscale(asv.table~div.depth$Mean.depth+div.depth$Realm.type+div.depth$Type,
                      distance = "euclidean")

#----------------------------------------------------------------------------------------------------#
#MGLMM####

asv.matrix<-as.matrix(asv.matrix)

mglmm<-manyglm(asv.matrix~div.depth$Mean.depth+div.depth$Realm.type+div.depth$Type, family = "poisson")
summary.manyglm(mglmm)


#----------------------------------------------------------------------------------------------------#
#GLMM####
#Generalized Mixed model with random term + Poisson distribution

#Set factors
div.depth$Depth.zone<-factor(div.depth$Depth.zone, 
                             levels = c("Intertidal","Shelf", "Bathyal", "Abyssal"))

fit.pois<-glmer(Observed~Depth.zone + (1 | Realm.type), family = "poisson", data = div.depth)

fit.pois.ext<-glmer(Observed~Depth.zone + Type + (1 | Realm.type), family = "poisson", data = div.depth)

fit.pois.cont.ext<-glmer(Observed~Mean.depth + Type + (1 | Realm.type), family = "poisson", data = div.depth)


fit.pois.cont<-glmer(Observed~Mean.depth + (1 | Realm.type), family = "poisson", data = div.depth)

#Check assumptions
#Residuals vs fitted

#Arrange factors
div.depth$Depth.zone<-factor(div.depth$Depth.zone, levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))
plot(residuals(fit.pois)~fitted(fit.pois), main="Residuals vs Fitted")
plot(residuals(fit.pois.cont)~fitted(fit.pois.cont), main="Residuals vs Fitted")

plot(residuals(fit.pois.ext)~fitted(fit.pois.ext), main="Residuals vs Fitted")
plot(residuals(fit.pois.cont.ext)~fitted(fit.pois.cont.ext), main="Residuals vs Fitted")

qqnorm(residuals(fit.pois.cont))
qqnorm(residuals(fit.pois.ext))
qqnorm(residuals(fit.pois.cont.ext))

summary(fit.pois)
car::Anova(fit.pois, test = "Chisq")
plot(allEffects(fit.pois))

plot(predictorEffects(fit.pois.ext))



summary(fit.pois.ext)
car::Anova(fit.pois.ext, test = "Chisq")
plot(allEffects(fit.pois.ext))


summary(fit.pois.cont)
car::Anova(fit.pois.cont, test = "Chisq")
plot(allEffects(fit.pois.cont))

summary(fit.pois.cont.ext)
car::Anova(fit.pois.cont.ext, test = "Chisq")
plot(allEffects(fit.pois.cont.ext))


#plot effects sequence incorrect - correct manually
fit.pois.ext@frame[["Depth.zone"]]<-factor(fit.pois.ext@frame[["Depth.zone"]], levels = c("Intertidal","Shelf", "Abyssal", "Bathyal"))

summary(fit.pois.ext)
car::Anova(fit.pois.ext, test = "Chisq")
plot(allEffects(fit.pois.ext))


#Compute least squares means of the fitted mode
lsmeans(fit.pois, "Depth.zone")

#Posthoc pairwise comparisons
lsmeans(fit.pois, list(pairwise~Depth.zone), adjust="tukey")

lsmeans(fit.pois.ext, list(pairwise~Type), adjust="tukey")

lsmeans(fit.pois.ext, list(pairwise~Depth.zone), adjust="tukey")

lsmeans(fit.pois.cont.ext, list(pairwise~Type), adjust="tukey")


#----------------------------------------------------------------------------------------------------#
#Stats:Alpha-div####

#Load data
div.depth<-read_excel("Nematoda_metanalysis_diversity.xlsx", 
                      sheet="rare.8548")

div.depth<-read_excel("Nematoda_metanalysis_diversity.xlsx", 
                      sheet="rare.8548.gen")

div.depth<-read_excel("Nematoda_metanalysis_diversity.xlsx", 
                      sheet="rare.8548.plain")

#Add factors
div.depth$Type<-as.factor(div.depth$Type)
div.depth$Realm.type<-as.factor(div.depth$Realm.type)
div.depth$Depth.zone<-as.factor(div.depth$Depth.zone)

#Test normality by TYPE
#Make subsets for samples n<3
abyssal<-subset(div.depth, Type=="AbyssalPlain")
hydrovent<-subset(div.depth, Type=="HydrothermalVent")
intertidal<-subset(div.depth, Type=="Intertidal")
mudvolcano<-subset(div.depth, Type=="MudVolcano")
seamount<-subset(div.depth, Type=="Seamount")
shelf<-subset(div.depth, Type=="Shelf")

#By Depth zone
abyssal<-subset(div.depth, Depth.zone=="Abyssal")
bathyal<-subset(div.depth, Depth.zone=="Bathyal")
intertidal<-subset(div.depth, Depth.zone=="Intertidal")
shelf<-subset(div.depth, Depth.zone=="Shelf")

#Index=Observed
shapiro.test(hydrovent$Observed)
shapiro.test(intertidal$Observed)
shapiro.test(mudvolcano$Observed)
shapiro.test(seamount$Observed)
shapiro.test(shelf$Observed)

shapiro.test(abyssal$Observed)
shapiro.test(bathyal$Observed)
shapiro.test(intertidal$Observed)
shapiro.test(shelf$Observed)

#Index=Shannon
shapiro.test(hydrovent$Shannon)
shapiro.test(intertidal$Shannon)
shapiro.test(mudvolcano$Shannon)
shapiro.test(seamount$Shannon)
shapiro.test(shelf$Shannon)

shapiro.test(abyssal$Shannon)
shapiro.test(bathyal$Shannon)
shapiro.test(intertidal$Shannon)
shapiro.test(shelf$Shannon)

#Index=Simpson
shapiro.test(hydrovent$Simpson)
shapiro.test(intertidal$Simpson)
shapiro.test(mudvolcano$Simpson)
shapiro.test(seamount$Simpson)
shapiro.test(shelf$Simpson)

shapiro.test(abyssal$Simpson)
shapiro.test(bathyal$Simpson)
shapiro.test(intertidal$Simpson)
shapiro.test(shelf$Simpson)

#Test homogeneity of Variances
leveneTest(div.depth$Observed, div.depth$Depth.zone)
leveneTest(div.depth$Shannon, div.depth$Depth.zone)
leveneTest(div.depth$Simpson, div.depth$Depth.zone)

#Make df excluding Antarctic, Bathyal, Canyon
alpha.div.sub<-rbind(abyssal,hydrovent, intertidal, mudvolcano,seamount, shelf)

#1way ANOVA by DEPTH ZONE
x<-aov(alpha.div.sub$Observed~alpha.div.sub$Depth.zone)
x<-aov(alpha.div.sub$Shannon~alpha.div.sub$Depth.zone)
x<-aov(alpha.div.sub$Simpson~alpha.div.sub$Depth.zone)

#ANOVA#
#1way ANOVA by TYPE
x<-aov(div.depth$Simpson~div.depth$Type)
x<-aov(div.depth$Simpson~div.depth$Type)
x<-aov(alpha.div.df$Observed~sample_data(phylo.filt.rare)$Type)

x<-aov(div.depth$Observed~div.depth$Depth.zone)
x<-aov(div.depth$Shannon~div.depth$Depth.zone)
x<-aov(div.depth$Simpson~div.depth$Depth.zone)

#ANOVA table
summary(x)

#Post-hoc
TukeyHSD(x)

#Kruskal-Wallis test
#Non-parametric alternative to 1way ANOVA >2 groups
kruskal.test(alpha.div.sub$Observed~alpha.div.sub$Type)

#Post-hoc test
dunnTest(div.depth$Observed, div.depth$Type, method = "bonferroni")
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
#Align####
#Use DECIPHER alignment

#Load sequences
asv.seqs<-readDNAStringSet("Nematoda_metanalysis_rm1_rare_sq13307.fna")
asv.seqs<-refseq(nematoda.rare.no.na)

#Generate chained guide tree
gT<-lapply(order(width(asv.seqs), decreasing=TRUE),
           function(x) {
             attr(x,"height")<-0
             attr(x,"label")<-names(asv.seqs)[x]
             attr(x,"members")<- 1L
             attr(x,"leaf")<-TRUE
             x
           })
attr(gT,"height")<-0.5
attr(gT,"members")<-length(asv.seqs)
class(gT)<-"dendrogram"

#Nucleotide sequences need to be in the same orientation
asv.seqs<-OrientNucleotides(asv.seqs)

#Perform the alignment
aligned<-AlignSeqs(asv.seqs,guideTree=gT,processors=8)
aligned<-AlignSeqs(asv.seqs,processors=8)

#Adjust alignment
adjust.aligned<-AdjustAlignment(aligned,processors=8)

#Write the alignment to a new FASTA file
writeXStringSet(adjust.aligned,file="Nematoda_metanalysis_rm1_rare_align_sq13307.fna")
#----------------------------------------------------------------------------------------------------#
#ModelTest####
#Select model with lowest AIC
#Load fasta
desmo.align<-read.dna(file = "SO239_desmoscolex_conn_57.fna", format = "fasta")

#Convert to phyDA7at format
desmo.phydat<-phyDat(acantho.align, type = "DNA", levels = NULL)

#Run modeltest
desmo.test<-modelTest(acantho.phydat, model = "all")
#In our case GTR+G+I
#----------------------------------------------------------------------------------------------------#
#Phylogenetic tree####
#Use FastTree externally
#Example code
#FastTree -gtr -nt -gamma Nematoda_metanalysis_rare_align_sq8548.fna > Nematoda_metanalysis_rare_align_sq8548.nwk

library(phangorn)

fit <- pml(treeNJ, data=nema.align)
fit.1<-update(fit, k=4, inv=0.2, bf=baseFreq(nema.align))

library(unix)
rlimit_as(1e19)

fitGTR <- optim.pml(fit.1, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0))

bs<-bootstrap.pml(fitGTR, bs=100, optNni=TRUE,
                  control = pml.control(trace = 0))


#----------------------------------------------------------------------------------------------------#
#Tree rooting####

#Read the tree made with FastTree using ape
tree.unrooted<-read.newick("Nematoda_metanalysis_rare_align_sq8548.nwk")

#Root the tree using phangorn
tree.rooted<-midpoint(tree.unrooted)
#----------------------------------------------------------------------------------------------------#
#Tree to phyloseq####

#Read tree
tree<-read.tree("Nematoda_metanalysis_treeGTR_rooted_mid_8548.nwk")

#Transform tree to be accessible to phyloseq
tree<-phy_tree(tree)
tree<-phy_tree(tree.unrooted)

phylo.rare<-merge_phyloseq(phylo.rare, tree)
phylo.rare.type<-merge_phyloseq(phylo.rare.type, tree)
phylo.rare.zone<-merge_phyloseq(phylo.rare.zone, tree)

#----------------------------------------------------------------------------------------------------#
#NNMDS####
unifracs<-GUniFrac(otu.tab=data,tree=tree_rooted,alpha=c(0, 0.5, 1))$unifracs

#Weighted UniFrac
dw<-unifracs[, , "d_1"]
#Unweighted UniFrac
du<-unifracs[, , "d_UW"]

#Make Metric MDS | k is the number of dim
#Unweighted
fit.uw<-cmdscale(du,eig=TRUE,k=2) 
fit.uw$points
#Weighted
fit.w<-cmdscale(dw,eig=TRUE,k=2) 
fit.w$points

#Save coordinates as dataframe
#Unweighted
mds.uw<-as.data.frame(fit.uw$points)
mds.uw
#Weighted
mds.w<-as.data.frame(fit.w$points)
mds.w

#Save coordinates externally to spreadsheet
write.xlsx2(mds.uw, "Nematoda_metanalysis_mds_coord.xlsx", sheetName = "unweighted.unifrac", 
            append = TRUE)

write.xlsx2(mds.w, "Nematoda_metanalysis_mds_coord.xlsx", sheetName = "weighted.unifrac", 
            append = TRUE)

#Add metadata in excel to colour groups by factors
#Load data  
mds.meta.uw<-read.xlsx("Nematoda_metanalysis_mds_coord.xlsx", sheetIndex = 1)
mds.meta.w<-read.xlsx("Nematoda_metanalysis_mds_coord.xlsx", sheetIndex = 2)

#Arrange factors
mds.meta.uw$Depth.zone<-factor(mds.meta.uw$Depth.zone, levels = c("Intertidal","Shelf", "Bathyal", "Abyssal"))
mds.meta.w$Depth.zone<-factor(mds.meta.w$Depth.zone, levels = c("Intertidal","Shelf", "Bathyal", "Abyssal"))


#Make plots
#Unweighted
#With samples labels
nmds.unweighted<-ggplot(mds.meta.uw, aes(x=Axis.1, y=Axis.2, colour=Type, fill=Type, 
                                         shape=Depth.zone, label=Realm.type))+
  geom_point(size=4)+
  theme(aspect.ratio=1)+ggtitle("PCoA Unweighted UniFrac")+
  scale_shape_manual(values = c(21,22,24,25))+
  geom_text_repel(max.overlaps = 10)
nmds.unweighted

#Unweighted
#Without samples labels
nmds.unweighted<-ggplot(mds.meta.uw, aes(x=Axis.1, y=Axis.2, colour=Type, fill=Type, 
                                         shape=Depth.zone))+
  geom_point(size=4)+
  theme(aspect.ratio=1)+ggtitle("PCoA Unweighted UniFrac")+
  scale_shape_manual(values = c(21,22,24,25))
nmds.unweighted


#Weighted
#With samples labels
nmds.weighted<-ggplot(mds.meta.w, aes(x=Axis.1, y=Axis.2, colour=Type, fill=Type, 
                                         shape=Depth.zone, label=Realm.type))+
  geom_point(size=4)+
  theme(aspect.ratio=1)+ggtitle("PCoA Weighted UniFrac")+
  scale_shape_manual(values = c(21,22,24,25))+
  geom_text_repel(max.overlaps = 10)
nmds.weighted

#Weighted
#Without samples labels
nmds.weighted<-ggplot(mds.meta.w, aes(x=Axis.1, y=Axis.2, colour=Type, fill=Type, 
                                      shape=Depth.zone))+
  geom_point(size=4)+
  theme(aspect.ratio=1)+ggtitle("PCoA Weighted UniFrac")+
  scale_shape_manual(values = c(21,22,24,25))
nmds.weighted


nmds.all<-ggarrange(nmds.unweighted, nmds.weighted, common.legend = TRUE, legend = "right", 
                    ncol = 1)
nmds.all






#PERMDISP on Unifrac distances
anova(betadisper(as.dist(dw),group = div.depth$Type,type="centroid"))
anova(betadisper(as.dist(dw),group = div.depth$Realm,type="centroid"))

anova(betadisper(as.dist(du),factor$Area,type="centroid"))

#PERMANOVA on Unifrac distances
adonis2(as.dist(dw)~div.depth$Type,permutations=9999)
adonis2(as.dist(du)~Area,data=factor,permutations=9999)

#Pairwise PERMANOVA
pairwise.perm.manova(resp=as.dist(dw),fact=div.depth$Type,p.method="bonferroni",nperm=9999)
pairwise.perm.manova(resp=as.dist(du),fact=div.depth$Type,p.method="bonferroni",nperm=9999)









#----------------------------------------------------------------------------------------------------#
#Beta diversity####

#Generate distance matrix
wunifrac.dist<-phyloseq::distance(phylo.rare, method="wunifrac")
unifrac.dist<-phyloseq::distance(phylo.rare, method="unifrac")

wunifrac.dist<-phyloseq::distance(phylo.rare.no.na, method="wunifrac")
unifrac.dist<-phyloseq::distance(phylo.rare.no.na, method="unifrac")

bray.dist<-phyloseq::distance(phylo.rare.gen, method="bray")

#Ordination
wunifrac.ord<-ordinate(phylo.rare, method="PCoA", distance=wunifrac.dist)
unifrac.ord<-ordinate(phylo.rare, method="PCoA", distance=unifrac.dist)
wunifrac.ord<-ordinate(phylo.rare.no.na, method="PCoA", distance=wunifrac.dist)
unifrac.ord<-ordinate(phylo.rare.no.na, method="PCoA", distance=unifrac.dist)
bray.ord<-ordinate(phylo.rare.gen, method="PCoA", distance=bray.dist)


#Plot PCoA####
#Add stat_ellipse() for ellipse
#Add label = "Mean.depth" to label depth
pcoa.wunifrac<-plot_ordination(phylo.rare, wunifrac.ord, shape="Type", 
                               color = "Type")+
  theme(aspect.ratio=1)+ geom_point(size=3)+ggtitle("PCoA Weighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
pcoa.wunifrac

pcoa.bray<-plot_ordination(phylo.rare.gen, bray.ord, shape="Depth.zone", 
                               color = "Type")+
  theme(aspect.ratio=1)+ geom_point(size=3)+ggtitle("PCoA Bray Genus")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
pcoa.bray

pcoa.wunifrac<-plot_ordination(phylo.rare.no.na, wunifrac.ord, shape="Depth.zone", 
                               color = "Type")+
  theme(aspect.ratio=1)+ geom_point(size=4)+ggtitle("PCoA Weighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
pcoa.wunifrac

pcoa.unifrac<-plot_ordination(phylo.rare.no.na, unifrac.ord, shape="Depth.zone", 
                               color = "Type")+
  theme(aspect.ratio=1)+ geom_point(size=4)+ggtitle("PCoA Unweighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
pcoa.unifrac

pcoa.wunifrac<-plot_ordination(phylo.rare, wunifrac.ord, shape="Depth.zone", 
                               color = "Type",label = "Sample")+
  theme(aspect.ratio=1)+ geom_point(size=4)+ggtitle("PCoA Weighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
pcoa.wunifrac

pcoa.unifrac<-plot_ordination(phylo.rare, unifrac.ord, shape="Depth.zone", 
                              color = "Type",label = "Sample")+
  theme(aspect.ratio=1)+ geom_point(size=4)+ggtitle("PCoA Unweighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))+
  geom_vline(xintercept = 0.2)
pcoa.unifrac


pcoa.wnifrac<-plot_ordination(phylo.rare.no.na, wnifrac.ord, shape="Depth.zone", 
                               color = "Type", label = "Sample")+
  theme(aspect.ratio=1)+ geom_point(size=4)+ggtitle("PCoA Unweighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
  #geom_polygon(aes(fill=Depth.zone)) + geom_point(size=3)
pcoa.unifrac

pcoa.unifrac<-plot_ordination(phylo.rare.no.na, unifrac.ord, shape="Depth.zone",
                              color = "Type")+
  theme(aspect.ratio=1)+ geom_point(size=3)+ggtitle("PCoA Unweighted UniFrac")+
  theme(legend.title = element_blank(), legend.text = element_text(size=12))+
  scale_shape_manual(values=c(0:14))
  pcoa.unifrac
  
#Add polygons
#geom_polygon(aes(fill=Depth.zone)) + geom_point(size=3)


#Arrange plots
pcoa.all<-ggarrange(pcoa.unifrac, pcoa.wunifrac, ncol = 2, 
          common.legend = TRUE, legend = "right",  widths = c(2,2,2), heights = c(2,2,2), labels = "AUTO")
pcoa.all

#Save as pdf
ggsave2(plot=pcoa.all,"Nematoda_metanalysis_pcoa_rare_sq11480.pdf",units = "in",
        width = 14,height = 8.5,scale = 0.95)
#----------------------------------------------------------------------------------------------------#
#Stats: Beta-div####
#Betadisper ~Levene test
anova(betadisper(d=unifrac.dist,group=sample_data(phylo.rare)$Type,type="centroid"))
anova(betadisper(d=unifrac.dist,group=sample_data(phylo.rare)$Depth.zone,type="centroid"))

#Permanova
adonis2(unifrac.dist ~ sample_data(phylo.rare)$Type)
adonis2(unifrac.dist ~ sample_data(phylo.rare)$Depth.zone)


pairwise.adonis(x=unifrac.dist, fact=mymetadata$Type,p.adjust.m = "BH")
#----------------------------------------------------------------------------------------------------#
#Clustering####
#Clustering is based on the whole distance whereas ordination represents parts of the distance
#the most it can with 2 dimensions

#Generate distance matrix
wunifrac.dist<-phyloseq::distance(phylo.rare.type, method="wunifrac")
wunifrac.dist<-phyloseq::distance(phylo.rare.zone, method="wunifrac")

wunifrac.dist<-phyloseq::distance(phylo.rare, method="wunifrac")

unifrac.dist<-phyloseq::distance(phylo.rare.type, method="unifrac")
unifrac.dist<-phyloseq::distance(phylo.rare.zone, method="unifrac")

#Define number of plots in window
par(mfrow=c(2, 4))
par(mfrow=c(1, 2))

#Set margins
#par(mar = c(bottom, left, top, right))
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Perform hierarchical clustering
clust<-hclust(wunifrac.dist, method = "ward.D2")
clust<-hclust(wunifrac.dist, method = "ward.D")
clust<-hclust(wunifrac.dist, method = "average")
clust<-hclust(wunifrac.dist, method = "single")
clust<-hclust(wunifrac.dist, method = "complete")
clust<-hclust(wunifrac.dist, method = "median")
clust<-hclust(wunifrac.dist, method = "centroid")
clust<-hclust(wunifrac.dist, method = "mcquitty")

clust<-hclust(unifrac.dist, method = "ward.D2")

#Plot clustering
plot(clust, hang=-1)

#Delinetate clusters with rectangle at specific similarity
#h=similarity, border=colours
rect.hclust(clust, h=0.9,border = 1:9)

#Weighted unifrac delivers more intuitive result
#----------------------------------------------------------------------------------------------------#
#CAP
#Constrained Analysis of Principal Coordinates

#Convert sample_data to data.frame
metadata<-as(sample_data(phylo.rare.type), "data.frame")

#Run cap
cap<-capscale(wunifrac.dist ~ Realm.type, data = metadata)
cap

#Test significance with anova
anova<-anova(cap, permutations = 999)
anova
#----------------------------------------------------------------------------------------------------#
#UpsetR####

#create pres-abs files from OTU tables
otu.table<-as.data.frame(otu_table(phylo.rare.type))
otu.table<-as.data.frame(otu_table(phylo.rare.zone))
otu.table<-as.data.frame(otu_table(phylo.rare.gen.type))
otu.table<-as.data.frame(otu_table(phylo.rare.gen.zone))

#Transpose for vegan
#TAXA: columns | SITE: row
otu.table<-as.data.frame(t(otu.table)) 

#Transform to pres-abs
otu.table.pa<-decostand(otu.table, method = "pa")

#Remove row names
rownames(otu.table.pa)<-NULL


#Change column name
base::colnames(otu.table.pa)[2]<-"Hydrothermal vent"
base::colnames(otu.table.pa)[3]<-"Mud volcano"


#Upset plot SHARED GENERA by DEPTH
nematoda.upset.plot.zone<-upset(otu.table.pa, nsets = 4, nintersects = NA, text.scale = 1.3,
                           sets = c("Abyssal", "Bathyal", "Shelf", "Intertidal"),
                           keep.order = TRUE, mb.ratio = c(0.75, 0.25))
nematoda.upset.plot.zone


#Upset Plot SHARED GENERA by TYPE
nematoda.upset.plot.type<-upset(otu.table.pa, nsets = 8, nintersects = NA, text.scale = 1.3,keep.order = TRUE)
nematoda.upset.plot.type

#Upset plot SHARED ASVs
nematoda.upset.plot.zone<-upset(otu.table.pa, nintersects = 20, text.scale=1.2, nsets = 10, 
                           sets = c("Abyssal", "Bathyal", "Shelf", "Intertidal"),
                           keep.order = TRUE)
nematoda.upset.plot.zone

ggarrange(nematoda.upset.plot.type, nematoda.upset.plot.zone)
#----------------------------------------------------------------------------------------------------#
#ses####
#ses.PD, ses.MNTD, ses.MPD
#Phylogenetic community structure indices with picante

#Load phylogenetic tree
tree<-read.tree("Nematoda_metanalysis_rare_tree_rooted_sq9088.nwk")
tree<-read.tree("Nematoda_metanalysis_treeGTR_rooted_mid_8548.nwk")

#Load community data | Labels MUST correspond to tree tips | Column:Species x Rows:Sites
comm<-as.data.frame(t(otu_table(phylo.rare)))

#PD is not statistically independent of species richness, it positively correlates 
#with species richness across samples. The function ses.pd compares observed PD to the values
#expected under various randomizations and allows a way to standardize for unequal richness across samples.
pd.ses<-ses.pd(samp=comm,tree=tree, null.model=c("taxa.labels","richness","frequency","sample.pool","phylogeny.pool","independentswap","trialswap"), include.root = TRUE, iterations = 100, 99)
pd.ses

pdGTR.ses<-ses.pd(samp=comm,tree=tree, null.model=c("taxa.labels","richness","frequency","sample.pool","phylogeny.pool","independentswap","trialswap"), include.root = TRUE, iterations = 100, 99)
pdGTR.ses

#Calculate SES Mean Pairwise Distance (=Nearest Relative Index [NRI]) unweighted | Requires distance matrix 
mpd.ses.unweighted<-ses.mpd(samp=comm,dis=cophenetic(tree),
                            null.model=c("taxa.labels","richness","frequency",
                                         "sample.pool","phylogeny.pool", 
                                         "independentswap", "trialswap"),abundance.weighted=FALSE,
                            runs=999,iterations=100)
mpd.ses.unweighted

mpdGTR.ses.unweighted<-ses.mpd(samp=comm,dis=cophenetic(tree),
                               null.model=c("taxa.labels","richness","frequency",
                                            "sample.pool","phylogeny.pool", 
                                            "independentswap", "trialswap"),abundance.weighted=FALSE,
                               runs=999,iterations=100)


#Calculate SES Mean Nearest Taxon Distance (=Nearest Taxon Index [NTI]) unweighted | Requires distance matrix
mntd.ses.unweighted<-ses.mntd(samp=comm,dis=cophenetic(tree),
                            null.model=c("taxa.labels","richness","frequency",
                                         "sample.pool","phylogeny.pool", 
                                         "independentswap", "trialswap"),abundance.weighted=FALSE,
                            runs=999,iterations=100)
mntd.ses.unweighted

mntdGTR.ses.unweighted<-ses.mntd(samp=comm,dis=cophenetic(tree),
                                 null.model=c("taxa.labels","richness","frequency",
                                              "sample.pool","phylogeny.pool", 
                                              "independentswap", "trialswap"),abundance.weighted=FALSE,
                                 runs=999,iterations=100)

#Write output to excel
write.xlsx2(pd.ses,"Nematoda_metanalysis_rare_8548_ses.xlsx", sheetName = "pd.ses")
write.xlsx2(mntd.ses.unweighted,"Nematoda_metanalysis_rare_8548_ses.xlsx", 
            sheetName = "mntd.ses", append = TRUE)
write.xlsx2(mpd.ses.unweighted,"Nematoda_metanalysis_rare_8548_ses.xlsx", 
            sheetName = "mpd.ses", append = TRUE)

write.xlsx2(pdGTR.ses,"Nematoda_metanalysis_rare_8548_ses.xlsx", sheetName = "pdGTR.ses")
write.xlsx2(mntdGTR.ses.unweighted,"Nematoda_metanalysis_rare_8548_ses.xlsx", 
            sheetName = "mntdGTR.ses", append = TRUE)
write.xlsx2(mpdGTR.ses.unweighted,"Nematoda_metanalysis_rare_8548_ses.xlsx", 
            sheetName = "mpdGTR.ses", append = TRUE)

#----------------------------------------------------------------------------------------------------#
#Stats: ses####

#Load data
ses.all<-read_xlsx("Nematoda_metanalysis_rare_8548_ses.xlsx",sheet = "ses.all")
ses.all<-read_xlsx("Nematoda_metanalysis_rare_8548_ses.xlsx",sheet = "sesGTR.all")

#Add factors
ses.all$Realm.type<-as.factor(ses.all$Realm.type)
ses.all$Type<-as.factor(ses.all$Type)
ses.all$Metric<-as.factor(ses.all$Metric)
ses.all$Depth.zone<-as.factor(ses.all$Depth.zone)

#Make subsets by TYPE
#Canyon n=1, Mudvolcano n=2
canyon.ses<-subset(ses.all, Type=="Canyon")
mudvolcano.ses<-subset(ses.all, Type=="MudVolcano")
plain.ses<-subset(ses.all,Type=="Plain")
hydrovent.ses<-subset(ses.all, Type=="HydrothermalVent")
seamount.ses<-subset(ses.all, Type=="Seamount")

#By Depth zone
abyssal.ses<-subset(ses.all, Depth.zone=="Abyssal")
bathyal.ses<-subset(ses.all, Depth.zone=="Bathyal")
intertidal.ses<-subset(ses.all, Depth.zone=="Intertidal")
shelf.ses<-subset(ses.all, Depth.zone=="Shelf")

#Row bind >3 observations to use tapply()
#Exclude: intertidal
data.comb<-rbind(abyssal.ses, bathyal.ses, shelf.ses)
data.comb2<-rbind(plain.ses, hydrovent.ses, seamount.ses)
data.comb3<-rbind(plain.ses, hydrovent.ses, seamount.ses, canyon.ses, mudvolcano.ses)

#Check normality with Shapiro-Wilks test | requires >3 observations
tapply(data.comb$ses, data.comb$Depth.zone:data.comb$Metric, shapiro.test)
tapply(data.comb2$ses, data.comb2$Type:data.comb2$Metric, shapiro.test)

#T-test
tapply(data.comb$ses, data.comb$Depth.zone:data.comb$Metric, t.test)
tapply(data.comb2$ses, data.comb2$Type:data.comb2$Metric, t.test)

#Non-parametric alternative to 1sample t-test
#Wilcoxon signed rank test
tapply(data.comb$ses, data.comb$Depth.zone:data.comb$Metric, wilcox.test)
tapply(data.comb2$ses, data.comb2$Type:data.comb2$Metric, wilcox.test)
tapply(intertidal.ses$ses, intertidal.ses$Metric, wilcox.test)


tapply(abyssal$ses, abyssal$Metric, t.test)
tapply(bathyal$ses, bathyal$Metric, t.test)
tapply(hydrovent$ses, hydrovent$Metric, t.test)
tapply(intertidal$ses, intertidal$Metric, t.test)
tapply(seamount$ses, seamount$Metric, t.test)
tapply(shelf$ses, shelf$Metric, t.test)

#Ses plot####
#Load data
data<-read_xlsx("Nematoda_metanalysis_rare_8548_ses.xlsx", sheet = "sesGTR.all")

#Add factors
data$Realm.type<-as.factor(data$Realm.type)
data$Type<-as.factor(data$Type)
data$Metric<-factor(data$Metric, levels = c("ses.PD", "ses.MNTD", "ses.MPD"))
data$Depth.zone<-factor(data$Depth.zone, levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))
data$Result.Depth.zone<-as.factor(data$Result.Depth.zone)
data$Result.Type<-as.factor(data$Result.Type)

ses.plot1<-ggplot(data,aes(x=Depth.zone,y=ses,shape=Result.Type, colour=Depth.zone))+
  geom_jitter(position=position_dodge(0.3),size=5)+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  scale_shape_manual(values=c(0,1,9))+
  geom_hline(yintercept=0,linetype="dashed",color="black")+
  facet_grid(.~Metric,scales="free_y")+
  ylab("Standard Effect Size")+
  theme(strip.text.x =element_text(face="bold",size=12))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_text(size=11))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank(),legend.position="right")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ses.plot1

#Greyscale

par(mar = c(4.1, 4.4, 4.1, 1.9))


ses.plot1<-ggplot(data,aes(x=Type,y=ses,shape=Type, colour=Result.Type))+
  geom_jitter(position=position_dodge(0.3),size=5)+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  scale_shape_manual(values=c(0,1,2,6,5))+
  scale_colour_grey(start=0.1, end=0.8)+
  geom_hline(yintercept=0,linetype="dashed",color="black")+
  facet_grid(.~Metric,scales="free_y")+
  ylab("Standard Effect Size")+
  theme(strip.text.x =element_text(face="bold",size=12))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_text(size=11))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank())+
  theme(legend.position = c(0.1, 0.2))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ses.plot1


ses.plot1<-ggplot(data,aes(x=Depth.zone,y=ses,shape=Depth.zone, colour=Result.Type))+
  geom_jitter(position=position_dodge(0.3),size=5)+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  scale_shape_manual(values=c(0,1,2,6,5))+
  scale_colour_grey(start=0.1, end=0.8)+
  geom_hline(yintercept=0,linetype="dashed",color="black")+
  facet_grid(.~Metric,scales="free_y")+
  ylab("Standard Effect Size")+
  theme(strip.text.x =element_text(face="bold",size=12))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_text(size=11))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank())+
  theme(legend.position = c(0.1, 0.2))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ses.plot1

ses.plot1<-ggplot(data,aes(x=Depth.zone,y=ses,shape=Depth.zone, colour=Result.Depth.zone))+
  geom_jitter(position=position_dodge(0.3),size=5)+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  scale_shape_manual(values=c(0,1,2,6,5))+
  scale_colour_grey(start=0.1, end=0.8)+
  geom_hline(yintercept=0,linetype="dashed",color="black")+
  facet_grid(.~Metric,scales="free_y")+
  ylab("Standard Effect Size")+
  theme(strip.text.x =element_text(face="bold",size=12))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_text(size=11))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank())+
  theme(legend.position = c(0.1, 0.2))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ses.plot1


ses.plot2<-ggplot(data,aes(x=Depth.zone,y=ses,shape=Result.Depth.zone))+
  geom_jitter(position=position_dodge(0.3),size=4)+
  theme_bw()+theme(panel.background=element_blank())+theme(legend.title=element_blank())+
  scale_shape_manual(values=c(0,9))+
  geom_hline(yintercept=0,linetype="dashed",color="black")+
  facet_grid(.~Metric,scales="free_y")+
  ylab("Value")+
  theme(axis.text.x=element_text(size=11, angle = 90, hjust = 1, vjust =0.5))+
  theme(axis.title.y=element_text(size=12))+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank(),legend.position="bottom")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(strip.text=element_text(face="bold",size=10))
ses.plot2

ggarrange(ses.plot1, ses.plot2,common.legend = TRUE, align = "hv", legend = "bottom")

#Save as pdf
ggsave2("Nematoda_metanalysis_ses.pdf",units = "in",width = 14,height = 8.5,scale = 1)

#----------------------------------------------------------------------------------------------------#
#Evolutionary Distinctiveness####

#Calculates evolutionary distinctiveness measures for a suite of species by
#(1) equal splits (Reddingand Mooers 2006)
#(2) fair proportions (Isaac et al., 2007)

#Read tree
tree.rooted<-read.newick("Nematoda_metanalysis_treeGTR_rooted_mid_8548.nwk")

#Calculate ED by (1)
equal.split<-evol.distinct(tree = tree, scale = FALSE, use.branch.lengths = TRUE,
                         type = "equal.splits")

fair.prop<-evol.distinct(tree = tree.rooted, scale = FALSE, use.branch.lengths = TRUE,
                           type = "fair.proportion")
#!!Much more computationally demanding!!
#Did not complete run

#Append seqtab to ED data to split into TYPE
#Transform otu table to pres-abs
otu.table<-as.data.frame(otu_table(phylo.rare.type))
otu.table.pa<-decostand(otu.table, method = "pa")
otu.table.pa<-as.data.frame(t(otu.table.pa))
ed.table<-cbind(equal.split, otu.table.pa)

otu.table<-as.data.frame(otu_table(phylo.rare.zone))
otu.table.pa<-decostand(otu.table, method = "pa")
otu.table.pa<-as.data.frame(t(otu.table.pa))
ed.table<-cbind(equal.split, otu.table.pa)

#Write results ED to file
write.xlsx2(ed.table,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_evol_distinct.xlsx", sheetName = "phylo.rare.zone.8548",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(ed.table,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_evol_distinct_GTR.xlsx", sheetName = "phylo.rare.zone.8548",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(ed.table,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_evol_distinct_GTR.xlsx", sheetName = "phylo.rare.type.8548",col.names = TRUE, row.names = TRUE, append = TRUE)

write.xlsx2(ed.table,file="F:/HTS/Nematoda_metanalysis/Nematoda_metanalysis_evol_distinct.xlsx", sheetName = "phylo.rare.type.8548",col.names = TRUE, row.names = TRUE, append = TRUE)

#!!Manually sort ASV list (A1) in excel


#----------------------------------------------------------------------------------------------------#
#Histogram ED####

#Load data
all.ed<-read.xlsx("Nematoda_metanalysis_evol_distinct.xlsx", sheetName = "phylo.rare.zone.8548")
all.ed<-read.xlsx("Nematoda_metanalysis_evol_distinct.xlsx", sheetName = "phylo.rare.type.8548")

all.ed<-read.xlsx("Nematoda_metanalysis_evol_distinct_GTR.xlsx", sheetName = "phylo.rare.type.8548")
all.ed<-read.xlsx("Nematoda_metanalysis_evol_distinct_GTR.xlsx", sheetName = "phylo.rare.zone.8548")


#Replace zeros with NA
all.ed[all.ed == 0] <- NA

#Re-arrange dataframe to long format
all.ed.long<-pivot_longer(all.ed, cols =4:8, names_to="Type", values_drop_na = TRUE)
all.ed.long<-pivot_longer(all.ed, cols =4:7, names_to="Zone", values_drop_na = TRUE)

#Add factors
all.ed.long$Type<-as.factor(all.ed.long$Type)
all.ed.long$Zone<-factor(all.ed.long$Zone, 
                         levels = c("Intertidal", "Shelf", "Bathyal", "Abyssal"))

#Calculate means
ed.means<-aggregate(x=all.ed.long$w, by=list(all.ed.long$Type), FUN=mean)
base::colnames(ed.means)[2]<-"mean"

ed.means<-aggregate(x=all.ed.long$w, by=list(all.ed.long$Zone), FUN=mean)
base::colnames(ed.means)[2]<-"mean"

#Calculate ED SD
ed.sd<-aggregate(x=all.ed.long$w, by=list(all.ed.long$Type), FUN=sd)
base::colnames(ed.sd)[2]<-"sd"

ed.sd<-aggregate(x=all.ed.long$w, by=list(all.ed.long$Zone), FUN=sd)
base::colnames(ed.sd)[2]<-"sd"

#Change column name to Type to match aes
colnames(ed.means)[1]<-"Type"
colnames(ed.means)[1]<-"Zone"

means.sd<-cbind(ed.means, ed.sd)
means.sd

#Historgram
ed.hist<-ggplot(all.ed.long, aes(x=w))+
  geom_histogram(binwidth = 0.01)+
  theme_bw()+theme(legend.position = "none")+
  facet_grid(Type~., scales = "free_y")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0),breaks = seq(from=0, to=1.6, by=0.05), limits = c(0,1.6)) +  geom_vline(data=ed.means, aes(xintercept=mean, color="red"),linetype="dashed")+
  theme(strip.text.y = element_text(size =10))+
  theme(axis.text.x = element_text(angle = -90))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  xlab("Evolutionary Distinctiveness")+ylab("Count")
ed.hist


ed.hist<-ggplot(all.ed.long, aes(x=w))+
  geom_histogram(binwidth = 0.01)+
  theme_bw()+theme(legend.position = "none")+
  facet_grid(Zone~., scales = "free_y")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0),breaks = seq(from=0, to=2, by=0.05), limits = c(0,2)) +
  theme(axis.text.x = element_text(angle = -90))+
  geom_vline(data=ed.means, aes(xintercept=mean, color="red"),linetype="dashed")+
  theme(strip.text.y = element_text(size =10))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  xlab("Evolutionary Distinctiveness")+ylab("Count")
ed.hist

#Save as pdf
ggsave2("Nematoda_metanalysis_ED_hylo_rare_sq11480.pdf",
        units = "in",width = 14,height = 8.5,scale = 1)

#----------------------------------------------------------------------------------------------------#
#Stast: ED####
#Normality check with Shapiro
tapply(X=all.ed.long$w, all.ed.long$Type, shapiro.test)
leveneTest(all.ed.long$w, group =all.ed.long$Type)

tapply(X=all.ed.long$w, all.ed.long$Zone, shapiro.test)
leveneTest(all.ed.long$w, group =all.ed.long$Zone)
#Non-normal for Shapiro + Levene

#Proceed with Kruskal-Wallis
kruskal.test(all.ed.long$w~all.ed.long$Type)
kruskal.test(all.ed.long$w~all.ed.long$Zone)

#Post-hoc Dunn test
#Correction method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. The false discovery rate is a less stringent condition than the family-wise error rate, so these methods are more powerful than the others.
dunnTest(all.ed.long$w, all.ed.long$Type, method = "bh")
dunnTest(all.ed.long$w, all.ed.long$Zone, method = "bh")
pairwise.wilcox.test(all.ed.long$w, all.ed.long$Zone,p.adjust.method = "bonf")

#----------------------------------------------------------------------------------------------------#
#CCA####
#Canonical Correspondence Analysis 
#Following https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html

#Required data
#1. Community composition (rows: samples | columns: species | A1 empty)
#2. Environmental data (rows: samples | columns: parameters | A1 empty)

#Load matrices
comm.mat<-read.table("Nematoda_metanalysis_cca_asv.txt", header = TRUE)
env.mat<-read.table("Nematoda_metanalysis_cca_env.txt", header = TRUE)

comm.mat.log<-decostand(comm.mat, "log")

ccamodel<-cca(comm.mat~.,env.mat)

summary(ccamodel)

anova.cca(ccamodel)











#DeSeq2####
#The package DESeq2 provides methods to test for differential expression 
#by use of negative binomial generalized linear models; the estimates of 
#dispersion and logarithmic fold changes incorporate data-driven prior distributions.

#Need to convert the Realm and Type column into a factor
sample_data(phylo.prune)$Type<-as.factor(sample_data(phylo.prune)$Type)
sample_data(phylo.prune)$Realm<-as.factor(sample_data(phylo.prune)$Realm)

#Convert the phyloseq object to a DESeqDataSet and run DESeq2
nema.deseq2<-phyloseq_to_deseq2(phylo.prune, ~Type)
nema.deseq2<-DESeq(nema.deseq2, sfType = "iterate")

#Outcome --> Size estimates did not converge
#----------------------------------------------------------------------------------------------------#