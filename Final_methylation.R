library(data.table)

inputA<-read.table("INTERGENIC_ISLANDS_COMMON_ALL_SAMPLES_PROMOTER_METHYLATION.txt", header=T)

#to convert mm9 to mm10, given file must be in "chrN:start-end"
chrN<-paste0("chr",inputA$chr)
mm9<-paste0(chrN,":",inputA$promoter_start,"-",inputA$promoter_end)

#exporting the mm9 ready csv file for coversion through liftover
write.table(data.frame(mm9)[1],"mm9.csv",row.names = FALSE)

#Conversion of mm9 to mm10 produced some errors.

#mm10 file produced in BED format
require(GenomicRanges)
require(rtracklayer)

#Loading mm10 converted file and error file

mm10 <- read.table("hglft_genome_7f10_1d4a50.bed",header=F,colClasses = "character")

error<-read.table("hglft_genome_7f10_1d4a50.err.txt")

#Comparing the error with the mm10 file to determine the location of coversion error and removing it from
#the data by creating a column called 'con-error'
index<-match(error$V1,mm9)
exp_data<-data.frame(chrN, inputA[2:13])
#colnames(exp_data)=c("chr",'start','end','stable_id','strand', 'TSS', 'promoter_start', 'promoter_end', 'br_methylation', 'lv_methylation', 'ts_methylation', 'sc_methylation', 'sd_methylation')
exp_data[,"con_error"]<-NA
exp_data$con_error[index]<-"error"
exp_data<-exp_data[-index,]

# Stripping the data in a table format from the mm10 BED file
a<-apply(mm10,1, function(x) strsplit(x, ":")[[1]][1])
b<-apply(mm10,1, function(x) strsplit(x, ":")[[1]][2])
b1<-lapply(b, function(x) as.numeric(strsplit(x, "-")[[1]][1]))
b2<-lapply(b, function(x) as.numeric(strsplit(x, "-")[[1]][2]))

#Cleaned mm10 file of the input
exp_data<-data.frame(a, exp_data[2:6], unlist(b1), unlist(b2), exp_data[9:14])
colnames(exp_data)=c("chr",'stable_id','start','end','strand', 'TSS', 'promoter_start', 'promoter_end', 'br_methylation', 'lv_methylation', 'ts_methylation', 'sc_methylation', 'sd_methylation','error')
remove(a,b,b1,b2)

#Importing the big Genome file in mm10 format
bgwg<-import("bigfile.bigWig", format = "bigwig")
#Converting bigwig to dataframe
bgwg<-data.frame(as.character(seqnames(bgwg)),as.numeric(start(ranges(bgwg))),as.numeric(end(ranges(bgwg))),(score(bgwg)),strand(bgwg), stringsAsFactors = FALSE)
colnames(bgwg)<-c('chr','start','end', 'score', 'strand')

#---- ANALYSIS ------
require(data.table)
require(plyr)
groupA<- data.table(exp_data)

#setting the table key
setkey(groupA, chr, promoter_start, promoter_end)
groupB<- data.table(bgwg)

#setting the table key
setkey(groupB, chr, start, end)

#overlapped files
over <- foverlaps(groupA, groupB, nomatch = 0)

#Calculating the number of times each row is repeated according to the stable_id
CPG_islands<-ddply(over,.(stable_id), nrow)
colnames(CPG_islands)<-c("stable_id", 'CPG_islands')
per_methylation<-ddply(over,.(stable_id), summarize, 'average'=mean(score))

data1<-merge(exp_data, CPG_islands , by="stable_id", all=TRUE)
data1<-merge(data1, per_methylation, by= 'stable_id', all=TRUE)

data1<-data1[order(data1$chr,data1$start),]
write.csv(data1, "output2.csv")
