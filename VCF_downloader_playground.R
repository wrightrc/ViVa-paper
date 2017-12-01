
#source("Natural_Variation_Shiny_APP/VCF_downloader.R")
require(BSgenome.Athaliana.TAIR.TAIR9)

# this section is mainly a playground for calling functions and testing workflows
#

#geneInfo <- getGeneInfo(IAAgenes)
#geneInfo <- geneInfo[-c(2,6,8,13,15,17,18,19,20,21,23,28,31,32,34,35,37,39,40,42,44,46,48,53),]

geneInfo <- getGeneInfo(TPLgenes)
#geneInfo <- geneInfo[c(3,6,8,12,14), ]

#run loadData on a single transcript_ID, useful for debuging
SNPData = loadData(geneInfo[2,], strains)
#diversity <- Nucleotide_diversity(SNPData)
plotData <- plotPi(SNPData)

syn_loci <- filtR(SNPData,split_var="names",col_name="Effect",value="synonymous_variant")
#NucDiv(syn_loci)

missense_loci <- filtR(SNPData,split_var="names",col_name="Effect",value="missense_variant")
#NucDiv(missense_loci)

tableData <- rbind(simplifySNP(missense_loci), simplifySNP(syn_loci))

#geneInfo <- getGeneInfo(genes)
#tableData = polymorphTable(geneInfo, strains)
#write.csv(tableData,"output\\table.csv")


#SNPData[SNPData$names=="3:9869957_C/T",]

data1 <- readVcf(file="data1.vcf", genome="Athaliana")
data1Granges <- rowRanges(data1)
myStringSet <- data1Granges[271]$ALT[[1]]

Athaliana<-BSgenome.Athaliana.TAIR.TAIR9
subseq(Athaliana$Chr1, start(data1Granges[271:272]), end(data1Granges[271:272]))

tair10 <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")
attributePages(tair10)
listAttributes(tair10, page="feature_page")

output <- getBM(attributes=c("ensembl_transcript_id", "chromosome_name", "start_position",
                             "end_position", "strand", "transcript_start", "transcript_end", "pfam", "pfam_start", "pfam_end",
                             "gene3d", "gene3d_start", "gene3d_end"
), filters="tair_locus", values=IAAgenes[1:5], mart=tair10)


# NUCLEOTIDE DIVERSITY CALCULATION AND PLOTTING WITH TIDYVCF
###=============================================================================
#get gene info for single gene
geneInfo <- getGeneInfo("AT2G38690",  useCache=FALSE)
#download VCF in tidy format
TIR1_VCF <- VCFByTranscript(geneInfo, strains)
#parse the EFF field
TIR1_VCF <- parseEFF(TIR1_VCF, geneInfo$transcript_ID)
#calculate Nucleotide diversity at each position
TIR1_VCF <- Nucleotide_diversity(TIR1_VCF)
#list all Effects:
unique(TIR1_VCF$Effect)
#make a data frame with only coding (missense and synonymous) variants
coding_variants <- TIR1_VCF$dat[TIR1_VCF$dat$Effect %in% c("missense_variant", "synonymous_variant"), ]
#extract uniuqe position and effect
unique_coding_variants <- unique(coding_variants[ , c("POS", "Effect",
                                                      "Amino_Acid_Change",
                                                      "Diversity") ])
#add codon number to unique_coding_variants
unique_coding_variants <-ddply(unique_coding_variants, .fun=codonNumberKernel,
                               .variables=c("POS", "Amino_Acid_Change"))

#plot the diversity
plot <- ggplot(unique_coding_variants, aes(x=Codon_Number,y=Diversity, colour=Effect)) +
  geom_point() +
  scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1),limits=c(0.0001, 1)) +
  #scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
  ylab("nucleotide diversity, log scale")
print(plot)
print(plyr::count(unique_coding_variants, "Effect"))


# POLYMORPH TABLE WITH TIDYVCF
###=============================================================================
# get gene info for multiple genes
geneInfo <- getGeneInfo(genes)
# set up empty table
effects <- c("5_prime_UTR_variant",
             "intron_variant",
             "3_prime_UTR_variant",
             "synonymous_variant",
             "missense_variant",
             "upstream_gene_variant")
tableData <- matrix(nrow=length(geneInfo$transcript_ID), ncol=length(effects) + 6)
rownames(tableData) <- geneInfo$transcript_ID
colnames(tableData) <- c(effects, "coding_total", "Pi_coding", "Pi_non_syn", "Pi_syn", "Pi_NS_Ratio", "Pi_transcript")
tableData <- data.frame(tableData)

#for each transcript
for (i in 1:length(geneInfo$transcript_ID)) {
 tidyVCF <- VCFByTranscript(geneInfo[i, ], strains)
 tidyVCF <- parseEFF(tidyVCF, geneInfo$transcript_ID[i])
 data <- Nucleotide_diversity(tidyVCF$dat)

 #fill in the first part of the table
 variant_counts <- plyr::count(data, "Effect")
 for (j in 1:length(effects)){
   if (effects[j] %in% variant_counts$Effect){
     tableData[i,j] <- variant_counts[variant_counts$Effect %in% effects[j], "freq"]
   } else {
     tableData[i,j] <- 0
   }
 }
 tableData[i, "coding_total"] <- tableData[i, "missense_variant"] + tableData[i, "synonymous_variant"]

 #nucleotide diversity sums:
 reducedData <- unique(data[, c("POS", "Effect", "Diversity")])
 AA_Length <- unique(data$Amino_Acid_Length)
 AA_Length <- as.numeric(AA_Length[!is.na(AA_Length)])

 tableData[i, "Pi_non_syn"] <- sum(reducedData[reducedData$Effect %in% "missense_variant", "Diversity"]) / (3*AA_Length)
 tableData[i, "Pi_syn"] <- sum(reducedData[reducedData$Effect %in% "synonymous_variant", "Diversity"]) / (3*AA_Length)
 tableData[i, "Pi_NS_Ratio"] <- tableData[i, "Pi_non_syn"] / tableData[i, "Pi_syn"]

 tableData[i, "Pi_coding"] <- sum(unique(reducedData[reducedData$Effect %in% c("synonymous_variant","missense_variant") , c("POS", "Diversity")])$Diversity) / (3*AA_Length)
 tableData[i, "Pi_transcript"] <- sum(unique(reducedData[, c("POS", "Diversity")])$Diversity) / geneInfo[i, "transcript_length"]

}


# Color coding test
###=============================================================================

geneInfo <- getGeneInfo("AT1G80490")



tidyVCF <- VCFByTranscript(geneInfo[1, ], strains)
# Parse the EFF field
data <- parseEFF(tidyVCF$dat, geneInfo$transcript_ID[1])

# calculate diversity
data <- Nucleotide_diversity(tidyVCF$dat)

data[which(data$Diversity >= 10**-4), ]




# for TIR1
geneInfo <- getGeneInfo("AT3G62980",  useCache=FALSE)
tidyVCF <- VCFByTranscript(geneInfo[1, ], strains)
data <- tidyVCF$dat[tidyVCF$dat$gt_GT != "0|0",]
# Parse the EFF field
data <- parseEFF(tidyVCF = data, Transcript_ID = geneInfo$transcript_ID[1])

# calculate diversity
data <- Nucleotide_diversity(data)

keyPOS <- unique(data[which(data$Diversity >= 0.5*max(data$Diversity)), "POS"])

#keyPOS

keydata <- data[data$POS %in% keyPOS, ]
#keydata
#keydata2 <- keydata[keydata$Effect %in% c("missense_variant"), ]

keydata_labeled <- label_bySNPs(keydata)

map_test(keydata_labeled)

keyGT <- ddply(keydata, .variables="Indiv", .fun=buildGT)

#-------------------------------------------------------------------


geneInfo <- getGeneInfo("AT1G80490")
tidyVCF <- VCFByTranscript(geneInfo[1, ], strains)
data <- tidyVCF$dat
# Parse the EFF field
data <- parseEFF(tidyVCF = data, Transcript_ID = geneInfo$transcript_ID[1])

# calculate diversity
data <- Nucleotide_diversity(data)

data <- data[data$gt_GT != "0|0", ]


effects <- unique(data$Effect)


keyPOS <- unique(data[which(data$Diversity >= 0.5*max(data$Diversity)), "POS"])

#keyPOS

keydata <- data[data$POS %in% keyPOS, ]
#keydata
#keydata2 <- keydata[keydata$Effect %in% c("missense_variant"), ]

keydata_labeled <- label_bySNPs(keydata)

map_test(keydata_labeled)

keyGT <- ddply(keydata, .variables="Indiv", .fun=buildGT)

