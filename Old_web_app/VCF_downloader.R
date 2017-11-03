###=============================================================================

### LIBRARIES ==================================================================

library(VariantAnnotation)
library(plyr)
library(stringr)
library(biomaRt)
library(ggplot2)
library(reshape2)
library(vcfR)
library(dplyr)

### CODE INPUTS ================================================================

strains <- c(88,108,139,159,265,350,351,403,410,424,428,430,470,476,484,504,506,531,544,
             546,628,630,680,681,685,687,728,742,763,765,766,768,772,801,853,854,867,868,
             870,915,932,991,992,997,1002,1006,1061,1062,1063,1066,1070,1158,1166,1254,
             1257,1313,1317,1552,1612,1622,1651,1652,1676,1684,1739,1741,1756,1757,1793,
             1797,1819,1820,1829,1834,1835,1851,1852,1853,1872,1890,1925,1942,1943,1954,
             2016,2017,2031,2053,2057,2081,2091,2106,2108,2141,2159,2166,2171,2191,2202,
             2212,2239,2240,2276,2278,2285,2286,2317,2370,2412,4779,4807,4826,4840,4857,
             4884,4900,4939,4958,5023,5104,5151,5165,5210,5236,5249,5253,5276,5279,5349,
             5353,5395,5486,5577,5644,5651,5717,5718,5720,5726,5741,5748,5757,5768,5772,
             5776,5779,5784,5798,5800,5811,5822,5830,5831,5832,5836,5837,5856,5860,5865,
             5867,5874,5890,5893,5907,5921,5950,5984,5993,6008,6009,6010,6011,6012,6013,
             6016,6017,6019,6020,6021,6022,6023,6024,6025,6030,6034,6035,6036,6038,6039,
             6040,6041,6042,6043,6046,6064,6069,6070,6071,6073,6074,6076,6077,6085,6086,
             6087,6088,6090,6091,6092,6094,6095,6096,6097,6098,6099,6100,6101,6102,6104,
             6105,6106,6107,6108,6109,6111,6112,6113,6114,6115,6118,6119,6122,6123,6124,
             6125,6126,6128,6131,6132,6133,6134,6136,6137,6138,6140,6141,6142,6145,6148,
             6149,6150,6151,6153,6154,6163,6166,6169,6172,6173,6174,6177,6180,6184,6188,
             6189,6191,6192,6193,6194,6195,6198,6201,6202,6203,6209,6210,6214,6216,6217,
             6218,6220,6221,6231,6235,6237,6238,6240,6241,6242,6243,6244,6252,6255,6258,
             6268,6276,6284,6296,6390,6396,6413,6424,6434,6445,6680,6739,6740,6744,6749,
             6750,6805,6806,6814,6830,6897,6898,6900,6901,6903,6904,6907,6908,6909,6911,
             6913,6915,6917,6918,6919,6920,6922,6923,6924,6926,6927,6929,6931,6932,6933,
             6938,6940,6943,6944,6945,6951,6956,6957,6958,6959,6960,6961,6963,6966,6967,
             6968,6969,6970,6971,6973,6974,6975,6976,6979,6981,6982,6984,6986,6987,6989,
             6990,6992,6997,7000,7002,7003,7008,7013,7014,7025,7026,7028,7031,7033,7036,
             7058,7061,7062,7063,7064,7067,7068,7071,7072,7077,7081,7092,7094,7096,7102,
             7103,7106,7107,7109,7111,7117,7119,7120,7125,7126,7127,7130,7133,7138,7143,
             7147,7158,7160,7161,7162,7163,7164,7165,7169,7177,7181,7183,7186,7192,7199,
             7202,7203,7207,7208,7209,7213,7217,7218,7223,7231,7236,7244,7248,7250,7255,
             7258,7268,7273,7276,7280,7282,7287,7288,7296,7298,7305,7306,7307,7314,7316,
             7319,7320,7322,7323,7327,7328,7332,7333,7337,7342,7343,7344,7346,7347,7349,
             7350,7353,7354,7356,7358,7359,7372,7373,7377,7378,7382,7383,7384,7387,7394,
             7396,7404,7411,7413,7415,7416,7417,7418,7419,7424,7427,7430,7460,7461,7471,
             7475,7477,7514,7515,7516,7517,7520,7521,7523,7525,7529,7530,7566,7568,7717,
             7757,7767,7917,7947,8037,8057,8077,8132,8171,8214,8222,8227,8230,8231,8233,
             8234,8235,8236,8237,8238,8239,8240,8241,8242,8243,8244,8246,8247,8249,8256,
             8258,8259,8264,8283,8284,8285,8290,8297,8306,8307,8311,8312,8326,8334,8335,
             8337,8343,8351,8354,8357,8365,8366,8369,8376,8386,8387,8419,8420,8422,8424,
             8426,8427,8464,8483,8699,8723,9027,9057,9058,9067,9069,9070,9075,9078,9079,
             9081,9084,9085,9089,9091,9095,9099,9100,9102,9103,9104,9105,9106,9111,9113,
             9114,9115,9121,9125,9128,9130,9131,9133,9134,9298,9312,9314,9321,9323,9332,
             9336,9339,9343,9352,9353,9363,9369,9370,9371,9380,9381,9382,9383,9386,9388,
             9390,9391,9392,9394,9395,9399,9402,9404,9405,9407,9408,9409,9412,9413,9416,
             9421,9427,9433,9436,9437,9442,9450,9451,9452,9453,9454,9455,9470,9471,9476,
             9481,9503,9506,9507,9508,9509,9510,9511,9512,9513,9514,9515,9517,9518,9519,
             9520,9521,9522,9523,9524,9525,9526,9527,9528,9529,9530,9531,9532,9533,9534,
             9535,9536,9537,9539,9540,9541,9542,9543,9544,9545,9546,9547,9548,9549,9550,
             9551,9552,9553,9554,9555,9556,9557,9558,9559,9560,9561,9562,9564,9565,9567,
             9568,9569,9571,9573,9574,9576,9577,9578,9579,9581,9582,9583,9584,9585,9586,
             9587,9588,9589,9590,9591,9592,9593,9594,9595,9596,9597,9598,9599,9600,9601,
             9602,9606,9607,9608,9609,9610,9611,9612,9613,9615,9616,9617,9619,9620,9621,
             9622,9624,9625,9626,9627,9628,9629,9630,9631,9632,9633,9634,9635,9636,9637,
             9638,9639,9640,9641,9642,9643,9644,9645,9646,9647,9648,9649,9651,9653,9655,
             9656,9657,9658,9659,9660,9661,9663,9664,9665,9666,9667,9668,9669,9670,9671,
             9672,9673,9676,9677,9678,9679,9680,9681,9682,9683,9684,9685,9686,9687,9689,
             9690,9691,9692,9693,9694,9695,9696,9697,9698,9699,9700,9701,9703,9704,9705,
             9706,9707,9708,9709,9710,9711,9712,9713,9714,9716,9717,9718,9719,9720,9721,
             9722,9723,9725,9726,9727,9728,9729,9730,9731,9732,9733,9735,9736,9737,9738,
             9739,9741,9743,9744,9745,9747,9748,9749,9754,9755,9756,9757,9758,9759,9761,
             9762,9764,9766,9768,9769,9770,9771,9772,9774,9775,9776,9777,9778,9779,9780,
             9781,9782,9783,9784,9785,9786,9787,9788,9789,9790,9791,9792,9793,9794,9795,
             9796,9797,9798,9799,9800,9801,9802,9803,9804,9805,9806,9807,9808,9809,9810,
             9811,9812,9813,9814,9815,9816,9817,9819,9820,9821,9822,9823,9824,9825,9826,
             9827,9828,9830,9831,9832,9833,9834,9835,9836,9837,9838,9839,9840,9841,9843,
             9844,9845,9846,9847,9848,9849,9850,9851,9852,9853,9854,9855,9856,9857,9858,
             9859,9860,9861,9862,9864,9866,9867,9868,9869,9870,9871,9873,9874,9875,9876,
             9877,9878,9879,9880,9881,9882,9883,9885,9886,9887,9888,9890,9891,9892,9894,
             9895,9897,9898,9899,9900,9901,9902,9903,9904,9905,9906,9908,9909,9910,9911,
             9912,9914,9915,9917,9918,9920,9921,9924,9925,9926,9927,9928,9929,9930,9932,
             9933,9935,9937,9938,9939,9941,9942,9943,9944,9945,9946,9947,9948,9949,9950,
             9951,9952,9953,9955,9956,9957,9958,9959,9960,9962,9963,9964,9965,9966,9968,
             9969,9970,9971,9972,9973,9974,9975,9976,9978,9979,9980,9981,9982,9983,9984,
             9985,9986,9987,9988,9990,9991,9993,9995,9996,9997,9998,9999,10001,10002,10004,
             10005,10006,10008,10009,10010,10011,10012,10013,10014,10015,10017,10018,10020,
             10022,10023,10027,14312,14313,14314,14315,14318,14319,15560,15591,15592,15593,
             18694,18696,19949,19950,19951)

effects <- c("5_prime_UTR_variant",
             "intron_variant",
             "3_prime_UTR_variant",
             "synonymous_variant",
             "missense_variant",
             "upstream_gene_variant")

genes <- c("AT3G62980",
           "AT4G03190",
           "AT3G26810",
           "AT1G12820",
           "AT4G24390",
           "AT5G49980",
           "AT2G39940")

IAAgenes <- c("AT4G14560",
           "AT3G23030",
           "AT1G04240",
           "AT5G43700",
           "AT1G15580",
           "AT1G52830",
           "AT3G23050",
           "AT2G22670",
           "AT5G65670",
           "AT1G04100",
           "AT4G28640",
           "AT1G04550",
           "AT2G33310",
           "AT4G14550",
           "AT1G80390",
           "AT3G04730",
           "AT1G04250",
           "AT1G51950",
           "AT3G15540",
           "AT2G46990",
           "AT3G16500",
           "AT4G29080",
           "AT5G25890",
           "AT4G32280",
           "AT3G62100",
           "AT3G17600",
           "AT2G01200",
           "AT5G57420",
           "AT1G15050")

TPLgenes <- c("AT1G15750",
              "AT1G80490",
              "AT3G16830",
              "AT5G27030",
              "AT3G15880")

### FUNCTIONS ====================================================================================

#' Make a string containing the desired genomic ranges
#' of the format "chrom:start-stop"
#' @param data a data.frame with columns `chromosome_name`, `transcript_start`, and `transcript_end`
#'
#' @return a vector of strings with the format "chrom:start-stop"
#' @export
#'
#' @examples
#'
makeRegionString <- function (data) {
  # format "chrom:start-stop"
  return(paste(c(as.character(data$chromosome_name),
                 ":",
                 as.character(data$transcript_start),
                 "-",
                 as.character(data$transcript_end)), collapse=""))
}

#' Download and save a .VCF file from 1001genomes.org.
#'
#' @param fName the file name the downloaded vcf will be saved as eg.
#'    "data1.vcf"
#' @param strainStr A character string of comma separated strain or
#'    "Ecotype ID"s used in the construction of the URL
#' @param regionStr A character string of the format "[chrom]:[start]-[stop]"
#'    where [chrom], [start], and [stop] are numbers. this string is used
#'    directly in the construction of the URL
#' @param download A logical, if TRUE (default) the file will be downloaded, if
#'     FALSE, the file will not be downloaded, the URL will still be returned
#'
#' @return The URL used for the download
#' @export
#'
#' @examples
downloadData <- function (fName, strainStr, regionStr,
                          download=TRUE) {
  url <- c("http://tools.1001genomes.org/api/v1/vcfsubset/strains/", strainStr, "/regions/",
           regionStr, "/type/snpeff/format/vcf")
  url <- paste(url, collapse="")
  if (download == TRUE) {
    download.file(url, fName)
  }
  return (url)
}




#' Download two vcf files from 1001genomes.org and merge their contents
#'
#' @param fName a character string of the file name to save the final .vcf file to
#' @param strainVect numeric vector list of the strains
#' @param regionStr A character string of the format "[chrom]:[start]-[stop]"
#'    where [chrom], [start], and [stop] are numbers. this string is used
#'    directly in the construction of the URL
#'
#' @return a vcfR object of teh combined vcf file
#' @export
#'
#' @examples \code{my.VCF <- downloadMerge("myFullVCF.vcf.gz", strains, "1:4368760-4371298")}
downloadMerge <- function (fName, strainVect, regionStr) {
  strains <- as.character(strainVect)
  # the URL to download the VCF for all 1135 strains is too long so we have to split it into two sets

  splitPoint <- round(length(strains) / 2) #where to split the strain list

  # first file
  strainString <- paste(as.character(strains[1:splitPoint]), collapse=",")
  fileName <- "tempData1.vcf"
  downloadData(fileName, strainString, regionStr)
  # second file
  strainString <- paste(as.character(strains[(splitPoint + 1):length(strains)]), collapse=",")
  fileName <- "tempData2.vcf"
  downloadData(fileName, strainString, regionStr)
  #
  #load the two temporary vcf files then delete the files
  data1 <- read.vcfR("tempData1.vcf", verbose=FALSE, convertNA=TRUE)
  if (file.exists("tempData1.vcf")) file.remove("tempData1.vcf")

  data2 <- read.vcfR("tempData2.vcf", verbose=FALSE, convertNA=TRUE)
  if (file.exists("tempData2.vcf")) file.remove("tempData2.vcf")

  #combine the two gt fields, and write a new combined vcf file.
  data1@gt <- cbind(data1@gt, data2@gt[,-1])
  write.vcf(data1, fName)

  return (data1)  # return the vcfR object of the combined .vcf
}


EFFParse <- function (data) {
  # when applied to a row with "EFF" field, appends the parsed effect fields to the row
  # for use with adply()
  #
  # Args:
  #   data: a single row of a dataframe with a column labeled "EFF" conaining a list of text
  #         strings of the following format:
  #           "Effect(Effect_Impact|Functional_Class|Codon_Change|Amino_Acid_Change|Amino_Acid_length|Gene_Name|Transcript_BioType|Gene_Coding|Transcript_ID|Exon_Rank|Genotype_Number)"
  #         data must contain columns "EFF", "GT"(numeric), "transcript_ID"
  #
  # Returns:
  #   original data row, with new columns containing parsed "EFF" field info appended to the end.
  #   parser selects the correct effect instance based on the "GT" and "transcript_ID" columns in data
  #
  # column names for the EFF fields after parsing. Note, this could be done by reading the VCF header
  EFFColNames = c("Effect", "Effect_Impact", "Functional_Class", "Codon_Change",
                  "Amino_Acid_Change", "Amino_Acid_Length", "Gene_Name", "Transcript_BioType",
                  "Gene_Coding", "Transcript_ID", "Exon_Rank", "Genotype_Number")
  effect <- data$EFF[[1]]
  effect <- data.frame(str_split(effect, pattern="\\(|\\||\\)", simplify=TRUE))
  colnames(effect) <- c(EFFColNames, "")
  effect <- effect[, 1:12, drop=FALSE] ##remove last column of empty strings ""
  effect <- effect[as.character(effect$Transcript_ID) == data$transcript_ID, , drop=FALSE]
  if (nrow(effect) == 0) {
    # if there are no rows left, exit function and return NULL
    return (NULL)
  }
  if ((as.character(data$GT) == "0") %in% TRUE) {
    # if the genotype we want is 0, set several fields to "REF"
    effect$"Effect" <- effect$"Effect_Impact" <- effect$"Functional_Class" <- effect$"Codon_Change" <- effect$"Amino_Acid_Change" <- "REF"
    effect$Genotype_Number <- 0
  }
  effect <- effect[effect$Genotype_Number == data$GT, , drop=FALSE]
  effect <- effect[1, , drop=FALSE]  ###use first effect, (temporary to get code to work)
  return (cbind(data, effect))
}

EFFParse2 <- function (tidyVCF, Transcript_ID){
  EFFColNames = c("Effect", "Effect_Impact", "Functional_Class", "Codon_Change",
                  "Amino_Acid_Change", "Amino_Acid_Length", "Gene_Name", "Transcript_BioType",
                  "Gene_Coding", "Transcript_ID", "Exon_Rank", "Genotype_Number")

  data <- tidyVCF$dat

  output <- tidyVCF
  output$dat <- ddply(data, "POS", .fun=EFFParse2Kernel, Transcript_ID, EFFColNames)

  return (output)

}

EFFParse2Kernel <- function (data, Transcript_ID, EFFColNames){
  if (length(unique(data$EFF)) > 1) {
    print("warning multiple effects found")
  }

  effect <- unique(data$EFF)
  # split by comma to generate a vector of different effects
  effect <- str_split(effect, pattern=",", simplify=TRUE)
  # split by "(", "|", and ")" to separate fields
  effect <- data.frame(str_split(effect, pattern="\\(|\\||\\)", simplify=TRUE), stringsAsFactors = FALSE)
  # remove last column of empty strings ""found after ")"
  effect <- effect[, 1:12, drop=FALSE]
  # add column names to effects
  colnames(effect) <- c(EFFColNames)
  # only keep effects that match the transcript ID
  effect <- effect[effect$Transcript_ID == Transcript_ID, ]

  if (nrow(effect) > 0){   # if there are some effects remaining:
    #create a "gt_GT" column in the effect dataframe that matches the format of the VCF$dat
    effect$gt_GT <- paste(effect$Genotype_Number, "|", effect$Genotype_Number, sep = "")
    # merge the effect df with the original data df by the gt_GT field
    output <- merge(data, effect, by="gt_GT", all.x=TRUE)
  }
  else{  #if there are no effects matching the transcript ID, return the data unaltered
    output <- data
  }

  return(output)
}




#' Parse the correct genotypes from a VRanges object
#' @description The `readVCFasVRanges` function for some reason always chooses the first genotype to use in creating the VRanges. This function matches the alternate allele with the correct genotype
#'
#' @param dataRow
#'
#' @return
#' @export
#'
#' @examples
GTParse <- function (dataRow) {

  # for use with adply, the dataframe dataRow comes from should contain
  Genotypes <- unique(dataRow)  #create a list of the unique genotypes present in this row
  Genotypes <-  Genotypes[which(Genotypes != ".")]  #remove "." (undetermined Genotype) from list
  if (length(Genotypes) == 0) {
    return (NULL)
  }
  GTNumeric <- as.numeric(str_split(Genotypes, pattern="\\|", simplify=TRUE)[,1])
  output <- cbind.data.frame("Genotype"=Genotypes, "Accession_count"=NA, "Accession_list"=NA, "GT"=GTNumeric)
  for (i in 1:length(Genotypes)) {
    # for each unique genotype present in the input row, create a output row with a count and list
    # of accessions with that genotype
    output$Accession_list[i] <- list(names(which(dataRow == Genotypes[i])))
    output$Accession_count[i] <- length(output$Accession_list[i][[1]])
  }
  return (output)
}


#' Load a
#'
#' @param geneInfo
#' @param strains
#'
#' @return
#' @export
#'
#' @examples
loadData <- function (geneInfo, strains) {
  # this function is intended to produce a useful dataframe with each row representing a unique
  # polymorphism, and containing a count and list of accessions that contain that polymorphism.
  # geneInfo table or DF must contain "transcript_ID" and "regionStrings" columns, and have only
  # a single row see getGeneInfo()
  transcript_ID <- as.character(geneInfo$transcript_ID)
  regionString <- as.character(geneInfo$regionString)
  strains <- as.character(strains)
  # the URL to download the VCF for all 1135 strains is too long so we have to split it into two sets
  splitPoint <- round(length(strains) / 2) #where to split the strain list
  # first file
  strainString <- paste(as.character(strains[1:splitPoint]), collapse=",")
  fileName <- "data1.vcf"
  downloadData(fileName, strainString, regionString)
  # second file
  strainString <- paste(as.character(strains[(splitPoint + 1):length(strains)]), collapse=",")
  fileName <- "data2.vcf"
  downloadData(fileName, strainString, regionString)
  #
  data1 <- readVcf(file="data1.vcf", genome="Athaliana")
  data2 <- readVcf(file="data2.vcf", genome="Athaliana")
  #file.remove("data1.vcf")     #comment out these lines for easy review of last file downloaded
  #file.remove("data2.vcf")

  GenotypeMatrix <- cbind(geno(data1)[[1]], geno(data2)[[1]])
  GT <- adply(GenotypeMatrix, .fun=GTParse, .margins=1, .id="names")

  data <- data1
  output <- cbind(data.frame(ranges(rowRanges(data))), mcols(rowRanges(data)))[c("names", "REF", "ALT")]
  output <- cbind(output, data.frame(info(data)["EFF"]))
  output <- merge(output, GT, by="names", all.GT=TRUE)
  output <- cbind.data.frame(output, "transcript_ID"=rep(transcript_ID,nrow(output)), "transcript_length"=rep(geneInfo$transcipt_length, nrow(output)))
  output <- adply(output, .fun = EFFParse, .margins = 1)
  return (output)
}


SNPDataByTranscript <- function (geneInfo, strains, tidy=TRUE){
  #download a full c
  #
  transcript_ID <- as.character(geneInfo$transcript_ID)
  regionString <- as.character(geneInfo$regionString)

  VCF.out <- downloadMerge("fullVCF.vcf.gz", strains, regionString)

  if (tidy == TRUE){
    VCF.out <- vcfR2tidy(VCF.out, single_frame = TRUE, info_fields = c("AC", "EFF"), format_fields = ("GT"))
    VCF.out$dat <- VCF.out$dat[!(is.na(VCF.out$dat$gt_GT)), ]
  }

  return (VCF.out)
}




#' Create a table of nucleotide diversity data
#'
#' @param geneInfo
#' @param strains
#'
#' @return
#' @export
#'
#' @examples
polymorphTable <- function (geneInfo, strains) {
    #create a table by loading data for each of a set of genes, and counting # of polymorphisms in each category
    # geneInfo table or DF must contain "transcript_ID" and "regionStrings" columns see getGeneInfo()
  tableData <- matrix(nrow=length(geneInfo$transcript_ID), ncol=length(effects) + 6)
  rownames(tableData) <- geneInfo$transcript_ID
  colnames(tableData) <- c(effects, "coding_total", "Pi_coding", "Pi_non_syn", "Pi_syn", "Pi_NS_Ratio", "Pi_transcript")
  tableData <- data.frame(tableData)
    #iterate through each cell of the table filling in counts of variant types, load new data for each column/transcript_ID
  for (i in 1:length(geneInfo$transcript_ID)) {
    data <- loadData(geneInfo[i, ], strains)
    for (j in 1:length(effects)) {
      tableData[i, j] <- sum(data[(data$Effect == effects[j]) %in% TRUE, "Accession_count"])
    }
    tableData$coding_total[i] <- tableData$synonymous_variant[i] + tableData$missense_variant[i]
    tableData[i, (j + 2):(j + 6)] <- Nucleotide_diversity(data)
  }
  tableData <- rbind(tableData, "MEANS"=colMeans(tableData))
  return (tableData)
}


#' Get gene information
#'
#' @description Get the start and end position of each transcript of a set of
#'  genes using `bioMart` package to query the TAIR10 database via plants.ensembl.org
#'
#' @param genes
#' @param firstOnly
#' @param inputType
#' @param useCache
#'
#' @return
#' @export
#'
#' @examples
getGeneInfo <- function (genes, firstOnly=TRUE, inputType="tair_locus", useCache=TRUE) {
  # input "genes" should be a character vector of tair IDs
  # use bioMart to find the start and end position and transcript IDs for a given set of genes
  # outputs a table containing "transcript_ID" and "regionString" columns required for other fuctions in this code

  retrievedInfo <- NULL
  genes2 <- genes
  if (useCache == TRUE){
    geneInfoCache <- read.table(file="geneInfoCache.txt", header=TRUE)
    retrievedInfo <- geneInfoCache[geneInfoCache$tair_locus %in% genes, ]
    genes2 <- genes[!(genes %in% geneInfoCache$tair_locus)] #remove genes present in the cache from genes list
    print("new genes:")
    print(genes2)   # list new genes, not found in cache
  }

  output <- NULL
  if (length(genes2) > 0){
    tair10 <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")
    output <- getBM(attributes=c("tair_locus", "tair_symbol","ensembl_transcript_id", "chromosome_name", "start_position",
                                  "end_position", "strand", "transcript_start", "transcript_end"
                                  ), filters=inputType,
                                  values=genes2, mart=tair10)
    # create a list of strings encoding the chromosome and start and end position of all transcript IDs to be analyzed
    output$regionString <- as.character(alply(output, .fun=makeRegionString, .margins=1, .expand=FALSE))
    names(output)[names(output) == "ensembl_transcript_id"] <- "transcript_ID"
    output$transcipt_length <- abs(output$transcript_end - output$transcript_start)

  }
  if (useCache == TRUE) {
    # append cache
    geneInfoCache <- unique(rbind(geneInfoCache, output))
    write.table(geneInfoCache, file="geneInfoCache.txt", row.names=FALSE)
  }


  output <- rbind(retrievedInfo, output)

  if (firstOnly == TRUE) {
    # if firstOnly is TRUE, only return transcript IDs containing ".1"
    output <- output[(grepl(".1", output$transcript_ID, fixed=TRUE)), ]
  }
  return (output)
}

#' Nucleotide diversity kernel calculation
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
NucDiv_Kernel <- function (data) {
  Nucleotide_Diversity <- (sum(data$Accession_count)**2 - sum(data$Accession_count**2)) / sum(data$Accession_count)**2
  Unique_Variants <- nrow(data)
  Readable_Sequences <- sum(data$Accession_count)
  Largest_set <- max(data$Accession_count)
  Remainder <- Readable_Sequences - Largest_set
  codonNumber <- NA
  if (TRUE %in% (substring(data$Amino_Acid_Change, 1, 1) == "p")) {
      AAChangeStr <- str_split(data$Amino_Acid_Change, pattern="/", simplify=TRUE)[, 1]
      codonNumber <- as.numeric(gsub(".*?([0-9]+).*", "\\1", AAChangeStr))
      codonNumber <- codonNumber[!is.na(codonNumber)][1]
  }
  output <- data.frame(codonNumber, Nucleotide_Diversity, Unique_Variants, Readable_Sequences, Largest_set, Remainder)
  return (output)
}

#' Calculate nucleotide diversity
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
NucDiv <- function (SNPData) {
  return (ddply(SNPData, .(names), .fun=NucDiv_Kernel))
}

#' Filter nucleotide diversity kernel
#'
#' @param data
#' @param col_name
#' @param value
#'
#' @return
#' @export
#'
#' @examples
filtR_Kernel <- function (data,col_name, value) {
  if(value %in% data[, colnames(data) == col_name]) {
    return (data)
  }
  return (NULL)
}

#' A filter function using ddply
#'
#' @param data
#' @param split_var
#' @param col_name
#' @param value
#'
#' @return
#' @export
#'
#' @examples
filtR <- function(data, split_var, col_name, value){
  #this function keeps all rows with the same split_var value if any row contains 'value' in the col_name column
  #example filtR(SNPData, "names","Effect", "missense_variant" )
  return (ddply(data, split_var, .fun=filtR_Kernel, col_name=col_name, value=value))
}

#' Calculate nucleotide diversity
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
Nucleotide_diversity <- function (SNPData) {

  missense_loci <- filtR(SNPData,split_var="names", col_name="Effect", value="missense_variant")
  Pi_non_syn <- sum(NucDiv(missense_loci)$Nucleotide_Diversity) / (as.numeric(as.character(SNPData$Amino_Acid_Length[1])) * 3)

  syn_loci <- filtR(SNPData, split_var="names", col_name="Effect", value="synonymous_variant")
  Pi_syn <- sum(NucDiv(syn_loci)$Nucleotide_Diversity) / (as.numeric(as.character(SNPData$Amino_Acid_Length[1])) * 3)

  coding_loci <- unique(rbind(syn_loci,missense_loci))
  Pi_coding <- sum(NucDiv(coding_loci)$Nucleotide_Diversity) / (as.numeric(as.character(SNPData$Amino_Acid_Length[1])) * 3)

  Pi_transcript <- sum(NucDiv(SNPData)$Nucleotide_Diversity) / (as.numeric(as.character(SNPData$transcript_length[1])))

  Pi_NS_Ratio <- Pi_non_syn / Pi_syn

  return (data.frame(Pi_coding, Pi_non_syn, Pi_syn, Pi_NS_Ratio, Pi_transcript))
}


#' Simplify the SNP table
#' @description By cutting out all columns except "names", "Genotype",
#' "Codon_Change", "Amino_Acid_Change", "Effect", and "Accession_count"
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
simplifySNP <- function (SNPData) {
  return (SNPData[, c("names",
                      "Genotype",
                      "Codon_Change",
                      "Amino_Acid_Change",
                      "Effect",
                      "Accession_count")])
}

#' Add codon numbering to the nucleotide diversity kernel
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
codonNumberKernel <- function (SNPData) {
  # add the codon number to a dataframe containing a "Amino_Acid_Change"
  # use with ddply()
  changeStr <- SNPData$Amino_Acid_Change[grepl( "p.", SNPData$Amino_Acid_Change)][1]
  codonNumber <- str_extract_all(str_extract_all(changeStr, "p.[A-z]{3}[0-9]*")[[1]], "[0-9]+")[[1]]
  rows <- nrow(SNPData)
  return (cbind(SNPData, "Codon_Number"=codonNumber))
}

#' Ploting function for the nucleotide diversity
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
plotPi <- function (SNPData) {
  AALength <- as.numeric(as.character(SNPData$Amino_Acid_Length[1]))
  plotData <- data.frame("locus"=1:AALength, "synonymous_diversity"=0,  "missense_diversity"=0)

  syn_loci <- filtR(SNPData,split_var="names",col_name="Effect",value="synonymous_variant")
  syn_diversity <- NucDiv(syn_loci)

  missense_loci <- filtR(SNPData,split_var="names",col_name="Effect",value="missense_variant")
  missense_diversity <- NucDiv(missense_loci)

  for (i in 1:nrow(syn_diversity)) {
    locus <- syn_diversity$codonNumber[i]
    plotData$synonymous_diversity[locus] <- plotData$synonymous_diversity[locus] + syn_diversity$Nucleotide_Diversity[i]
  }

  for (i in 1:nrow(missense_diversity)) {
    locus <- missense_diversity$codonNumber[i]
    plotData$missense_diversity[locus] <- plotData$missense_diversity[locus] + missense_diversity$Nucleotide_Diversity[i]
  }

  plotData$synonymous_diversity <- plotData$synonymous_diversity / 3
  plotData$missense_diversity <- plotData$missense_diversity / 3

  plotDataMelt <-  melt(plotData, id.vars="locus")
  plot <- ggplot(plotDataMelt, aes(x=locus,y=value, colour=variable)) +
    geom_point() +
    scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1),limits=c(0.0001, 1)) +
    scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
    ylab("nucleotide diversity, log scale")

  print(plot)
  return (plot)
}


