plotCodingDiv <- function(uniqueCodingVars){

  effectClasses <- readRDS(system.file("data", "effect_classes.rds", package="r1001genomes"))

  classColors <- data.frame("color" = brewer.pal(n = 5, name = "RdYlBu")[5:1],
                            "Class" = c("Synonymous", "Non_Coding", "Splice",
                                        "Missense", "Nonsense"),
                            "labels" = c("Synonymous", "Non_Coding", "Splice",
                                         "Missense", "Nonsense"), stringsAsFactors=FALSE)

  uniqueCodingVars <- dplyr::left_join(uniqueCodingVars, effectClasses, by = "Effect")
  uniqueCodingVars <- dplyr::left_join(uniqueCodingVars, classColors, by = "Class")
  classes <- classColors$Class %in% unique(uniqueCodingVars$Class)
  #plot the diversity
  plot <- ggplot(uniqueCodingVars, aes(x=Codon_Number,y=Diversity, colour=color, shape = Effect)) +
    geom_point(size = 4, alpha = .9, position = position_jitter(height = 0.2)) +
    scale_y_log10(breaks=c(0.001, 0.01, 0.1),limits=c(0.001, 1)) +
    #scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
    ylab("nucleotide diversity, log scale") + theme_few(base_size = 18) +
    scale_color_identity("Class", breaks = classColors$color[classes],
                         labels = classColors$labels[classes],
                         guide = "legend")
  return(plot)
}
