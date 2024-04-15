## Allele frequency spectra for artifical hybrid.

# Load required package for plots
library(ggplot2)

# Read in files.
arenosa <- read.csv('newtxts/arenosa_632.txt',sep = '\t')
lyrata <- read.csv('newtxts/lyrata_272_with_some_hybrids.txt',sep = '\t')

# Select rows present in both samples.
merged <- merge(arenosa,lyrata, by="POS")

# Make mean allele frequency column.
merged$mean_AF <- (merged$AF.x + merged$AF.y)/2

# Select rows over 0.1 allele frequency.
merged$meanover1 <- ifelse(merged$mean_AF > 0.1,merged$mean_AF, NA)

# Plot histogram.
ggplot(merged,aes(x=meanover1))+
  geom_histogram(bins=100,alpha=0.5,color='cyan4',fill='cyan4')+
  xlab('Allele frequency')+
  ylab('Count')+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color='black'),
    plot.title = element_text(size=10)
  )+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  scale_x_continuous(expand = expansion(mult = c(0,0.05)))+
  ggtitle('Number of bins = 100')


## Allele frequency spectra for selected cohort.

# Read in filtered VCF of selected cohort.
vcf <- read.csv('freqedit3.csv')

# Select hybrid individuals.
vcf$meanhybrids <- rowMeans(vcf[c("OCH.03tl","OCH.04tl","OCH.06tl",
                                  "OCH.07tl","OCH.08tl","HAB.01tl",
                                  "HAB.02tl", "HAB.03tl","FRE.03tl",
                                  "FRE.04tl","FRE.07tl")],na.rm = TRUE)

# Select rows over 0.1 allele frequency.

vcf$hover1 <- ifelse(vcf$meanhybrids>0.1,vcf$meanhybrids, NA)

# Plot histogram.
ggplot(vcf,aes(x=hover1))+
  geom_histogram(bins=40,alpha=0.5,color='cyan4',fill='cyan4')+
  xlab('Allele frequency')+
  ylab('Count')+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color='black'),
    plot.title = element_text(size=10)
  )+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  scale_x_continuous(expand = expansion(mult = c(0,0.05)))+
  ggtitle('Number of bins = 40')
