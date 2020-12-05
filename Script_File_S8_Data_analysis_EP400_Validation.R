## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5 substrates validation
##
## Script name: Script_File_S8_Data_analysis_EP400_Validation.R
##
## Purpose of script: Clean up MaxQuant data and analyse
##
## Author: Florian Weiland
##
## Date Created: 2020-05-07
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Libraries ----

library(vsn)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(extrafont)
library(matrixStats)

## Folders and filenames ----

data.dir   <- "./Data/"
output.dir <- "./Output_R/"

meanssd.file   <- "MeansSDPlot_EP400.png"
qq.file.1      <- "QQ_plot_EP400_1.png"
vsn.tmt.file.1 <- "VSN_intensities_TMT_EP400_1.png" 
box.file.1     <- "Boxplots_EP400_1.png"
stats.file     <- "Statistics_EP400.csv"
data.file.1    <- "Data_EP400_1.csv"

## Set theme for plots

theme_set(theme_classic(base_family = "Arial"))
pal <- wes_palette("Darjeeling1")

## Read in data ----

message("Reading data for EP400")

data.mq <- read.csv(
  paste0(data.dir, "evidence.txt"),
  sep = "\t",
  stringsAsFactors = FALSE
)

peptides.df <- read.csv(
  paste0(data.dir, "peptides.txt"),
  sep = "\t",
  stringsAsFactors = FALSE
)

peptides.df <- peptides.df[which(peptides.df$Gene.names == "EP400"), ]

## Data clean-up ----
## Remove Reverse hits, potential contaminants and PEP >= 0.05

message("Clean up data")

rem.rev <- which(data.mq$Reverse == "+")
rem.cont <- which(data.mq$Potential.contaminant == "+")
rem.pep <- which(data.mq$PEP >= 0.05)

data <- data.mq[-c(rem.rev, rem.cont, rem.pep), ]

## Last intensity column is TMT131, not used in these experiments.

intensity.cols <- grep(
  "Reporter.intensity.corrected.[1-9]$", 
  colnames(data)
)

## Extract EP400 data

data <- data[grep("EP400\\.1", data$Experiment), ]
data <- data[which(data$Gene.names == "EP400"), ]

## Extract corrected reporter ion intensities and normalize using vsn ----

data.int <- as.matrix(data[, intensity.cols])

##Remove Zeroes from intensities and replace with NA

data.int[data.int == 0] <- NA

## VSN ----

message("VSN transform")

data.model <- vsn2(data.int)
data.vsn = predict(data.model, newdata = data.int)
colnames(data.vsn) <- paste0("VSN.", colnames(data.vsn))

##Check if standard deviation approx. constant

meanSdPlot(data.vsn)
ggsave(
  paste0(output.dir, meanssd.file),
  dpi = 1200
)

data.1 <- cbind.data.frame(data, data.vsn)

## QQ-Plots ----
## EP400

message("QQ plots EP400")

data.1.int <- data.1[, grep("^VSN", colnames(data.1))]

data.qq <- melt(data.1.int)
data.qq <- data.qq[-which(is.na(data.qq$value) == TRUE), ]
colnames(data.qq) <- c("TMT.label", "VSN.Intensity")

ggplot(data = data.qq, aes(sample = VSN.Intensity)) +
  stat_qq() + 
  stat_qq_line() +
  ggtitle("Q-Q Plot VSN transformed data")

ggsave(
  paste0(output.dir, qq.file.1),
  dpi = 1200
)

## VSN intensities ----
## EP400

message("VSN transformed TMT intensities plot EP400")

data.qq$Group <- data.qq$TMT.label
data.qq$Group <- sub("VSN.Reporter.intensity.corrected.", "", data.qq$TMT.label)
data.qq$Group <- sub("1|2|3", "FLAG+H2O2", data.qq$Group)
data.qq$Group <- sub("4|5|6", "WT+H2O2", data.qq$Group)
data.qq$Group <- sub("7|8|9", "K42R+H2O2", data.qq$Group)
data.qq$Group <- factor(data.qq$Group, levels = c("FLAG+H2O2", "WT+H2O2", "K42R+H2O2"))

data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.1", "126", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.2", "127N", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.3", "127C", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.4", "128N", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.5", "128C", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.6", "129N", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.7", "129C", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.8", "130N", data.qq$TMT.label)
data.qq$TMT.label <- sub("VSN.Reporter.intensity.corrected.9", "130C", data.qq$TMT.label)

data.qq$TMT.label <- factor(
  data.qq$TMT.label,
  levels = c(
    "126",
    "127N",
    "127C",
    "128N",
    "128C",
    "129N",
    "129C",
    "130N",
    "130C"
  ))

ggplot(
  data.qq,
  mapping = aes(x = TMT.label, y = VSN.Intensity)
) +
  geom_boxplot(
    lwd = 0.8,
    fill = c(rep(pal[1],3), rep(pal[5],3), rep(pal[3],3))
  ) +
  theme(
    axis.text.x = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(size = 14, family = "Arial"),
    axis.title.y = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(family = "Roboto Mono", size = 14)
  ) +
  labs(
    x = "TMT channel",
    y = "VSN Intensity",
    title = "VSN intensities of peptides"
  ) +
  annotate(
    geom = "text",
    x = 2 + 3 * (0:2),
    y = 29,
    label = c(
      expression(FLAG+H[2]*O[2]),
      expression(WT+H[2]*O[2]),
      expression(K^42*R+H[2]*O[2])
    ),
    size = 6
  )

ggsave(
  paste0(output.dir, vsn.tmt.file.1),
  dpi = 1200
)

## Extract peptides of interest ----
## Need to go by amino acid number of S of interest to account for missed cleavages

message("Prepare data")

## Position of amino acid of interest in amino acid sequence. These positions are isoform specific!
## MaxqQuant reports amino acid numbers for isoform 2. Serine 729 in isoform 1 is serine 693 in isoform 2.
serine.pos <- 693

peptide.spans <- data.frame(
  Sequence = peptides.df$Sequence,
  Start.Position = peptides.df$Start.position,
  End.Position = peptides.df$End.position
)

for (cur.row in 1:nrow(peptide.spans)) {
  
  aa.numbers.temp <- seq(
    peptide.spans[cur.row, "Start.Position"],
    peptide.spans[cur.row, "End.Position"],
    1
  )
  
  if (serine.pos %in% aa.numbers.temp) {
    
    peptide.spans[cur.row, "Contains.Serine.of.Interest"] <- TRUE
    
  } else {
    
    peptide.spans[cur.row, "Contains.Serine.of.Interest"] <- FALSE
    
  }
  
}

## Statistical testing ----

## Extract peptides containing serine of interest

pep.soi <- peptide.spans[which(peptide.spans$Contains.Serine.of.Interest == TRUE), ]

## Sum intensity data of all peptides containing Serine of interest, build dataframes for t-test

pep.soi.1.rows <- which(data.1$Sequence %in% pep.soi$Sequence)

## Build dataframe for t-test

data.1.quant <- as.data.frame(
  cbind(
    as.character(data.1[pep.soi.1.rows, "Modified.sequence"]),
    as.character(data.1[pep.soi.1.rows, "Modifications"]),
    as.character(data.1[pep.soi.1.rows, "Phospho..STY..Probabilities"]),
    data.1.int[pep.soi.1.rows, ]
  ),
  stringsAsFactors = FALSE
)

colnames(data.1.quant) <- c(
  "Modified.Sequence",
  "Modifications",
  "Phospho.Probabilities",
  "FLAG+CDKL5_WT+H2O2.1",
  "FLAG+CDKL5_WT+H2O2.2",
  "FLAG+CDKL5_WT+H2O2.3",
  "EP400+CDKL5_WT+H2O2.1",
  "EP400+CDKL5_WT+H2O2.2",
  "EP400+CDKL5_WT+H2O2.3",
  "EP400+CDKL5_K42R+H2O2.1",
  "EP400+CDKL5_K42R+H2O2.2",
  "EP400+CDKL5_K42R+H2O2.3"
)

p.pep.soi.1.rows <- which(data.1.quant$Modifications == "Phospho (STY)")
data.1.quant <- data.1.quant[p.pep.soi.1.rows , ]

## Average intensity values (median) ----

message("Average intensities")

data.1.av <- as.matrix(data.1.quant[, 4:12])
data.1.av <- colMedians(data.1.av, na.rm = TRUE)

names(data.1.av) <- c(
  "FLAG+CDKL5_WT+H2O2.1",
  "FLAG+CDKL5_WT+H2O2.2",
  "FLAG+CDKL5_WT+H2O2.3",
  "EP400+CDKL5_WT+H2O2.1",
  "EP400+CDKL5_WT+H2O2.2",
  "EP400+CDKL5_WT+H2O2.3",
  "EP400+CDKL5_K42R+H2O2.1",
  "EP400+CDKL5_K42R+H2O2.2",
  "EP400+CDKL5_K42R+H2O2.3"
)

## t-test ----
## var.equal = TRUE set due to vsn normalization

message("Statistical testing")

flag.1   <- grep("FLAG\\+CDKL5_WT\\+H2O2", names(data.1.av))
wt.1   <- grep("EP400\\+CDKL5_WT\\+H2O2", names(data.1.av))
k42r.1 <- grep("EP400\\+CDKL5_K42R\\+H2O2", names(data.1.av))

WT.K42R.1.p <- t.test(
  as.numeric(data.1.av[wt.1]),
  as.numeric(data.1.av[k42r.1]),
  alternative = "two.sided",
  var.equal = TRUE,
  na.action = na.omit
)$p.value

WT.FLAG.1.p <- t.test(
  as.numeric(data.1.av[wt.1]),
  as.numeric(data.1.av[flag.1]),
  alternative = "two.sided",
  var.equal = TRUE,
  na.action = na.omit
)$p.value

## Multiple testing, "standard" p-value threshold adjusted by Bonferroni Correction. 
## Two t-tests: α = 0.05/2 = 0.025
## All significant

## Visualization ----

message("Boxplots EP400")

peptide.1 <- as.character(unique(data.1.quant$Modified.Sequence))
peptide.1 <- substring(peptide.1, 2, nchar(peptide.1) - 1)

data.1.gg <- as.data.frame(data.1.av)
data.1.gg$Group <- rownames(data.1.gg)
data.1.gg$Group <- sub(".\\d{1}$", "", data.1.gg$Group)

data.1.gg$Group <- factor(
  data.1.gg$Group,
  levels = c("FLAG+CDKL5_WT+H2O2","EP400+CDKL5_WT+H2O2", "EP400+CDKL5_K42R+H2O2")
)

colnames(data.1.gg)[1] <- "VSN.Intensity"

ggplot(
  data = data.1.gg, 
  mapping = aes(x = Group, y = VSN.Intensity)
) +
  geom_boxplot(fill = c(pal[1], pal[5], pal[3]), size = 1, fatten = 1.5) + 
  geom_point(size = 3) +
  labs(title = paste("EP400 - phospho-Serine 729(693) -", peptide.1)) +
  ylab("VSN Intensity") +
  theme(
    axis.text.x = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(family = "Roboto Mono", size = 14),
  ) +
  scale_x_discrete(
    labels = c(
      expression(phantom(0)^phantom(0)*FLAG+H[2]*O[2]),
      expression(phantom(0)^phantom(0)*WT+H[2]*O[2]),
      expression(K^42*R~+H[2]*O[2])
    )
  )

ggsave(
  paste0(output.dir, box.file.1),
  dpi = 1200
)

## Assemble statistics dataframe ----
## Fold change

message("Prepare statistics dataframe")

fc.wt.k42r.1 <- round(
  2^(abs(mean(as.numeric(data.1.av[wt.1]), na.rm = TRUE) -
           mean(as.numeric(data.1.av[k42r.1]), na.rm = TRUE))),
  2
)

fc.wt.flag.1 <- round(
  2^(abs(mean(as.numeric(data.1.av[wt.1]), na.rm = TRUE) -
           mean(as.numeric(data.1.av[flag.1]), na.rm = TRUE))),
  2
)

stats.df <- data.frame(
  Sample = "EP400",
  Sequence = peptide.1,
  p.value.WT.FLAG = WT.FLAG.1.p,
  p.value.WT.K42R = WT.K42R.1.p,
  Fold.change.WT.FLAG = fc.wt.flag.1,
  Fold.change.WT.K42.R = fc.wt.k42r.1
)

message("Write statistics dataframe")

write.csv(
  stats.df, 
  paste0(output.dir, stats.file),
  row.names = FALSE
)

write.csv(
  data.1.quant,
  paste0(output.dir, data.file.1),
  row.names = FALSE
)
