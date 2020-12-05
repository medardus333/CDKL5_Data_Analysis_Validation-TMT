## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5 substrates validation
##
## Script name: Script_File_S5_Wrapper_Validation.R
##
## Purpose of script: Wrapper function to call data anlysis scripts used for
## statistical analysis of TMT XIC data
## 
## Instructions:
## 1) Place "evidence.txt" and "peptides.txt" into ./Data/
## 2) Install missing libraries, version numbers used during original analysis
##    can be found in "Text_File_S2_Session_Info_Validation.txt"
## 3) Run this script
##
## Author: Florian Weiland
##
## Date Created: 2020-05-07
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Check for missing libraries

libs <- c(
  "vsn",
  "ggplot2",
  "reshape2",
  "wesanderson",
  "extrafont",
  "matrixStats"
)

libs.install <- which(libs %in% rownames(installed.packages()) == FALSE)

if (length(rlibs.install) >= 1) {

  message(paste0("\nPlease install ", libs[rlibs.install]))

  message("Version numbers used is stated in \"Text_File_S2_Session_Info_Validation.txt\"")

}

## Directories

if (dir.exists("./Output_R/") == FALSE) {

dir.create("./Output_R/")

}

## Source scripts

## ggplot throws warnings, but they can be ignored. Caused by annotation of boxplots.
## The warning states: is.na() applied to non-(list or vector) of type 'expression'

source("Script_File_S6_Data_analysis_ELOA_Validation.R")
source("Script_File_S7_Data_analysis_TTDN1_Validation.R")
source("Script_File_S8_Data_analysis_EP400_Validation.R")

sink("Text_File_S2_Session_Info_Validation.txt")
sessionInfo()
sink()


