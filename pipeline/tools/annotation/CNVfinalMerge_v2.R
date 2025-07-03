#### FINAL CNV FILE ####

library(dplyr)
library(optparse)


option_list=list(
  make_option(c('-d','--mix'),type="character", help="mix"),
  make_option(c('-a','--annot'),type="character", help="annot file"))

opt <- parse_args(OptionParser(option_list = option_list))


mixer=opt$d
annotsv=opt$a



df_mixer <- read.csv( mixer, fill=TRUE, sep= '\t') # .combined.txt
print("MIXER... READ")

df_annot <- read.csv(annotsv, fill=TRUE, sep= '\t') # .combinedAnnotated.tsv
print("ANNOTSV OUTPUT... READ")

# df_mixer_to_annot <- df_mixer[,c("SAMPLE", "HOMO_HETEROZYGOUS", "CHR", "START", "STOP", "NUM_OF_PROGRAMS" )]

drop_columns <- c("AnnotSV.ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TADcoordinates", "ENCODEexperiments",  "GHid_not_elite", "ACMG" )

df_annot_corrected <-  df_annot[, ! names(df_annot) %in% drop_columns, drop = FALSE]

df_intermediate <- left_join( df_mixer, df_annot_corrected, 
	by=c("CHR"="SV_chrom","START"="SV_start", "STOP"="SV_end", "CNV_TYPE"="SV_type" )) %>% distinct()

unique(df_mixer$CHR)
unique(df_annot_corrected$SV_chrom)

colnames(df_intermediate["SV_length"])
which(colnames(df_intermediate)=="Annotation_mode")
columnas<-colnames(df_intermediate)
columnas[1:30]

df_final <- df_intermediate[, c(135:138, 35 , 1:30, 31:135)] # For the new version 

split<-df_final %>% filter(Annotation_mode=="split")
write.table(split, "TFM_Patri_200625.final_split.txt", sep="\t", row.names=FALSE,  quote= FALSE )

full<-df_final %>% filter(Annotation_mode=="full")
write.table(full, "TFM_Patri_200625.final_full.txt", sep="\t", row.names=FALSE,  quote= FALSE )



