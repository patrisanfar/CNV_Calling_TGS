rm(list=ls())
library(optparse)


#**************#
#   Arguments  #
#**************#
option_list=list(
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Project name."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"))

opt<-parse_args(OptionParser(option_list=option_list))

output_dir <- opt$outputdir
bedFile <- opt$bed
projectname <- opt$name


output_dir <- "/mnt/beegfs/home/pasanfar/results/mixer/"
bedFile <- "/mnt/beegfs/home/pasanfar/panel/P08_AREPA_capture_targets_v2.bed"
projectname <- "TFMPatri"
setwd(output_dir)

#=====================#
#                     #
#   #=============#   #
#   # IMPORT DATA #   # 
#   #=============#   #
#                     #
#=====================#

#*****************#
# Import Bed File #
#*****************#
bedFile_data <- read.delim(bedFile, header = FALSE, stringsAsFactors = FALSE)

if(length(colnames(bedFile_data))<4){
  cat("ERROR: bed file without NAME field")
  quit(status = 1)} else if (length(colnames(bedFile_data))>4){
  bedFile_data = bedFile_data[,1:4]}
  
colnames(bedFile_data) <- c("CHR", "START", "STOP", "GENE_NAME")
bedFile_data$AVG_POS <- (bedFile_data$START + bedFile_data$STOP) / 2 # Create a column with the medium point of the exons
bedFile_data$CHR <- gsub("chr", "", bedFile_data$CHR)

# Create exon names. It changes the order in case there are overlapping genes.
bedFile_data <- bedFile_data[order(bedFile_data$CHR, bedFile_data$GENE_NAME, bedFile_data$START),]
j <- 1
for (i in 1:nrow(bedFile_data)){
  bedFile_data$TARGET[i] <- paste(bedFile_data$GENE_NAME[i], j, sep = "_")
  j <- j + 1
  if (bedFile_data$GENE_NAME[i] != bedFile_data$GENE_NAME[i+1] & i != nrow(bedFile_data)){j <- 1}
}

bedFile_data <- bedFile_data[order(bedFile_data$CHR, bedFile_data$START),]

# Create total_data. It contains the output of the programs
total_data <- data.frame(stringsAsFactors = FALSE)



### Import Decon data
decon_output <- paste(output_dir, 'decon', '.results.txt', sep = '')
tryCatch(
  {
    decon_data <- read.delim(decon_output, stringsAsFactors = FALSE,sep=",") # Get the information from decon output
    decon_data$PROGRAM <- "DC"
    df_parcial <- data.frame(decon_data$Sample,
                             gsub("chr", "", decon_data$Chromosome),
                             decon_data$Start,
                             decon_data$End,
                             tolower(decon_data$CNV.type),  # para que sea 'duplication' o 'deletion' en minúsculas
                             decon_data$PROGRAM ,
                             decon_data$Reads.ratio,
                             NA, NA, NA, NA,
                             decon_data$BF,
                             NA, NA, NA, NA, NA, NA,
                             stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM","DC_ratio", "CN_ratio", "C2_ratio", "PM_ratio","CL_ratio", "DC_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC","CL_potentialAF","CL_loglikelihood","CL_qvalue")
    
    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no decon output"), 
  warning=function(e) print("There is no decon output")
)
  
### Import Convading data 
convading_output <- paste(output_dir, 'CoNVaDING', '.results.txt', sep = '')
tryCatch(
  {
    convading_data <- read.delim(convading_output, stringsAsFactors = FALSE)
    convading_data$PROGRAM <- "CN"
    
    df_parcial <- data.frame(convading_data$SAMPLE_NAME,
                             convading_data$CHR,
                             convading_data$START,
                             convading_data$STOP,
                             tolower(convading_data$ABBERATION),
                             convading_data$PROGRAM,
                             NA,
                             convading_data$AUTO_RATIO,
                             NA, NA, NA,
                             NA,
                             convading_data$AUTO_ZSCORE,
                             convading_data[[grep("shapiro", colnames(convading_data), ignore.case = TRUE)]],
                             NA, NA, NA, NA, 
                             stringsAsFactors = FALSE)
    
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM","DC_ratio", "CN_ratio","C2_ratio", "PM_ratio","CL_ratio", "DC_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC","CL_potentialAF","CL_loglikelihood","CL_qvalue")
    
    total_data <- rbind(total_data, df_parcial)
  },
  error=function(e) print("There is no CoNVaDING output"),
  warning=function(e) print("There is no CoNVaDING output")
)

### Import Codex2 data
codex2_output <- paste(output_dir, 'CODEX2', '.results.txt', sep = '')
tryCatch(
  {
    codex2_data <- read.delim(codex2_output, stringsAsFactors = FALSE,sep = ",")
    codex2_data$chr <- gsub("chr", "", codex2_data$chr)
    codex2_data$cnv <- toupper(codex2_data$cnv)
    codex2_data$PROGRAM <- "C2"
    
    df_parcial <- data.frame(codex2_data$sample_name, codex2_data$chr, codex2_data$st_bp, codex2_data$ed_bp,
                             codex2_data$cnv, codex2_data$PROGRAM, NA, NA, codex2_data$copy_no/2, NA, NA, NA, NA, NA, codex2_data$mBIC, NA, NA, NA, stringsAsFactors = FALSE) # mBIC: modified bayesian information criterion
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM","DC_ratio", "CN_ratio", "C2_ratio", "PM_ratio","CL_ratio", "DC_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC","CL_potentialAF","CL_loglikelihood","CL_qvalue")
    
    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no CODEX2 output"), 
  warning=function(e) print("There is no CODEX2 output")
)

### Import Panelcn.MOPS  
panelMops_output <- paste(output_dir, 'panelcn.MOPS', '.results.txt', sep = '')
tryCatch(
  {
    panelMops_data <- read.delim(panelMops_output, stringsAsFactors = FALSE)
    panelMops_data$PROGRAM <- "PM"
    
    df_parcial <- data.frame(panelMops_data$Sample, panelMops_data$Chr, panelMops_data$Start, panelMops_data$End,
                             panelMops_data$CN, panelMops_data$PROGRAM, NA, NA, NA, panelMops_data$RC.norm/panelMops_data$medRC.norm, NA, NA, NA, NA, NA, NA, NA, NA, stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM","DC_ratio", "CN_ratio", "C2_ratio", "PM_ratio","CL_ratio", "DC_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC","CL_potentialAF","CL_loglikelihood","CL_qvalue")
    
    total_data <- rbind(total_data, df_parcial)

  },
  error=function(e) print("There is no Panelcn.MOPS output"), 
  warning=function(e) print("There is no Panelcn.MOPS output")
)

### Import ClinCNV  
ClinCNV_output <- paste(output_dir, 'ClinCNV', '.results.txt', sep = '')
tryCatch(
  {
    ClinCNV_data <- read.delim(ClinCNV_output, stringsAsFactors = FALSE,sep=",")
    ClinCNV_data$X.chr<-gsub("chr", "", ClinCNV_data$X.chr)
    ClinCNV_data$PROGRAM <- "CL"
    ClinCNV_data$CNV_type<-ifelse(ClinCNV_data$CN_change>=2,"DUP","DEL") # cuidao ahi
    df_parcial <- data.frame(ClinCNV_data$filename, ClinCNV_data$X.chr, ClinCNV_data$start, ClinCNV_data$end,
                             ClinCNV_data$CNV_type, ClinCNV_data$PROGRAM, NA, NA, NA, NA, ClinCNV_data$CN_change/2, NA, NA, NA, NA,ClinCNV_data$potential_AF,ClinCNV_data$loglikelihood,ClinCNV_data$qvalue, stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM","DC_ratio", "CN_ratio", "C2_ratio", "PM_ratio","CL_ratio", "DC_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC","CL_potentialAF","CL_loglikelihood","CL_qvalue")
    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no ClinCNV output"), 
  warning=function(e) print("There is no ClinCNV output")
)

total_data <- data.frame(lapply(total_data, function(x) if (is.factor(x)) as.character(x) else {x}), stringsAsFactors=FALSE)

total_data$SAMPLE_NAME <- gsub("\\.recal","", total_data$SAMPLE_NAME, perl = TRUE)
total_data$SAMPLE_NAME <- gsub("\\_cnvs.tsv","", total_data$SAMPLE_NAME, perl = TRUE)
total_data$SAMPLE_NAME <- gsub("_.*","", total_data$SAMPLE_NAME, perl = TRUE)
total_data$SAMPLE_NAME <- gsub("\\.","\\-", total_data$SAMPLE_NAME, perl = TRUE)

total_data$CNV_TYPE <- toupper(total_data$CNV_TYPE)
total_data$CNV_TYPE <- gsub("DELETION", "DEL", total_data$CNV_TYPE)
total_data$CNV_TYPE <- gsub("DUPLICATION", "DUP", total_data$CNV_TYPE)

#===================#
#                   #
#   #===========#   #
#   # PROCESING #   # 
#   #===========#   #
#                   #
#===================#

#*****************#
# Exon assignment #
#*****************#
# basicamente va linea por linea del bed y añade las cnvs que se han detectado en esos exones

total_exons <- do.call("rbind", apply(bedFile_data, 1, function(x){
  # para cada fila del bed ejecuta una función y va guardando en un dataframe con rbind
  
  positionRef = as.integer(x[5]) # punto medio del exon
  chrRef = x[1] # cromosoma
  exonName = x[6] # nombre del exon
  
  total_data_WE = total_data[total_data$CHR == chrRef & positionRef > total_data$START & positionRef < total_data$STOP,]
  # se filtran las líneas de total data que contienen todos los CNVs detectados por todas las herramientas, 
  # y selecciona solo aquellos CNVs que solapan con el exón actual (x) del BED.
  
  if (nrow(total_data_WE) != 0){ # si hay cnvs con esa condición
    total_data_WE$EXON <- exonName # se añade una nueva columna con el nombre del exon
  } 
  return(total_data_WE) # devuelve cnvs etiquetados con el exón correspondiente
}))

# SE HA CREADO TOTAL EXONS QUE CONTIENE LAS CNVS ASOCIADAS A EXONES

total_exons$HOMOZYGOUS <- NA
# crea una nueva columna en total exons y la inicializa con NA
total_exons$HOMOZYGOUS[grep("HOM", total_exons$CNV_TYPE)] <- "HOM" # Mark homozygous CNVs (CN and PM)
# busca filas donde el tipo de CNV contiene HOM, y en esas filas pone HOM en la columna homocygous
total_exons$CNV_TYPE <- gsub("HOM_", "", total_exons$CNV_TYPE) # Convert HOM-DUP/DEL to DUP/DEL
# limpia la columna CNV TYPE quitando HOM
total_exons$EXON_MIX <- paste(total_exons$SAMPLE_NAME, total_exons$EXON, total_exons$CNV_TYPE, sep = "_") # Create a key for comparison
# crea una clave única para cada muestra, exón y tipo de CNV, ejemplo IM167_PARK7_2_DEL
total_exons$GENE <- gsub("_[0-9]*$", "", total_exons$EXON)
# extrae el nombre del gen desde el nombre del exon: PARK7_2 -> PARK7
total_exons$GENE_MIX <- paste(total_exons$SAMPLE_NAME, total_exons$GENE, total_exons$CNV_TYPE, sep = "_")
# crea una clave única para cada muestra, gene y tipo de cnv, ejemplo IM167_PARK7_DEL

total_exons$PROGRAM_MIX <- paste(total_exons$PROGRAM, total_exons$EXON_MIX, sep = "_")
# crea una nueva columna llamada PROGRAM MIX, que une nombre herramienta con clave exon_mix
total_exons <- total_exons[!duplicated(total_exons$PROGRAM_MIX),] # Delete rows with duplicated $PROGRAM_MIX
# elimina las filas duplicadas de program mix
total_exons$EXON_POS <- as.numeric(gsub("^.*_", "", total_exons$EXON))
# extraer la posición del exón dentro del gen EXON = "PARK7_2" → "2" → EXON_POS = 2


#*****************#
# Output creation #
#*****************#

# Función que agrupa números consecutivos en intervalos, ejemplo: exones consecutivos
interval.maker <- function(numbers, intput_type){ # define dos argumentos: números y tipo de input (string o vector)
  # Transform a list of numbers to intervals
  if (intput_type == "str"){num_vector <- as.numeric(strsplit(numbers, split = ", ")[[1]])
  }else if (intput_type == "vec") {num_vector <- numbers
  }else{return(print("intput_type not detected"))}
  # se transforma la entrada a un vector de números, si el input es string lo hace, si ya es vector no hace nada
  num_vector <- unique(num_vector)
  # elimina duplicados
  if (length(num_vector) <= 1){return(numbers)} # si solo hay un número lo devuelve directamente
  num_vector <- sort(num_vector) # se ordena el vector de menor a mayor para detectar consecutivos
  output_interval <- c("") # inicializa vector vacío donde se irán añadiendo los intervalos
  for (i in 2:length(num_vector)){ # el bucle recorre los números y detecta 
    if (num_vector[i] - 1 == num_vector[i - 1]){ # si el número actual es consecutivo al anterior
      if (tail(output_interval, 1) == "-"){ # verifica si estamos dentro de intervalo abierto, ya se añadió "-"
        if (i == length(num_vector)){ # si es el último número
          output_interval <- c(output_interval, num_vector[i]) # cierra el intervalo añadiendo el número final
        }
      }else{ # si no está dentro de un intervalo, lo inicia
        if (i == length(num_vector)){ # si es el último número
          output_interval <- c(output_interval, num_vector[i - 1], "-", num_vector[i]) # se añade el rango completo
        }else{ # si no es el último número
          output_interval <- c(output_interval, num_vector[i - 1], "-") # se añade el número inicial y se espera a cerrarlo
        }
      }
    }else{ # si no son consecutivos
      if (i == length(num_vector)){ # si era el último número
        output_interval <- c(output_interval, num_vector[i - 1], ", ", num_vector[i]) # se añade el número anterior y este separados por coma
      }else{ # si no era el último número
        output_interval <- c(output_interval, num_vector[i - 1], ", ") # se añade el número anterior seguido por comas y ver que pasa con los siguientes
      }
    }
  }
  output_interval <- paste(output_interval, collapse = "") # une todo el contenido en un solo string
  if (intput_type == "vec") {output_interval <- strsplit(output_interval, split = ", ")[[1]]}
  # si la entrada fue un vector numérico, devuelve la salida como vector de strings individuales
  return(output_interval)
}

total_list <- split(total_exons, total_exons$GENE_MIX)
# divide el data frame totalexon en una lista
  ## cada elemento de la lista es un gene_mix y contiene todas las filas correspondientes (exon_mix: exones afectados en ese gen)
    ## ejemplo: elemento de lista IM001_PARK7_DEL, contiene "IM001_PARK7_1_DEL" "IM001_PARK7_2_DEL" "IM001_PARK7_3_DEL" "IM001_PARK7_4_DEL"

data_out <- data.frame(stringsAsFactors = FALSE)
# se crea un data frame vacío que se rellenará con los resultados

row_idx <- 1
# inicializa el índice de la fila para ir rellenando data out
index_iter<-1
for (i in total_list){
  print((index_iter/length(total_list)) *100 )
  intervals <- interval.maker(unique(i$EXON_POS), "vec") # resume las posiciones de los exones afectados
  
  for (j in intervals){ # recorre cada intervalo
    min <- as.numeric(gsub("-.*$", "", j)) # extrae número mínimo del intervalo
    max <- as.numeric(gsub("^.*-", "", j)) # extrae número máximo del intervalo
    
    CNV = i[i$EXON_POS >= min & i$EXON_POS <= max,] # selecciona los exones dentro de ese intervalo (ej. PARK7_1, PARK7_2, PARK7_3).
    
    # se guardan los campos generales de la muestra: nombre, gen y tipo cnv (desglosar el nombre de gene_mix)
    data_out[row_idx, "SAMPLE"] <- CNV$SAMPLE_NAME[1]
    data_out[row_idx, "GENE"] <- CNV$GENE[1]
    data_out[row_idx, "CNV_TYPE"] <- CNV$CNV_TYPE[1]
    
    # selecciona exones donde DC_ratio es anómalo
    homo_exons_DC<- gsub("^.*_", "", CNV[CNV$DC_ratio < 0.1 | CNV$DC_ratio > 1.75, "EXON"]) # se queda con el número del exon
    homo_exons_DC <- paste(homo_exons_DC[!is.na(homo_exons_DC)], collapse = ",") # valores distintos de NA los pega
    if (homo_exons_DC != ""){homo_exons_DC <- paste("DC(", homo_exons_DC,") ", sep = "")} # ej: "DC(2,3) "
    
    # igual que el de antes
    homo_exons_C2<- gsub("^.*_", "", CNV[CNV$C2_ratio < 0.1 | CNV$C2_ratio > 1.75, "EXON"])
    homo_exons_C2 <- paste(homo_exons_C2[!is.na(homo_exons_C2)], collapse = ",")
    if (homo_exons_C2 != ""){homo_exons_C2 <- paste("C2(", homo_exons_C2,") ", sep = "")}
    
    # filtra cnvs marcadas como homocigotas
    homo_exons_CN<- gsub("^.*_", "", CNV[CNV$PROGRAM == "CN" & CNV$HOMOZYGOUS == "HOM", "EXON"])
    homo_exons_CN <- paste(homo_exons_CN[!is.na(homo_exons_CN)], collapse = ",")
    if (homo_exons_CN != ""){homo_exons_CN <- paste("CN(", homo_exons_CN,") ", sep = "")}
    
    # igual
    homo_exons_PM<- gsub("^.*_", "", CNV[CNV$PROGRAM == "PM" & CNV$HOMOZYGOUS == "HOM", "EXON"])
    homo_exons_PM <- paste(homo_exons_PM[!is.na(homo_exons_PM)], collapse = ",")
    if (homo_exons_PM != ""){homo_exons_PM <- paste("PM(", homo_exons_PM,") ", sep = "")}
    
    # igual
    homo_exons_CL<- gsub("^.*_", "", CNV[CNV$CL_ratio <= 0 | CNV$CL_ratio >= 2, "EXON"])
    homo_exons_CL <- paste(homo_exons_CL[!is.na(homo_exons_CL)], collapse = ",")
    if (homo_exons_CL != ""){homo_exons_CL <- paste("CL(", homo_exons_CL,") ", sep = "")}
    
    # une todos los resultados de homocigosis detectados por las herramientas
    homo_exons = paste("HOM: ", homo_exons_DC, homo_exons_CN, homo_exons_C2, homo_exons_PM, sep = "")
    if (homo_exons == "HOM: "){data_out[row_idx, "HOMO_HETEROZYGOUS"] <- "HET" # si no se encontraron marca fila como hetero
    }else{data_out[row_idx, "HOMO_HETEROZYGOUS"] <- homo_exons} # si se encontró, se escribe la columna de hoocigosis
    
    # asignan a cada CNV su cromosoma (CHR), la posición inicial (START) y final (STOP) del intervalo genómico afectado
    data_out[row_idx, "CHR"] <- CNV$CHR[1]
    data_out[row_idx, "START"] <- CNV$START[which.min(CNV$START)]
    data_out[row_idx, "STOP"] <- CNV$STOP[which.max(CNV$STOP)]
    
    data_out[row_idx, "NUM_OF_TARGETS"] <- length(unique(gsub("^.*_","",CNV$EXON))) # cuantos exones únicos hay
    data_out[row_idx, "TARGETS"] <- paste(unique(gsub("^.*_","",CNV$EXON)), collapse = ", ") # cuáles son los exones
    data_out[row_idx, "DC_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "DC",]$EXON), collapse = ", ")
    data_out[row_idx, "CN_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "CN",]$EXON), collapse = ", ")
    data_out[row_idx, "C2_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "C2",]$EXON), collapse = ", ")
    data_out[row_idx, "PM_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "PM",]$EXON), collapse = ", ")
    data_out[row_idx, "CL_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "CL",]$EXON), collapse = ", ")
   
    data_out[row_idx, "MAX_PROGRAMS_SAME_TARGETS"] <- max(as.vector(table(CNV$EXON))) # cuántas veces se ha detectado el mismo exón por varios programas
    data_out[row_idx, "MAX_PROGRAMS_ANY_TARGET"] <- length(unique(CNV$PROGRAM))

    # Obtener cuántos programas detectaron cada exón
    programs_per_exon <- table(CNV$EXON)
    
    # Cuál fue el máximo número de programas detectando el mismo exón
    max_program_count <- max(programs_per_exon)
    
    # Qué exones tienen exactamente ese número de programas
    exons_with_max_programs <- names(programs_per_exon[programs_per_exon == max_program_count])
    
    # Extraer solo el número del exón (de FCGBP_12 → 12)
    targets_with_max_programs <- gsub("^.*_", "", exons_with_max_programs)
    
    # Convertir a intervalo legible
    if (length(targets_with_max_programs) > 0) {
      data_out[row_idx, "TARGETS_WITH_MAX_PROGRAMS"] <- paste(interval.maker(as.numeric(targets_with_max_programs), "vec"), collapse = ", ")
      
    } else {
      data_out[row_idx, "TARGETS_WITH_MAX_PROGRAMS"] <- NA
    }
    
    # ============================
    # Programas que detectan cada exón
    # ============================

    
    if (length(unique(CNV$EXON)) > 0) {
      programs_by_target_list <- lapply(sort(unique(CNV$EXON)), function(exon_name) {
        programs <- sort(unique(CNV[CNV$EXON == exon_name, "PROGRAM"]))
        exon_label <- paste0(exon_name, ":", paste(programs, collapse = ","))
        return(exon_label)
      })
      data_out[row_idx, "PROGRAMS_BY_TARGET"] <- paste(programs_by_target_list, collapse = "; ")
    } else {
      data_out[row_idx, "PROGRAMS_BY_TARGET"] <- NA
    }
    # promedio de los ratios de cobertura de cada herramienta
    data_out[row_idx, "DC_ratio"] <- mean(CNV[CNV$PROGRAM == "DC",]$DC_ratio)
    data_out[row_idx, "CN_ratio"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_ratio)
    data_out[row_idx, "C2_ratio"] <- mean(CNV[CNV$PROGRAM == "C2",]$C2_ratio)
    data_out[row_idx, "PM_ratio"] <- mean(CNV[CNV$PROGRAM == "PM",]$PM_ratio)
    data_out[row_idx, "CL_ratio"] <- mean(CNV[CNV$PROGRAM == "CL",]$CL_ratio)
    
    # valores estadísticos calculados por cada programa
    data_out[row_idx, "DC_BF"] <- mean(CNV[CNV$PROGRAM == "DC",]$DC_BF)
    data_out[row_idx, "CN_ZSCORE"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_ZSCORE)
    data_out[row_idx, "CN_SHAPIRO.WILK"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_SHAPIRO.WILK)
    data_out[row_idx, "C2_mBIC"] <- mean(CNV[CNV$PROGRAM == "C2",]$C2_mBIC)
    data_out[row_idx, "CL_potentialAF"] <- mean(CNV[CNV$PROGRAM == "CL",]$CL_potentialAF)
    data_out[row_idx, "CL_loglikelihood"] <- mean(CNV[CNV$PROGRAM == "CL",]$CL_loglikelihood)
    data_out[row_idx, "CL_qvalue"] <- mean(CNV[CNV$PROGRAM == "CL",]$CL_qvalue)

    row_idx <- row_idx + 1 # aumenta en 1 el índice de la fila
  }
  index_iter<-index_iter+1
}

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
# se crea una copia para realizar modificaciones sin modificar el original
data_out2 <- data_out

data_out2$TARGETS <- sapply(data_out2$TARGETS, function(x) interval.maker(x, "str"))
data_out2$DC_TARGETS <- sapply(data_out2$DC_TARGETS, function(x) interval.maker(x, "str"))
data_out2$CN_TARGETS <- sapply(data_out2$CN_TARGETS, function(x) interval.maker(x, "str")) 
data_out2$C2_TARGETS <- sapply(data_out2$C2_TARGETS, function(x) interval.maker(x, "str"))
data_out2$PM_TARGETS <- sapply(data_out2$PM_TARGETS, function(x) interval.maker(x, "str"))
data_out2$CL_TARGETS <- sapply(data_out2$CL_TARGETS, function(x) interval.maker(x, "str"))

# Estas líneas aplican la función interval.maker con el modo "str" a varias columnas del data_out2, 
# que contienen listas de exones como "1, 2, 3, 5, 6", 
# para convertir esas listas en intervalos compactos como "1-3, 5-6"


# Order the output
if(nrow(data_out2)>0){ # verifica si el data frame tiene alguna fila (si se detectaron CNVs)
  data_out2 <- data_out2[order(data_out2$NUM_OF_PROGRAMS, decreasing = TRUE),] # ordena de forma descendente el número de programas que han detectado cnvs
  data_out2[data_out2 == ""] <- NA # en celdas vacías se pone NA
}

# Escribe una cabecera de ayuda al principio del archivo .combined.txt, explicando cada columna de data_out2. 


# Final table, tabla final data out
write.table(data_out2, file = paste(projectname, ".combined_200625.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)

# To Annotate table
# SE PREPARA UN ARCHIVO LIMPIO DE CNVS PARA ANOTAR
data_out_toAnnotate <- data.frame(data_out2$CHR, data_out2$START, data_out2$STOP, data_out2$CNV_TYPE, stringsAsFactors = FALSE)
data_out_toAnnotate <-unique(data_out_toAnnotate)
write.table(data_out_toAnnotate, file = paste(projectname, "_combined_v2_toAnnotate.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


