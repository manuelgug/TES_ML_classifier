

library(data.table)
library(dplyr)
library(dcifer)
library(reshape2)
library(tidyr)
library(purrr)
library(stringr)


site <- "Tete"


#select data type betweem "TRAINING_DATA" or "REAL_DATA"
DATA_TYPE = "TRAINING_DATA"


if (DATA_TYPE == "TRAINING_DATA") {
  
  # for the training data:
  PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_",site,".RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_",site,".RDS"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
} else if (DATA_TYPE == "REAL_DATA") {
  
  #the actual data
  PAIRS_GENOMIC <- read.csv(paste0("genomic_updated_",site,".csv"))
  PAIRS_METADATA <- read.csv(paste0("metadata_updated_",site,".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)

  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(read_counts = reads)  
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(NIDA = sampleID) # format genomic file
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>%
    separate(NIDA, into = c("NIDA", "run"), sep = "__", remove = TRUE)
  PAIRS_GENOMIC <- inner_join(PAIRS_METADATA, PAIRS_GENOMIC, by = "NIDA")
  
  PAIRS_METADATA <- PAIRS_METADATA %>% select(PairsID, NIDA, time_point) # format metadata file
  PAIRS_METADATA <- PAIRS_METADATA %>%
    pivot_wider(names_from = time_point, values_from = NIDA, names_prefix = "NIDA") %>%
    rename(NIDA1 = NIDAD0, NIDA2 = NIDADx)
  
} else {
  
  print("Incorrect data type. Options are 'TRAINING_DATA' and 'REAL_DATA'.")
  
}




####### DCIFER'S IBD #######------------------  

dsmp <- formatDat(PAIRS_GENOMIC, svar = "NIDA", lvar = "locus", avar = "allele")

#use already calculated coi instead?
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)

afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 

dres0 <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, 
                alpha = 0.05, nr = 1e3)   

gc()

suppressWarnings({
  dres0_long <- melt(dres0)
})
dres0_long$value <- ifelse(dres0_long$Var1 == dres0_long$Var2, 1, dres0_long$value) # put 1 if the sample is compared with itself

#need extra foramtting for real data after passing through dcifer because the nidas...
if (DATA_TYPE == "REAL_DATA"){
  
  dres0_long$Var1 <- as.character(dres0_long$Var1)
  dres0_long$Var2 <- as.character(dres0_long$Var2)
  dres0_long <- dres0_long %>%
    mutate(Var1 = if_else(!str_detect(Var1, "\\."), paste0(Var1, ".0"), Var1),
           Var2 = if_else(!str_detect(Var2, "\\."), paste0(Var2, ".0"), Var2))  # If no ".", add ".0"
}

dres0_long <- dres0_long[dres0_long$Var3 == "estimate" & !is.na(dres0_long$value),]
dres0_long <- dres0_long %>% select(-Var3)
colnames(dres0_long) <- c("infection1", "infection2", "IBD_estimate")


# First match: D0 with infection1 and Dx with infection2
match1 <- PAIRS_METADATA %>%
  left_join(dres0_long, by = c("NIDA1" = "infection1", "NIDA2" = "infection2")) %>%
  select(PairsID, NIDA1, NIDA2, IBD_estimate)

match1 <- match1[!is.na(match1$IBD_estimate),]


# Second match: Dx with infection1 and D0 with infection2
match2 <- PAIRS_METADATA %>%
  left_join(dres0_long, by = c("NIDA2" = "infection1", "NIDA1" = "infection2")) %>%
  select(PairsID, NIDA1, NIDA2, IBD_estimate)

match2 <- match2[!is.na(match2$IBD_estimate),]

dres0_long_final<- rbind(match1, match2)
dres0_long_final <- distinct(dres0_long_final)

dres0_long_final_summarized <- dres0_long_final %>%
  select(PairsID, IBD_estimate) %>%
  arrange(PairsID)

dres0_long_final_summarized <- dres0_long_final_summarized[complete.cases(dres0_long_final_summarized),]

# MERGE WITH METADATA
dres0_long_final_summarized <- merge(dres0_long_final_summarized, PAIRS_METADATA, by = "PairsID")

if (DATA_TYPE == "REAL_DATA"){
  
  metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  metadata_updated$post_effective_coi_med <- round(metadata_updated$post_effective_coi_med)
  
  metadata_updated_wide <- metadata_updated %>%
    pivot_wider(
      id_cols = PairsID, 
      names_from = time_point, 
      values_from = c(NIDA, post_effective_coi_med), 
      names_glue = "{.value}_{time_point}"
    )
 
  metadata_updated_wide$eCOI_pairs <- paste0(metadata_updated_wide$post_effective_coi_med_D0, "__", metadata_updated_wide$post_effective_coi_med_Dx)
   
  dres0_long_final_summarized <- merge(dres0_long_final_summarized, metadata_updated_wide[c("PairsID", "eCOI_pairs")], by = "PairsID")
  
}


write.csv(dres0_long_final_summarized, paste0(site,"_", DATA_TYPE,".csv"), row.names = F)

