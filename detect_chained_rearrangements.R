library(tidyverse)
library(circlize)

######dbGAP vcf#######
library(vcfR)
library(tidyverse)
setwd("pathtovcffiles")

vcf_files <- read.table("vcf_files.txt")
vcf_file <- vcf_files[,1]
ES_all <- rbind()
chr_search <- c("22","11","16","2") %>% paste(collapse = "|")
for (i in 1:length(vcf_file)) {
  file <- vcf_file[i]
  vcf <- read.vcfR(file)
  df <- vcfR2tidy(vcf)
  sv <- df$fix
  ES_genes <- sv %>% filter(str_detect(CHROM, chr_search))
  ES_genes["sample"] <- vcf_file[i]
  ES_all <- rbind(ES_all, ES_genes)
}

sample_canonical <- ES_all %>%  
  mutate(gene1 = case_when(CHROM == "22" & POS > 29664257 & POS < 29696511 ~ "EWSR1",
                           CHROM == "16" & POS > 31191460 & POS < 31202926 ~ "FUS",
                           CHROM == "11" & POS > 128563967 & POS < 128683162 ~ "FLI1",
                           CHROM == "21" & POS > 39751958 & POS < 39870401 ~ "ERG",
                           CHROM == "2" & POS > 219845809 & POS < 219849906 ~ "FEV",
                           CHROM == "7" & POS > 13930854 & POS < 14029291 ~ "ETV1",
                           CHROM == "17" & POS > 41605214 & POS < 41623267 ~ "ETV4")) %>% 
  mutate(gene2 = case_when(CHR2 == "22" & POS2 > 29664257 & POS2 < 29696511 ~ "EWSR1",
                           CHR2 == "16" & POS2 > 31191460 & POS2 < 31202926 ~ "FUS",
                           CHR2 == "11" & POS2 > 128563967 & POS2 < 128683162 ~ "FLI1",
                           CHR2 == "21" & POS2 > 39751958 & POS2 < 39870401 ~ "ERG",
                           CHR2 == "2" & POS2 > 219845809 & POS2 < 219849906 ~ "FEV",
                           CHR2 == "7" & POS2 > 13930854 & POS2 < 14029291 ~ "ETV1",
                           CHR2 == "17" & POS2 > 41605214 & POS2 < 41623267 ~ "ETV4")) %>% 
  filter(!(is.na(gene1)) | !(is.na(gene2))) %>% select(c(2:15,21,25,"gene1","gene2","sample")) %>% filter(SVTYPE=="BND")
vcf_canonical <- unique(sample_canonical$sample)

#filter germline variants
vcf_keep <- rbind()
for (i in 1:length(vcf_canonical)) {
  file <- vcf_canonical[i]
  vcf <- read.vcfR(file)
  df <- vcfR2tidy(vcf)
  sv <- df$fix
  sv["sample"] <- vcf_canonical[i]
  vcf_keep <- rbind(vcf_keep, sv)
}

vcf_keep <- vcf_keep %>% select(c(2:15,21,25,"sample"))

vcf_canonical_pre <- paste0(sub(".vcf$",".pre.vcf", vcf_canonical))
vcf_pre <- rbind()
for (i in 1:length(vcf_canonical_pre)) {
  file <- vcf_canonical_pre[i]
  vcf <- read.vcfR(file)
  df <- vcfR2tidy(vcf)
  sv <- df$fix
  sv["sample"] <- vcf_canonical_pre[i]
  vcf_pre <- rbind(vcf_pre, sv)
}

vcf_pre <- vcf_pre %>% select(c(2:15,21,25,"sample"))

vcf_diff <- setdiff(vcf_keep$POS, vcf_pre$POS)
vcf_keep <- vcf_keep %>% filter(POS %in% vcf_diff)
vcf_ES <- vcf_keep %>%  
  mutate(gene1 = case_when(CHROM == "22" & POS > 29664257 & POS < 29696511 ~ "EWSR1",
                           CHROM == "16" & POS > 31191460 & POS < 31202926 ~ "FUS",
                           CHROM == "11" & POS > 128563967 & POS < 128683162 ~ "FLI1",
                           CHROM == "21" & POS > 39751958 & POS < 39870401 ~ "ERG",
                           CHROM == "2" & POS > 219845809 & POS < 219849906 ~ "FEV",
                           CHROM == "7" & POS > 13930854 & POS < 14029291 ~ "ETV1",
                           CHROM == "17" & POS > 41605214 & POS < 41623267 ~ "ETV4")) %>% 
  mutate(gene2 = case_when(CHR2 == "22" & POS2 > 29664257 & POS2 < 29696511 ~ "EWSR1",
                           CHR2 == "16" & POS2 > 31191460 & POS2 < 31202926 ~ "FUS",
                           CHR2 == "11" & POS2 > 128563967 & POS2 < 128683162 ~ "FLI1",
                           CHR2 == "21" & POS2 > 39751958 & POS2 < 39870401 ~ "ERG",
                           CHR2 == "2" & POS2 > 219845809 & POS2 < 219849906 ~ "FEV",
                           CHR2 == "7" & POS2 > 13930854 & POS2 < 14029291 ~ "ETV1",
                           CHR2 == "17" & POS2 > 41605214 & POS2 < 41623267 ~ "ETV4")) %>% 
  filter(!is.na(gene1) | !is.na(gene2)) %>% select(-c(gene1,gene2)) %>% 
  filter(SVTYPE == "BND")

sample_names <- unique(vcf_ES$sample)

for (sample_name in sample_names) {
  # Remove the ".vcf" extension
  sample_name_no_ext <- gsub(".vcf", "", sample_name)
  
  filtered_data <- vcf_keep %>%
    filter(sample == sample_name) %>% 
    select(CHROM, POS, CHR2, POS2, END, SVTYPE)
  
  vcf_ES_data <- vcf_ES[vcf_ES$"sample" == sample_name, ] %>% 
    select(CHROM, POS, CHR2, POS2, END, SVTYPE)
  
  combined <- rbind(filtered_data, vcf_ES_data) %>% mutate(identifier = paste0("row",row_number()))
  
  assigned_variable_name <- paste(sample_name_no_ext, sep="")
  assign(assigned_variable_name, combined)
}


#####chained rearrangements#####
#find chain rearrangements
# Assuming your dataframe is named 'sv_data'
# Create a subset of the data where CHROM is "22" and POS is either 29668451 or 29693850
result_list <- list()

for (i in seq_along(sample_names)) {
  sv_data_name <- paste0(gsub(".vcf", "", sample_names[i]))
  sv_data <- get(sv_data_name)
  
  start_rows <- rbind(subset(sv_data, CHROM == "22" & POS > 29664257 & POS < 29696511 & SVTYPE == "BND"),
                      subset(sv_data, CHROM == "16" & POS > 31191460 & POS < 31202926 & SVTYPE == "BND"))
  # Initialize an empty list to store the chains
  current_result_list <- list()
  
  # Loop through each start row
  for (j in seq_len(nrow(start_rows))) {
    # Initialize the current chain with the current start row
    current_chain <- start_rows[j, , drop = FALSE]
    
    # Determine the initial current position based on SVTYPE
    if (current_chain$SVTYPE == "BND") {
      current_pos <- current_chain$POS2
      current_chrom <- current_chain$CHR2
    } else {
      current_pos <- current_chain$END
      current_chrom <- current_chain$CHROM  # When SVTYPE is not "BND"
    }
    
    # Iterate through the remaining rows
    while (TRUE) {
      # Find rows where POS is within 10 mb of the current position
      next_row <- sv_data[
        abs(sv_data$POS - current_pos) <= 10000000 & 
          sv_data$CHROM == current_chrom &
          !(sv_data$identifier %in% current_chain$identifier)
        ,
      ]
      
      # If no matching rows are found, break the loop
      if (nrow(next_row) == 0) {
        break
      }
      
      # Add the matching rows to the current chain
      current_chain <- rbind(current_chain, next_row[1, , drop = FALSE])  # Use only the first row
      
      # Update the current position for the next iteration based on SVTYPE
      if (next_row$SVTYPE[1] == "BND") {
        current_pos <- unique(next_row$POS2)[1]
        current_chrom <- next_row$CHR2[1]
      } else {
        current_pos <- unique(next_row$END)[1]
        current_chrom <- next_row$CHROM[1]  # When SVTYPE is not "BND"
      }
    }
    
    # Append the current chain to the list if it's not empty
    if (nrow(current_chain) > 0) {
      current_result_list[[j]] <- current_chain
    }
  }
  
  # Create a named list element for the current sample
  named_list_element <- list(sample_name = sv_data_name, result_list = current_result_list)
  
  # Append the named list element for the current sample to the overall result list
  result_list[[i]] <- named_list_element
  
  # # Combine the list of chains for the current sample into a single dataframe
  # result_df <- do.call(rbind, current_result_list)
  # 
  # #  Assign the result to a variable with a name based on the sample
  # assigned_variable_name <- paste(sv_data_name, "result", sep = "_")
  # assign(assigned_variable_name, result_df)
  
}
result_df <- do.call(rbind, current_result_list)
# Loop through each sample's result_list
for (i in seq_along(result_list)) {
  sample_name <- result_list[[i]]$sample_name
  sample_result_list <- result_list[[i]]$result_list
  
  # Combine all chains into a single dataframe (assuming they have the same structure)
  combined_df <- do.call(rbind, sample_result_list)
  combined_df <- combined_df %>% distinct(CHROM, POS, CHR2, POS2, END)
  # Create a unique dataframe name based on the sample name
  df_name <- paste0(sample_name, "_combined_df")
  
  # Assign the combined dataframe to the global environment
  assign(df_name, combined_df)
}

# result_df <- SRR3538477_combined_df %>% 
# result_df <- SRR3589489_combined_df %>% 
# result_df <- SRR3765892_combined_df %>% 
# result_df <- SRR3784385_combined_df %>% 
result_df <- SRR3635106_combined_df %>%  
  mutate(gene1 = case_when(CHROM == "22" & POS > 29664257 & POS < 29696511 ~ "EWSR1",
                           CHROM == "16" & POS > 31191460 & POS < 31202926 ~ "FUS",
                           CHROM == "11" & POS > 128563967 & POS < 128683162 ~ "FLI1",
                           CHROM == "21" & POS > 39751958 & POS < 39870401 ~ "ERG",
                           CHROM == "2" & POS > 219845809 & POS < 219849906 ~ "FEV",
                           CHROM == "7" & POS > 13930854 & POS < 14029291 ~ "ETV1",
                           CHROM == "17" & POS > 41605214 & POS < 41623267 ~ "ETV4")) %>% 
  mutate(gene2 = case_when(CHR2 == "22" & POS2 > 29664257 & POS2 < 29696511 ~ "EWSR1",
                           CHR2 == "16" & POS2 > 31191460 & POS2 < 31202926 ~ "FUS",
                           CHR2 == "11" & POS2 > 128563967 & POS2 < 128683162 ~ "FLI1",
                           CHR2 == "21" & POS2 > 39751958 & POS2 < 39870401 ~ "ERG",
                           CHR2 == "2" & POS2 > 219845809 & POS2 < 219849906 ~ "FEV",
                           CHR2 == "7" & POS2 > 13930854 & POS2 < 14029291 ~ "ETV1",
                           CHR2 == "17" & POS2 > 41605214 & POS2 < 41623267 ~ "ETV4")) 
#SRR3635106; biospecimen_respository_subjectd_id SJDES006; BioSample: SAMN04621745; SRA: SRS1466425
#1538537	SJDES006	Female	8	Rib	"Lung, Liver, BM"	VAC/IE	White/Europe	EWS-FLI1 type 1	DOD
#SRR3589489; SJDES005; SAMN04621790
#1538536	SJDES005	Male	13	Gluteus muscle	Lung	VAC/IE	White/Europe	EWS-FLI1 type 2	NED
#SRR3765892; SJDES003-R; SAMN04621725
#1538534	SJDES003-R	Male	0.42	Femur	BM	VAC/IE	White/Europe	EWS-FLI1 type 1	NED
#SRR3784385; SJDES013; SAMN04621689
#1538544	SJDES013	Male	11	Vertebral body + rib	Lung	VAC/IE	White/Europe	EWS-FLI1 10-5	DOD
#SRR3538477; SJDES004; SAMN04621806

library(circlize)
bed1 <- result_df %>% 
  rename(chr = CHROM, start = POS) %>% 
  mutate(chr = paste0("chr",chr), end = start) %>% 
  select(chr, start, end)
bed2 <- result_df %>% 
  mutate(chr = ifelse(is.na(CHR2), CHROM, CHR2)) %>% 
  mutate(start = ifelse(is.na(POS2), END, POS2), end = start) %>% 
  mutate(chr = paste0("chr",chr)) %>% 
  select(chr, start, end)

circos.initializeWithIdeogram(plotType = c("labels", "ideogram"))
circos.initializeWithIdeogram()

circos.genomicLink(bed1, bed2)
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright"), 
                   border = NA, lwd = 2.5)
circos.clear()
#########################################