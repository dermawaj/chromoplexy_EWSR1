############################################################
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
###EWSR1::FLI1 multiSV oncoprint####
#oncoprint CDC73, FGFR1, EZH2, BCOR, COP1, PTEN, TP53, STAG2, CDKN2A
#overall heatmap#

ES_oncoprint <- read.csv("path/ES_FLI1_oncoprint.csv") 
ES_oncoprint <- ES_oncoprint %>% mutate(multiSV = ifelse(multiSV == 0, "No", "Yes"))
ES_oncoprint <- ES_oncoprint %>% 
  select(SAMPLE_ID, multiSV, TP53, CDKN2A, STAG2, ERF, EZH2, PTEN, BCOR, FGFR1, CREBBP, BCOR,
         CREBBP, ZFHX3, KRAS, RAD21, TERT,  MYC, DDR2, NF2,
         chr8p, chr8q)
ES_oncoprint1 <- ES_oncoprint %>% filter(multiSV == "No") 
ES_oncoprint2 <- ES_oncoprint %>% filter(multiSV == "Yes")

ha1 <- HeatmapAnnotation(multiSV = ES_oncoprint1[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))
ha2 <- HeatmapAnnotation(multiSV = ES_oncoprint2[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))

col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_insertion"="#8c510a",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             "upstream"="#542788",
             "splicing"="#542788",
             "IntragenicDeletion"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  IntragenicDeletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["IntragenicDeletion"], col = NA))
  },
  upstream = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["upstream"], col = NA))
  },
  nonframeshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_insertion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_insertion",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "upstream",
                                  "Amplification",
                                  "IntragenicDeletion",
                                  "Deletion"),
                           labels = c("Missense Mutation",
                                      "In-frame insertion",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Promoter mutation",
                                      "Amplification",
                                      "Intragenic Deletion",
                                      "Deletion"),
                           title_position = "topleft",
                           ncol = 2)
#multiSV No
ES_mut_mat1 <- t(ES_oncoprint1[3:ncol(ES_oncoprint1)])
colnames(ES_mut_mat1) <- ES_oncoprint1[,1]
ES_mut_mat1[is.na(ES_mut_mat1)] <- ""

ES_mut_mat_ht1 <- oncoPrint(ES_mut_mat1,
                            column_title = NULL, cluster_column_slices = F, 
                            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                            remove_empty_rows = F, remove_empty_columns = F,
                            column_split = factor(ES_oncoprint1[[2]]),
                            # column_order = 1:ncol(ES_mut_mat),
                            alter_fun = alter_fun,
                            col = col_mut, 
                            show_heatmap_legend = F,
                            heatmap_legend_param = heatmap_legend_mut,
                            height = unit(0.5 * nrow(ES_mut_mat1), "cm"),
                            width = unit(0.1 * ncol(ES_mut_mat1), "cm"),
                            bottom_annotation = ha1,
                            top_annotation = NULL,
                            # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                            #                                                                  height = unit(1.0, "cm"),
                            #                                                                  ylim = c(0, 10))),
                            right_annotation = NULL,
                            # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                            #                                                                width = unit(1, "cm"),
                            #                                                                ylim = c(0,1))),
                            show_pct = TRUE, show_row_names = F,
                            row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                            #row_title = "Alterations",
                            row_title_gp = gpar(fontsize = 10),
                            row_order = 1:nrow(ES_mut_mat1),
                            row_title_rot = 0, row_title_side = "left"
)
#multiSV yes
ES_mut_mat2 <- t(ES_oncoprint2[3:ncol(ES_oncoprint2)])
colnames(ES_mut_mat2) <- ES_oncoprint2[,1]
ES_mut_mat2[is.na(ES_mut_mat2)] <- ""

ES_mut_mat_ht2 <- oncoPrint(ES_mut_mat2,
                            column_title = NULL, cluster_column_slices = F, 
                            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                            remove_empty_rows = F, remove_empty_columns = F,
                            column_split = factor(ES_oncoprint2[[2]]),
                            # column_order = 1:ncol(ES_mut_mat),
                            alter_fun = alter_fun,
                            col = col_mut, 
                            show_heatmap_legend = T,
                            heatmap_legend_param = heatmap_legend_mut,
                            height = unit(0.5 * nrow(ES_mut_mat2), "cm"),
                            width = unit(0.1 * ncol(ES_mut_mat2), "cm"),
                            bottom_annotation = ha2,
                            top_annotation = NULL,
                            # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                            #                                                                  height = unit(1.0, "cm"),
                            #                                                                  ylim = c(0, 10))),
                            right_annotation = NULL,
                            # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                            #                                                                width = unit(1, "cm"),
                            #                                                                ylim = c(0,1))),
                            show_pct = TRUE, show_row_names = T,
                            row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                            #row_title = "Alterations",
                            row_title_gp = gpar(fontsize = 10),
                            row_order = 1:nrow(ES_mut_mat2),
                            row_title_rot = 0, row_title_side = "left"
)

ht_list <- ES_mut_mat_ht1 + ES_mut_mat_ht2
draw(ht_list, heatmap_legend_side = "bottom", column_title = "Ewing Sarcoma (EWSR1::FLI1)")
############################################################
###EWSR1::WT1 multiSV oncoprint####
#overall heatmap#

DSRCT_oncoprint <- read.csv("path/DSRCT_oncoprint.csv") 
DSRCT_oncoprint <- DSRCT_oncoprint %>% mutate(multiSV = ifelse(multiSV == 0, "No", "Yes"))
DSRCT_oncoprint <- DSRCT_oncoprint %>% 
  select(SAMPLE_ID, multiSV, ARID1A, TERT, STAG2, CRLF2, 
         PTEN, EZH2, TP53,
         chr1q, chr5p, chr5q, chr16p, chr16q, chr20p, chr20q, chr21q)
DSRCT_oncoprint1 <- DSRCT_oncoprint %>% filter(multiSV == "No") 
DSRCT_oncoprint2 <- DSRCT_oncoprint %>% filter(multiSV == "Yes")

ha1 <- HeatmapAnnotation(multiSV = DSRCT_oncoprint1[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))
ha2 <- HeatmapAnnotation(multiSV = DSRCT_oncoprint2[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))

col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_insertion"="#8c510a",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             "upstream"="#542788",
             "splicing"="#542788",
             "IntragenicDeletion"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  IntragenicDeletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["IntragenicDeletion"], col = NA))
  },
  upstream = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["upstream"], col = NA))
  },
  nonframeshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_insertion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_insertion",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "upstream",
                                  "Amplification",
                                  "IntragenicDeletion",
                                  "Deletion"),
                           labels = c("Missense Mutation",
                                      "In-frame insertion",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Promoter mutation",
                                      "Amplification",
                                      "Intragenic Deletion",
                                      "Deletion"),
                           title_position = "topleft",
                           ncol = 2)
#multiSV No
DSRCT_mut_mat1 <- t(DSRCT_oncoprint1[3:ncol(DSRCT_oncoprint1)])
colnames(DSRCT_mut_mat1) <- DSRCT_oncoprint1[,1]
DSRCT_mut_mat1[is.na(DSRCT_mut_mat1)] <- ""

DSRCT_mut_mat_ht1 <- oncoPrint(DSRCT_mut_mat1,
                               column_title = NULL, cluster_column_slices = F, 
                               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                               remove_empty_rows = F, remove_empty_columns = F,
                               column_split = factor(DSRCT_oncoprint1[[2]]),
                               # column_order = 1:ncol(DSRCT_mut_mat),
                               alter_fun = alter_fun,
                               col = col_mut, 
                               show_heatmap_legend = F,
                               heatmap_legend_param = heatmap_legend_mut,
                               height = unit(0.5 * nrow(DSRCT_mut_mat1), "cm"),
                               width = unit(0.1 * ncol(DSRCT_mut_mat1), "cm"),
                               bottom_annotation = ha1,
                               top_annotation = NULL,
                               # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                               #                                                                  height = unit(1.0, "cm"),
                               #                                                                  ylim = c(0, 10))),
                               right_annotation = NULL,
                               # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                               #                                                                width = unit(1, "cm"),
                               #                                                                ylim = c(0,1))),
                               show_pct = TRUE, show_row_names = F,
                               row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                               #row_title = "Alterations",
                               row_title_gp = gpar(fontsize = 10),
                               row_order = 1:nrow(DSRCT_mut_mat1),
                               row_title_rot = 0, row_title_side = "left"
)
#multiSV yes
DSRCT_mut_mat2 <- t(DSRCT_oncoprint2[3:ncol(DSRCT_oncoprint2)])
colnames(DSRCT_mut_mat2) <- DSRCT_oncoprint2[,1]
DSRCT_mut_mat2[is.na(DSRCT_mut_mat2)] <- ""

DSRCT_mut_mat_ht2 <- oncoPrint(DSRCT_mut_mat2,
                               column_title = NULL, cluster_column_slices = F, 
                               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                               remove_empty_rows = F, remove_empty_columns = F,
                               column_split = factor(DSRCT_oncoprint2[[2]]),
                               # column_order = 1:ncol(DSRCT_mut_mat),
                               alter_fun = alter_fun,
                               col = col_mut, 
                               show_heatmap_legend = T,
                               heatmap_legend_param = heatmap_legend_mut,
                               height = unit(0.5 * nrow(DSRCT_mut_mat2), "cm"),
                               width = unit(0.1 * ncol(DSRCT_mut_mat2), "cm"),
                               bottom_annotation = ha2,
                               top_annotation = NULL,
                               # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                               #                                                                  height = unit(1.0, "cm"),
                               #                                                                  ylim = c(0, 10))),
                               right_annotation = NULL,
                               # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                               #                                                                width = unit(1, "cm"),
                               #                                                                ylim = c(0,1))),
                               show_pct = TRUE, show_row_names = T,
                               row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                               #row_title = "Alterations",
                               row_title_gp = gpar(fontsize = 10),
                               row_order = 1:nrow(DSRCT_mut_mat2),
                               row_title_rot = 0, row_title_side = "left"
)

ht_list <- DSRCT_mut_mat_ht1 + DSRCT_mut_mat_ht2
draw(ht_list, heatmap_legend_side = "bottom", column_title = "Desmoplastic Small Round Cell Tumor (EWSR1::WT1)")

###############################################################################
###EWSR1::ERG multiSV oncoprint####
#overall heatmap#

ES_oncoprint <- read.csv("path/ES_ERG_oncoprint.csv") 
ES_oncoprint <- ES_oncoprint %>% mutate(multiSV = ifelse(multiSV == 0, "No", "Yes"))
ES_oncoprint1 <- ES_oncoprint %>% filter(multiSV == "No") 
ES_oncoprint2 <- ES_oncoprint %>% filter(multiSV == "Yes")

ha1 <- HeatmapAnnotation(multiSV = ES_oncoprint1[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))
ha2 <- HeatmapAnnotation(multiSV = ES_oncoprint2[[2]],
                         col = list(multiSV = c("No" = "#e0ecf4", "Yes" = "#8856a7")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           multiSV = list(title = "Chromoplexy",
                                          at = c("No", "Yes"),
                                          labels = c("No", "Yes"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))

col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_insertion"="#8c510a",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             "upstream"="#542788",
             "splicing"="#542788",
             "IntragenicDeletion"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  IntragenicDeletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["IntragenicDeletion"], col = NA))
  },
  upstream = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["upstream"], col = NA))
  },
  nonframeshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_insertion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_insertion",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "upstream",
                                  "Amplification",
                                  "IntragenicDeletion",
                                  "Deletion"),
                           labels = c("Missense Mutation",
                                      "In-frame insertion",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Promoter mutation",
                                      "Amplification",
                                      "Intragenic Deletion",
                                      "Deletion"),
                           title_position = "topleft",
                           ncol = 2)
#multiSV No
ES_mut_mat1 <- t(ES_oncoprint1[3:ncol(ES_oncoprint1)])
colnames(ES_mut_mat1) <- ES_oncoprint1[,1]
ES_mut_mat1[is.na(ES_mut_mat1)] <- ""

ES_mut_mat_ht1 <- oncoPrint(ES_mut_mat1,
                            column_title = NULL, cluster_column_slices = F, 
                            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                            remove_empty_rows = F, remove_empty_columns = F,
                            column_split = factor(ES_oncoprint1[[2]]),
                            #column_order = 1:ncol(ES_mut_mat1),
                            alter_fun = alter_fun,
                            col = col_mut, 
                            show_heatmap_legend = F,
                            heatmap_legend_param = heatmap_legend_mut,
                            height = unit(0.5 * nrow(ES_mut_mat1), "cm"),
                            width = unit(0.2 * ncol(ES_mut_mat1), "cm"),
                            bottom_annotation = ha1,
                            top_annotation = NULL,
                            # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                            #                                                                  height = unit(1.0, "cm"),
                            #                                                                  ylim = c(0, 10))),
                            right_annotation = NULL,
                            # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                            #                                                                width = unit(1, "cm"),
                            #                                                                ylim = c(0,1))),
                            show_pct = TRUE, show_row_names = F,
                            row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                            #row_title = "Alterations",
                            row_title_gp = gpar(fontsize = 10),
                            row_order = 1:nrow(ES_mut_mat1),
                            row_title_rot = 0, row_title_side = "left"
)
#multiSV yes
ES_mut_mat2 <- t(ES_oncoprint2[3:ncol(ES_oncoprint2)])
colnames(ES_mut_mat2) <- ES_oncoprint2[,1]
ES_mut_mat2[is.na(ES_mut_mat2)] <- ""

ES_mut_mat_ht2 <- oncoPrint(ES_mut_mat2,
                            column_title = NULL, cluster_column_slices = F, 
                            column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                            remove_empty_rows = F, remove_empty_columns = F,
                            column_split = factor(ES_oncoprint2[[2]]),
                            #column_order = 1:ncol(ES_mut_mat2),
                            alter_fun = alter_fun,
                            col = col_mut, 
                            show_heatmap_legend = T,
                            heatmap_legend_param = heatmap_legend_mut,
                            height = unit(0.5 * nrow(ES_mut_mat2), "cm"),
                            width = unit(0.2 * ncol(ES_mut_mat2), "cm"),
                            bottom_annotation = ha2,
                            top_annotation = NULL,
                            # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                            #                                                                  height = unit(1.0, "cm"),
                            #                                                                  ylim = c(0, 10))),
                            right_annotation = NULL,
                            # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                            #                                                                width = unit(1, "cm"),
                            #                                                                ylim = c(0,1))),
                            show_pct = TRUE, show_row_names = T,
                            row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                            #row_title = "Alterations",
                            row_title_gp = gpar(fontsize = 10),
                            row_order = 1:nrow(ES_mut_mat2),
                            row_title_rot = 0, row_title_side = "left"
)

ht_list <- ES_mut_mat_ht1 + ES_mut_mat_ht2
draw(ht_list, heatmap_legend_side = "bottom", column_title = "Ewing Sarcoma (EWSR1::ERG)")
###############################################################################
###DSRCT sequential multiSV oncoprint####

#overall heatmap#

DSRCT_oncoprint <- read.csv("path/DSRCT_sequential_oncoprint.csv") 
DSRCT_oncoprint <- DSRCT_oncoprint %>%
  mutate(Sample_type = factor(Sample_type, levels = c("Primary","Local Recurrence","Metastasis"))) %>% 
  dplyr::select(-Tumor_number)

DSRCT_longitudinal_Tumor_Purity <- read.delim("path/DSRCT_longitudinal_Tumor_Purity.txt")

DSRCT_longitudinal_BAF <- read.csv("Gpath/DSRCT_longitudinal_BAF.csv")
BAF <- DSRCT_longitudinal_BAF %>% 
  #filter(EWSR1_partner == "noncanonical") %>% 
  mutate(BAF_both = BAF_both *2) %>% dplyr::rename(BAF = BAF_both)

DSRCT_oncoprint <- DSRCT_oncoprint %>% 
  left_join(BAF %>% dplyr::select(SAMPLE_ID, BAF)) %>% distinct(SAMPLE_ID, .keep_all = T) %>% left_join(Tumor_Purity) 
write.csv(DSRCT_oncoprint, "path/DSRCT_sequential_oncoprint_.csv", row.names = F) 
DSRCT_oncoprint <- read.csv("path/DSRCT_sequential_oncoprint_.csv") 

col_CN = colorRamp2(c(1, 0.5, 0), c("#a50f15", "#fb6a4a", "#fff5f0"))  

heatmap_legend_CN <- list(title = "Breakpoint Allele Frequency", title_position = "topleft",
                          direction = "horizontal", legend_width = unit(3, "cm"))

BAF_mat <- DSRCT_oncoprint %>% 
  mutate(BAF_c = BAF_c/(Tumor_Purity/100)) %>% 
  mutate(BAF_nc = BAF_nc/(Tumor_Purity/100)) %>% 
  dplyr::select(SAMPLE_ID, BAF_c, BAF_nc) 

BAF_ht <- HeatmapAnnotation(BAF_nc = BAF_mat[[3]],
                            BAF_c = BAF_mat[[2]],
                            col = list(BAF_nc = col_CN,
                                       BAF_c = col_CN),
                            show_annotation_name = T,
                            annotation_name_rot = 0,
                            annotation_name_gp = gpar(fontsize = 10),
                            annotation_legend_param = list(title = "BAF (normalized)",
                                                           direction = "horizontal"))

DSRCT_oncoprint <- DSRCT_oncoprint %>% dplyr::select(-c(BAF_c, BAF_nc, Tumor_Purity)) 

ha1 <- HeatmapAnnotation(Sample = DSRCT_oncoprint[[2]],
                         col = list(Sample = c("Primary" = "#f0f0f0", "Local Recurrence" = "#bdbdbd", "Metastasis" = "#636363")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           Sample = list(title = "Sample Type",
                                         at = c("Primary", "Local Recurrence","Metastasis"),
                                         labels = c("Primary","Local Recurrence", "Metastasis"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))


col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_insertion"="#8c510a",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             "upstream"="#542788",
             "splicing"="#542788",
             "IntragenicDeletion"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac",
             "Translocation" = "#88419d",
             "DEL" = "#045a8d"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  IntragenicDeletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["IntragenicDeletion"], col = NA))
  },
  upstream = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["upstream"], col = NA))
  },
  nonframeshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_insertion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  },
  Translocation = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col_mut["Translocation"], col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col_mut["DEL"], col = NA))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_insertion",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "upstream",
                                  "Amplification",
                                  "IntragenicDeletion",
                                  "Deletion",
                                  "Translocation",
                                  "DEL"),
                           labels = c("Missense Mutation",
                                      "In-frame insertion",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Promoter mutation",
                                      "Amplification",
                                      "Intragenic Deletion",
                                      "Deletion",
                                      "Translocation",
                                      "Intrachromosomal Deletion"),
                           title_position = "topleft",
                           ncol = 3)

DSRCT_mut_mat <- t(DSRCT_oncoprint[4:ncol(DSRCT_oncoprint)])
colnames(DSRCT_mut_mat) <- DSRCT_oncoprint[,1]
DSRCT_mut_mat[is.na(DSRCT_mut_mat)] <- ""

DSRCT_mut_mat_ht <- oncoPrint(DSRCT_mut_mat,
                              column_title = NULL, cluster_column_slices = F, 
                              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                              remove_empty_rows = F, remove_empty_columns = F,
                              column_split = factor(DSRCT_oncoprint[[3]]),
                              column_order = 1:ncol(DSRCT_mut_mat),
                              alter_fun = alter_fun,
                              col = col_mut, 
                              show_heatmap_legend = T,
                              heatmap_legend_param = heatmap_legend_mut,
                              height = unit(0.3 * nrow(DSRCT_mut_mat), "cm"),
                              width = unit(0.5 * ncol(DSRCT_mut_mat), "cm"),
                              bottom_annotation = ha1,
                              top_annotation = BAF_ht,
                              # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                              #                                                                  height = unit(1.0, "cm"),
                              #                                                                  ylim = c(0, 10))),
                              right_annotation = NULL,
                              # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                              #                                                                width = unit(1, "cm"),
                              #                                                                ylim = c(0,1))),
                              show_pct = TRUE, show_row_names = T,
                              row_names_side = "right", row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                              #row_title = "Alterations",
                              row_title_gp = gpar(fontsize = 10),
                              row_order = 1:nrow(DSRCT_mut_mat),
                              row_title_rot = 0, row_title_side = "left"
)


ht_list <- DSRCT_mut_mat_ht
draw(ht_list, heatmap_legend_side = "bottom", column_title = "DSRCT (sequential)")
###############################################################################
###ES sequential multiSV oncoprint####
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


#overall heatmap#

ES_oncoprint <- read.csv("path/ES_sequential_oncoprint.csv") 
ES_oncoprint <- ES_oncoprint %>%
  mutate(Sample_type = factor(Sample_type, levels = c("Primary","Local Recurrence","Metastasis"))) %>% 
  dplyr::select(-Tumor_number)

ES_longitudinal_Tumor_Purity <- read.delim("path/ES_longitudinal_Tumor_Purity.txt")

ES_longitudinal_BAF <- read.csv("path/ES_longitudinal_BAF.csv")
BAF <- ES_longitudinal_BAF %>% 
  filter(EWSR1_partner == "noncanonical") %>% 
  mutate(BAF_both = BAF_both *2) %>% dplyr::rename(BAF = BAF_both)

ES_oncoprint <- ES_oncoprint %>% 
  left_join(BAF %>% dplyr::select(SAMPLE_ID, BAF)) %>% distinct(SAMPLE_ID, .keep_all = T) %>% left_join(Tumor_Purity) 

write.csv(ES_oncoprint, "path/ES_sequential_oncoprint_.csv", row.names = F) 
ES_oncoprint <- read.csv("path/ES_sequential_oncoprint_.csv") 

col_CN = colorRamp2(c(1, 0.5, 0), c("#a50f15", "#fb6a4a", "#fff5f0"))  

heatmap_legend_CN <- list(title = "Breakpoint Allele Frequency", title_position = "topleft",
                          direction = "horizontal", legend_width = unit(3, "cm"))

BAF_mat <- ES_oncoprint %>% 
  mutate(BAF_c = BAF_c/(Tumor_Purity/100)) %>% 
  mutate(BAF_nc = BAF_nc/(Tumor_Purity/100)) %>% 
  dplyr::select(SAMPLE_ID, BAF_c, BAF_nc) 

BAF_ht <- HeatmapAnnotation(BAF_nc = BAF_mat[[3]],
                            BAF_c = BAF_mat[[2]],
                            col = list(BAF_nc = col_CN,
                                       BAF_c = col_CN),
                            show_annotation_name = T,
                            annotation_name_rot = 0,
                            annotation_name_gp = gpar(fontsize = 10),
                            annotation_legend_param = list(title = "BAF (normalized)",
                                                           direction = "horizontal"))

ES_oncoprint <- ES_oncoprint %>% dplyr::select(-c(BAF_c, BAF_nc, Tumor_Purity)) 

ha1 <- HeatmapAnnotation(Sample = ES_oncoprint[[2]],
                         col = list(Sample = c("Primary" = "#f0f0f0", "Local Recurrence" = "#bdbdbd", "Metastasis" = "#636363")
                         ),
                         show_annotation_name = F,
                         annotation_legend_param = list(
                           Sample = list(title = "Sample Type",
                                         at = c("Primary","Local Recurrence","Metastasis"),
                                         labels = c("Primary","Local Recurrence","Metastasis"))
                         ),
                         annotation_name_side = "right",
                         annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 10))


col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_insertion"="#8c510a",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             "upstream"="#542788",
             "splicing"="#542788",
             "IntragenicDeletion"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac",
             "Translocation" = "#88419d",
             "DEL" = "#045a8d"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  IntragenicDeletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["IntragenicDeletion"], col = NA))
  },
  upstream = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.6, 
              gp = gpar(fill = col_mut["upstream"], col = NA))
  },
  nonframeshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_insertion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  },
  Translocation = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col_mut["Translocation"], col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col_mut["DEL"], col = NA))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_insertion",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "upstream",
                                  "Amplification",
                                  "IntragenicDeletion",
                                  "Deletion",
                                  "Translocation",
                                  "DEL"),
                           labels = c("Missense Mutation",
                                      "In-frame insertion",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Promoter mutation",
                                      "Amplification",
                                      "Intragenic Deletion",
                                      "Deletion",
                                      "Translocation",
                                      "Intrachromosomal Deletion"),
                           title_position = "topleft",
                           ncol = 3)

ES_mut_mat <- t(ES_oncoprint[4:ncol(ES_oncoprint)])
colnames(ES_mut_mat) <- ES_oncoprint[,1]
ES_mut_mat[is.na(ES_mut_mat)] <- ""

ES_mut_mat_ht <- oncoPrint(ES_mut_mat,
                           column_title = NULL, cluster_column_slices = F, 
                           column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                           remove_empty_rows = F, remove_empty_columns = F,
                           column_split = factor(ES_oncoprint[[3]]),
                           column_order = 1:ncol(ES_mut_mat),
                           alter_fun = alter_fun,
                           col = col_mut, 
                           show_heatmap_legend = T,
                           heatmap_legend_param = heatmap_legend_mut,
                           height = unit(0.3 * nrow(ES_mut_mat), "cm"),
                           width = unit(0.5 * ncol(ES_mut_mat), "cm"),
                           bottom_annotation = ha1,
                           top_annotation = BAF_ht,
                           # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                           #                                                                  height = unit(1.0, "cm"),
                           #                                                                  ylim = c(0, 10))),
                           right_annotation = NULL,
                           # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                           #                                                                width = unit(1, "cm"),
                           #                                                                ylim = c(0,1))),
                           show_pct = TRUE, show_row_names = T,
                           row_names_side = "right", row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                           #row_title = "Alterations",
                           row_title_gp = gpar(fontsize = 10),
                           row_order = 1:nrow(ES_mut_mat),
                           row_title_rot = 0, row_title_side = "left"
)


ht_list <- ES_mut_mat_ht
draw(ht_list, heatmap_legend_side = "bottom", column_title = "Ewing sarcoma (sequential)")
############################################################