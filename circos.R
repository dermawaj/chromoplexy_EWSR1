#####circos map######
library(circlize)

cytoband = read.cytoband()
cytoband_df = cytoband$df
chromosome = cytoband$chromosome

#########WT1########
bed1 <- read.delim("path/bed1_WT1.txt")
bed2 <- read.delim("path/bed2_WT1.txt")

circos.par(start.degree = 0)
xrange = c(cytoband$chr.len)
normal_chr_index = c(1:10,12:21,23:24)
zoomed_chr_index1 = c(11)
zoomed_chr_index2 = c(22)

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[1:10] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index1] * 0.2 / sum(xrange[zoomed_chr_index1]),
                 xrange[12:21] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index2] * 0.5 / sum(xrange[zoomed_chr_index1]),
                 xrange[23:24] / sum(xrange[normal_chr_index])) 

circos.initializeWithIdeogram(sector.width = sector.width, plotType = c("labels", "ideogram"))

circos.genomicLink(bed1, bed2, #border = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright", hue = "green"),
                   col = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright"),
                   lwd = c(rep(1,13),rep(18,18),rep(1,13)))
text(0.85, 0.65, "EWSR1", cex = 1, col = "#1b7837")
text(-1, -0.15, "WT1", cex = 1, col = "#2171b5")

############################
#####FLI1#########
bed1 <- read.delim("path/bed1_FLI1.txt")
bed2 <- read.delim("path/bed2_FLI1.txt")

xrange = c(cytoband$chr.len)
normal_chr_index = c(1:10,12:21,23:24)
zoomed_chr_index1 = c(11)
zoomed_chr_index2 = c(22)

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[1:10] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index1] * 0.2 / sum(xrange[zoomed_chr_index1]),
                 xrange[12:21] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index2] * 0.5 / sum(xrange[zoomed_chr_index1]),
                 xrange[23:24] / sum(xrange[normal_chr_index])) 

circos.initializeWithIdeogram(sector.width = sector.width, plotType = c("labels", "ideogram"))

circos.genomicLink(bed1, bed2, border = NA,
                   col = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright"),
                   lwd = c(rep(1,5),rep(34,34),rep(1,27)))
text(0.85, 0.65, "EWSR1", cex = 1, col = "#1b7837")
text(-0.9, 0.45, "FLI1", cex = 0.8, col = "#1b9e77")

########################################################
########ERG#######
bed1 <- read.delim("path/bed1_ERG.txt")
bed2 <- read.delim("path/bed2_ERG.txt")

xrange = c(cytoband$chr.len)
normal_chr_index = c(1:20,23:24)
zoomed_chr_index1 = c(21)
zoomed_chr_index2 = c(22)

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[1:20] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index1] * 0.3 / sum(xrange[zoomed_chr_index1]),
                 xrange[zoomed_chr_index2] * 0.3 / sum(xrange[zoomed_chr_index2]),
                 xrange[23:24] / sum(xrange[normal_chr_index])) 

circos.initializeWithIdeogram(sector.width = sector.width, plotType = c("labels", "ideogram"))
circos.genomicLink(bed1, bed2, border = NA,
                   col = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright", hue = "orange"),
                   lwd = c(rep(1,2),rep(5,5),rep(1,6),rep(6,6),rep(1,5)))
text(0.8, 0.7, "EWSR1", cex = 1, col = "#1b7837")
text(-0.05, 1.0, "ERG", cex = 1, col = "#984ea3")

circos.clear()

#######TMPRSS::ERG#######
bed1 <- read.delim("path/bed1_TMPRSS2.txt")
bed2 <- read.delim("path/bed2_TMPRSS2.txt")

xrange = c(cytoband$chr.len)
normal_chr_index = c(1:20,22:24)
zoomed_chr_index1 = c(21)

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[1:20] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index1] * 0.3 / sum(xrange[zoomed_chr_index1]),
                 xrange[22:24] / sum(xrange[normal_chr_index])) 

circos.initializeWithIdeogram(sector.width = sector.width, plotType = c("labels", "ideogram"))
circos.genomicLink(bed1, bed2, border = NA,
                   col = rand_color(nrow(bed1), transparency = 0.5, friendly = T, luminosity = "bright"),
                   lwd = c(rep(1,11),rep(17,17),rep(1,11)))
text(0.9, 0.55, "TMPRSS2", cex = 1, col = "#1b7837")
text(0.75, 0.65, "ERG", cex = 1, col = "#984ea3")

circos.clear()


############################################################
######circos all genes involved######
library(tidyverse)
library(circlize)
chromoplexy_genes <- read.csv("path/chromoplexy_genes.csv")

chromoplexy_genes_FLI1 <- chromoplexy_genes %>% 
  filter(Fusion == "EWSR1::FLI1") %>% 
  select(Case, Genes) %>% 
  mutate(Partner = lead(Genes)) %>% 
  slice(which(row_number() %% 2 == 1)) %>% 
  group_by(Genes, Partner) %>% tally() %>% 
  dplyr::rename(to = Genes, from = Partner, value = n) %>% select(from, to, value) %>% 
  arrange(desc(value))
chromoplexy_genes_FLI1$value <- c(4,2,rep(1,16))

chordDiagram(chromoplexy_genes_FLI1)

chordDiagram(chromoplexy_genes_FLI1, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chromoplexy_genes_FLI1))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

circos.clear()

chromoplexy_genes_ERG <- chromoplexy_genes %>% 
  filter(Fusion == "EWSR1::ERG") %>% 
  select(Case, Genes) %>% 
  mutate(Partner = lead(Genes)) %>% 
  slice(which(row_number() %% 2 == 1)) %>% 
  group_by(Genes, Partner) %>% tally() %>% 
  dplyr::rename(to = Genes, from = Partner, value = n) %>% select(from, to, value) %>% 
  arrange(desc(value))
chromoplexy_genes_ERG$value <- c(4,rep(1,7))

chordDiagram(chromoplexy_genes_ERG)
chordDiagram(chromoplexy_genes_ERG, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chromoplexy_genes_FLI1))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

circos.clear()


chromoplexy_genes_WT1 <- chromoplexy_genes %>% 
  filter(Fusion == "EWSR1::WT1") %>% 
  select(Case, Genes) %>% 
  mutate(Partner = lead(Genes)) %>% 
  slice(which(row_number() %% 2 == 1)) %>% 
  group_by(Genes, Partner) %>% tally() %>% 
  dplyr::rename(to = Genes, from = Partner, value = n) %>% select(from, to, value) %>% 
  arrange(desc(value))
chromoplexy_genes_WT1$value <- c(4,rep(1,18))

chordDiagram(chromoplexy_genes_WT1)
chordDiagram(chromoplexy_genes_WT1, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chromoplexy_genes_FLI1))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important


circos.clear()
#############################################################