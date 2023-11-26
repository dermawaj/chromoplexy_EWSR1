
#####BAF#####
library(ggplot2)
library(tidyverse)

EWSR1_multi_SV <- read.csv("path/EWSR1_multi_SV.csv")

BAF <- read.csv("path/BAF.csv")
BAF <- BAF %>% left_join(EWSR1_multi_SV %>% dplyr::select(Mnumber, SAMPLE_ID)) 

BAF <- BAF %>% mutate(chromoplexy = ifelse(chromoplexy == 1, "Yes", "No")) %>% 
  mutate(BAF_both = BAF_both *2) %>% 
  mutate(EWSR1_partner = factor(EWSR1_partner, levels = c("noncanonical", "canonical"))) %>% 
  # filter(fusion == "EWSR1::WT1") %>% 
  filter(fusion %in% c("EWSR1::FLI1","EWSR1::ERG"))

BAF_c <- BAF %>% filter(chromoplexy == "Yes")

plot2 <- BAF_c %>% 
  ggplot(aes(x = BAF_both, color = EWSR1_partner)) +
  geom_density() +
  scale_color_manual(values = c("#4daf4a","#e41a1c")) +
  theme_bw()


ks_result <- ks.test(BAF_c$BAF_both[BAF_c$EWSR1_partner=="noncanonical"],BAF_c$BAF_both[BAF_c$EWSR1_partner=="canonical"])
ks_text <- paste("K-S statistic", signif(ks_result$statistic, 4), "\n",
                 "p-value", signif(ks_result$p.value, 4))

density_group1 <- density(BAF_c$BAF_both[BAF_c$EWSR1_partner=="noncanonical"])
density_group2 <- density(BAF_c$BAF_both[BAF_c$EWSR1_partner=="canonical"])

density_group1$x[which.max(density_group1$y)]
density_group2$x[which.max(density_group2$y)]


plot_KS <- plot2 +
  annotate("text",
           label = ks_text,
           x = max(density_group1$x),
           y = max(density_group1$y),
           hjust = 1, vjust = 1,
           parse = FALSE, color = "black") +
  labs(title = "Ewing Sarcoma Chromoplexy Cases", x = "BAF", y = "Density") #+
labs(title = "Desmoplastic Small Round Cell Tumor Chromoplexy Cases", x = "BAF", y = "Density")

plot_KS 

#############################################################
####intrachromosomal breakpoint distance (breakpoint adjacency)######
library(tidyverse)

chromoplexy_intrachromosomal <- read.csv("path/chromoplexy_intrachromosomal.csv")

chromoplexy_intrachromosomal <- chromoplexy_intrachromosomal %>% 
  mutate(Chr = paste0("chr",Chr)) %>% 
  group_by(Case, Chr) %>% 
  mutate(dist = ((Pos - lag(Pos))/1000) + 0.0001) %>% 
  filter(!is.na(dist)) %>% filter(Chr %in% c("chr11","chr21","chr22"))

options(scipen = 999)
chromoplexy_intrachromosomal %>% ggplot(aes(x = dist, col = Chr)) +
  geom_density(linewidth = 2) +
  facet_wrap(~Chr, ncol = 1) +
  scale_x_continuous(trans = "log2",
                     labels = scales::number_format(accuracy = 0.001,
                                                    decimal.mark = '.')) +
  labs(y = "Density", title = "Intrachromosomal Breakpoint Distance (kb)", x = "") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.background =  element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
#############################################################
######chromoplexy types######
library(tidyverse)
library(ggsci)

#SV types
chromoplexy_types <- read.csv("path/chromoplexy_types.csv")

chromoplexy_types <- chromoplexy_types %>% group_by(Fusion,Type) %>% tally(name = "Count") 

chromoplexy_types %>% ggplot(aes(x=Fusion, y=Count, fill=Type)) +
  geom_bar(position="fill", stat = "identity") +
  scale_fill_npg(name = "", labels = c("Deletion","Inversion","Translocation")) +
  labs(title = "EWSR1-associated Chromoplectic Event Types", y = "Frequency") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

#interchromosomal vs intrachromosomal
chromoplexy_types <- read.csv("path/chromoplexy_types.csv")

chromoplexy_types <- chromoplexy_types %>% 
  mutate(Type = ifelse(Type == "TRA", "Interchromosomal", "Intrachromosomal")) %>% 
  group_by(Fusion,Type) %>% tally(name = "Count") 

chromoplexy_types %>% ggplot(aes(x=Fusion, y=Count, fill=Type)) +
  geom_bar(position="fill", stat = "identity", width = 0.5) +
  scale_fill_manual(name = "", labels = c("Interchromosomal", "Intrachromosomal"),
                    values = c("#4daf4a","#984ea3")) +
  labs(title = "Inter- vs. Intrachromosomal Events", y = "Frequency") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
#############################################################
#####ggbio#########
library(ggbio)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SummarizedExperiment)

data(genesymbol,package="biovizBase")

# columns(Homo.sapiens)
# wh <- genesymbol[c("EWSR1")]
# wh <- range(wh, ignore.strand = TRUE)
# p.txdb <- autoplot(Homo.sapiens, which = wh, 
#                    columns = c("GENENAME","REFSEQ"), names.expr = "GENENAME::REFSEQ",
#                    color="darkblue", fill="darkblue") + 
#   theme_bw() + theme_alignment(grid=F, border=FALSE)
# tracks(EWSR1 = p.txdb)
# 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# columns(txdb)


gr <- GRanges("chr22",IRanges(29664257,29700000))
gr
p2 <- autoplot(txdb, which = gr, color="darkblue", fill="darkblue", names.expr="gene_id:::tx_name") + 
  scale_x_sequnit("kb") + theme_bw() + theme_alignment(grid=FALSE, border=FALSE)

#EWSR1 transcript NM_013986 uc003aev.3 transcript variant 1 ENST000000397938.2
EWSR1.FLI1_exons <- read.csv("path/EWSR1-FLI1_exons.csv")
EWSR1.WT1_exons <- read.csv("path/EWSR1-WT1_exon.csv")
EWSR1.ERG_exons <- read.csv("path/EWSR1-ERG_exons.csv")

# gr.fli1 <- makeGRangesFromDataFrame(EWSR1.FLI1_exons)
# mcols(gr.fli1) <- DataFrame(score=EWSR1.FLI1_exons$score)
# gr.fli1 <- makeGRangesFromDataFrame(EWSR1.FLI1_exons,
#                          keep.extra.columns = F, ignore.strand = T,
#                          seqinfo = NULL,
#                          seqnames.field = c("seqnames"),
#                          start.field = "start",
#                          end.field = "end")
# gr.b <- GRanges(seqnames=EWSR1.FLI1_exons$seqnames, 
#             IRanges(start=EWSR1.FLI1_exons$start,width=(rep(1, nrow(EWSR1.FLI1_exons)))),
#             ignore.strand = T,
#             value=EWSR1.FLI1_exons$Partner,
#             score=rnorm(nrow(EWSR1.FLI1_exons),1,0.001))
gr.b <- GRanges(seqnames=EWSR1.FLI1_exons$seqnames, 
                IRanges(start=EWSR1.FLI1_exons$start,width=(rep(1, nrow(EWSR1.FLI1_exons)))),
                ignore.strand = T,
                value=factor(EWSR1.FLI1_exons$Partner, levels = c("FLI1","non-FLI1")),
                score=EWSR1.FLI1_exons$score)
gr.b

gr.c <- GRanges(seqnames=EWSR1.WT1_exons$seqnames, 
                IRanges(start=EWSR1.WT1_exons$start,width=(rep(1, nrow(EWSR1.WT1_exons)))),
                ignore.strand = T,
                value=factor(EWSR1.WT1_exons$Partner, levels = c("WT1","non-WT1")),
                score=EWSR1.WT1_exons$score)
gr.c

gr.d <- GRanges(seqnames=EWSR1.ERG_exons$seqnames, 
                IRanges(start=EWSR1.ERG_exons$start,width=(rep(1, nrow(EWSR1.ERG_exons)))),
                ignore.strand = T,
                value=EWSR1.ERG_exons$Partner,
                score=EWSR1.ERG_exons$score)
gr.d

# p3 <- autoplot(gr.fli1, geom="point", aes(y=score), position=position_jitter(width = 0.00, height = 0.01)) + 
#   theme_bw() + theme_alignment(grid=T, border=FALSE)
p4 <- autoplot(gr.b, geom = "point", aes(y=score, color = value, shape = value)) + 
  scale_color_manual(values = c("darkgreen","purple")) +
  scale_x_sequnit("kb") + theme_bw() + theme_alignment(grid=F, border=FALSE) + theme(legend.position = "left")
p4 <- autoplot(gr.c, geom = "point", aes(y=score, color = value, shape = value)) + 
  scale_color_manual(values = c("orange","purple")) +
  scale_x_sequnit("kb") + theme_bw() + theme_alignment(grid=F, border=FALSE) + theme(legend.position = "left")
p4 <- autoplot(gr.d, geom = "point", aes(y=score, color = value, shape = value)) + 
  scale_color_manual(values = c("blue","purple")) +
  scale_x_sequnit("kb") + theme_bw() + theme_alignment(grid=F, border=FALSE) + theme(legend.position = "left")

tracks(p4,p2,heights = c(0.5,2)) + 
  scale_x_sequnit("kb")

############################################################