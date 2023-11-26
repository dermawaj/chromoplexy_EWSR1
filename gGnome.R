######gGnome#######
library(vcfR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gGnome)
library(data.table)
#load Junctions into gGnome from a variety of input formats from common junction / SV callers, 
#including SvaBa, DELLY, Novobreak using the function jJ.
# import DELLY
delly <- gGnome::jJ("~/PATH/.vcf")

head(delly)
delly$grl[1:2]
seqnames_grl <- delly$grl@unlistData@seqnames

#gGraph
gg <- gG(junctions = delly)

## we use gTrack to plot the gTrack associated with this gGraph
## the second argument to gTrack plot is a string or GRanges representing the
## window to plot, the links argument enables drawing of junctions from GRangesList
win_all = GRanges(seqnames=as.character(1:22),
                  ranges=IRanges(start=rep(1,22),
                                 end=seqlengths)
)

plot(gg$gt, win)
plot(gg$gt, win_all)
plot(gg$gt, '22', links = delly$grl)

#genome "breaks"
seqinfo(gg)
seqinfo_df <- data.frame(seqnames = seqinfo(gg)@seqnames,seqlengths = seqinfo(gg)@seqlengths)

### returns nodes data.table
nodes <- gg$nodes$dt[1:nrow(gg$nodes$dt)]

nodes %>% distinct(seqnames)

### returns edges data.table
edges <- gg$edges$dt[1:nrow(gg$edges$dt)]

### returns junctions data.table
junctions <- gg$edges$junctions$dt[1:nrow(gg$edges$junctions$dt)]

#gGraph
gg <- gG(junctions = delly)
#gwalk
gg$set(gr.labelfield = 'node.id')

#example case
hg_seqlengths()
win = GRanges(seqnames=c("11","21","22"),
              ranges=IRanges(start=c(1,1,1),
                             end=c(135006516,48129895,51304566)))

## define a simple window on chromosome 22
win = GRanges('22:1-51304566')

win = GRanges(seqnames=c("1","2","13","16","19","21","22"),
              ranges=IRanges(start=c(1,1,1,1,1,1,1),
                             end=c(249250621,243199373,115169878,90354753,59128983,48129895,51304566)))

seqlengths <- unname(hg_seqlengths())[1:22]

win_all = GRanges(seqnames=as.character(1:22),
                  ranges=IRanges(start=rep(1,22),
                                 end=seqlengths)
)

plot(gg$gt, win)

gTrack::plot.gTrack(c(gencode,gg$gt),win)
#EWSR1 22:29684827 node 914
#FLI1 11:128660148 node 302
#EWSR1 22:29688706 node 918
#PRR5-ARHGAP8 22:45165105 node 924
#EWSR1 22:29684530 node 914
#NCRNA00159 21:33419628 node 871

p2 = gg$paths(302,914)
p2$dt
p2$mark(col = 'green')
p1 = gg$paths(917,924)
p1$dt
p1$mark(col = 'blue')
p3 = gg$paths(871,914)
p3$dt
p3$mark(col = 'purple')
p4 = gg$paths(279,917)
p4$dt
p4$mark(col = 'darkgreen')

plot(gg$gt, win)

# subset ALT edges
## enumerate ALT edges classes
table(gg$edges[type == 'ALT']$dt$class)
## subset TRA-like edges
edges_tra_like <- gg$edges$dt[class == 'DEL-like']

#gtrack
## track.gencode, pulls hg19 GENCODE by default, modulates the stacking and font sizes
gencode = track.gencode(stack.gap = 1e5, cex.label = 0.8, height = 20)
gencode = track.gencode(build = 'hg19', genes = c("EWSR1"))

#add nodes and junctions
## randomly sample 4 width 1 GRanges representing SNV
snv = gr.sample(win, 4, wid = 1)

## disjoin with gr= argument breaks the graph at these SNV
gg1d = gg1$copy$disjoin(gr = snv)

## plot results
plot(gg1d$gt, win)

## new edges are specified as data.table with n1, n2, n1.side and n2.side
gg1d$add(edges = data.table(n1 = 3, n1.side = 'left', n2 = 7, n2.side = 'right'))

## plot 
plot(gg1d$gt, win)

## connect syntax specifies edges as pairs of "signed" node ids
## this means that the edge leaves the <right> side of the first signed node
## and enters the <left> side of the second signed node

## thus here we create an edge leaving the right side of 5 and entering the right side of 8
## (i.e. the left side of -8)
gg1d$connect(5, -8)

## this connects the right side of 3 and the left side of 9
gg1d$connect(3, 9)

## plot 
plot(gg1d$gt, win)

#edge clusters
#Clusters of quasi-reciprocal ALT edge can reveal “cycles” or “chains” of rearrangements in genome graphs. 
#Such patterns have been dubbed “chromoplexy” (multi-way balanced rearrangements) 
#and also linked to TIC (templated insertion chains). 
#We enable detection of these edge clusters with the gGraph method $eclusters. 
gg$eclusters()
gg$eclusters(thresh = 50)
## paths are labeled by a "p" prefix, and cycles labeled by a "c" prefix
## here we see a multi-junction cluster p16 with 8 edges
sort(table(gg$edges$dt$ecluster))
## we can mark these edges (and their associated nodes) with a special color
gg$edges[ecluster == 38]$mark(col = 'gray')
gg$edges[ecluster == 38]$nodes$mark(col = 'gray')
gg$edges[ecluster == 42]$mark(col = 'darkgreen')
gg$edges[ecluster == 42]$nodes$mark(col = 'darkgreen')
gg$edges[ecluster == 5]$mark(col = 'green')
gg$edges[ecluster == 5]$nodes$mark(col = 'green')
gg$edges[ecluster == 6]$mark(col = 'blue')
gg$edges[ecluster == 6]$nodes$mark(col = 'blue')
gg$edges[ecluster == 41]$mark(col = 'purple')
gg$edges[ecluster == 41]$nodes$mark(col = 'purple')
## here, the edges and nodes of the cluster that we have discovered
## are highlighted in purple 
plot(c(gg$gt), gg$edges[ecluster == 42]$nodes$gr %>% streduce(1e5))
plot(c(gg$gt), gg$edges[ecluster == 39]$nodes$gr %>% streduce(1e5))

#walk decomposition
gg.sub = gg[, ecluster == 17]
## walk decompositiion of this small subgraph
## generates 28 possible linear alleles
walks = gg.sub$walks()
## we can choose the longest walk (most nodes traversed)
walks = walks[which.max(length)]
## and plot with the walk track up top (with nodes already marked purple from before)
plot(c(gg$gt, walks$gtrack(name = "Longest walk")), walks$footprint+1e4)

#gwalk de novo
## retrieve signed node ids associated with a gWalk
nid = p5$snode.id

## retrieve signed edge ids associated with a gWalk
eid = p5$sedge.id

## you can instantiate a gWalk from signed node ids
gW(snode.id = nid, graph = p5$graph)
##    walk.id name length     wid circular
## 1:       1    1      9 6643813    FALSE
##                                                                                                                                                              gr
## 1: 2:202339385-203161868+ -> 2:203161869-205496567+ -> 3:177097274-177097491+ -> ... -> 21:16266971-16500987- -> 21:15890640-16266970- -> 21:14779256-15890639-
## or from signed edge ids
gW(sedge.id = eid, graph = p5$graph)

#protein fusions
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
# we are looking for any fusions connecting the genes EWSR1 using the genes= argument
fus <- fusions(gg, gff, genes = c('EWSR1'))
## fusions will output many "near duplicates" which just represent various combinations
## of near equivalent transcripts, we can filter these down using gWalk operations
ufus = fus[in.frame == TRUE][!duplicated(genes)]

#chromoplexy
gg <- chromoplexy(gg, max.cn = NA, mark = T, mark.col = 'purple')
gg <- chromoplexy(gg)
gg$gt
gg_gr <- as.data.table(gg$gr)
gg$edgesdt

gg$edges[which(chromoplexy>0)]


Sys.setenv(DEFAULT_GENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")

hg_seqlengths()

plot(gg$gt, gg$edges[which(chromoplexy>0)]$shadow %>% streduce(5e6)); title("Chromoplexy in EWSR1::ERG")
plot(gg$gt)

#simple
gg <- simple(gg)

############################################################