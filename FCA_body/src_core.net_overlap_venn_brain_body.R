
library(VennDiagram)
## core net
x1=read.table('~/Documents/fly.brain_coexpr.net_2022.02/2022_02_fly.brain.core_coexpr.net/07.detect_core/brain_scRNA-seq_n15c0.005/cut43_net.txt',header=T);
x2=read.table('./brain_snRNA-seq_n15c0.005/cut12_net.txt',header=T);
head(x1)
head(x2)
dim(x1) #2140
dim(x2) #1894

# #overlap edge
edges1=paste(x1$from,x1$to) #brain
edges2=paste(x2$from,x2$to) #body
sum(edges1 %in% edges2) #572
overlap.edges=edges1[edges1 %in% edges2] #all RpS, "DnaJ-1 Hsp83"

gene1=unique(c(x1[,1],x1[,2]))
gene2=unique(c(x2[,1],x2[,2]))

(k1=length(gene1)) #205, brain
(k2=length(gene2)) #206, body
overlap.genes=gene2[gene2 %in% gene1]

# for overlapped edges
overlap <- calculate.overlap(x = list('Brain' = edges1, 'Body' = edges2))
area_overlap <- sapply(overlap, length)

g2 <- draw.pairwise.venn(area1 = area_overlap[1],
                        area2 = area_overlap[2],
                        cross.area = area_overlap[3],
                        category = c("Brain", "Body"),
                        ind = FALSE,
                        inverted = FALSE,
                        scaled = TRUE,
                        ext.text = TRUE,
                        lwd = 1, 
                        ext.line.lwd = 1,
                        ext.dist = -0.15,
                        ext.length = 0.9,
                        ext.pos = -4,
                        fill = c("darkorchid1",'cornflowerblue'),
                        cex = 6,
                        cat.cex = 6,
                        cat.col = c("black", "black"),
                        cat.just = list(c(1, -1), c(0, -1)),
                        cat.dist = c(-0.09, -0.09),#vertical set name position
                        cat.pos = c(0,0), #vertical set name position
                        rotation.degree = 0)
plot.new()
grid.draw(g2)

# for overlapped genes
overlap <- calculate.overlap(x = list('Brain' = gene1, 'Body' = gene2))
area_overlap <- sapply(overlap, length)
 g1 <- draw.pairwise.venn(area1 = area_overlap[1],
                          area2 = area_overlap[2],
                          cross.area = area_overlap[3],
                          category = c("Brain", "Body"),
                          ind = FALSE,
                          inverted = TRUE,
                          scaled = TRUE,
                          ext.text = TRUE,
                          lwd = 1, 
                          ext.line.lwd = 1,
                          ext.dist = -0.15,
                          ext.length = 0.9,
                          ext.pos = -4,
                          fill = c( "darkorchid1",'cornflowerblue'),
                          cex = 6,
                          cat.cex = 6,
                          cat.col = c("black", "black"),
                          cat.just = list(c(1, -3), c(0, -3)),
                          cat.dist = c(-0.16, -0.16),
                          cat.pos = c(10,10),
                          rotation.degree = 0)
  plot.new()
  grid.draw(g1)
  

jpeg("overlap.genes.jpg", width = 1200, height = 1200)
plot.new()
title(main = "Venn diagram of core genes", cex.main = 6, cex.sub = 5, line = -4, outer = TRUE)
grid.draw(g1)
dev.off()


jpeg("overlap.edges.jpg", width = 1200, height = 1200)
plot.new()
title(main = "Venn diagram of core edges", cex.main = 6, cex.sub = 5, line = -4, outer = TRUE)
grid.draw(g2)
dev.off()
