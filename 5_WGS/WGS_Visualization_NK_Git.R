###vcfR
library(vcfR)
library(ape)
#library(VariantAnnotation)
library(dplyr)
library(CMplot)
library(tidyverse)
library(data.table)

setwd("~/Yale University Dropbox/Kaiyuan Tang/Chenlab_Share/_shared_manuscript_files/Project_NKGOF/SourceData/WGS/VCF/")

#read in files
OE <- fread("grep -v '^##' ./Donor0958_OR7A10OE_Unique_Indel_Variant.vcf") %>%
  rename(Chromosome = 'CHROM')
Stop <- fread("grep -v '^##' ./Donor0958_OR7A10Stop_Unique_Indel_Variant.vcf") %>%
  rename(Chromosome = "CHROM")

#bining the variants
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
chrom_info <- fread("hg38.chrom.sizes", col.names = c("Chromosome", "Length")) %>%
  slice_head(n = 24) %>%
  mutate(Chromosome = factor(Chromosome, levels = chrom_order)) %>%
  arrange(Chromosome)

segment_size <- 1e7
chromosome = rep(chrom_info$Chromosome, ceiling(chrom_info$Length / segment_size))
start = c()
end = c()
for (i in 1:nrow(chrom_info)){
  start = c(start, seq(0, chrom_info$Length[i], by = segment_size))
  sequence = seq(segment_size, chrom_info$Length[i], by = segment_size)
  end = c(end, sequence, chrom_info$Length[i])
}
segment = data.frame(chromosome, start, end)

#NK variant freq
OE_seg <- segment %>% mutate(count = 0)
for (i in 1:nrow(OE)) {
  chr <- OE$Chromosome[i]
  pos <- OE$POS[i]
  depth <- as.numeric(str_split(str_split(OE$`Donor0958_OR7A10OE`[i], ":")[[1]][2], ",")[[1]][2])
  print(chr)
  OE_seg[chromosome == chr & start < pos & end > pos, "count"] <- OE_seg[chromosome == chr & start < pos & end > pos, "count"] + 1
}

Stop_seg <- segment %>% mutate(count = 0)
for (i in 1:nrow(Stop)) {
  chr <- Stop$Chromosome[i]
  pos <- Stop$POS[i]
  depth <- as.numeric(str_split(str_split(Stop$`Donor0958_OR7A10Stop`[i], ":")[[1]][2], ",")[[1]][2])
  print(chr)
  Stop_seg[chromosome == chr & start < pos & end > pos, "count"] <- Stop_seg[chromosome == chr & start < pos & end > pos, "count"] + 1
}

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("StructuralVariantAnnotation")
library(StructuralVariantAnnotation)

OE_SV <- VariantAnnotation::readVcf("Donor0958_OR7A10OE_Unique_SV1000.vcf")
Stop_SV <- VariantAnnotation::readVcf("Donor0958_OR7A10Stop_Unique_SV1000.vcf")

#filter out those breakpoints separated by only 50bp apart
OE_bpgr <- breakpointRanges(OE_SV)
OE_bedpe <- breakpointgr2bedpe(OE_bpgr)
OE_bedpe_filter <- OE_bedpe %>%
  filter(!((abs(start2 - start1) < 50) & (chrom1 == chrom2))) %>%
  dplyr::select(chrom1, start1, end1, chrom2, start2, end2)
OE_bedpe_filter

Stop_bpgr <- breakpointRanges(Stop_SV)
Stop_bedpe <- breakpointgr2bedpe(Stop_bpgr)
Stop_bedpe_filter <- Stop_bedpe %>%
  filter(!((abs(start2 - start1) < 50) & (chrom1 == chrom2))) %>%
  dplyr::select(chrom1, start1, end1, chrom2, start2, end2)
Stop_bedpe_filter

###Circos plot
library(RCircos)
GRCh38 <- fread("GRCh38_cytoBandIdeogram.txt") %>% #downloaded from UCSC table browser
  dplyr::rename(Chromosome = "#chrom", ChromStart = "chromStart", ChromEnd = "chromEnd", Band = "name", Stain = "gieStain") %>%
  dplyr::slice(1:864) %>%
  filter(Chromosome != "chrM")

chrX_Y <- GRCh38 %>%
  filter(Chromosome == "chrX" | Chromosome == "chrY")

GRCh38 <- GRCh38 %>% 
  filter(Chromosome != "chrX" & Chromosome != "chrY" & Chromosome != "chrUn_GL000195v1") %>%
  bind_rows(chrX_Y) %>%
  mutate(Chromosome = as.factor(Chromosome))

GRCh38_final <- GRCh38 %>%
  mutate(Stain = "gneg")
#data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL
cyto.info <- GRCh38
tracks.inside <- 2
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside) 
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$track.background <- "wheat"
rcircos.params$text.size <- 4
rcircos.params$track.height <- 0.1
rcircos.params$chrom.width <- 0.1
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Plot.Parameters();

#plot histogram 
out.file <- "CircosPlot_Donor0958_OR7A10OE.pdf";
pdf(file=out.file, height=10, width=10, compress=T)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
#RCircos.Gene.Name.Plot(gene_data, name.col=4, track.num = 2, side = "in")

rcircos.params$track.background <- "lightblue"
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.Histogram.Plot(OE_seg, data.col = 4, track.num = 1, side = "in")
RCircos.Link.Plot(OE_bedpe_filter, track.num = 2, TRUE)
dev.off()

#plot histogram 
out.file <- "CircosPlot_Donor0958_OR7A10Stop.pdf";
pdf(file=out.file, height=10, width=10, compress=T)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
#RCircos.Gene.Name.Plot(gene_data, name.col=4, track.num = 2, side = "in")

rcircos.params$track.background <- "grey"
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.Histogram.Plot(Stop_seg, data.col = 4, track.num = 1, side = "in")
RCircos.Link.Plot(Stop_bedpe_filter, track.num = 2, TRUE)
dev.off()
