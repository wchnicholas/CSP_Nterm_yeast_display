#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(stringr)
require(cowplot)

plot_enrich_heatmap <- function(norm_enrich_table, WTresibox){
  textsize <- 7
  p <-  ggplot() +
          geom_tile(data=norm_enrich_table,aes(x=resi,y=aa,fill=log10(norm_affinity))) +
          labs("Dose (mg)") +
          scale_fill_gradientn(colours=c("blue", "white", "red"),
                limits=c(-5.5,1.5),
                values=rescale(c(-5.5, 0, 1.5)),
                breaks=c(-5,-4,-3,-2,-1,0,1),
                labels=c('-5','-4','-3','-2','-1','0','1'),
                guide="colorbar",
                na.value="grey") +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=textsize,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=textsize,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 1,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.5, barheight = 6, title="Relative\naffinity")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='gray', size=0.5) +
          xlab("Position") +
          ylab("Amino acid")
  ggsave('graph/norm_affinity_heatmap.png',p,width=7, height=2.2, dpi=1200)
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W'))
WTresibox  <- read_tsv('data/WT_heatmap.tsv')
count_table  <- read_tsv('result/nterm_CSP_sub_count.tsv')
frac_R1G1 <- (0.3829*470455+0.3945*470408)/(470455+470408)
frac_R1G2 <- (0.0205*470455+0.0214*470408)/(470455+470408)
frac_R1G3 <- (0.0270*470455+0.0257*470408)/(470455+470408)
frac_R2G1 <- frac_R1G1*0.4878
frac_R2G2 <- frac_R1G2*0.4939
frac_R2G3 <- frac_R1G3*0.0747
enrich_table <- count_table %>% 
                  spread(Sample, Count) %>%
                  mutate(norm_Input=(Input+1)/sum(.$Input)) %>%
                  mutate(norm_R1G1=(R1G1+1)/sum(.$R1G1)*frac_R1G1) %>%
                  mutate(norm_R1G2=(R1G2+1)/sum(.$R1G2)*frac_R1G2) %>%
                  mutate(norm_R1G3=(R1G3+1)/sum(.$R1G3)*frac_R1G3) %>%
                  mutate(norm_R2G1=(R2G1+1)/sum(.$R2G1)*frac_R2G1) %>%
                  mutate(norm_R2G2=(R2G2+1)/sum(.$R2G2)*frac_R2G2) %>%
                  mutate(norm_R2G3=(R2G3+1)/sum(.$R2G3)*frac_R2G3) %>%
                  mutate(affinity=norm_R2G1/norm_R2G2)
WT_affinity   <- filter(enrich_table, Mut=='WT')$affinity
norm_enrich_table <- enrich_table %>%
                       filter(Mut!='WT') %>%
                       mutate(norm_affinity=affinity/WT_affinity) %>%
                       mutate(Pos=str_sub(Mut,2,-2)) %>%
                       mutate(resi=str_sub(Mut,1,-2)) %>%
                       mutate(aa=str_sub(Mut,-1,-1)) %>%
                       filter(aa %in% aa_level) %>%
                       mutate(aa=factor(aa,levels=aa_level)) %>%
                       mutate(Pos=factor(Pos,levels=as.character(seq(26,97)))) %>%
                       arrange(Pos) %>%
                       mutate(resi=factor(resi,levels=unique(resi))) %>%
                       mutate(Pos=as.numeric(as.character(Pos))) %>%
                       select(Mut, resi, Pos, aa, norm_affinity)
plot_enrich_heatmap(norm_enrich_table, WTresibox)
