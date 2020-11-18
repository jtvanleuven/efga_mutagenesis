#libraries
{
  library(stringr)
  library(ggplot2)
  library(seqinr)
  library(reshape2)
  library(cowplot)
}

ref <- read.fasta("data/efgA_6c0z.fasta")

seqs <- read.table("data/efgA_all.blc", stringsAsFactors = F)
names <- seqs$V1[which(str_detect(seqs$V1,">"))]
seqs <- seqs$V1[which(!str_detect(seqs$V1,">"))]
seqs <- seqs[!nchar(seqs) < length(names)]
#aling_length <- length()
seqs.tab <- data.frame(matrix(, nrow=length(seqs), ncol=length(names)), stringsAsFactors = F)   
names.s <- str_replace_all(names,">", "")
names.s <- str_split(names.s,"/",simplify = T)[,1]
names(seqs.tab) <- names.s
for(i in 1:length(seqs)){
  seqs.tab[i,] <- unlist(str_split(seqs[i],""))
}
seqs.tab.t <- data.frame(t(seqs.tab), stringsAsFactors = F)
row.names(seqs.tab.t) <- names
ref_index <- which(str_detect(names,"6c0z"))
row.names(seqs.tab.t)[ref_index]
seqs.tab.t[ref_index,]
ref_numbs <- which(str_detect(seqs.tab.t[ref_index,], "[A-Z]"))
seqs.tab.t.short <- seqs.tab.t[,ref_numbs]
names(seqs.tab.t.short) <- 1:ncol(seqs.tab.t.short)
tmp <- seqs.tab.t.short
tmp$name <- row.names(tmp)
seqs.tab.wide <- melt(tmp, id.vars = c("name"))
seqs.tab.wide[seqs.tab.wide$value=="-",]$value <- NA
ggplot(seqs.tab.wide, aes(x=variable, y=name, fill=value, color=value)) +
  geom_tile()
ggsave(filename = "plots/aa_colors.pdf", width = 6, height= 20, units = "in")

#conservation score from jalview
#This is an automatically calculated quantitative alignment annotation which measures the number of conserved physico-chemical properties conserved for each column of the alignment. Its calculation is based on the one used in the AMAS method of multiple sequence alignment analysis :
#Livingstone C.D. and Barton G.J. (1993), Protein Sequence Alignments: A Strategy for the Hierarchical Analysis of Residue Conservation. CABIOS Vol. 9 No. 6 (745-756)).
efga_cons <- read.csv("data/efga_only_cons.csv", header = F, row.names = 1, stringsAsFactors = F)
all_cons <- read.csv("data/efga_all_cons.csv", header = F, row.names = 1, stringsAsFactors = F)
all.cons.short <- all_cons[,ref_numbs]
efga.cons.short <- efga_cons[,ref_numbs]
names(all.cons.short) <- as.character(1:length(ref_numbs))
names(efga.cons.short) <- as.character(1:length(ref_numbs))
row.names(efga.cons.short) <- "EfgA Conservation"
seqs.tab.t.short <- rbind(efga.cons.short, all.cons.short, seqs.tab.t.short)

plot.cons <- seqs.tab.t.short[1:2,]
plot.cons$measure <- row.names(plot.cons)
plot.cons.w <- melt(plot.cons, id.vars = c("measure"))
names(plot.cons.w) <- c("measure", "site", "conservation")
plot.cons.w$conservation <- as.numeric(plot.cons.w$conservation)
plot.cons.w$site <- as.numeric(plot.cons.w$site)

p.cons <- ggplot(plot.cons.w, aes(x=site, y=conservation, color=measure, group=measure)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~measure, nrow = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), legend.position = "none")
p.cons

p.cons.all <- ggplot(plot.cons.w[plot.cons.w$measure=="Conservation",], aes(x=site, y=conservation, color=measure, group=measure)) +
  geom_line(color="grey40", size=1.1) +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
  ylab("DUF336 Conservation")
p.cons.all

p.cons.efga <- ggplot(plot.cons.w[plot.cons.w$measure=="EfgA Conservation",], aes(x=site, y=conservation, color=measure, group=measure)) +
  geom_line(color="grey40", size=1.1) +
  theme_classic() +
  theme(text = element_text(size=15),axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "none")+
  ylab("EfgA Conservation")
p.cons.efga



foldx <- read.csv("data/efga_mut_ddgfold_sd_ddgbind_sd.csv", sep="\t")
names(foldx) <- c("Mutation", "ddG_fold", "foldErr", "ddG_bind", "binderr")
foldx$site <- as.numeric(str_extract(foldx$Mutation, "\\d+"))

p.fold <- ggplot(foldx[,c("ddG_fold","site")], aes(x=site, y=ddG_fold)) +
  geom_hline(yintercept = 0, color="grey") +
  geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab(expression(Delta~Delta~G[fold])) +
  xlab("") +
  theme(text = element_text(size=15), axis.title.x = element_blank(), axis.text.x = element_blank())+
  xlim(c(0,144))
p.fold

p.bind <- ggplot(foldx[,c("ddG_bind","site")], aes(x=site, y=ddG_bind)) +
  geom_hline(yintercept = 0, color="grey") +
  geom_jitter(height = 0.1, width=0.3, alpha=0.3) +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=2, color="red") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ylab(expression(Delta~Delta~G[bind])) +
  xlab("") +
  theme(text = element_text(size=15))+
  xlim(c(0,144))
p.bind

plot_grid(p.cons.all, p.cons.efga, p.fold, p.bind, nrow=4, align = 'v')
ggsave(filename = "plots/combo.pdf", width = 6.5, height= 9, units = "in")

write.csv(seqs.tab.t.short, "results/indexed_consv.csv", quote = F, row.names = T) 
write.csv(foldx, "results/foldx.csv", quote = F, row.names = F)

tmp <- foldx[,c("ddG_fold", "ddG_bind", "site")]
foldx.agg <- tmp %>% group_by(site) %>%
  summarise_each(funs(.[which.max(abs(.))]))
})
names(foldx.agg) <- c("site", "ddG_fold_max", "ddG_bind_max")
write.csv(foldx.agg, "results/foldx_maxvals.csv", quote = F, row.names = F)
