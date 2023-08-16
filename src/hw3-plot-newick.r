# Run with:
# Rscript hw3-plot-newick.r tree.txt tip-labels.txt tree.pdf
library('ape')
library('RColorBrewer')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('Need 3 args.')

# load tree
newtree <- read.tree(args[1])
cat('\nLoaded tree:\n')
print(newtree)

cat('\nStructure of tree object:\n')
print(str(newtree))

# load tip labels and colors
tip.labels <- read.table(args[2], sep='\t',head=F,row=1,check=F,comment='')

tip.labels <- tip.labels[newtree$tip.label,]
cat('\nLoaded tip label table:\n')
print(tip.labels)

cols <- as.character(tip.labels[,2])
tip.labels <- rownames(tip.labels)
cat('\nTip colors are:\n')
print(cols)
cat('\nTip labels are:\n')
print(tip.labels)

# plot tree with labeled nodes (boostrap) and tips (phylum)
cat('Plotting tree...\n')
pdf(args[3],width=8,height=8)
par(oma=c(0,0,0,0), mar=c(1,1,1,1))
plot.phylo(newtree,show.tip.label=FALSE,type='fan')
tiplabels(tip.labels,col=cols, frame='none', cex=.5)
dev.off()

