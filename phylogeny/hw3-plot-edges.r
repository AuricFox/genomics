# Run with:
# Rscript hw3-plot-edges.r edges.txt tip-labels.txt tree.pdf
# or 
# Rscript hw3-plot-edges.r edges.txt tip-labels.txt bootstrap.txt tree-bootstrap.pdf
library('ape')
library('RColorBrewer')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3 && length(args) != 4) stop('Need 3 or 4 args.')

# make a translator between two edges matrices
# that have the same tip _labels_ but different
# tip IDs and internal nodes IDs.
# Translator output will be a 2-column matrix
# with the node ID of the first tree in col 1
# and node ID in second tree in col 2
# tiplabels1 order corresponds to tip IDs 1...ntips
# in edges1
# tiplabels2 order corresponds to tip IDs 1...ntips
# in edges2
translate.internal.nodes <- function(edges1, edges2, tiplabels1, tiplabels2){
	res <- matrix(0,nr=length(unique(edges1[,1])),ncol=2)
	res[,1] <- sort(unique(edges1[,1]))

	# first assign parent IDs of tips
	for(i in 1:nrow(edges1)){
		if(!(edges1[i,2] %in% edges1[,1])){
			# this is a tip, it's never a parent
			# find it in the edges2 mat and
			# set the correct parent in res
			tip.id1 <- edges1[i,2]
			tiplabel1 <- tiplabels1[tip.id1]
			tip.id2 <- which(tiplabels2 == tiplabel1)
			parent1.id <- edges1[i,1]
			parent2.id <- edges2[edges2[,2] == tip.id2,1]
			res[res[,1] == parent1.id,2] <- parent2.id
		}
	}

	# now for each translated pair, look up their
	# parents and assign those;
	# repeat while there are any unassigned

	while(any(res[,2] == 0)){
	# for(j in 1){
		for(i in 1:nrow(res)){
			if(res[i,2] != 0){
				parent1.id <- edges1[edges1[,2] == res[i,1],1]
				parent2.id <- edges2[edges2[,2] == res[i,2],1]
				res[res[,1] == parent1.id,2] <- parent2.id
			}
		}
	}
	return(res)
}

helper <- function(edge,edge.length,tip.label,node.label,root.ix=NULL){
# assumes tips are the first indices (i.e. 1, 2, 3, ...)
# and that those indices correspond to tip.label string vector
# If root is NULL, then this is the initial call to the function.
# Find the root node as the internal node that is never a child
# Otherwise, this is a recursive call; print postorder traversal 
# of subtree starting with root.ix.

	terminal.semicolon <- FALSE
	if(is.null(root.ix)){
		terminal.semicolon <- TRUE
		root.ix <- setdiff(edge[,1], edge[,2])
	}
	if(!any(edge[,1] == root.ix)){
		# Base case: If root is a tip, print its name
		row.ix <- which(edge[,2] == root.ix)
		ret.val <- paste(tip.label[root.ix],':',edge.length[row.ix],sep='')
	} else {
		# Recursion: get trees for children, then add this node
		# like this: (child1tree,child2tree,...):edgelength.to.parent
		ret.val <- '('
		for(row.ix in which(edge[,1] == root.ix)){
			subtree.val <- helper(edge, edge.length, tip.label, root.ix=edge[row.ix,2])
			delimiter <- ',' 
			if(ret.val == '(') delimiter <- '' # only add delimiter if not first child
			ret.val <- paste(ret.val,delimiter,subtree.val,sep='')
		}
		ret.val <- paste(ret.val,')',sep='')


		# test whether this is a child node and has a length to parent
		if(any(edge[,2] == root.ix)){
			ret.val <- paste(ret.val,':',edge.length[which(edge[,2] == root.ix)],sep='')
		}
	}
	if(terminal.semicolon) ret.val <- paste(ret.val,';',sep='')
	return(ret.val)
}

# used to convert tree-like list to actual "phylo" object
# because I couldn't manage to cast a list directly as a phylo object
# due to limited documentation of package
make.phylo.from.treelike.list <- function(t){
	t <- read.tree(text=helper(t$edge,t$edge.length,t$tip.label))
	return(t)
}



# load edges
edges <- as.matrix(read.table(args[1],sep='\t',head=F,check=F,comment=''))
cat('\nLoaded edge matrix:\n')
print(edges)

# load tip labels and colors
tip.labels <- read.table(args[2], sep='\t',head=F,row=1,check=F,comment='')
cat('\nLoaded tip label table:\n')
print(tip.labels)

cols <- as.character(tip.labels[,2])
tip.labels <- rownames(tip.labels)
names(cols) <- tip.labels
cat('\nTip colors are:\n')
print(cols)
cat('\nTip labels are:\n')
print(tip.labels)

# construct new tree
newtree <- list(edge=as.matrix(edges[,1:2]),
                edge.length=edges[,3],
                tip.label=tip.labels,
                Nnode=as.integer(max(edges[,1]) - min(edges[,1]) + 1))


class(newtree) <- 'phylo'
raw.edges <- as.matrix(edges)
newtree <- make.phylo.from.treelike.list(newtree)

cat("\nOriginal tip label order:\n")
print(tip.labels)

cat("\nTip labels order after converting to tree object:\n")
print(newtree$tip.label)

cat("\nOriginal edge matrix:\n")
print(raw.edges)

cat("\nEdges matrix after converting to tree object:\n")
print(cbind(newtree$edge, newtree$edge.length))

# bootstrap support
boots <- NULL
if(length(args) == 4) {
	raw.boots <- read.delim(args[3],sep='\t',head=F)
	cat('\nBootstrap values are:\n')
	print(raw.boots)
	# convert original bootstrap order to match new order in tree,
	# using the fact that tip label strings are preserved after
	# writing/reading newick
	mapping <- translate.internal.nodes(raw.edges, newtree$edge, tip.labels, newtree$tip.label)
	cat('\nMapping of original internal node ID to tree object internal node ID: \n')
	print(mapping)
	boots <- numeric(max(raw.boots[,1]))
	for(i in 1:nrow(mapping)){
		boot.value <- raw.boots[raw.boots[,1] == mapping[i,1],2]
		boots[mapping[i,2]] <- boot.value
	}
	# truncate boots vector to start with first internal node
	boots <- boots[min(raw.boots[,1]):max(raw.boots[,1])]
	cat('\nFinal bootstrap values after translation to new internal node indices:')
	print(boots)
}


# plot tree with labeled nodes (boostrap) and tips (phylum)
cat('\nPlotting tree...\n')
pdf(args[length(args)],width=8,height=8)
par(oma=c(0,0,0,0), mar=c(1,1,1,1))
plot.phylo(newtree,show.tip.label=FALSE,type='fan')
if(!is.null(boots)){
    nodelabels(col='black',frame='none',pie=boots,width=10,height=10,cex=.3)
}
tiplabels(newtree$tip.label,col=cols[newtree$tip.label], frame='none', cex=.5)
dev.off()

