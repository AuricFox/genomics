# Example input and outputs for HW3

# Your program should run with:
<programname> example.fna

# Your program should output these files:
genetic-distances.txt
edges.txt
tree.tre
bootstrap.txt

# Plots were generated with provided scripts:
Rscript ../hw3-plot-newick.r tree.tre example-tip-labels.txt tree-newick.pdf
Rscript ../hw3-plot-edges.r edges.txt example-tip-labels.txt tree.pdf
Rscript ../hw3-plot-edges.r edges.txt example-tip-labels.txt bootstrap.txt tree-bootstrap.pdf

# Note that tree-newick.pdf and tree.pdf are identical
# This is because the same tree is represented in
# the Newick format tree.tre
# and the edges format edges.txt
