## CSCI 5481

Samuel Kitzerow, kitze012

### Running Program

Main Program:
The default file input is ./input/Homework4-seqs-with-primers.fna
```
python .\main.py
python .\main.py filename
```

Phylogeny Program:  
```
python3 phylogeny.py input_file
python3 phylogeny.py hw3.fna
```

R Program:  
```
Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt tree.pdf
Rscript hw3-plot-newick.r tree.tre hw3-tip-labels.txt tree-newick.pdf

Rscript hw3-plot-edges.r ./output/edges.txt hw3-tip-labels.txt ./output/tree.pdf
Rscript hw3-plot-newick.r ./output/tree.tre hw3-tip-labels.txt ./output/tree-newick.pdf
```


### Special Note

The phylogney program was tested on VOLE3D. The module skbio was throwing unresolved errors
so the version was down graded to 0.5.6. `pip install scikit-bio==0.5.6`. The main program was
tested on windows with the EC portion omitted.
