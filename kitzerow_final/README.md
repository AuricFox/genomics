## Final Project
**Genome Sequencing of the SARS-CoV-2 Spike Protein using the de Bruijn Graph**   
Jonathan Dekraker<sup>1</sup>, José Solórzano<sup>2</sup>, Samuel Kitzerow<sup>3</sup>   
1. Department of Computer Science.
2. Department of Plant Pathology.
3. Department of Computer Science.


### Run Main Program

Use default file (sars_spike_protein_reads.fastq):  
```
python .\main.py
```

User entered parameters:
```
python .\main.py -l			# Input list indices [start:end]
python .\main.py -k			# Input kmer range
python .\main.py -kl		# Input both
```

### Run Alignment

```
python .\alignment.py
```

### Run Trimming

```
java -jar trimmomatic-0.39.jar SE -phred33 ../input/sars_spike_protein_raw_reads.fastq ../output.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Running Cleanup
Allows trimming, discards bad reads, and alpha value is 0.1.
```
./lighter -r ../output.fastq -k 17 5000000 0.1 -t 10 -trim -discard
```

edges.txt:
	- File that contains all directed edges. Created by printing set from main.py

eulerianCycle.py:
	- Computes and prints eulerian circuit from spike_protein_directed_graph.txt

eulerianPath.py:
	- Computes and prints eulerian path from spike_protein_directed_graph.txt

main.py:
	- Only change is the create_directed_graph function which creates spike_protein_directed_graph.txt

spike_protein_directed_graph.txt:
	- A directed graph where the graph is given in the form of an adjacency list (from edges.txt)

## Function Runtimes

The system becomes sluggish when k-mers greater than 5 are used. The addition of progress bars helps judge the completion of 
the program; however, the image function that creates the directed graph is slow. Upon investigation it was determined that
`nx.draw(fig, pos, **options)` and `if(save_fig): plt.savefig(file, dpi=500)` takes a considerable amount of time (with the later 
taking the longest to complete). 
