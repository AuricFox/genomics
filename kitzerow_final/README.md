## Final Project



### Run Trimming

```
java -jar trimmomatic-0.39.jar SE -phred33 ../input/sars_spike_protein_raw_reads.fastq ../output.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Running Cleanup
Allows trimming, discards bad reads, and alpha value is 0.1.
```
./lighter -r ../output.fastq -k 17 5000000 0.1 -t 10 -trim -discard
```
