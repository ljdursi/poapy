# A Simple Partial Order Alignment implementation

This is a simple python implementation of a Partial Order Alignment for MSA,
based on

[Multiple sequence alignment using partial order graphs](http://bioinformatics.oxfordjournals.org/content/18/3/452.short) (2002) by Lee, Grasso, and Sharlow 
 and
[Generating consensus sequences from partial order ...](http://bioinformatics.oxfordjournals.org/content/19/8/999.short) (2003) by Lee

for education/demonstration purposes.

Some examples are provided; they are run as

```
$ ./poa.py examples/example4.fa --html ex4.html
seq1		------------------------------TGTACNT-GTTTGTGA-GGCTA
seq2		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq3		---A-GTTCCTGCTGCGTTTGCTGGACTTATGTACTT-GTTTGTGA-GGCAA
seq4		A--A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTTGGTTTGTGANGGCAA
seq5		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTNA-GGCAA
seq6		---A-GTTCCTGCTGCGTTTGCT-----------------------------
seq7		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTT----------
seq8		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq9		---A-GTTNCTGNTGNGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq10		-------------------------------GTACNT-GTTTGTGA-GGCTA
seq11		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq12		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq13		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq14		---A-GTTCCTGCTGCTTTTGCTGGACTGATGTACTT-GATTGTGA-GGCAA
seq15		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq16		---A-GTTCCTGCTGCGCTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq17		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGC-GGCAA
seq18		-A--G--TCCTGC-GCGTTTGC-GGACGGATGTACTT-G-TTGTGA-GGCAA
seq19		------------------------------------------------GCAA
seq20		-----------------------------------------------GGCAA
seq21		--------------------------CTGATGTACTT-G-TTGTGGAGGCAA
seq22		---A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
seq23		-----GTTCTGCCTGCGTTTGCTGAACTGATGTACTT-G-TTAGTA-AGCAA
seq24		--C--GTTACTGCGGGGTTTGCTGGACTCATG-ACTTTG-TTNGTA-GGCAA

Consensus 0	A--A-GTTCCTGCTGCGTTTGCTGGACTGATGTACTT-GTTTGTGA-GGCAA
```

which produces an HTML file ex4.html, which uses [viz.js](http://visjs.org/)
to give an interactive output of small alignments; the output from above looks
like:

![Example 4 Output](imgs/screenshot.png  "Example 4 Output")
