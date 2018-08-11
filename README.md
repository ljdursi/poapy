# A Simple Partial Order Alignment implementation

[![Build Status](https://travis-ci.org/ljdursi/poapy.svg?branch=master)](https://travis-ci.org/ljdursi/poapy)

This is a simple python implementation of a Partial Order Alignment for MSA,
based on

[Multiple sequence alignment using partial order graphs](http://bioinformatics.oxfordjournals.org/content/18/3/452.short) (2002) by Lee, Grasso, and Sharlow 
 and
[Generating consensus sequences from partial order ...](http://bioinformatics.oxfordjournals.org/content/19/8/999.short) (2003) by Lee

for education/demonstration purposes. A short description of the method can be found on [this article](http://simpsonlab.github.io/2015/05/01/understanding-poa/) on the SimpsonLab group blog.

Some examples are provided; they are run as

```
$ ./poa.py examples/example4.fa --html ex4.html
seq1            -------------------------------TGTACNTGT-TTGTGA-GG-CTA
seq2            ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq3            ---AGTTCCTGCTGCG--TTTGCTGGACTTATGTACTTGT-TTGTGA-GG-CAA
seq4            --AAGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGTGTTGTGANGG-CAA
seq5            ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGT-ANGG-CAA
seq6            ---AGTTCCTGCTGCG--TTTGCT------------------------------
seq7            ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TT-----------
seq8            ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq9            ---AGTTNCTGNTGNG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq10           --------------------------------GTACNTGT-TTGTGA-GG-CTA
seq11           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq12           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq13           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq14           ---AGTTCCTGCTGCT--TTTGCTGGACTGATGTACTTGA-TTGTGA-GG-CAA
seq15           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq16           ---AGTTCCTGCTGCG--CTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq17           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGC-GG-CAA
seq18           ---AGT-CCTGC-GCG--TTTGC-GGACGGATGTACTTG--TTGTGA-GG-CAA
seq19           -------------------------------------------------G-CAA
seq20           ------------------------------------------------GG-CAA
seq21           ---------------------------CTGATGTACTTG--TTGTGA-GGGCAA
seq22           ---AGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
seq23           -G----TTCTGCCTGCG-TTTGCTGAACTGATGTACTTGT-TAGT-A-AG-CAA
seq24           C---GTTACTGC-GGG-GTTTGCTGGACTCATG-ACTTTTGTNGTAG--G-CAA
Consensus0      --AAGTTCCTGCTGCG--TTTGCTGGACTGATGTACTTGT-TTGTGA-GG-CAA
```

which produces an HTML file ex4.html, which uses [viz.js](http://visjs.org/)
to give an interactive output of small alignments; the output from above looks
like:

![Example 4 Output](imgs/screenshot.png  "Example 4 Output")

**NOTE**: This version has been modified from the earlier version to remove
affine gap penalties (which I hadn't originally been doing correctly!).
The correct affine gap calculation introduces enough complication to somewhat
obscure the important part of this example code, which is the graph-string
alignment, so what remains is a simple linear gap model.
