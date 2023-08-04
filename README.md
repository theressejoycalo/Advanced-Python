# Advanced Programming in Python 
## Project: DNA Sequence Assembly

<div align="justify">
  
### Objective & General Information
A nucleotide is the basic building block of DNA sequences. The four nucleobases that can appear in a DNA sequence are A,C,G and T. Thus, a DNA sequence is a long chain (string) of nucleotides such as (the following sequence will be referred to it as SQ1)<br />
<br />
                  TTAATTACTCACTGGCTAATTACTCACTGGGTCACTACGCACTG<br />
<br />
In order to construct a DNA sequence of a given DNA, typically the sequence is multiplied and segmented into different lengths. These smaller segments of the DNA undergo a chemical process in order to know the order of nucleotides in each segment. The following sequences are segments of SQ1:<br />
<br />
TTAATTA ATTACTC ACTCAC TCACTGGCTAA CTAATTACTCACTGG CTGGGT GGGTCACT CACTACGCACTG<br />
<br />
These DNA sequences of the smaller segments typically overlap. Thus, they have to be aligned or connected together in order to know the correct sequence. <br />
<br />
Without a reference DNA sequence, there is no algorithm that will get the correct sequence for all the possible cases of segments. However, we are going to add some assumptions under which we guarantee to find the correct DNA sequence from the given segments.

<br />

Main Objective In this project, you will implement a program that given overlapping sequences of small DNA segments assembles these segments and returns the whole DNA sequence.


</div>
