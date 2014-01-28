Organellar Genomes of White Spruce (*Picea glauca*)
===================================================

Shaun D Jackman, Anthony Raymond, Ben Vandervalk, Hamid Mohamadi,
René Warren, Stephen Pleasance, Robin Coope, Macaire MS Yuen,
Christopher Keeling, Carol Ritland, Jean Bousquet, Alvin Yanchuk,
Kermit Ritland, John MacKay, Steven JM Jones, Joerg C Bohlmann,
İnanç Birol

Abstract
========

The genome sequences of the plastid and mitochondrion of white spruce
(*Picea glauca*) are assembled from whole genome Illumina sequencing
data using ABySS and aligned to the Norway spruce (*Picea abies*)
using BWA-MEM. The putative mitochondrial sequences are classified
using *k*-means clustering in R. The plastid genome is 120 kbp, and
the putative mitochondrial genome is 6 Mbp.

Introduction
============

The SMarTForests project published the draft genome sequence of the
20 gigabase [white spruce (*Picea glauca*) genome][whitespruce], seven
times the size of the human genome, sequenced using the Illumina HiSeq
and MiSeq sequencing platforms. Whole genome sequencing data contains
reads originating from both the nuclear and organellar genomes.
Whereas one copy of the diploid nuclear genome is found in each cell,
hundreds of organelles are present, and thus hundreds of copies of the
organellar genomes. This abundance results in an overrepresentation of
the organellar genomes in whole genome sequencing data.

Assembling a single lane of whole genome sequencing data using
[ABySS][abyss] yields an assembly composed of organellar sequences and
nuclear repeat elements. The assembled sequences that originate from
the organellar genomes are separated from those of nuclear origin by
classifying the sequences using their length, depth of coverage and GC
content. The organellar genomes of white spruce are compared to those
of [Norway spruce (*Picea abies*)][norwayspruce].

Methods
=======

Plastid
-------

+ Merged the overlapping paired-end reads (ABySS-mergepairs)
+ Assembled with [ABySS][abyss] 1.3.7
+ Separated putative plastid by length and depth of coverage
+ Scaffolded with mate-pair reads
+ Aligned to the Norway spruce plastid (124 kbp) using [BWA-MEM][bwamem]
+ All 117 annotated genes are covered, 114 full and 3 partial

Mitochondrion
-------------

+ Filled the gap between the paired-end reads using a Bloom filter de Bruijn Graph (ABySS-connectpairs)
+ Assembled with [ABySS][abyss] 1.3.7
+ Separated putative mitochondrial sequence by length, depth of coverage and GC content
+ Scaffolded with mate-pair reads
+ Aligned to the Norway spruce putative mitochondrion (5.5 Mbp) using [BWA-MEM][bwamem]

Results
=======

Sequencing/assembly metrics     |Plastid         |Mitochondrion
------------------------------- |--------------- |-------------
Number of lanes                 |1 MiSeq lane    |1 HiSeq lane
Number of read pairs            |4.7 million     |133 million
Read length                     |300 bp          |150 bp
Number of merged reads          |3.0 million     |1.4 million
Median merged read length       |492 bp          |465 bp
Number of assembled reads       |21 thousand     |377 thousand
Proportion of organellar reads  |1/140 or 0.7%   |1/350 or 0.3%
Depth of coverage               |80x             |30x
Assembled genome size           |125 kbp         |6.0 Mbp
Number of contigs†              |6 contigs       |223 contigs
Contig N50†                     |70 kbp          |39 kbp
Number of scaffolds             |1 scaffold      |78 scaffolds
Scaffold N50                    |125 kbp         |157 kbp
Largest scaffold                |125 kbp         |519 kbp
Identity to Norway spruce       |99.2%           |98.3%
Coverage of Norway spruce       |98.8%           |59.6%

† Permitting gaps less than 500 bp

Conclusion
==========

+ One lane of MiSeq data assembles the 124 kbp plastid genome
+ One lane of HiSeq data assembles the putative 6 Mbp mitochondrial
  genome
+ The assembly is composed of organellar sequences as well as nuclear
  repeat elements
+ The organellar sequences are separated from the nuclear sequences by
  classifying the sequences using their length, depth of coverage and
  GC content
+ These white spruce sequences are aligned to the
  [Norway spruce (*Picea abies*)][norwayspruce] complete plastid
  genome (NC_021456) and putative mitochondrial sequences with 99.2%
  and 98.3% identity, respectively.

References
==========

+ [ABySS: a parallel assembler for short read sequence data][abyss]
+ [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM][bwamem]
+ [Assembling the 20Gb white spruce (*Picea glauca*) genome from whole-genome shotgun sequencing data][whitespruce]
+ [The Norway spruce genome sequence and conifer genome evolution][norwayspruce]
+ [CGAP: a new comprehensive platform for the comparative analysis of chloroplast genomes][cgap]
+ [Genomic Clues to the Ancestral Flowering Plant][amborellaperspective]
+ [The Amborella Genome and the Evolution of Flowering Plants][amborellanuc]
+ [Horizontal Transfer of Entire Genomes via Mitochondrial Fusion in the Angiosperm Amborella][amborellamt]
+ [Comparative chloroplast genomics reveals the evolution of *Pinaceae* genera and subfamilies][pinaceae]
+ [Automatic annotation of organellar genomes with DOGMA][dogma]
+ [OrganellarGenomeDRAW (OGDRAW): a tool for the easy generation of high-quality custom graphical maps of plastid and mitochondrial genomes][ogdraw]

[abyss]: http://genome.cshlp.org/content/19/6/1117
[bwamem]: http://arxiv.org/pdf/1303.3997.pdf
[whitespruce]: http://bioinformatics.oxfordjournals.org/content/29/12/1492
[norwayspruce]: http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html
[cgap]: http://www.biomedcentral.com/1471-2105/14/95/abstract
[amborellaperspective]: http://www.sciencemag.org/content/342/6165/1456
[amborellanuc]: http://www.sciencemag.org/content/342/6165/1241089
[amborellamt]: http://www.sciencemag.org/content/342/6165/1468
[pinaceae]: http://gbe.oxfordjournals.org/content/2/504
[dogma]: http://bioinformatics.oxfordjournals.org/content/20/17/3252
[ogdraw]: http://nar.oxfordjournals.org/content/41/W1/W575
