Organellar Genomes of White Spruce (*Picea glauca*): Assembly and Annotation
============================================================================

Shaun D Jackman, Anthony Raymond, Ben Vandervalk, Hamid Mohamadi,
René Warren, Stephen Pleasance, Robin Coope, Macaire MS Yuen,
Christopher Keeling, Carol Ritland, Jean Bousquet, Alvin Yanchuk,
Kermit Ritland, John MacKay, Steven JM Jones, Joerg C Bohlmann,
İnanç Birol

Abstract
========

The genome sequences of the plastid and mitochondrion of white spruce
(*Picea glauca*) are assembled from whole genome Illumina sequencing
data using ABySS. Whole genome sequencing data contains reads from
both the nuclear and organellar genomes. Reads of the organellar
genomes are abundant, because each cell contains hundreds of
mitochondria and plastids. One lane of MiSeq data assembles the 124
kbp plastid genome in a single scaffold, and one lane of HiSeq data
assembles a putative 6 Mbp mitochondrial genome. The raw assembly is
expected to be composed of organellar sequence as well as nuclear
repeat elements. The organellar sequences are separated from the
assembly by classifying the sequences using their length, depth of
coverage and GC content. The genes and repeats of the plastid and
mitochondrial genomes are annotated using MAKER-P.

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

### Assembly

The overlapping paired-end reads were merged using ABySS-mergepairs
(distributed with ABySS). These merged reads were assembled using
[ABySS][abyss] 1.3.7. Contigs that are putatively derived from the
plastid were separated by length and depth of coverage using
thresholds chosen by inspection (see supplementary figure 1). These
putative plastid contigs were scaffolded using ABySS-scaffold
(distributed with ABySS).

### Annotation

The assembled plastid genome was annotated using [DOGMA][dogma] with
default parameters (shown in supplementary material).

### Comparative genomics

The assembled plastid genome was aligned to the
[Norway spruce][norwayspruce] complete plastid genome ([NC_021456][])
using [BWA-MEM][bwamem]. Coverage and identity of these alignments
were calculated using the script `bam-identity` (see supplementary
materials). The two genomes were compared using [QUAST][quast] to
confirm the presence of the annotated genes of the Norway spruce
plastid in the white spruce plastid.

[NC_021456]: http://www.ncbi.nlm.nih.gov/nuccore/NC_021456

Mitochondrion
-------------

+ Filled the gap between the paired-end reads using a Bloom filter de Bruijn Graph (ABySS-connectpairs)
+ Assembled with [ABySS][abyss] 1.3.7
+ Separated putative mitochondrial sequence by length, depth of coverage and GC content
+ Scaffolded with mate-pair reads
+ Aligned to the Norway spruce putative mitochondrion (5.5 Mbp) using [BWA-MEM][bwamem]

The mitochondrial genome was annotated using [MAKER][maker]
(parameters shown in supplementary material). The Cycas taitungensis
mitochondrion (NC_010303.1) was used for protein homology evidence,
aligned using [BLAST][blast]. The tRNA genes were discovered using
[tRNAscan][trnascan]. Repeats were identified using
[RepeatMasker][repeatmasker].

Results
=======

**Table 1**: Sequencing, assembly and alignment metrics of the white
spruce organellar genomes

Metric                          |Plastid         |Mitochondrion
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

Plastid
-------

### Assembly

The plastid genome was assembled into a single circular scaffold with
seven gaps. The assembly metrics are shown in table 1.

### Annotation

The plastid genome has 92 protein coding (mRNA) genes, 41 transfer RNA
(tRNA) genes and 4 ribosomal RNA (rRNA) genes, shown in figure 1.

![Figure 1: Plastid genes](plastid-annotation.png)

**Figure 1**: The annotated plastid genome, which was annotated using
[DOGMA][dogma] and plotted using [OGDRAW][ogdraw].

### Comparative genomics

The genomes of the white spruce plastid and Norway spruce plastid show
perfect synteny with no structural rearrangements. All 117 genes (114
in full and 3 partial) of the Norway spruce plastid genome
([NC_021456][]) are present in the white spruce plastid genome.

Mitochondrion
-------------

The mitochondrial genome contains 60 protein coding genes and 27
transfer RNA (tRNA) genes, shown in figure 2.

![Figure 2: Mitochondrial genes](mt-annotation.png)

**Figure 2**: The annotated mitochondrial genome, plotted using [OGDRAW][ogdraw].

Simple repeats and the LTR Copia and Gypsy are the most common
repeats found in the mitochondrial genome, shown in figure 3.

![Figure 3: Mitochondrial repeats](mt-repeats.png)

**Figure 3**: Repetitive sequence of the mitochondrial genome

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
  genome ([NC_021456][]) and putative mitochondrial sequences with
  99.2% and 98.3% identity, respectively.

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
+ [MAKER-P: a tool-kit for the rapid creation, management, and quality control of plant genome annotations][maker]
+ [Basic Local Alignment Search Tool][blast]
+ [tRNAscan-SE: A Program for Improved Detection of Transfer RNA Genes in Genomic Sequence][trnascan]
+ [Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-3.0. 1996-2010 http://www.repeatmasker.org][repeatmasker]
+ [QUAST: quality assessment tool for genome assemblies][quast]

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
[maker]: http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144
[blast]: http://www.sciencedirect.com/science/article/pii/S0022283605803602
[trnascan]: http://nar.oxfordjournals.org/content/25/5/0955
[repeatmasker]: http://www.repeatmasker.org/
[quast]: http://bioinformatics.oxfordjournals.org/content/29/8/1072

Supplementary material
======================

DOGMA parameters
----------------

The following parameters were used for the annotation of the plastid
using [DOGMA][dogma].

Parameter                                        | Value
-------------------------------------------------|-------------
Genome type                                      | Chloroplast
Gapped alignment?                                | Yes
Genetic Code for Blastx                          | 11 Plant plastid
Percent identity cutoff for protein coding genes | 60
Percent identity cutoff for RNAs                 | 80
E-value                                          | 1e-5
Number of blast hits to return                   | 5

Mitochondrion MAKER parameters
------------------------------

```
genome=pg29mt-concat.fa
organism_type=eukaryotic
protein=NC_010303.faa
model_org=all
rmlib=PICEAGLAUCA_rpt2.0.fa
repeat_protein=/usr/local/Cellar/maker/2.31/libexec/data/te_proteins.fasta
protein2genome=1
trna=1
est_forward=1
single_exon=1
```
