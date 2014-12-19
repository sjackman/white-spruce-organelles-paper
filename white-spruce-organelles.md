---
title: 'Organellar Genomes of White Spruce (*Picea glauca*): Assembly and Annotation'
author:
  Shaun D Jackman, Anthony Raymond, Ben Vandervalk, Hamid Mohamadi,
  Rene Warren, Stephen Pleasance, Robin Coope, Macaire MS Yuen,
  Christopher Keeling, Carol Ritland, Jean Bousquet, Alvin Yanchuk,
  Kermit Ritland, John MacKay, Steven JM Jones, Joerg C Bohlmann,
  Inanc Birol
leftrunninghead: Jackman et al.
keywords: organelle, genome, genome sequencing, genome sequence assembly,
  white spruce, Picea glauca
abstract:
  The genome sequences of the plastid and mitochondrion of white spruce
  (*Picea glauca*) are assembled from whole genome Illumina sequencing
  data using ABySS. Whole genome sequencing data contains reads from
  both the nuclear and organellar genomes. Reads of the organellar
  genomes are abundant, because each cell contains hundreds of
  mitochondria and plastids. One lane of MiSeq data assembles the 123
  kbp plastid genome in a single contig, and one lane of HiSeq data
  assembles a putative 5.9 Mbp mitochondrial genome. The raw assembly is
  expected to be composed of organellar sequence as well as nuclear
  repeat elements. The organellar sequences are separated from the
  assembly by classifying the sequences using their length, depth of
  coverage and GC content. The genes and repeats of the plastid and
  mitochondrial genomes are annotated using MAKER-P.
---

Introduction
================================================================================

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
================================================================================

The software used in this analysis, their versions and the digital
object identifiers (DOI) of their respective publications are listed
in supplementary Table S1.

Plastid
------------------------------------------------------------

The overlapping paired-end reads were merged using ABySS-mergepairs.
These merged reads were assembled using [ABySS][abyss]. Contigs that
are putatively derived from the plastid were separated by length and
depth of coverage using thresholds chosen by inspection (see
supplementary Figure S1). These putative plastid contigs were
assembled into scaffolds using ABySS-scaffold. The assembled plastid
genome was initially annotated using [DOGMA][dogma], but DOGMA is an
interactive web application, which is not convenient for an automated
pipeline. We instead used [MAKER-P][maker] for annotation, which is
intended for automated pipelines, and used the [Norway
spruce][norwayspruce] complete plastid genome ([NC_021456][]) for
both protein-coding and non-coding gene homology evidence. The
parameters of MAKER are show in supplementary Table S2. The
inverted repeat was identified using [MUMmer][], shown in
supplementary Figure S3.

The assembled plastid genome was aligned to the Norway spruce plastid
using [BWA-MEM][bwamem]. Coverage and identity of these alignments
were calculated using the script `bam-identity` (see supplementary
materials). The two genomes were compared using [QUAST][quast] to
confirm the presence of the annotated genes of the Norway spruce
plastid in the white spruce plastid.

[NC_021456]: http://www.ncbi.nlm.nih.gov/nuccore/NC_021456

Mitochondrion
------------------------------------------------------------

ABySS-konnector was used to fill the gap between the paired-end reads
of a single lane of Illumina HiSeq sequencing of a paired-end library.
These connected paired-end reads were assembled using [ABySS][abyss].
Putative mitochondrial sequences were separated from the assembly by
their length, depth of coverage and GC content using k-means
clustering in R (see supplementary Figure S2). These putative
mitochondrial contigs were then assembled into scaffolds using
ABySS-scaffold with a single lane of Illumina HiSeq sequencing of a
mate-pair library.

The mitochondrial genome was annotated using [MAKER-P][maker]
(parameters shown in supplementary Table S3). The proteins of all
green plants (viridiplantae) with complete mitochondrial genome
sequences in NCBI GenBank, 51 species, were used for protein homology
evidence and aligned using [BLAST][blast] and [Exonerate][exonerate].
The [prince sago palm (Cycas taitungensis)
mitochondrion][ctaitungensis] ([NC_010303][]) is the closest related
species, being the only gymnosperm with a complete mitochondrial
genome. Transfer RNA (tRNA) were annotated using [tRNAscan][trnascan].
Ribosomal RNA (rRNA) were annotated using [Barrnap][barrnap]. Repeats
were identified using [RepeatMasker][repeatmasker] and RepeatModeler.

[NC_010303]: http://www.ncbi.nlm.nih.gov/nuccore/NC_010303

The putative mitochondrial sequences of white spruce were aligned
to the putative mitochondrial sequences of the
[Norway spruce][norwayspruce] using [BWA-MEM][bwamem]. Coverage and
identity of these alignments were calculated using the script
`bam-identity` (see supplementary materials).

Results
================================================================================

The assembly and annotation metrics are summarized in [Table 1][].

<a name="table-1"></a>

[Table 1]: #table-1

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
Assembled genome size           |123,266 bp      |5.92 Mbp
Number of contigs               |1 contig        |223 contigs x
Contig N50                      |123 kbp         |39 kbp x
Number of scaffolds             |1 scaffold      |61 scaffolds
Scaffold N50                    |123 kbp         |287 kbp
Largest scaffold                |123 kbp         |598 kbp
GC content                      |38.8%           |44.7%
Protein coding genes (mRNA)     |74              |54
Transfer RNA genes (tRNA)       |36              |23
Ribosomal RNA genes (rRNA)      |4               |4
Coding genes containing introns |4               |11
Introns in coding genes         |5               |15
tRNA genes containing introns   |5               |0
Identity to Norway spruce       |99.2%           |98.3%
Coverage of Norway spruce       |99.7%           |59.6%
x Permitting gaps less than 500 bp | |

Plastid
------------------------------------------------------------

The plastid genome was assembled into a single circular contig of
123,266 bp. The assembly metrics are shown in [Table 1][]. The plastid
genome contains 114 genes: 74 protein coding (mRNA) genes, 36 transfer RNA
(tRNA) genes and 4 ribosomal RNA (rRNA) genes, shown in [Figure 1][].

All protein-coding genes are single copy, except *psbI* and *ycf12*, which have
two copies each. All tRNA genes are single copy, except *trnH-GUG*, *trnI-CAU*,
*trnS-GCU* and *trnT-GGU*, which have two copies each. All rRNA genes are
single copy.

The protein-coding genes
*atpF*, *petB*, *petD*, *rpl2*, *rpl16*, *rpoC1* and *rps12*
each contain one intron, and *ycf3* contains two introns.
The tRNA genes
*trnA-UGC*, *trnG-GCC*, *trnI-GAU*, *trnK-UUU*, *trnL-UAA* and *trnV-UAC*
each contain one intron.
The rRNA genes are not spliced.

The first and smallest exons of the genes *petB*, *petD* and *rpl16* are 6, 8
and 9 bp respectively. These genes likely belong to
[polycistronic transcripts][] of their respective protein complexes, but the
short size of their initial exons make them difficult to annotate all the same.
The initial exons of these genes were added to their annotations manually.

The gene *rps12* of a plastid genome is typically [trans-spliced][], which makes
it difficult to annotate using [MAKER][]. It is composed of three exons and one
cis-spliced intron. It required manually editing the gene annotation to
incorporate trans-splicing in the gene model.

Each copy of the inverted repeat (IR) is 445 bp in size, much smaller than most
plants, but typical of *Pinaceae* ([Lin, 2010][]). Unlike most inverted repeats,
which are typically identical, the two copies differ by a single base. The IR
contains a single gene, the tRNA *trnI-CAU*.

All 114 genes of the Norway spruce plastid genome ([NC_021456][]) are present
in the white spruce plastid genome. The genomes of the white spruce plastid and
Norway spruce plastid show perfect gene synteny with no structural
rearrangements.

<a name="figure-1"></a>

[Figure 1]: #figure-1

![The annotated plastid genome, which was annotated using
[MAKER-P][maker] and plotted using [OGDRAW][ogdraw].
](figure/plastid-annotation.png)

Mitochondrion
------------------------------------------------------------

The mitochondrial genome was assembled into 61 scaffolds (223 contigs,
permitting gaps less than 500 bp) with a scaffold N50 of 287 kbp
(contig N50 of 39 kbp). The largest scaffold is 598 kbp. The assembly
metrics are shown in [Table 1][].

The mitochondrial genome contains 54 protein coding (mRNA) genes, 23
transfer RNA (tRNA) genes and 4 ribosomal RNA (rRNA) genes, shown in
[Figure 2][]. The coding genes compose 50 kbp (<1%) of the genome,
shown in [Figure 3][].

The protein-coding genes
*atp8*, *cox1*, *matR*, *nad7*, *rpl10*, *rps1* and *rps4*
each contain one intron, and
*ccmFn*, *nad4*, *nad5* and *rps3-2*
each contain two introns.

The putative mitochondrial sequences of white spruce and Norway spruce
show high sequence similarity, over 98% nucleotide identity, but only
60% of the Norway spruce putative mitochondrial sequences are covered
by alignments of the white spruce sequences.

Repeats compose 400 kbp (~7%) of the mitochondrial genome. Simple
repeats, the LINE Jockey and the LTR Copia and Gypsy are the most
common repeats, shown in [Figure 4][].

<a name="figure-2"></a>

[Figure 2]: #figure-2

![The annotated mitochondrial genome, which was annotated using
[MAKER-P][maker] and plotted using [OGDRAW][ogdraw].
](figure/mt-annotation.png)

<a name="figure-3"></a>

[Figure 3]: #figure-3

![The sizes of the mitochondrial genes, grouped by family.
](figure/mt-genes.png)

<a name="figure-4"></a>

[Figure 4]: #figure-4

![Repetitive sequence of the mitochondrial genome.
](figure/mt-repeats.png)

Conclusion
================================================================================

One lane of MiSeq sequencing of whole genome DNA is sufficient to
assemble the 123 kbp complete plastid genome, and one lane of HiSeq
sequencing of whole genome DNA is sufficient to assemble a draft 5.9 Mbp
mitochondrial genome of white spruce. Scaffold contiguity is improved
with additional mate-pair library sequencing. The resulting assembly
of whole genome sequencing data is composed of organellar sequences as
well as high-copy-number nuclear repeat elements. The mitochondrial
sequences are separated using a k-means classifier based on their
length, depth of coverage and GC content of the sequences.

The white spruce plastid genome shows no structural rearrangements
when compared with Norway spruce. The mitochondrial genome in contrast
shows much structural rearrangement, though more work is needed to
determine what is due to the draft nature of these mitochondrial
assemblies and what is true structural rearrangement.

The protein coding gene content of the mitochondrial genome is quite
sparse, with 54 protein coding genes in 5.9 Mbp, in comparison to the
plastid genome, with 74 protein coding genes in 123 kbp. 7% of the
mitochondrial genome is composed of repeats, and <1% is composed of
coding genes. A significant portion, over 90%, of the unusually large
size of the white spruce mitochondrial genome is yet unexplained.

Acknowledgements
================================================================================

Shaun Jackman would like to thank his supervisors Inanc Birol and
Joerg Bohlmann for their guidance in the preparation of this
manuscript, and Carson Holt for being exceedingly responsive and
helpful in tweaking MAKER.

References
================================================================================

| [ABySS: a parallel assembler for short read sequence data][abyss]
| [Horizontal Transfer of Entire Genomes via Mitochondrial Fusion in the Angiosperm Amborella][amborellamt]
| [The Amborella Genome and the Evolution of Flowering Plants][amborellanuc]
| [Genomic Clues to the Ancestral Flowering Plant][amborellaperspective]
| [barrnap 0.4.2 - rapid ribosomal RNA prediction][barrnap]
| [Basic Local Alignment Search Tool][blast]
| [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM][bwamem]
| [CGAP: a new comprehensive platform for the comparative analysis of chloroplast genomes][cgap]
| [The Mitochondrial Genome of the Gymnosperm Cycas taitungensis Contains a Novel Family of Short Interspersed Elements, Bpu Sequences, and Abundant RNA Editing Sites][ctaitungensis]
| [Automatic annotation of organellar genomes with DOGMA][dogma]
| [Automated generation of heuristics for biological sequence comparison][exonerate]
| [MAKER-P: a tool-kit for the rapid creation, management, and quality control of plant genome annotations][maker]
| [The Norway spruce genome sequence and conifer genome evolution][norwayspruce]
| [OrganellarGenomeDRAW (OGDRAW): a tool for the easy generation of high-quality custom graphical maps of plastid and mitochondrial genomes][ogdraw]
| [Comparative chloroplast genomics reveals the evolution of *Pinaceae* genera and subfamilies][pinaceae]
| [QUAST: quality assessment tool for genome assemblies][quast]
| [Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-3.0. 1996-2010 http://www.repeatmasker.org][repeatmasker]
| [The Sequence Alignment/Map format and SAMtools][samtools]
| [tRNAscan-SE: A Program for Improved Detection of Transfer RNA Genes in Genomic Sequence][trnascan]
| [Assembling the 20Gb white spruce (*Picea glauca*) genome from whole-genome shotgun sequencing data][whitespruce]
| [Gremme, G., Steinbiss, S., & Kurtz, S. (2013)][GenomeTools]
  GenomeTools: a comprehensive software library for efficient processing of structured genome annotations.
  *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, 10(3), 645-656.
| [Kurtz, S., Phillippy, A., Delcher, A. L., Smoot, M., Shumway, M., Antonescu, C., & Salzberg, S. L. (2004)][MUMmer]
  Versatile and open software for comparing large genomes.
  *Genome biology*, 5(2), R12.
| [Lin, C. P., Huang, J. P., Wu, C. S., Hsu, C. Y., & Chaw, S. M. (2010)][Lin, 2010]
  Comparative chloroplast genomics reveals the evolution of Pinaceae genera and subfamilies.
  *Genome biology and evolution*, 2, 504-517.
| [Barkan, A. (1988).][polycistronic transcripts]
  Proteins encoded by a complex chloroplast transcription unit are each translated from both monocistronic and polycistronic mRNAs.
  *The EMBO journal*, 7(9), 2637.
| [Hildebrand, M., Hallick, R. B., Passavant, C. W., & Bourque, D. P. (1988)][trans-spliced]
  Trans-splicing in chloroplasts: the rps 12 loci of Nicotiana tabacum.
  *Proceedings of the National Academy of Sciences*, 85(2), 372-376.

[abyss]: http://genome.cshlp.org/content/19/6/1117
[amborellamt]: http://www.sciencemag.org/content/342/6165/1468
[amborellanuc]: http://www.sciencemag.org/content/342/6165/1241089
[amborellaperspective]: http://www.sciencemag.org/content/342/6165/1456
[barrnap]: http://www.vicbioinformatics.com/software.barrnap.shtml
[blast]: http://www.sciencedirect.com/science/article/pii/S0022283605803602
[bwamem]: http://arxiv.org/pdf/1303.3997.pdf
[cgap]: http://www.biomedcentral.com/1471-2105/14/95/abstract
[ctaitungensis]: http://mbe.oxfordjournals.org/content/25/3/603.short
[dogma]: http://bioinformatics.oxfordjournals.org/content/20/17/3252
[exonerate]: http://www.biomedcentral.com/1471-2105/6/31
[GenomeTools]: http://www.computer.org/csdl/trans/tb/2013/03/ttb2013030645-abs.html
[Lin, 2010]: http://gbe.oxfordjournals.org/content/2/504
[maker]: http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144
[MUMmer]: http://genomebiology.com/content/5/2/R12
[norwayspruce]: http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html
[ogdraw]: http://nar.oxfordjournals.org/content/41/W1/W575
[pinaceae]: http://gbe.oxfordjournals.org/content/2/504
[polycistronic transcripts]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC457051/
[quast]: http://bioinformatics.oxfordjournals.org/content/29/8/1072
[repeatmasker]: http://www.repeatmasker.org/
[samtools]: http://bioinformatics.oxfordjournals.org/content/25/16/2078
[trans-spliced]: http://www.pnas.org/content/85/2/372.short
[trnascan]: http://nar.oxfordjournals.org/content/25/5/0955
[whitespruce]: http://bioinformatics.oxfordjournals.org/content/29/12/1492
