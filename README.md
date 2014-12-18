Shaun D Jackman, Anthony Raymond, Ben Vandervalk, Hamid Mohamadi, René
Warren, Stephen Pleasance, Robin Coope, Macaire MS Yuen, Christopher
Keeling, Carol Ritland, Jean Bousquet, Alvin Yanchuk, Kermit Ritland,
John MacKay, Steven JM Jones, Joerg C Bohlmann, İnanç Birol

Abstract
========

The genome sequences of the plastid and mitochondrion of white spruce
(*Picea glauca*) are assembled from whole genome Illumina sequencing
data using ABySS. Whole genome sequencing data contains reads from both
the nuclear and organellar genomes. Reads of the organellar genomes are
abundant, because each cell contains hundreds of mitochondria and
plastids. One lane of MiSeq data assembles the 123 kbp plastid genome in
a single contig, and one lane of HiSeq data assembles a putative 5.9 Mbp
mitochondrial genome. The raw assembly is expected to be composed of
organellar sequence as well as nuclear repeat elements. The organellar
sequences are separated from the assembly by classifying the sequences
using their length, depth of coverage and GC content. The genes and
repeats of the plastid and mitochondrial genomes are annotated using
MAKER-P.

Introduction
============

The SMarTForests project published the draft genome sequence of the 20
gigabase [white spruce (*Picea glauca*)
genome](http://bioinformatics.oxfordjournals.org/content/29/12/1492),
seven times the size of the human genome, sequenced using the Illumina
HiSeq and MiSeq sequencing platforms. Whole genome sequencing data
contains reads originating from both the nuclear and organellar genomes.
Whereas one copy of the diploid nuclear genome is found in each cell,
hundreds of organelles are present, and thus hundreds of copies of the
organellar genomes. This abundance results in an overrepresentation of
the organellar genomes in whole genome sequencing data.

Assembling a single lane of whole genome sequencing data using
[ABySS](http://genome.cshlp.org/content/19/6/1117) yields an assembly
composed of organellar sequences and nuclear repeat elements. The
assembled sequences that originate from the organellar genomes are
separated from those of nuclear origin by classifying the sequences
using their length, depth of coverage and GC content. The organellar
genomes of white spruce are compared to those of [Norway spruce (*Picea
abies*)](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html).

Methods
=======

The software used in this analysis, their versions and the digital
object identifiers (DOI) of their respective publications are listed in
supplementary [Table S1](#table-s1-software).

Plastid
-------

The overlapping paired-end reads were merged using ABySS-mergepairs.
These merged reads were assembled using
[ABySS](http://genome.cshlp.org/content/19/6/1117). Contigs that are
putatively derived from the plastid were separated by length and depth
of coverage using thresholds chosen by inspection (see supplementary
[Figure S1](#figure-s1-classify-plastid-sequences)). These putative
plastid contigs were assembled into scaffolds using ABySS-scaffold. The
assembled plastid genome was initially annotated using
[DOGMA](http://bioinformatics.oxfordjournals.org/content/20/17/3252),
but DOGMA is an interactive web application, which is not convenient for
an automated pipeline. We instead used
[MAKER-P](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144)
for annotation, which is intended for automated pipelines, and used the
[Norway
spruce](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html)
complete plastid genome
([NC\_021456](http://www.ncbi.nlm.nih.gov/nuccore/NC_021456)) for both
protein-coding and non-coding gene homology evidence. The parameters of
MAKER are show in supplementary [Table
S2](#table-s2-plastid-maker-parameters). The inverted repeat was
identified using [MUMmer](http://genomebiology.com/content/5/2/R12),
shown in supplementary [Figure
S3](#figure-s3-mummer-alignment-of-the-plastid).

The assembled plastid genome was aligned to the Norway spruce plastid
using [BWA-MEM](http://arxiv.org/pdf/1303.3997.pdf). Coverage and
identity of these alignments were calculated using the script
`bam-identity` (see supplementary materials). The two genomes were
compared using
[QUAST](http://bioinformatics.oxfordjournals.org/content/29/8/1072) to
confirm the presence of the annotated genes of the Norway spruce plastid
in the white spruce plastid.

Mitochondrion
-------------

ABySS-konnector was used to fill the gap between the paired-end reads of
a single lane of Illumina HiSeq sequencing of a paired-end library.
These connected paired-end reads were assembled using
[ABySS](http://genome.cshlp.org/content/19/6/1117). Putative
mitochondrial sequences were separated from the assembly by their
length, depth of coverage and GC content using k-means clustering in R
(see supplementary [Figure
S2](#figure-s2-classify-mitochondrial-sequences)). These putative
mitochondrial contigs were then assembled into scaffolds using
ABySS-scaffold with a single lane of Illumina HiSeq sequencing of a
mate-pair library.

The mitochondrial genome was annotated using
[MAKER-P](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144)
(parameters shown in supplementary [Table
S3](#table-s3-mitochondrion-maker-parameters)). The proteins of all
green plants (viridiplantae) with complete mitochondrial genome
sequences in NCBI GenBank, 51 species, were used for protein homology
evidence and aligned using
[BLAST](http://www.sciencedirect.com/science/article/pii/S0022283605803602)
and [Exonerate](http://www.biomedcentral.com/1471-2105/6/31). The
[prince sago palm (Cycas taitungensis)
mitochondrion](http://mbe.oxfordjournals.org/content/25/3/603.short)
([NC\_010303](http://www.ncbi.nlm.nih.gov/nuccore/NC_010303)) is the
closest related species, being the only gymnosperm with a complete
mitochondrial genome. Transfer RNA (tRNA) were annotated using
[tRNAscan](http://nar.oxfordjournals.org/content/25/5/0955). Ribosomal
RNA (rRNA) were annotated using
[Barrnap](http://www.vicbioinformatics.com/software.barrnap.shtml).
Repeats were identified using
[RepeatMasker](http://www.repeatmasker.org/) and RepeatModeler.

The putative mitochondrial sequences of white spruce were aligned to the
putative mitochondrial sequences of the [Norway
spruce](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html)
using [BWA-MEM](http://arxiv.org/pdf/1303.3997.pdf). Coverage and
identity of these alignments were calculated using the script
`bam-identity` (see supplementary materials).

Results
=======

The assembly and annotation metrics are summarized in [Table
1](#table-1).

<a name="table-1"></a>

**Table 1**: Sequencing, assembly and alignment metrics of the white
spruce organellar genomes

<table>
<thead>
<tr class="header">
<th align="left">Metric</th>
<th align="left">Plastid</th>
<th align="left">Mitochondrion</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Number of lanes</td>
<td align="left">1 MiSeq lane</td>
<td align="left">1 HiSeq lane</td>
</tr>
<tr class="even">
<td align="left">Number of read pairs</td>
<td align="left">4.7 million</td>
<td align="left">133 million</td>
</tr>
<tr class="odd">
<td align="left">Read length</td>
<td align="left">300 bp</td>
<td align="left">150 bp</td>
</tr>
<tr class="even">
<td align="left">Number of merged reads</td>
<td align="left">3.0 million</td>
<td align="left">1.4 million</td>
</tr>
<tr class="odd">
<td align="left">Median merged read length</td>
<td align="left">492 bp</td>
<td align="left">465 bp</td>
</tr>
<tr class="even">
<td align="left">Number of assembled reads</td>
<td align="left">21 thousand</td>
<td align="left">377 thousand</td>
</tr>
<tr class="odd">
<td align="left">Proportion of organellar reads</td>
<td align="left">1/140 or 0.7%</td>
<td align="left">1/350 or 0.3%</td>
</tr>
<tr class="even">
<td align="left">Depth of coverage</td>
<td align="left">80x</td>
<td align="left">30x</td>
</tr>
<tr class="odd">
<td align="left">Assembled genome size</td>
<td align="left">123,266 bp</td>
<td align="left">5.92 Mbp</td>
</tr>
<tr class="even">
<td align="left">Number of contigs</td>
<td align="left">1 contig</td>
<td align="left">223 contigs †</td>
</tr>
<tr class="odd">
<td align="left">Contig N50</td>
<td align="left">123 kbp</td>
<td align="left">39 kbp †</td>
</tr>
<tr class="even">
<td align="left">Number of scaffolds</td>
<td align="left">1 scaffold</td>
<td align="left">61 scaffolds</td>
</tr>
<tr class="odd">
<td align="left">Scaffold N50</td>
<td align="left">123 kbp</td>
<td align="left">287 kbp</td>
</tr>
<tr class="even">
<td align="left">Largest scaffold</td>
<td align="left">123 kbp</td>
<td align="left">598 kbp</td>
</tr>
<tr class="odd">
<td align="left">GC content</td>
<td align="left">38.8%</td>
<td align="left">44.7%</td>
</tr>
<tr class="even">
<td align="left">Protein coding genes (mRNA)</td>
<td align="left">74</td>
<td align="left">54</td>
</tr>
<tr class="odd">
<td align="left">Transfer RNA genes (tRNA)</td>
<td align="left">36</td>
<td align="left">23</td>
</tr>
<tr class="even">
<td align="left">Ribosomal RNA genes (rRNA)</td>
<td align="left">4</td>
<td align="left">4</td>
</tr>
<tr class="odd">
<td align="left">Coding genes containing introns</td>
<td align="left">4</td>
<td align="left">11</td>
</tr>
<tr class="even">
<td align="left">Introns in coding genes</td>
<td align="left">5</td>
<td align="left">15</td>
</tr>
<tr class="odd">
<td align="left">tRNA genes containing introns</td>
<td align="left">5</td>
<td align="left">0</td>
</tr>
<tr class="even">
<td align="left">Identity to Norway spruce</td>
<td align="left">99.2%</td>
<td align="left">98.3%</td>
</tr>
<tr class="odd">
<td align="left">Coverage of Norway spruce</td>
<td align="left">99.7%</td>
<td align="left">59.6%</td>
</tr>
</tbody>
</table>

† Permitting gaps less than 500 bp

Plastid
-------

The plastid genome was assembled into a single circular contig of
123,266 bp. The assembly metrics are shown in [Table 1](#table-1). The
plastid genome contains 114 genes: 74 protein coding (mRNA) genes, 36
transfer RNA (tRNA) genes and 4 ribosomal RNA (rRNA) genes, shown in
[Figure 1](#figure-1).

All protein-coding genes are single copy, except *psbI* and *ycf12*,
which have two copies each. All tRNA genes are single copy, except
*trnH-GUG*, *trnI-CAU*, *trnS-GCU* and *trnT-GGU*, which have two copies
each. All rRNA genes are single copy.

The protein-coding genes *atpF*, *petB*, *petD*, *rpl2*, *rpl16*,
*rpoC1* and *rps12* each contain one intron, and *ycf3* contains two
introns. The tRNA genes *trnA-UGC*, *trnG-GCC*, *trnI-GAU*, *trnK-UUU*,
*trnL-UAA* and *trnV-UAC* each contain one intron. The rRNA genes are
not spliced.

The first and smallest exons of the genes *petB*, *petD* and *rpl16* are
6, 8 and 9 bp respectively. These genes likely belong to [polycistronic
transcripts](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC457051/) of
their respective protein complexes, but the short size of their initial
exons make them difficult to annotate all the same. The initial exons of
these genes were added to their annotations manually.

The gene *rps12* of a plastid genome is typically
[trans-spliced](http://www.pnas.org/content/85/2/372.short), which makes
it difficult to annotate using
[MAKER](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144).
It is composed of three exons and one cis-spliced intron. It required
manually editing the gene annotation to incorporate trans-splicing in
the gene model.

Each copy of the inverted repeat (IR) is 445 bp in size, much smaller
than most plants, but typical of *Pinaceae* ([Lin,
2010](http://gbe.oxfordjournals.org/content/2/504)). Unlike most
inverted repeats, which are typically identical, the two copies differ
by a single base. The IR contains a single gene, the tRNA *trnI-CAU*.

All 114 genes of the Norway spruce plastid genome
([NC\_021456](http://www.ncbi.nlm.nih.gov/nuccore/NC_021456)) are
present in the white spruce plastid genome. The genomes of the white
spruce plastid and Norway spruce plastid show perfect gene synteny with
no structural rearrangements.

<a name="figure-1"></a>

![The annotated plastid genome, which was annotated using
[MAKER-P](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144)
and plotted using
[OGDRAW](http://nar.oxfordjournals.org/content/41/W1/W575).](figure/plastid-annotation.png)

Mitochondrion
-------------

The mitochondrial genome was assembled into 61 scaffolds (223 contigs,
permitting gaps less than 500 bp) with a scaffold N50 of 287 kbp (contig
N50 of 39 kbp). The largest scaffold is 598 kbp. The assembly metrics
are shown in [Table 1](#table-1).

The mitochondrial genome contains 54 protein coding (mRNA) genes, 23
transfer RNA (tRNA) genes and 4 ribosomal RNA (rRNA) genes, shown in
[Figure 2](#figure-2). The coding genes compose 50 kbp (\<1%) of the
genome, shown in [Figure 3](#figure-3).

The protein-coding genes *atp8*, *cox1*, *matR*, *nad7*, *rpl10*, *rps1*
and *rps4* each contain one intron, and *ccmFn*, *nad4*, *nad5* and
*rps3-2* each contain two introns.

The putative mitochondrial sequences of white spruce and Norway spruce
show high sequence similarity, over 98% nucleotide identity, but only
60% of the Norway spruce putative mitochondrial sequences are covered by
alignments of the white spruce sequences.

Repeats compose 400 kbp (~7%) of the mitochondrial genome. Simple
repeats, the LINE Jockey and the LTR Copia and Gypsy are the most common
repeats, shown in [Figure 4](#figure-4).

<a name="figure-2"></a>

![The annotated mitochondrial genome, which was annotated using
[MAKER-P](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144)
and plotted using
[OGDRAW](http://nar.oxfordjournals.org/content/41/W1/W575).](figure/mt-annotation.png)

<a name="figure-3"></a>

![The sizes of the mitochondrial genes, grouped by
family.](figure/mt-genes.png)

<a name="figure-4"></a>

![Repetitive sequence of the mitochondrial
genome.](figure/mt-repeats.png)

Conclusion
==========

One lane of MiSeq sequencing of whole genome DNA is sufficient to
assemble the 123 kbp complete plastid genome, and one lane of HiSeq
sequencing of whole genome DNA is sufficient to assemble a draft 5.9 Mbp
mitochondrial genome of white spruce. Scaffold contiguity is improved
with additional mate-pair library sequencing. The resulting assembly of
whole genome sequencing data is composed of organellar sequences as well
as high-copy-number nuclear repeat elements. The mitochondrial sequences
are separated using a k-means classifier based on their length, depth of
coverage and GC content of the sequences.

The white spruce plastid genome shows no structural rearrangements when
compared with Norway spruce. The mitochondrial genome in contrast shows
much structural rearrangement, though more work is needed to determine
what is due to the draft nature of these mitochondrial assemblies and
what is true structural rearrangement.

The protein coding gene content of the mitochondrial genome is quite
sparse, with 54 protein coding genes in 5.9 Mbp, in comparison to the
plastid genome, with 74 protein coding genes in 123 kbp. 7% of the
mitochondrial genome is composed of repeats, and \<1% is composed of
coding genes. A significant portion, over 90%, of the unusually large
size of the white spruce mitochondrial genome is yet unexplained.

Acknowledgements
================

Shaun Jackman would like to thank his supervisors Inanc Birol and Joerg
Bohlmann for their guidance in the preparation of this manuscript, and
Carson Holt for being exceedingly responsive and helpful in tweaking
MAKER.

References
==========

-   [ABySS: a parallel assembler for short read sequence
    data](http://genome.cshlp.org/content/19/6/1117)
-   [Horizontal Transfer of Entire Genomes via Mitochondrial Fusion in
    the Angiosperm
    Amborella](http://www.sciencemag.org/content/342/6165/1468)
-   [The Amborella Genome and the Evolution of Flowering
    Plants](http://www.sciencemag.org/content/342/6165/1241089)
-   [Genomic Clues to the Ancestral Flowering
    Plant](http://www.sciencemag.org/content/342/6165/1456)
-   [barrnap 0.4.2 - rapid ribosomal RNA
    prediction](http://www.vicbioinformatics.com/software.barrnap.shtml)
-   [Basic Local Alignment Search
    Tool](http://www.sciencedirect.com/science/article/pii/S0022283605803602)
-   [Aligning sequence reads, clone sequences and assembly contigs with
    BWA-MEM](http://arxiv.org/pdf/1303.3997.pdf)
-   [CGAP: a new comprehensive platform for the comparative analysis of
    chloroplast
    genomes](http://www.biomedcentral.com/1471-2105/14/95/abstract)
-   [The Mitochondrial Genome of the Gymnosperm Cycas taitungensis
    Contains a Novel Family of Short Interspersed Elements, Bpu
    Sequences, and Abundant RNA Editing
    Sites](http://mbe.oxfordjournals.org/content/25/3/603.short)
-   [Automatic annotation of organellar genomes with
    DOGMA](http://bioinformatics.oxfordjournals.org/content/20/17/3252)
-   [Automated generation of heuristics for biological sequence
    comparison](http://www.biomedcentral.com/1471-2105/6/31)
-   [MAKER-P: a tool-kit for the rapid creation, management, and quality
    control of plant genome
    annotations](http://www.plantphysiol.org/content/early/2013/12/06/pp.113.230144)
-   [The Norway spruce genome sequence and conifer genome
    evolution](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature12211.html)
-   [OrganellarGenomeDRAW (OGDRAW): a tool for the easy generation of
    high-quality custom graphical maps of plastid and mitochondrial
    genomes](http://nar.oxfordjournals.org/content/41/W1/W575)
-   [Comparative chloroplast genomics reveals the evolution of
    *Pinaceae* genera and
    subfamilies](http://gbe.oxfordjournals.org/content/2/504)
-   [QUAST: quality assessment tool for genome
    assemblies](http://bioinformatics.oxfordjournals.org/content/29/8/1072)
-   [Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-3.0. 1996-2010
    http://www.repeatmasker.org](http://www.repeatmasker.org/)
-   [The Sequence Alignment/Map format and
    SAMtools](http://bioinformatics.oxfordjournals.org/content/25/16/2078)
-   [tRNAscan-SE: A Program for Improved Detection of Transfer RNA Genes
    in Genomic
    Sequence](http://nar.oxfordjournals.org/content/25/5/0955)
-   [Assembling the 20Gb white spruce (*Picea glauca*) genome from
    whole-genome shotgun sequencing
    data](http://bioinformatics.oxfordjournals.org/content/29/12/1492)

[Gremme, G., Steinbiss, S., & Kurtz, S.
(2013)](http://www.computer.org/csdl/trans/tb/2013/03/ttb2013030645-abs.html)
GenomeTools: a comprehensive software library for efficient processing
of structured genome annotations. *IEEE/ACM Transactions on
Computational Biology and Bioinformatics*, 10(3), 645-656.

[Kurtz, S., Phillippy, A., Delcher, A. L., Smoot, M., Shumway, M.,
Antonescu, C., & Salzberg, S. L.
(2004)](http://genomebiology.com/content/5/2/R12) Versatile and open
software for comparing large genomes. *Genome biology*, 5(2), R12.

[Lin, C. P., Huang, J. P., Wu, C. S., Hsu, C. Y., & Chaw, S. M.
(2010)](http://gbe.oxfordjournals.org/content/2/504) Comparative
chloroplast genomics reveals the evolution of Pinaceae genera and
subfamilies. *Genome biology and evolution*, 2, 504-517.

[Barkan, A.
(1988).](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC457051/) Proteins
encoded by a complex chloroplast transcription unit are each translated
from both monocistronic and polycistronic mRNAs. *The EMBO journal*,
7(9), 2637.

[Hildebrand, M., Hallick, R. B., Passavant, C. W., & Bourque, D. P.
(1988)][trans-splicing] Trans-splicing in chloroplasts: the rps 12 loci
of Nicotiana tabacum. *Proceedings of the National Academy of Sciences*,
85(2), 372-376.

------------------------------------------------------------------------

Supplementary material
======================

Figure S1: Classify plastid sequences
-------------------------------------

![Figure S1](figure/plastid-classify.png)

**Figure S1**: Six plastid sequences were separated by length and depth
of coverage using thresholds chosen by inspection

Figure S2: Classify mitochondrial sequences
-------------------------------------------

![Figure S2](figure/mt-classify.png)

**Figure S2**: Mitochondrial sequences were separated by length, depth
of coverage and GC content using k-means clustering in R

Figure S3: MUMmer alignment of the plastid
------------------------------------------

![Figure S3](figure/plastid-mummer.png)

**Figure S3**: MUMmer was used to identify the inverted repeat of the
plastid

Table S1: Software
------------------

<table>
<thead>
<tr class="header">
<th align="left">Software</th>
<th align="left">Version</th>
<th align="left">DOI</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">ABySS</td>
<td align="left">1.3.7</td>
<td align="left"><a href="http://dx.doi.org/10.1101/gr.089532.108">10.1101/gr.089532.108</a></td>
</tr>
<tr class="even">
<td align="left">BLAST</td>
<td align="left">2.2.29</td>
<td align="left"><a href="http://dx.doi.org/10.1016/S0022-2836(05)80360-2">10.1016/S0022-2836(05)80360-2</a></td>
</tr>
<tr class="odd">
<td align="left">BWA</td>
<td align="left">0.7.8</td>
<td align="left"><a href="http://dx.doi.org/10.1093/bioinformatics/btp324">10.1093/bioinformatics/btp324</a></td>
</tr>
<tr class="even">
<td align="left">Barrnap</td>
<td align="left">0.4.2</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">DOGMA</td>
<td align="left">NA</td>
<td align="left"><a href="http://dx.doi.org/10.1093/bioinformatics/bth352">10.1093/bioinformatics/bth352</a></td>
</tr>
<tr class="even">
<td align="left">Exonerate</td>
<td align="left">2.2.0</td>
<td align="left"><a href="http://dx.doi.org/10.1186/1471-2105-6-31">10.1186/1471-2105-6-31</a></td>
</tr>
<tr class="odd">
<td align="left">GenomeTools</td>
<td align="left">1.5.3</td>
<td align="left"><a href="http://dx.doi.org/10.1109/TCBB.2013.68">10.1109/TCBB.2013.68</a></td>
</tr>
<tr class="even">
<td align="left">HMMER</td>
<td align="left">3.1b1</td>
<td align="left"><a href="http://dx.doi.org/10.1371/journal.pcbi.1002195">10.1371/journal.pcbi.1002195</a></td>
</tr>
<tr class="odd">
<td align="left">MAKER-P</td>
<td align="left">2.31.4</td>
<td align="left"><a href="http://dx.doi.org/10.1104/pp.113.230144">10.1104/pp.113.230144</a></td>
</tr>
<tr class="even">
<td align="left">MUMmer</td>
<td align="left">3.23</td>
<td align="left"><a href="http://dx.doi.org/10.1186/gb-2004-5-2-r12">10.1186/gb-2004-5-2-r12</a></td>
</tr>
<tr class="odd">
<td align="left">QUAST</td>
<td align="left">2.3</td>
<td align="left"><a href="http://dx.doi.org/10.1093/bioinformatics/btt086">10.1093/bioinformatics/btt086</a></td>
</tr>
<tr class="even">
<td align="left">RECON</td>
<td align="left">1.0.7</td>
<td align="left"><a href="http://dx.doi.org/10.1101/gr.88502">10.1101/gr.88502</a></td>
</tr>
<tr class="odd">
<td align="left">RMBlast</td>
<td align="left">2.2.28</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RepeatMasker</td>
<td align="left">4.0.5</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RepeatModeler</td>
<td align="left">1.0.7</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RepeatScout</td>
<td align="left">1.0.5</td>
<td align="left"><a href="http://dx.doi.org/10.1093/bioinformatics/bti1018">10.1093/bioinformatics/bti1018</a></td>
</tr>
<tr class="odd">
<td align="left">SAMtools</td>
<td align="left">0.1.19</td>
<td align="left"><a href="http://dx.doi.org/10.1093/bioinformatics/btp352">10.1093/bioinformatics/btp352</a></td>
</tr>
<tr class="even">
<td align="left">TRF</td>
<td align="left">4.07b</td>
<td align="left"><a href="http://dx.doi.org/10.1093/nar/27.2.573">10.1093/nar/27.2.573</a></td>
</tr>
<tr class="odd">
<td align="left">tRNAscan-SE</td>
<td align="left">1.23</td>
<td align="left"><a href="http://dx.doi.org/10.1093/nar/25.5.0955">10.1093/nar/25.5.0955</a></td>
</tr>
</tbody>
</table>

Table S2: Plastid MAKER parameters
----------------------------------

    #-----Genome (these are always required)
    genome=pg29-plastid.fa #genome sequence (fasta file or fasta embeded in GFF3 file)
    organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

    #-----EST Evidence (for best results provide a file for at least one)
    est=NC_021456.frn #set of ESTs or assembled mRNA-seq in fasta format

    #-----Protein Homology Evidence (for best results provide a file for at least one)
    protein=cds_aa.fa #protein sequence file in fasta format (i.e. from mutiple oransisms)

    #-----Repeat Masking (leave values blank to skip repeat masking)
    model_org=
    repeat_protein=

    #-----Gene Prediction
    est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no

    #-----MAKER Behavior Options
    est_forward=1 #map names and attributes forward from EST evidence, 1 = yes, 0 = no
    single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
    single_length=50 #min length required for single exon ESTs if 'single_exon is enabled'

Table S3: Mitochondrion MAKER parameters
----------------------------------------

    genome=pg29mt-concat.fa
    organism_type=eukaryotic
    protein=cds_aa.fa
    model_org=picea
    rmlib=rmlib.fa
    repeat_protein=/usr/local/opt/maker/libexec/data/te_proteins.fasta
    protein2genome=1
    trna=1
    other_gff=pg29mt-concat.rrna.gff
    est_forward=1
    single_exon=1
