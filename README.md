# Organellar Genomes of White Spruce (*Picea glauca*): Assembly and Annotation

Shaun D Jackman, Rene Warren, Ewan Gibb, Ben Vandervalk, Hamid Mohamadi, Justin
Chu, Anthony Raymond, Stephen Pleasance, Robin Coope, Mark R Wildung, Carol
Ritland, Steven JM Jones, Joerg C Bohlmann, Inanc Birol

# Abstract

The genome sequences of the plastid and mitochondrion of white spruce (*Picea
glauca*) are assembled from whole genome Illumina sequencing data using ABySS.
Whole genome sequencing data contains reads from both the nuclear and organellar
genomes. Reads of the organellar genomes are abundant, because each cell
contains hundreds of mitochondria and plastids. One lane of MiSeq data assembles
the 123 kbp plastid genome, and one lane of HiSeq data assembles the 5.9 Mbp
draft mitochondrial genome. The coding genes, ribosomal RNA and transfer RNA of
both genomes are annotated. The transcript abundance of the mitochondrial genes
is quantified in three developmental tissues and five mature tissues. C-to-U RNA
editing is observed in the majority of mitochondrial genes and modifies ACG
codons to create cryptic AUG start codons in four genes. NCBI GenBank contains a
single mitochondrial genome of a gymnosperm, and this work contributes the first
coniferous mitochondrial genome.

Introduction
============

Most plant cells contain two types of organelles that harbour their own genomes,
mitochondria and plastids. In the *Pinaceae*, mitochondrial genomes are
inherited maternally, and plastid genomes are inherited paternally (Whittle &
Johnston 2002).

Complete plastid genomes of the gymnosperms Norway spruce (*Picea abies*)
(Nystedt et al. 2013), *Podocarpus lambertii* (Nascimento Vieira, Faoro,
Rogalski, et al. 2014), *Taxus chinensis* var. *mairei* (Zhang et al. 2014) and
four *Juniperus* species (Guo et al. 2014) have recently been published in NCBI
Genbank (Benson et al. 2014). These projects used a variety of strategies for
isolating cpDNA, using physical separation methods in the lab or computationally
separating cpDNA sequences from nuclear sequences, sequencing and assembly,
shown in Table 1.

The *Picea abies* genome used 454 GS FLX Titanium sequencing and Sanger
sequencing of PCR amplicons for finishing, BLAST (Altschul et al. 1990) to
isolate the cpDNA reads, and the software Newbler to assemble the reads. The *P.
lambertii* genome assembly isolated the cpDNA using the saline Percoll gradient
protocol of Nascimento Vieira, Faoro, Freitas Fraga, et al. (2014), Illumina
MiSeq sequencing and the software Newbler to assemble the reads. The *Juniperus
bermudiana* genome assembly used long-range PCR to amplify the plastid DNA, a
combination of Illumina GAII and Sanger sequencing, and the software Geneious to
assemble the reads using *C. japonica* as a reference genome. The other three
*Juniperus* genome assemblies used Illumina MiSeq sequencing and the software
Velvet (Zerbino & Birney 2008) to assemble the reads. The *T. chinensis* genome
assembly used whole-genome Illumina HiSeq 2000 sequencing, BLAT (Kent 2002) to
isolate the cpDNA reads and SOAPdenovo (Luo et al. 2012) to assemble the
isolated cpDNA reads. All of these projects used DOGMA (Wyman et al. 2004) to
annotate the assembly.

<table>
<caption>Table 1: Methods of cpDNA separation, sequencing and assembly of complete plastid genomes of gymnosperms published. <sup>a</sup> Finished with PCR and Sanger sequencing</caption>
<thead>
<tr class="header">
<th align="left">Species</th>
<th align="left">cpDNA Separation</th>
<th align="left">Sequencing</th>
<th align="left">Assembler software</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><em>Picea abies</em></td>
<td align="left">BLAST <em>in silico</em></td>
<td align="left">454 GS FLX Titanium <sup>a</sup></td>
<td align="left">Newbler</td>
</tr>
<tr class="even">
<td align="left"><em>Podocarpus lambertii</em></td>
<td align="left">Saline Percoll gradient</td>
<td align="left">Illumina MiSeq</td>
<td align="left">Newbler</td>
</tr>
<tr class="odd">
<td align="left"><em>Juniperus bermudiana</em></td>
<td align="left">Longer-range PCR</td>
<td align="left">Illumina GAII <sup>a</sup></td>
<td align="left">Geneious</td>
</tr>
<tr class="even">
<td align="left">Other <em>Juniperus</em></td>
<td align="left">Unspecified</td>
<td align="left">Illumina MiSeq</td>
<td align="left">Velvet</td>
</tr>
<tr class="odd">
<td align="left"><em>Taxus chinensis</em></td>
<td align="left">BLAT <em>in silico</em></td>
<td align="left">Illumina HiSeq 2000</td>
<td align="left">SOAPdenovo</td>
</tr>
</tbody>
</table>

Only one complete mitochondrial genome of a gymnosperm has been published,
*Cycas taitungensis* (Chaw et al. 2008), while complete mitochondrial genome
sequences of the angiosperms *Brassica maritima* (Grewe et al. 2014), *Brassica
oleracea* (*ibid.*), *Capsicum annuum* (Jo et al. 2014), *Eruca sativa* (Wang et
al. 2014), *Helianthus tuberosus* (Bock et al. 2014), *Raphanus sativus* (Jeong
et al. 2014), *Rhazya stricta* (Park et al. 2014) and *Vaccinium macrocarpon*
(Fajardo et al. 2014) have been published in NCBI Genbank. Six of these projects
gave details of the sample preparation, sequencing, assembly and annotation
strategy. Three projects enriched organellar DNA using varying laboratory
methods (Keren et al. 2009; Kim et al. 2007; Chen et al. 2011), and the
remainder used total genomic DNA. Three projects used Illumina HiSeq 2000
sequencing and Velvet for assembly, and three projects used Roche 454 GS-FLX
sequencing and Newbler for assembly. Most projects used an aligner such as BLAST
(Altschul et al. 1990) to isolate sequences with similarity to known
mitochondrial sequence, either before or after assembly. Two projects used
Mitofy (Alverson et al. 2010) to annotate the genome, and the remainder used a
collection of tools such as BLAST, tRNAscan-SE (Lowe & Eddy 1997) and ORF Finder
to annotate genes. Plant mitochondrial genomes can substantially vary in size,
with some of the largest mitochondrial genomes reported for the basal angiosperm
*Amborella trichopoda* (3.9 Mbp; Rice et al. 2013) and the two Silene species
*S. noctiflora* and *S. conica* (6.7 Mbp and 11.3 Mbp, respectively; Sloan et
al. 2012).

The SMarTForests project has recently published a set of stepwise improved
assemblies of the 20 gigabase white spruce (*Picea glauca*) genome (Birol et al.
2013; Warren, Keeling, et al. 2015), a gymnosperm genome seven times the size of
the human genome, sequenced using the Illumina HiSeq and MiSeq sequencing
platforms. The whole genome sequencing data contained reads originating from
both the nuclear and organellar genomes. Whereas one copy of the diploid nuclear
genome is found in each cell, hundreds of organelles are present, and thus
hundreds of copies of the organellar genomes. This abundance results in an
overrepresentation of the organellar genomes in whole genome sequencing data.

Assembling a single lane of white spruce whole genome sequencing data using the
software ABySS (Simpson et al. 2009) yielded an assembly composed of organellar
sequences and nuclear repeat elements. The assembled sequences that originate
from the organellar genomes were separated from those of nuclear origin by
classifying the sequences using their length, depth of coverage and GC content.
The plastid genome of white spruce is compared to that of Norway spruce (*Picea
abies*) (Nystedt et al. 2013), and the mitochondrial genome of white spruce is
compared to that of prince sago palm (*Cycas taitungensis*) (Chaw et al. 2008).

Analysis of cpDNA is useful in reconstructing phylogenies of plants (Wu et al.
2007), in determining the origin of an expanding population (Aizawa et al. 2012)
and in determining when distinct lineages of a species resulted from multiple
colonization events (Jardón-Barbolla et al. 2011). These contrasting inheritance
schemes of plastids and mitochondria can be useful in the characterization of
species expanding their range. In the case of two previously allopatric species
now found in sympatry, the mitochondrial DNA (mtDNA) is contributed by the
resident species, whereas introgression of the plastid genome into the expanding
species is limited, since pollen is more readily dispersed than seeds (Du et al.
2011). Differential gene flow of cpDNA and mtDNA due to different methods of
inheritance and dispersion results in new assemblages of organellar genomes and
an increase of genetic diversity after expansion from a refugium (Gerardi et al.
2010).

Material and Methods
====================

DNA, RNA and software materials
-------------------------------

Genomic DNA was collected from the apical shoot tissues of a single interior
white spruce tree, clone PG29, and sequencing libraries constructed as described
in Birol et al. (2013). Because the original intention of this sequencing
project was to assemble the nuclear genome, an organelle exclusion method was
used to preferentially extract nuclear DNA. Sequencing reads from both
organellar genomes were present in sufficient depth however to assemble their
genomes.

RNA was extracted from eight samples, three developmental stages and five mature
tissues: megagametophyte, embryo, seedling, young buds, xylem, mature needles,
flushing buds and bark, described in Warren, Keeling, et al. (2015). These
samples were sequenced with the Illumina HiSeq 2000 (Warren, Keeling, et al.
2015). The RNA-seq data was used to quantify the transcript abundance of the
annotated mitochondrial genes using the software Salmon (Patro et al. 2014).

The software used in this analysis and their versions are listed in
supplementary Table S1. All software tools were installed using Homebrew
(<http://brew.sh>).

Methods used to assemble the plastid genome
-------------------------------------------

A single lane of Illumina MiSeq paired-end sequencing (SRR525215) was used to
assemble the plastid genome. Paired-end sequencing usually leaves a gap of
unsequenced nucleotides in the middle of the DNA fragment. Because 300 bp
paired-end reads were sequenced from a library of 500 bp DNA fragments, the
reads are expected to overlap by 100 bp. These overlapping paired-end reads were
merged using ABySS-mergepairs, a component of the software ABySS (Simpson et al.
2009). These merged reads were assembled using ABySS. Contigs that are
putatively derived from the plastid were separated by length and depth of
coverage using thresholds chosen by inspection of a scatter plot (see
supplementary Figure S1). These putative plastid contigs were assembled into
scaffolds using ABySS-scaffold.

We ran the gap-filling application Sealer (Paulino et al., in review; options
`-v -j 12 -b 30G -B 300 -F 700` with `-k` from 18 to 108 with step size 6) on
the ABySS assembly of the plastid genome, closing 5 of the remaining 7 gaps,
with a resulting assembly consisting of two large (~50 and ~70 kbp) scaftigs.
Given the small size of the plastid genome, we opted to manually finish the
assembly using the software Consed 20.0 (Gordon & Green 2013). We loaded the
resulting gap-filled assembly into Consed and imported Pacific Biosciences
(PacBio) sequencing data (SRR2148116 and SRR2148117), 9204 reads 500 bp and
larger, into the assembly and aligned them to the plastid genome using
cross\_match (Green 1999) from within Consed. For each scaftig end, 6 PacBio
reads were pulled out and assembled using the mini-assembly feature in Consed.
Cross\_match alignments of the resulting contigs to the plastid assembly were
used to merge the two scaftigs and confirm that the complete circular genome
sequence was obtained. In a subsequent step, 7,742 Illumina HiSeq reads were
imported and aligned to the assembly using Consed. These reads were selected
from the library of 133 million reads used to assemble the mitochondrion on the
basis of alignment to our draft plastid genome using BWA 0.7.5a (Li 2013),
focusing on regions that would benefit from read import by restricting our
search to regions with ambiguity and regions covered by PacBio reads
exclusively. The subset of Illumina reads were selected using samtools 0.1.18,
mini-assembled with Phrap (Green 1999) and the resulting contigs re-merged to
correct bases in gaps filled only by PacBio, namely one gap and sequence at
edges confirming the circular topology. The starting base was chosen using the
Norway spruce plastid genome sequence (NC\_021456, Nystedt et al. 2013). Our
assembly was further polished using the Genome Analysis Toolkit (GATK)
2.8-1-g932cd3a FastaAlternateReferenceMaker (McKenna et al. 2010).

The assembled plastid genome was initially annotated using DOGMA (Wyman et al.
2004). Being an interactive web application, it is not convenient for automated
annotation. The software MAKER (Campbell et al. 2014) is not interactive and is
designed for automated annotation, and we used it to annotate the white spruce
plastid using the Norway spruce plastid genome (NC\_021456, Nystedt et al. 2013)
for both protein-coding and non-coding gene homology evidence. The parameters of
MAKER are show in supplementary Table S2. The inverted repeat was identified
using MUMmer (Kurtz et al. 2004), shown in supplementary Figure S3.

The assembled plastid genome was aligned to the Norway spruce plastid using
BWA-MEM (Li 2013). The two genomes were compared using QUAST (Gurevich et al.
2013) to confirm the presence of the annotated genes of the Norway spruce
plastid in the white spruce plastid.

Methods used to assemble the mitochondrial genome
-------------------------------------------------

ABySS-Konnector (Vandervalk et al. 2014) was used to fill the gap between the
paired-end reads of a single lane of Illumina HiSeq 2000 paired end-sequencing
(SRR525196). These connected paired-end reads were assembled using ABySS.
Putative mitochondrial sequences were separated from nuclear sequences by their
length, depth of coverage and GC content using k-means clustering in R (see
supplementary Figure S2). The putative mitochondrial contigs were then assembled
into scaffolds using ABySS-scaffold with a single lane of Illumina HiSeq
sequencing of a mate-pair library.

The ABySS assembly of the white spruce mitochondrial genome resulted in 71
scaffolds. We ran the gap-filling application Sealer attempting to close the
gaps between every combination of two scaffolds. This approach closed 10 gaps
and yielded 61 scaffolds, which we used as input to the LINKS scaffolder 1.1
(Warren, Vandervalk, et al. 2015) (options `-k 15 -t 1 -l 3 -r 0.4`, 19
iterations with `-d` from 500 to 6000 with step size 250) in conjunction with
long PacBio reads, further decreasing the number of scaffolds to 58. The
Konnector pseudoreads were aligned to the 58 LINKS scaffolds with BWA 0.7.5a
(`bwa mem -a multimap`), and we created links between two scaffolds when reads
aligned within 1000 bp of the edges of any two scaffolds. We modified LINKS to
read the resulting SAM alignment file and link scaffolds satisfying this
criteria (options `LINKS-sam -e 0.9 -a 0.5`), bringing the final number of
scaffolds to 38. We confirmed the merges using mate-pair reads. The white spruce
mate-pair libraries used for confirmation are presented in Birol et al. (2013)
and available from DNAnexus (SRP014489
<http://sra.dnanexus.com/studies/SRP014489>). In brief, mate-pair reads from
three fragment size libraries (5, 8 and 12 kbp) were aligned to the 38-scaffold
assembly with BWA-MEM 0.7.10-r789 and the resulting alignments parsed with a
PERL script. A summary of this validation is presented in supplemental Table S4.
Automated gap-closing was performed with Sealer 1.0 (options
`-j 12 -B 1000 -F 700 -P10 -k96 -k80`) using Bloom filters built from the entire
white spruce PG29 read data set (Warren, Keeling, et al. 2015) and closed 55 of
the 182 total gaps (30.2%). We polished the gap-filled assembly using GATK, as
described for the plastid genome.

The assembled scaffolds were aligned to the NCBI nucleotide (nt) database using
BLAST to check for hits to mitochondrial genomes and to screen for
contamination.

The mitochondrial genome was annotated using MAKER (parameters shown in
supplementary Table S3) and Prokka (Seemann 2014), and the two sets of
annotations were merged using BEDTools (Quinlan & Hall 2010) and GenomeTools
(Gremme et al. 2013), selecting the MAKER annotation when the two tools had
overlapping annotations. The proteins of all green plants (*Viridiplantae*) with
complete mitochondrial genome sequences in NCBI GenBank (Benson et al. 2014),
142 species, were used for protein homology evidence, the most closely related
of which is the prince sago palm (*Cycas taitungensis*) (NC\_010303 Chaw et al.
2008), being the only gymnosperm with a complete mitochondrial genome. Transfer
RNA (tRNA) were annotated using ARAGORN (Laslett & Canback 2004). Ribosomal RNA
(rRNA) were annotated using RNAmmer (Lagesen et al. 2007). Prokka uses Prodigal
(Hyatt et al. 2010) to annotate open reading frames. Repeats were identified
using RepeatMasker (Smit et al. 1996) and RepeatModeler.

The RNA-seq reads were aligned to the annotated mitochondrial genes using
BWA-MEM and variants were called using samtools and bcftools requiring a minimum
genotype quality of 50 to identify possible sites of C-to-U RNA editing.

Results
=======

The white spruce plastid genome
-------------------------------

The assembly and annotation metrics for the plastid and mitochondrial genomes
are summarized in Table 2. The plastid genome was assembled into a single
circular contig of 123,266 bp containing 114 identified genes: 74 protein coding
(mRNA) genes, 36 transfer RNA (tRNA) genes and 4 ribosomal RNA (rRNA) genes,
shown in Figure 1.

<table>
<caption>Table 2: Sequencing, assembly and annotation metrics of the white spruce organellar genomes. The number of distinct genes are shown in parentheses.</caption>
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
<td align="left">4.9 million</td>
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
<td align="left">5.94 Mbp</td>
</tr>
<tr class="even">
<td align="left">Number of contigs</td>
<td align="left">1 contig</td>
<td align="left">130 contigs</td>
</tr>
<tr class="odd">
<td align="left">Contig N50</td>
<td align="left">123 kbp</td>
<td align="left">102 kbp</td>
</tr>
<tr class="even">
<td align="left">Number of scaffolds</td>
<td align="left">1 scaffold</td>
<td align="left">36 scaffolds</td>
</tr>
<tr class="odd">
<td align="left">Scaffold N50</td>
<td align="left">123 kbp</td>
<td align="left">369 kbp</td>
</tr>
<tr class="even">
<td align="left">Largest scaffold</td>
<td align="left">123 kbp</td>
<td align="left">1222 kbp</td>
</tr>
<tr class="odd">
<td align="left">GC content</td>
<td align="left">38.8%</td>
<td align="left">44.7%</td>
</tr>
<tr class="even">
<td align="left">Number of genes without ORFs</td>
<td align="left">114 (108)</td>
<td align="left">143 (74)</td>
</tr>
<tr class="odd">
<td align="left">Protein coding genes (mRNA)</td>
<td align="left">74 (72)</td>
<td align="left">106 (51)</td>
</tr>
<tr class="even">
<td align="left">Ribosomal RNA genes (rRNA)</td>
<td align="left">4 (4)</td>
<td align="left">8 (3)</td>
</tr>
<tr class="odd">
<td align="left">Transfer RNA genes (tRNA)</td>
<td align="left">36 (32)</td>
<td align="left">29 (20)</td>
</tr>
<tr class="even">
<td align="left">Open reading frames (ORF) <span class="math inline"> ≥ </span> 300 bp</td>
<td align="left">NA</td>
<td align="left">1065</td>
</tr>
<tr class="odd">
<td align="left">Coding genes containing introns</td>
<td align="left">8</td>
<td align="left">5</td>
</tr>
<tr class="even">
<td align="left">Introns in coding genes</td>
<td align="left">9</td>
<td align="left">7</td>
</tr>
<tr class="odd">
<td align="left">tRNA genes containing introns</td>
<td align="left">6</td>
<td align="left">0</td>
</tr>
</tbody>
</table>

All protein-coding genes are single copy, except *psbI* and *ycf12*, which have
two copies each. All tRNA genes are single copy, except *trnH-GUG*, *trnI-CAU*,
*trnS-GCU* and *trnT-GGU*, which have two copies each. All rRNA genes are single
copy.

The protein-coding genes *atpF*, *petB*, *petD*, *rpl2*, *rpl16*, *rpoC1* and
*rps12* each contain one intron, and *ycf3* contains two introns. The tRNA genes
*trnA-UGC*, *trnG-GCC*, *trnI-GAU*, *trnK-UUU*, *trnL-UAA* and *trnV-UAC* each
contain one intron. Of the 15 introns, 11 are determined to be group II introns
by RNAweasel (Lang et al. 2007).

The first and smallest exons of the genes *petB*, *petD* and *rpl16* are 6, 8
and 9 bp respectively. These genes likely belong to polycistronic transcripts
(Barkan 1988) of their respective protein complexes, but the short size of their
initial exons make them difficult to annotate all the same. The initial exons of
these genes were added to their annotations manually.

The gene *rps12* of a plastid genome is typically trans-spliced (Hildebrand et
al. 1988), which makes it difficult to annotate using MAKER. It is composed of
three exons and one cis-spliced intron. It required manually editing the gene
annotation to incorporate trans-splicing in the gene model.

Each copy of the inverted repeat (IR) is 445 bp in size, much smaller than most
plants, but typical of *Pinaceae* (Lin et al. 2010). Unlike most inverted
repeats, which are typically identical, the two copies differ by a single base.
The IR contains a single gene, the tRNA *trnI-CAU*.

All 114 genes of the Norway spruce plastid genome (Nystedt et al. 2013) are
present in the white spruce plastid genome in perfect synteny. Alignment of the
white spruce genome to the Norway spruce genome using BWA-MEM (Li 2013) reveal
no large-scale structural rearrangements. Alignments of the white spruce plastid
genome cover 99.7% of the Norway spruce plastid genome, and the sequence
identity in aligned regions is 99.2%.

![Figure 1: The complete plastid genome of white spruce, annotated using MAKER
and plotted using OrganellarGenomeDRAW (Lohse et al.
2007).](figure/plastid-annotation.png)

The white spruce mitochondrial genome
-------------------------------------

The mitochondrial genome was assembled into 38 scaffolds (132 contigs) with a
scaffold N50 of 369 kbp (contig N50 of 102 kbp). The largest scaffold is 1222
kbp (Table 2). The scaffolds were aligned to the NCBI nucleotide (nt) database
using BLAST. Of the 38 scaffolds, 26 scaffolds align to mitochondrial genomes, 3
small scaffolds (&lt;10 kbp) align to *Picea glauca* mRNA clones and BAC
sequences, 7 small scaffolds (&lt;10 kbp) had no significant hits, and 2 small
scaffolds (&lt;5 kbp) align to cloning vectors. These last two scaffolds were
removed from the assembly.

The mitochondrial genome contains 106 protein coding (mRNA) genes, 29 transfer
RNA (tRNA) genes and 8 ribosomal RNA (rRNA) genes. The 106 protein-coding genes
(51 distinct genes) compose 75 kbp (1.3%) of the genome. The 29 tRNA genes are
found in 20 distinct species for 15 amino acids. The relative order of the genes
on the scaffolds and gene size is shown in Figure 2. The size of each gene
family is shown in Figure 3. The precise position of each gene on its scaffold
is shown in supplementary Figure S4.

All tRNA genes are single copy, except *trnD-GUC* which has 3 copies, *trnM-CAU*
which has 7 copies, and *trnY-GUA* which has 2 copies. The rRNA gene *rrn5* has
4 copies, *rrn18* has 3 copies, and *rrn26* has 1 copy.

A large number of open reading frames are identified: 6265 of at least 90 bp,
composing 1.4 Mbp, and 1065 of at least 300 bp, composing 413 kbp. These open
reading frames do not have sufficient sequence similarity to the genes of the
*Viridiplantae* mitochondria used for protein homology evidence to be annotated
by either MAKER or Prokka.

A total of 7 introns are found in 5 distinct protein-coding genes. The
protein-coding genes *nad2*, *nad5*, and *nad7* each contain one intron, and
*nad4* and *rps3* each contain two introns. All introns are determined to be
group II introns by RNAweasel (Lang et al. 2007).

Repeats compose 390 kbp (6.6%) of the mitochondrial genome. Simple repeats and
LTR Copia, ERV1 and Gypsy are the most common repeats, shown in Figure 4.

All 39 protein coding genes and 3 rRNA genes of the *Cycas taitungensis*
mitochondrion are seen in white spruce. Of the 22 tRNA genes of *Cycas
taitungensis*, 13 are found in white spruce, and 8 tRNA genes are seen in white
spruce that are not seen in *Cycas taitungensis*.

![Figure 2: The relative order of the genes on the scaffolds, and the size of
each gene. Each box is proportional to the size of the gene including introns,
except that genes smaller than 200 bp are shown as 200 bp. The space between
genes is not to scale. An asterisk indicates that the gene name is truncated.
Only scaffolds that have annotated genes are shown.](figure/mt-gene-order.png)

![Figure 3: The gene content of the white spruce mitochondrial genome, grouped
by gene family. Each box is proportional to the size of the gene including
introns. The colour of each gene is unique within its gene
family.](figure/mt-genes.png)

![Figure 4: The repetitive sequence content of the white spruce mitochondrial
genome, annotated using RepeatMasker and RepeatModeler.](figure/mt-repeats.png)

The transcriptome of the white spruce mitochondrial genome
----------------------------------------------------------

The transcript abundance of the mitochondrial coding genes with known function
is shown in Figure 5. The transcript abundance of the mitochondrial coding genes
including open read frames is shown in Figure 6. Of the samples analyzed, the
transcriptomes of megagametophyte and embryo have the highest abundance of
coding mitochondrial genes and cluster together.

Of the 106 coding genes with known function, 60 are expressed in at least one of
the mature tissues, 29 are expressed in one of the developing tissues but not in
a mature tissue, and 17 are not found to be expressed. Of the 6265 ORFs at least
90 bp, 427 (7%) are expressed in at least one of the mature tissues, 2809 (45%)
are expressed in one of the developing tissues but not in a mature tissue, and
3029 (48%) are not found to be expressed. A gene with an abundance of at least
ten transcripts per million as quantified by Salmon is considered to be
expressed. These results are shown in Table 3.

<table>
<caption>Table 3: Number of expressed protein-coding genes and open reading frames tabulated by developmental stage.</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="right">Both</th>
<th align="right">Mature only</th>
<th align="right">Developing only</th>
<th align="right">Neither</th>
<th align="right">Sum</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">CDS</td>
<td align="right">60</td>
<td align="right">0</td>
<td align="right">29</td>
<td align="right">17</td>
<td align="right">106</td>
</tr>
<tr class="even">
<td align="left">ORF</td>
<td align="right">411</td>
<td align="right">16</td>
<td align="right">2809</td>
<td align="right">3029</td>
<td align="right">6265</td>
</tr>
<tr class="odd">
<td align="left">Sum</td>
<td align="right">471</td>
<td align="right">16</td>
<td align="right">2838</td>
<td align="right">3046</td>
<td align="right">6371</td>
</tr>
</tbody>
</table>

Possible C-to-U RNA editing, positions where the genome sequence shows C but the
RNA-seq reads shows T, is observed in 68 of 106 coding genes shown in
supplementary Table S5, with the most highly edited gene, *nad3*, seeing 9 edits
per 100 bp. It can be difficult to distinguish RNA editing events from genomic
SNV and miscalled variants caused by misaligned reads. We note however that 91%
(1601 of 1751) of the variants called from the RNA-seq data are C-to-T variants
shown in supplementary Table S6, which indicates that a large fraction of these
variants are due to C-to-U RNA editing. C-to-U RNA editing can create new start
and stop codons, but it is not able to destroy existing start and stop codons.
Editing of the ACG codon to AUG to create a cryptic start codon is frequently
seen in organellar genomes (Neckermann et al. 1994). Four genes have cryptic ACG
start codons and corroborating C-to-U RNA editing evidence in the RNA-seq data:
*mttB*, *nad1*, *rps3* and *rps4*.

![Figure 5: A heatmap of the transcript abundance of mitochondrial protein
coding genes. Each column is a tissue sample. Each row is a gene. Each cell
represents the transcript abundance of one gene in one sample. The colour scale
is log<sub>10</sub>(TPM), where TPM is transcripts per million as measured by
Salmon.](figure/mt-cds-heatmap.png)

![Figure 6: A heatmap of the transcript abundance of mitochondrial protein
coding genes, including open reading frames. Each column is a tissue sample.
Each row is a gene. Each cell represents the transcript abundance of one gene in
one sample. The colour scale is log<sub>10</sub>(TPM), where TPM is transcripts
per million as measured by Salmon.](figure/mt-cds-orf-heatmap.png)

Conclusion
==========

One lane of MiSeq sequencing of whole genome DNA is sufficient to assemble the
123 kbp plastid genome, and one lane of HiSeq sequencing of whole genome DNA is
sufficient to assemble the 5.9 Mbp mitochondrial genome of white spruce.
Additional Illumina and PacBio sequencing is used to improved scaffold
contiguity and to close scaffold gaps, after which the plastid genome is
assembled in a single contig and the largest mitochondrial scaffold is 1.2 Mbp.

The white spruce plastid genome shows no structural rearrangements when compared
with Norway spruce, and all genes of the Norway spruce (*Picea abies*) plastid
are present in the white spruce plastid. All genes of the prince sago palm
(*Cycas taitungensis*) mitochondrion are present in the white spruce
mitochondrion.

The protein coding gene content of the mitochondrial genome is quite sparse,
with 106 protein coding genes in 5.9 Mbp, in comparison to the plastid genome,
with 74 protein coding genes in 123 kbp. Nearly 7% of the mitochondrial genome
is composed of repeats, and roughly 1% is composed of coding genes. A
significant portion, over 90%, of the unusually large size of the white spruce
mitochondrial genome is yet unexplained.

Acknowledgements
================

This work was supported by Genome Canada, Genome British Columbia and Genome
Quebec as part of the SMarTForests Project (www.smartforests.ca). We thank
Carson Holt for his help with the MAKER analysis and Martin Krzywinski for his
help with Figure 2.

References
==========

Aizawa M, Kim Z-S, Yoshimaru H. 2012. Phylogeography of the korean pine (pinus
koraiensis) in northeast asia: Inferences from organelle gene sequences. Journal
of plant research. 125:713–723.

Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment
search tool. Journal of molecular biology. 215:403–410.

Alverson AJ et al. 2010. Insights into the evolution of mitochondrial genome
size from complete sequences of citrullus lanatus and cucurbita pepo
(cucurbitaceae). Molecular biology and evolution. 27:1436–1448.

Barkan A. 1988. Proteins encoded by a complex chloroplast transcription unit are
each translated from both monocistronic and polycistronic mRNAs. The EMBO
journal. 7:2637.

Benson DA et al. 2014. GenBank. Nucleic Acids Research. 43:D30–D35.
[doi: 10.1093/nar/gku1216](http://doi.org/10.1093/nar/gku1216).

Birol I et al. 2013. Assembling the 20 gb white spruce (picea glauca) genome
from whole-genome shotgun sequencing data. Bioinformatics. btt178.

Bock DG, Kane NC, Ebert DP, Rieseberg LH. 2014. Genome skimming reveals the
origin of the jerusalem artichoke tuber crop species: Neither from jerusalem nor
an artichoke. New Phytologist. 201:1021–1030.

Campbell MS et al. 2014. MAKER-p: A tool kit for the rapid creation, management,
and quality control of plant genome annotations. Plant physiology. 164:513–524.

Chaw S-M et al. 2008. The mitochondrial genome of the gymnosperm cycas
taitungensis contains a novel family of short interspersed elements, bpu
sequences, and abundant rNA editing sites. Molecular biology and evolution.
25:603–615.

Chen J et al. 2011. Substoichiometrically different mitotypes coexist in
mitochondrial genomes of brassica napus l. PLoS One. 6:e17662.

Du FK et al. 2011. Direction and extent of organelle dNA introgression between
two spruce species in the qinghai-tibetan plateau. New Phytologist.
192:1024–1033.

Fajardo D et al. 2014. The american cranberry mitochondrial genome reveals the
presence of selenocysteine (tRNA-sec and sECIS) insertion machinery in land
plants. Gene. 536:336–343.

Gerardi S, JARAMILLO-CORREA JP, Beaulieu J, Bousquet J. 2010. From glacial
refugia to modern populations: New assemblages of organelle genomes generated by
differential cytoplasmic gene flow in transcontinental black spruce. Molecular
ecology. 19:5265–5280.

Gordon D, Green P. 2013. Consed: A graphical editor for next-generation
sequencing. Bioinformatics. 29:2936–2937.
[doi: 10.1093/bioinformatics/btt515](http://doi.org/10.1093/bioinformatics/btt515).

Green P. 1999. Documentation for phrap and cross-match. University of
Washington, Seattle. <http://www.phrap.org/phredphrapconsed.html>.

Gremme G, Steinbiss S, Kurtz S. 2013. GenomeTools: A comprehensive software
library for efficient processing of structured genome annotations. Computational
Biology and Bioinformatics, IEEE/ACM Transactions on. 10:645–656.

Grewe F et al. 2014. Comparative analysis of 11 brassicales mitochondrial
genomes and the mitochondrial transcriptome of&lt; i&gt; brassica
oleracea&lt;/i&gt;. Mitochondrion.

Guo W et al. 2014. Predominant and substoichiometric isomers of the plastid
genome coexist within juniperus plants and have shifted multiple times during
cupressophyte evolution. Genome biology and evolution. 6:580–590.

Gurevich A, Saveliev V, Vyahhi N, Tesler G. 2013. QUAST: Quality assessment tool
for genome assemblies. Bioinformatics. 29:1072–1075.

Hildebrand M, Hallick RB, Passavant CW, Bourque DP. 1988. Trans-splicing in
chloroplasts: The rps 12 loci of nicotiana tabacum. Proceedings of the National
Academy of Sciences. 85:372–376.

Hyatt D et al. 2010. Prodigal: Prokaryotic gene recognition and translation
initiation site identification. BMC Bioinformatics. 11:119.
[doi: 10.1186/1471-2105-11-119](http://doi.org/10.1186/1471-2105-11-119).

Jardón-Barbolla L, Delgado-Valerio P, Geada-López G, Vázquez-Lobo A, Piñero D.
2011. Phylogeography of pinus subsection australes in the caribbean basin.
Annals of botany. 107:229–241.

Jeong Y-M et al. 2014. The complete mitochondrial genome of cultivated radish
wK10039 (raphanus sativus l.). Mitochondrial DNA. 1–2.

Jo YD, Choi Y, Kim D-H, Kim B-D, Kang B-C. 2014. Extensive structural variations
between mitochondrial genomes of cMS and normal peppers (capsicum annuum l.)
revealed by complete nucleotide sequencing. BMC genomics. 15:561.

Kent WJ. 2002. BLAT—the bLAST-like alignment tool. Genome research. 12:656–664.

Keren I et al. 2009. AtnMat2, a nuclear-encoded maturase required for splicing
of group-iI introns in arabidopsis mitochondria. Rna. 15:2299–2311.

Kim DH, Kang JG, Kim B-D. 2007. Isolation and characterization of the
cytoplasmic male sterility-associated orf456 gene of chili pepper (capsicum
annuum l.). Plant molecular biology. 63:519–532.

Kurtz S et al. 2004. Versatile and open software for comparing large genomes.
Genome biology. 5:R12.

Lagesen K et al. 2007. RNAmmer: Consistent and rapid annotation of ribosomal rNA
genes. Nucleic acids research. 35:3100–3108.

Lang BF, Laforest M-J, Burger G. 2007. Mitochondrial introns: A critical view.
Trends in Genetics. 23:119–125.
[doi: 10.1016/j.tig.2007.01.006](http://doi.org/10.1016/j.tig.2007.01.006).

Laslett D, Canback B. 2004. ARAGORN, a program to detect tRNA genes and tmRNA
genes in nucleotide sequences. Nucleic Acids Research. 32:11–16.

Li H. 2013. Aligning sequence reads, clone sequences and assembly contigs with
bWA-mEM. arXiv preprint arXiv:1303.3997.

Lin C-P, Huang J-P, Wu C-S, Hsu C-Y, Chaw S-M. 2010. Comparative chloroplast
genomics reveals the evolution of pinaceae genera and subfamilies. Genome
biology and evolution. 2:504–517.

Lohse M, Drechsel O, Bock R. 2007. OrganellarGenomeDRAW (oGDRAW): A tool for the
easy generation of high-quality custom graphical maps of plastid and
mitochondrial genomes. Current genetics. 52:267–274.

Lowe TM, Eddy SR. 1997. TRNAscan-sE: A program for improved detection of
transfer rNA genes in genomic sequence. Nucleic acids research. 25:0955–964.

Luo R et al. 2012. SOAPdenovo2: An empirically improved memory-efficient
short-read de novo assembler. Gigascience. 1:18.

McKenna A et al. 2010. The genome analysis toolkit: A mapReduce framework for
analyzing next-generation dNA sequencing data. Genome research. 20:1297–1303.
[doi: 10.1101/gr.107524.110](http://doi.org/10.1101/gr.107524.110).

Nascimento Vieira L do et al. 2014. An improved protocol for intact chloroplasts
and cpDNA isolation in conifers. PloS one. 9:e84792.

Nascimento Vieira L do et al. 2014. The complete chloroplast genome sequence of
podocarpus lambertii: Genome structure, evolutionary aspects, gene content and
sSR detection. PloS one. 9:e90618.

Neckermann K, Zeltz P, Igloi GL, Kössel H, Maier RM. 1994. The role of RNA
editing in conservation of start codons in chloroplast genomes. Gene.
146:177–182.
[doi: 10.1016/0378-1119(94)90290-9](http://doi.org/10.1016/0378-1119(94)90290-9).

Nystedt B et al. 2013. The norway spruce genome sequence and conifer genome
evolution. Nature. 497:579–584.

Park S et al. 2014. Complete sequences of organelle genomes from the medicinal
plant rhazya stricta (apocynaceae) and contrasting patterns of mitochondrial
genome evolution across asterids. BMC genomics. 15:405.

Patro R, Mount SM, Kingsford C. 2014. Sailfish enables alignment-free isoform
quantification from rNA-seq reads using lightweight algorithms. Nature
biotechnology. 32:462–464.
[doi: 10.1038/nbt.2862](http://doi.org/10.1038/nbt.2862).

Quinlan AR, Hall IM. 2010. BEDTools: A flexible suite of utilities for comparing
genomic features. Bioinformatics. 26:841–842.
[doi: 10.1093/bioinformatics/btq033](http://doi.org/10.1093/bioinformatics/btq033).

Rice DW et al. 2013. Horizontal transfer of entire genomes via mitochondrial
fusion in the angiosperm amborella. Science. 342:1468–1473.

Seemann T. 2014. Prokka: Rapid prokaryotic genome annotation. Bioinformatics.
btu153.

Simpson JT et al. 2009. ABySS: A parallel assembler for short read sequence
data. Genome research. 19:1117–1123.

Sloan DB et al. 2012. Rapid evolution of enormous, multichromosomal genomes in
flowering plant mitochondria with exceptionally high mutation rates. PLoS
biology. 10:e1001241.

Smit A, Hubley R, Green P. 1996. RepeatMasker open-3.0.
<http://www.repeatmasker.org>.

Vandervalk BP et al. 2014. Konnector: Connecting paired-end reads using a bloom
filter de bruijn graph. In: 2014 IEEE international conference on bioinformatics
and biomedicine (BIBM). Institute of Electrical & Electronics Engineers (IEEE).
[doi: 10.1109/bibm.2014.6999126](http://doi.org/10.1109/bibm.2014.6999126).

Wang Y et al. 2014. Complete mitochondrial genome of eruca sativa mill.(Garden
rocket). PloS one. 9:e105748.

Warren RL et al. 2015. Improved white spruce (picea glauca) genome assemblies
and annotation of large gene families of conifer terpenoid and phenolic defense
metabolism. The Plant Journal.
[doi: 10.1111/tpj.12886](http://doi.org/10.1111/tpj.12886).

Warren RL, Vandervalk BP, Jones SJ, Birol I. 2015. LINKS: Scaffolding genome
assemblies with kilobase-long nanopore reads. bioRxiv. 016519.
[doi: 10.1101/016519](http://doi.org/10.1101/016519).

Whittle C-A, Johnston MO. 2002. Male-driven evolution of mitochondrial and
chloroplastidial dNA sequences in plants. Molecular biology and evolution.
19:938–949.

Wu C-S, Wang Y-N, Liu S-M, Chaw S-M. 2007. Chloroplast genome (cpDNA) of cycas
taitungensis and 56 cp protein-coding genes of gnetum parvifolium: Insights into
cpDNA evolution and phylogeny of extant seed plants. Molecular biology and
evolution. 24:1366–1379.

Wyman SK, Jansen RK, Boore JL. 2004. Automatic annotation of organellar genomes
with dOGMA. Bioinformatics. 20:3252–3255.

Zerbino DR, Birney E. 2008. Velvet: Algorithms for de novo short read assembly
using de bruijn graphs. Genome research. 18:821–829.

Zhang Y et al. 2014. The complete chloroplast genome sequence of taxus chinensis
var. mairei (taxaceae): Loss of an inverted repeat region and comparative
analysis with related species. Gene. 540:201–209.
