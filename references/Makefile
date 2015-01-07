# Fetch relevant literature using Entrez Direct
# Written by Shaun Jackman http://sjackman.ca

all: \
	gymnosperm-cp.2014.doi.tsv \
	gymnosperm-cp.2014.organism.tsv \
	spermatophyta-mt.2014.doi.tsv \
	spermatophyta-mt.2014.organism.tsv \
	viridiplantae-mt.organism.tsv

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

%.query.xml: %.query
	esearch -db nuccore -query "`<$<`" >$@

%.2014.query.xml: %.query.xml
	efilter -datetype CreateDate -mindate 2014 -maxdate 2014 <$< >$@

%.docsum.xml: %.query.xml
	efetch -format docsum <$< >$@

%.xml: %.query.xml
	efetch -format xml <$< >$@

%.organism.orig.tsv: %.docsum.xml
	 xtract -pattern DocumentSummary -element Organism <$< |sort -u >$@

%.organism.tsv: %.organism.orig.tsv
	sed 's/ subsp.*//' $< |uniq >$@

%.taxid.tsv: %.docsum.xml
	 xtract -pattern DocumentSummary -element TaxId <$< |sort -nu >$@

%.title.tsv: %.docsum.xml
	 xtract -pattern DocumentSummary -element Title <$< |sort -u >$@

%.doi.tsv: %.xml
	 xtract -pattern Seq-entry -element DOI <$< |sort -u >$@

# Convert a list of species names to phylogenetic tree
# Web app: http://phylot.biobyte.de/
%.organism.tree: %.organism.tsv
	bin/bp_taxonomy2tree.pl -e $< >$@

# Render a phylogenetic tree
#%.organism.tree.svg: %.organism.tree:
# brew cask install figtree
# open -a FigTree
