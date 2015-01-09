all: index.html \
	README.md \
	white-spruce-organelles.pdf \
	white-spruce-organelles-supp.html \
	white-spruce-organelles-supp.pdf

clean:
	rm -f \
		index.html \
		README.md \
		white-spruce-organelles.pdf \
		white-spruce-organelles.tex \
		white-spruce-organelles-supp.html \
		white-spruce-organelles-supp.pdf \
		white-spruce-organelles-supp.tex

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

# Dependencies

white-spruce-organelles.pdf: \
	figure/mt-annotation.png \
	figure/mt-classify.png \
	figure/mt-genes.png \
	figure/mt-repeats.png \
	figure/plastid-annotation.png \
	figure/plastid-classify.png \
	figure/plastid-mummer.png

# Rules

# Render HTML from Markdown
index.html: white-spruce-organelles.md
	pandoc -s --bibliography=white-spruce-organelles.bib --csl=gbe.csl -o $@ $<

# Render strict Markdown from Pandoc Markdown
README.md: white-spruce-organelles.md white-spruce-organelles.bib gbe.csl readme.markdown_strict
	pandoc --template=readme --bibliography=white-spruce-organelles.bib --csl=gbe.csl -t markdown_strict --columns=80 -o $@ $<

# Render docx from Markdown
%.docx: %.md %.bib gbe.csl
	pandoc --bibliography=$*.bib --csl=gbe.csl -o $@ $<

# Render the TeX from the Markdown
%.orig.tex: %.md %.bib gbe.csl gbe.latex
	pandoc --template=gbe --bibliography=$*.bib --csl=gbe.csl -o $@ $<

# Munge the TeX
%.tex: %.orig.tex
	sed -e 's/\\begin{longtable}/\\begin{table}[!b]\\centering\\begin{tabular}/' \
		-e 's/\\end{longtable}/\\end{tabular}\\end{table}/' \
		-e 's/\\endhead//' $< >$@

# Render the manuscript PDF
white-spruce-organelles.pdf: %.pdf: %.tex gbe/gbe.cls
	TEXINPUTS=gbe: pdflatex -interaction=batchmode $<

# Render the supplementary material HTML
white-spruce-organelles-supp.html: %.html: %.md
	pandoc -s -o $@ $<

# Render the supplementary material TeX
white-spruce-organelles-supp.tex: %.tex: %.md
	pandoc -s -o $@ $<

# Render the supplementary material PDF
white-spruce-organelles-supp.pdf: %.pdf: %.tex
	pdflatex -interaction=batchmode $<

# Download the Genome Biology and Evolution LaTeX template
gbe_tex_template.zip:
	curl -O http://www.oxfordjournals.org/our_journals/gbe/for_authors/gbe_tex_template.zip

# Unzip the LaTeX template
gbe/gbe.cls: gbe_tex_template.zip
	unzip -od gbe_tex_template $<
	unzip gbe_tex_template/GBE_TEX_Samples.zip
	rm -rf gbe gbe_tex_template
	mv gbe_sample gbe
	touch $@

# Download the GBE citation style language (CSL)
gbe.csl:
	curl -o $@ https://raw.githubusercontent.com/citation-style-language/styles/master/genome-biology-and-evolution.csl
