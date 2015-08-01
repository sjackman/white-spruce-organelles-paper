all: index.html \
	README.md \
	white-spruce-organelles.pdf \
	white-spruce-organelles-supp.html \
	white-spruce-organelles-supp.pdf

docx: \
	white-spruce-organelles.docx \
	white-spruce-organelles-supp.docx

clean:
	rm -f \
		index.html \
		README.md \
		white-spruce-organelles.pdf \
		white-spruce-organelles.tex \
		white-spruce-organelles-supp.html \
		white-spruce-organelles-supp.pdf \
		white-spruce-organelles-supp.tex

.PHONY: all clean docx
.DELETE_ON_ERROR:
.SECONDARY:

# Dependencies

white-spruce-organelles.pdf: \
	figure/plastid-annotation.png \
	figure/mt-gene-order.png \
	figure/mt-genes.png \
	figure/mt-repeats.png \
	figure/mt-cds-heatmap.png \
	figure/mt-cds-orf-heatmap.png

# Rules

# Render HTML from Markdown
index.html: white-spruce-organelles.md
	pandoc -s --filter=pandoc-citeproc -o $@ $<

# Render strict Markdown from Pandoc Markdown
README.md: white-spruce-organelles.md readme.markdown_strict
	pandoc --template=readme -t markdown_strict --columns=80 --filter=pandoc-citeproc -o $@ $<

# Render docx from Markdown
%.docx: %.md
	pandoc --filter=pandoc-citeproc -o $@ $<

# Render the TeX from the Markdown
%.orig.tex: %.md gbe.latex
	pandoc --template=gbe --filter=pandoc-citeproc -o $@ $<

# Munge the TeX
%.tex: %.orig.tex
	sed -e '/\\begin{longtable}/{h;s/.*/\\begin{table*}[!bt]/;}' \
		-e '/}\\tabularnewline/{s/\\tabularnewline//;p;g;}' \
		-e 's/\\begin{longtable}/\\centering\\begin{tabular}/' \
		-e '/\\endfirsthead/,/\\endhead/d' \
		-e 's/\\end{longtable}/\\end{tabular}\\end{table*}/' \
		-e 's/{figure}/{figure*}/' \
		$< >$@

# Render the manuscript PDF
white-spruce-organelles.pdf: %.pdf: %.tex gbe/gbe.cls
	TEXINPUTS=gbe: pdflatex -interaction=batchmode $<

# Render the supplementary material HTML
white-spruce-organelles-supp.html: %.html: %.md
	pandoc -s --toc --filter=pandoc-citeproc -o $@ $<

# Render the supplementary material TeX
white-spruce-organelles-supp.tex: %.tex: %.md
	pandoc -s --toc --filter=pandoc-citeproc -o $@ $<

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
