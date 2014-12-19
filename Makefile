all: index.html README.md white-spruce-organelles.pdf

clean:
	rm -f index.html README.md white-spruce-organelles.pdf white-spruce-organelles.tex

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
	pandoc -s -o $@ $<

# Render strict Markdown from Pandoc Markdown
README.md: white-spruce-organelles.md
	pandoc -t markdown_strict -o $@ $<

# Render the TeX from the Markdown
%.orig.tex: %.md gbe.latex
	pandoc --template=gbe -o $@ $<

# Munge the TeX
%.tex: %.orig.tex
	sed -e 's/\\begin{longtable}/\\begin{table}[!b]\\centering\\begin{tabular}/' \
		-e 's/\\end{longtable}/\\end{tabular}\\end{table}/' \
		-e 's/\\endhead//' $< >$@

# Render the PDF from the TeX
%.pdf: %.tex gbe/gbe.cls
	TEXINPUTS=.:gbe: pdflatex -interaction=batchmode $<

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
