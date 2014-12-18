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

index.html: white-spruce-organelles.md
	pandoc -s -o $@ $<

README.md: white-spruce-organelles.md
	pandoc -t markdown_strict -o $@ $<

%.tex: %.md
	pandoc -s -o $@ $<

%.pdf: %.tex
	pdflatex $<
