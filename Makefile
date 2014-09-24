all: white-spruce-organelles.html white-spruce-organelles.pdf

clean:
	rm -f white-spruce-organelles.html white-spruce-organelles.pdf

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

white-spruce-organelles.html: README.md
	pandoc -s -o $@ $<

white-spruce-organelles.pdf: README.md
	pandoc -o $@ $<
