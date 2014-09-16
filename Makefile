all: white-spruce-organelles.html white-spruce-organelles.pdf

clean:
	rm -f white-spruce-organelles.html white-spruce-organelles.pdf

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

white-spruce-organelles.html: README.md
	pandoc -s -o $@ $<

white-spruce-organelles.pdf: README.md
	pandoc -o $@ $<
