all:

sync:
	rsync -a -v -z --exclude .DS_Store --exclude "*~" --exclude doxygen --exclude Makefile --delete --delete-excluded ./ carpet@www.carpetcode.org:www.carpetcode.org/htdocs

.PHONY: all sync
