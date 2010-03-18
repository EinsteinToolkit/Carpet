all:

sync:
	rsync -a -v -z --exclude .DS_Store --exclude "*~" --exclude doxygen --exclude Makefile --delete --delete-excluded ./ carpet@carpetcode.dyndns.org:www.carpetcode.org/htdocs
	rsync -a -v -z --exclude .DS_Store --exclude "*~" --exclude doxygen --exclude Makefile --delete --delete-excluded -e ssh --rsync-path=/home/schnette/bin/rsync ./ carpet@www.carpetcode.org:www.carpetcode.org/htdocs

.PHONY: all sync
