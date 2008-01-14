all:

sync:
	rsync -a -v -z --exclude .DS_Store --exclude "*~" --exclude doxygen --exclude Makefile --delete --delete-excluded -e ssh --rsync-path=/home/schnette/bin/rsync ./ carpet@www.carpetcode.org:www.carpetcode.org/htdocs

.PSEUDO: all sync
