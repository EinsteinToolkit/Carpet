all:

sync:
	rsync -a -v -z --exclude "*~" --exclude doxygen --delete -e ssh --rsync-path=/home/schnette/bin/rsync ./ carpet@www.carpetcode.org:www.carpetcode.org/htdocs

.PSEUDO: all sync
