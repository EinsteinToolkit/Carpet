all:

sync:
	rsync -a -v -z --exclude "*~" --exclude doxygen --delete -e ssh ./ carpet@www.carpetcode.org:www.carpetcode.org/htdocs

.PSEUDO: all sync
