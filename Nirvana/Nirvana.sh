#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
set -x                          # Output commands
set -e                          # Abort on errors



################################################################################
# Build
################################################################################

echo "BEGIN MESSAGE"
echo "Building Nirvana..."
echo "END MESSAGE"

# Set locations
THORN=Nirvana
NAME=Nirvana
SRCDIR=$(dirname $0)
INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
NIRVANA_DIR=${INSTALL_DIR}

(
    exec >&2                    # Redirect stdout to stderr
    set -x                      # Output commands
    set -e                      # Abort on errors
    cd ${SCRATCH_BUILD}
    if [ -e ${DONE_FILE} -a ${DONE_FILE} -nt ${SRCDIR}/dist/${NAME}.tar.gz \
                         -a ${DONE_FILE} -nt ${SRCDIR}/Nirvana.sh ]
    then
        echo "Nirvana: The enclosed Nirvana library has already been built; doing nothing"
    else
        echo "Nirvana: Building enclosed Nirvana library"
        
        # Should we use gtar or tar?
        TAR=$(gtar --help > /dev/null 2> /dev/null && echo gtar || echo tar)
        
        # Set up environment
        unset LIBS
        
        echo "Nirvana: Preparing directory structure..."
        mkdir external done 2> /dev/null || true
        rm -rf ${INSTALL_DIR}
        mkdir ${INSTALL_DIR}
        
        echo "Nirvana: Unpacking archive..."
        pushd ${INSTALL_DIR}
        ${TAR} xzf ${SRCDIR}/dist/${NAME}.tar.gz
        
        echo "Nirvana: Building..."
        cd ${NAME}
        for dir in $SYS_INC_DIRS; do
            CPPFLAGS="$CPPFLAGS -I$dir"
        done
        ${CXX} ${CPPFLAGS} ${CXXFLAGS} -c *.cc $(for dir in ${HDF5_INC_DIRS}; do echo -I${dir}; done)
        ${AR} ${ARFLAGS} libNirvana.a *.o
	if [ ${USE_RANLIB} = 'yes' ]; then
	    ${RANLIB} ${RANLIBFLAGS} libNirvana.a
        fi
        popd
        
        date > ${DONE_FILE}
        echo "Nirvana: Done."
    fi
)
    
if (( $? )); then
    echo 'BEGIN ERROR'
    echo 'Error while building Nirvana. Aborting.'
    echo 'END ERROR'
    exit 1
fi



# Set options
NIRVANA_INC_DIRS="${NIRVANA_DIR}/Nirvana"
NIRVANA_LIB_DIRS="${NIRVANA_DIR}/Nirvana"
NIRVANA_LIBS='Nirvana'



################################################################################
# Configure Cactus
################################################################################

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "HAVE_NIRVANA     = 1"
echo "NIRVANA_DIR      = ${NIRVANA_DIR}"
echo "NIRVANA_INC_DIRS = ${NIRVANA_INC_DIRS}"
echo "NIRVANA_LIB_DIRS = ${NIRVANA_LIB_DIRS}"
echo "NIRVANA_LIBS     = ${NIRVANA_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(NIRVANA_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(NIRVANA_LIB_DIRS)'
echo 'LIBRARY           $(NIRVANA_LIBS)'
