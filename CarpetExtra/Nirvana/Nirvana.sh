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
NAME=Nirvana
SRCDIR=$(dirname $0)
INSTALL_DIR=${SCRATCH_BUILD}
NIRVANA_DIR=${INSTALL_DIR}/build-${NAME}/${NAME}

# Set up environment
unset LIBS
    
(
    exec >&2                    # Redirect stdout to stderr
    set -x                      # Output commands
    set -e                      # Abort on errors
    cd ${INSTALL_DIR}
    if [ -e done-${NAME} -a done-${NAME} -nt ${SRCDIR}/dist/${NAME}.tar.gz \
                         -a done-${NAME} -nt ${SRCDIR}/Nirvana.sh ]
    then
        echo "Nirvana: The enclosed Nirvana library has already been built; doing nothing"
    else
        echo "Nirvana: Building enclosed Nirvana library"
        
        # Should we use gtar or tar?
        TAR=$(gtar --help > /dev/null 2> /dev/null && echo gtar || echo tar)
        
        echo "Nirvana: Unpacking archive..."
        rm -rf build-${NAME}
        mkdir build-${NAME}
        pushd build-${NAME}
        ${TAR} xzf ${SRCDIR}/dist/${NAME}.tar.gz
        popd
        
        echo "Nirvana: Building..."
        pushd build-${NAME}/${NAME}
        ${CXX} ${CPPFLAGS} ${CXXFLAGS} -c *.cc $(for dir in ${HDF5_INC_DIRS}; do echo -I${dir}; done)
        ${AR} ${ARFLAGS} libNirvana.a *.o
	if [ ${USE_RANLIB} = 'yes' ]; then
	    ${RANLIB} ${RANLIBFLAGS} libNirvana.a
        fi
        popd
        
        echo 'done' > done-${NAME}
        echo "Nirvana: Done."
    fi
)
    
if (( $? )); then
    echo 'BEGIN ERROR'
    echo 'Error while building Nirvana.  Aborting.'
    echo 'END ERROR'
    exit 1
fi



# Set options
NIRVANA_INC_DIRS="${NIRVANA_DIR}"
NIRVANA_LIB_DIRS="${NIRVANA_DIR}"
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
