NPROCESS = 4
SITE_DIR = ../../site_scons

all : .PHONY

.PHONY :
	scons --site-dir=${SITE_DIR} -j${NPROCESS}

tree :
	scons --site-dir=${SITE_DIR} --tree=all

clean :
	rm -rf ./build
	rm estuary.log
