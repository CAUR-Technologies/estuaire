#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#

.SUFFIXES: templates/%.mako
SHELL = /bin/sh

TEMPLATES_DIR:= templates
TEMPLATES_SRC:= $(wildcard ${TEMPLATES_DIR}/*hpp.mako)
TEMPLATES_OBJ:= $(TEMPLATES_SRC:.mako=)

SRCDIR:= ./src
BUILDDIR:= ./build
RENDERER:= mako-render

CC:= g++
CXXFLAGS:= -O3 -g -Wall
INCLUDE_DIR:= -I${BUILDDIR}/include

all: templates python_package

test : src/test.cpp templates
	${CC} src/test.cpp ${INCLUDE_DIR} ${CXXFLAGS} -o ${BUILDDIR}/$@
	${BUILDDIR}/$@ > test.txt

clean:
	rm -rf ./build

templates : build/include/nvect.hpp build/include/narray.hpp build/include/solver.hpp build/include/raytrace.hpp build/solver.pyx build/raytrace.pyx

python_package: build/include/narray.hpp build/include/solver.hpp build/include/raytrace.hpp eikonal/* templates/solver.hpp.mako
	python setup.py build


build/include/%.hpp : templates/%.hpp.mako
	@if test ! -d build/include ; then mkdir -p build/include; fi
	$(RENDERER) $< > $@

build/%.pyx : templates/%.pyx.mako src/eikonal/%.pxd 
	@if test ! -d build/eikonal ; then mkdir -p build/eikonal; fi
	@touch build/eikonal/__init__.py
	@cp src/eikonal/*.pxd build/eikonal/
	$(RENDERER) $< > $@

