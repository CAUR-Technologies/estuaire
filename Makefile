#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
# @Copyright (C) 2010
# All rights reserved.
#

SHELL := /bin/sh

# ---- Paths & tools ---------------------------------------------------------

VENV        := .venv
PY          := $(VENV)/bin/python
PIP         := $(VENV)/bin/pip
RENDERER    := $(VENV)/bin/mako-render

SRCDIR      := ./src
BUILDDIR    := ./build
INCLUDE_DIR := -I$(BUILDDIR)/include

TEMPLATES_DIR := templates
TEMPLATES_SRC := $(wildcard $(TEMPLATES_DIR)/*hpp.mako)
TEMPLATES_OBJ := $(TEMPLATES_SRC:.mako=)

CC       := g++
CXXFLAGS := -O3 -g -Wall

# ---- Default ---------------------------------------------------------------

.PHONY: all
all: venv templates wheel

# ---- Virtualenv & tooling --------------------------------------------------

.PHONY: venv
venv:
	@if [ ! -d "$(VENV)" ]; then \
		python3 -m venv "$(VENV)"; \
	fi
	@$(PIP) install --upgrade pip
	# Tools used by Makefile recipes (not runtime deps)
	@$(PIP) install build Mako

# Optional: install cibuildwheel when needed
.PHONY: venv-cibw
venv-cibw: venv
	@$(PIP) install cibuildwheel

# ---- C++ sanity test (optional) -------------------------------------------

.PHONY: test
test: venv templates
	@mkdir -p $(BUILDDIR)
	$(CC) $(SRCDIR)/test.cpp $(INCLUDE_DIR) $(CXXFLAGS) -o $(BUILDDIR)/$@
	$(BUILDDIR)/$@ > test.txt

# ---- Clean -----------------------------------------------------------------

.PHONY: clean
clean:
	rm -rf ./build

.PHONY: distclean
distclean: clean
	rm -rf $(VENV) ./dist ./*.egg-info

# ---- Template rendering & generated sources --------------------------------

.PHONY: templates
templates: $(BUILDDIR)/include/nvect.hpp \
           $(BUILDDIR)/include/narray.hpp \
           $(BUILDDIR)/include/solver.hpp \
           $(BUILDDIR)/include/raytrace.hpp \
           $(BUILDDIR)/solver.pyx \
           $(BUILDDIR)/raytrace.pyx

$(BUILDDIR)/include/%.hpp: $(TEMPLATES_DIR)/%.hpp.mako | venv
	@mkdir -p $(BUILDDIR)/include
	$(RENDERER) $< > $@

$(BUILDDIR)/%.pyx: $(TEMPLATES_DIR)/%.pyx.mako $(SRCDIR)/eikonal/%.pxd | venv
	@mkdir -p $(BUILDDIR)/eikonal
	@touch $(BUILDDIR)/eikonal/__init__.py
	@cp $(SRCDIR)/eikonal/*.pxd $(BUILDDIR)/eikonal/
	$(RENDERER) $< > $@

# ---- Python package build/install (PEP 517 isolation) ---------------------

# Builds both wheel and sdist with build-isolation (respects pyproject.toml)
.PHONY: build wheel sdist
build: wheel sdist

wheel: venv templates
	$(PY) -m build --wheel

sdist: venv templates
	$(PY) -m build --sdist

# Local editable install for dev (keeps isolation of your global Python)
.PHONY: dev
dev: venv
	$(PIP) install -e .

# Install the wheel you just built (clean runtime install in the venv)
.PHONY: install
install: venv wheel
	$(PIP) install dist/*.whl

# ---- Manylinux/macOS/Windows wheels (optional) -----------------------------

.PHONY: cibuildwheel
cibuildwheel: venv-cibw
	$(PY) -m cibuildwheel --output-dir dist
