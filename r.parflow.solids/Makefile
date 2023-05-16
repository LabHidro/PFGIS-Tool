# fix this relative to include/
# or use absolute path to the GRASS source code
MODULE_TOPDIR = ../../..

PGM = r.parflow.solids

EXTRA_CLEAN_DIRS=doxygenhtml

LIBES = $(GISLIB) $(RASTERLIB)
DEPENDENCIES = $(GISDEP) $(RASTERDEP)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd

man:
	doxygen $(MODULE_TOPDIR)/Doxyfile
