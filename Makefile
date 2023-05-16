
MODULE_TOPDIR = ../..

SUBDIRS1 = \
	r.parflow \
	r.parflow.subsurfacedepth \
	r.parflow.solids \
	r.parflow.writepfb

SUBDIRS = pflib $(SUBDIRS1)

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsubdirs


$(SUBDIRS1): pflib
