# $Id: Makefile 512 2006-11-21 15:22:16Z goerke $
include common.mk

.PHONY: all doc test dist distclean message example

relpathheaders=$(addprefix $(INCDIR)/, $(headers))
relpathlibraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfiles))

# how to build the stuff in the lib directory:
SUBDIRS = lib
.PHONY: $(SUBDIRS)

new:
	$(MAKE) -C $(LIBSRCDIRNAME)
	$(MAKE) clean -C $(SUBDIRS)
	$(MAKE) -C $(SUBDIRS)
 
clean:
	$(MAKE) $@ -C $(SUBDIRS)
	$(MAKE) $@ -C $(LIBSRCDIRNAME)
    
#-----------------End of general section---------------------------------------

