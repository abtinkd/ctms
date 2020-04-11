#
# Makefile for non-Microsoft compilers
#

all: MakeAll

MakeAll:
	$(MAKE) -C ctms
	$(MAKE) -C converter
	$(MAKE) -C extractor
	$(MAKE) -C predictor
	$(MAKE) -C analyzer

clean:	
	$(MAKE) clean -C ctms
	$(MAKE) clean -C converter
	$(MAKE) clean -C extractor
	$(MAKE) clean -C predictor
	$(MAKE) clean -C analyzer
