##-----------------------------------------------------------------------------
## Top level SPECC Makefile
##-----------------------------------------------------------------------------

## BASEDIR should be set by sourcing the appropriate env config file. 
#BASEDIR := $(shell pwd)
BUILD_BASE = $(BASEDIR)/build

## GIO Modules
GENERICIO_MODULE  := $(GIO_DIR)

ifndef SPECC_PLATFORM
  $(error SPECC environment not set, please source an appropriate env script)
endif

##-----------------------------------------------------------------------------
## 		T A R G E T S
##-----------------------------------------------------------------------------

default: sim

## Target to create build-directory
$(BUILD_BASE):
	mkdir -p $(BUILD_BASE)
	mkdir -p $(BUILD_BASE)/common
	mkdir -p $(BUILD_BASE)/analysis
	mkdir -p $(BUILD_BASE)/initializer
	mkdir -p $(BUILD_BASE)/initializer/tests
	mkdir -p $(BUILD_BASE)/wavefunction
	mkdir -p $(BUILD_BASE)/simulation/
	mkdir -p $(BUILD_BASE)/simulation/tests
	mkdir -p $(BASEDIR)/bin

## Target for main sim
sim: $(BUILD_BASE)
	cd wavefunction && $(MAKE)
	cd common && $(MAKE)
	cd initializer && $(MAKE)
	cd simulation && $(MAKE)
	cd analysis && $(MAKE)

## Clean	
.PHONY: clean
clean:
	cd wavefunction && $(MAKE) clean
	cd common && $(MAKE) clean
	cd initializer && $(MAKE) clean
	cd simulation && $(MAKE) clean
	cd analysis && $(MAKE) clean

.PHONY: depclean
depclean:
	cd wavefunction && $(MAKE) depclean
	cd common && $(MAKE) depclean
	cd initializer && $(MAKE) depclean
	cd simulation && $(MAKE) depclean
	cd analysis && $(MAKE) depclean

clean-all: clean depclean
