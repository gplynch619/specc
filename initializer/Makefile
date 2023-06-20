CPP := mpicxx 
CPPFLAGS := -O3 -std=c++11 -DFFTW3 -fopenmp -MMD -MP -DVELTEST

BUILD_BASE := $(BASEDIR)/build
BUILD_DIR := $(BUILD_BASE)/initializer
BIN_DIR := $(BASEDIR)/bin
SRCDIR := .

PROJECT = $(BIN_DIR)/initread $(BIN_DIR)/test_CC $(BIN_DIR)/initDensity $(BIN_DIR)/growth_factor
########### SOURCE FILES ###################
SUBDIRS := $(shell ls -F | grep "\/" )

SOURCE_FILES := $(wildcard *.cxx)
CORE_OBJECTS := $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)
SOURCE_FILES += $(foreach dir, $(SUBDIRS), $(wildcard $(dir)*.cxx))

###########################################
########### OBJ AND DEP FILES ###################

OBJECTS := $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

############ HEADERS ######################
FFTW3_DIR = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))

HEADERS = -I$(SRCDIR) -I$(BASEDIR)/common -I$(SWFFT_DIR) -I$(BASEDIR)/wavefunction -I$(FFTW3_DIR)/include

#LIBS = $(FFTW3_DIR)/libfftw3.a $(FFTW3_DIR)/libgsl.a $(SWFFT_DIR)/build/distribution.o 
LIBS = $(GSL_DIR) $(SWFFT_DIR)/build/distribution.o 
LIBS += $(BUILD_BASE)/wavefunction/Wavefunction.cxx.o 

LDFLAGS= -L$(FFTW3_DIR) -L$(SWFFT_DIR)/build -lfftw3 -lm

########## GENERICIO ######################
ifndef GENERICIO_MODULE
 GENERICIO_MODULE := $(GIO_DIR)
endif

ifdef GENERICIO_MODULE
 $(info [INFO]: using GenericIO from $(GIO_DIR))
 include $(GENERICIO_MODULE)/include.mk
 CPPFLAGS += ${GENERICIO_INCLUDES}
 HEADERS += ${GENERICIO_INCLUDES}
 HEADERS += -I$(GIO_DIR)
 LIBS += ${GIO_LIBS} 
endif
###########################################

DEPENDENCY_OPTIONS = -MM

#Create .d files

$(BUILD_DIR)/%.cxx.d : %.cxx
	@$(CPP) $(DEPENDENCY_OPTIONS) $(HEADERS) $< -MT "$(BUILD_DIR)/$*.o $(BUILD_DIR)/$*.d" -MF $(BUILD_DIR)/$*.cxx.d 
	@echo "[ $(CPP) ] $<"
all : $(DEPS) $(PROJECT)

#linking
$(BIN_DIR)/initread: $(BUILD_DIR)/tests/initread_test.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/initread_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"
	
$(BIN_DIR)/test_CC: $(BUILD_DIR)/tests/cosmoclass_test.cxx.o $(CORE_OBJECTS)
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/cosmoclass_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/growth_factor: $(BUILD_DIR)/tests/growth_factor.cxx.o $(CORE_OBJECTS)
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/growth_factor.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/initDensity: $(BUILD_DIR)/tests/density_test.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/density_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

ifneq "$(strip $(DEPS))" ""
-include $(DEPS)
endif

$(BUILD_DIR)/%.cxx.o: %.cxx
	@$(CPP) $(CPPFLAGS) $(HEADERS) -c $< -o $@
	@echo "[ $(CPP) ] $@"

.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(PROJECT) $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(DEPS)

.PHONY: test
test:
	echo $(SWFFT_DIR)
	echo $(SOURCE_FILES)
	echo $(OBJECTS)

clean-all: clean depclean
