CPP := mpicxx
CPPFLAGS := -O3 -std=c++11 -DFFTW3 -fopenmp -MMD -MP 

BUILD_BASE := $(BASEDIR)/build
BUILD_DIR := $(BUILD_BASE)/analysis
BIN_DIR := $(BASEDIR)/bin
SRCDIR := .

PROJECT = $(BIN_DIR)/massProfile
PROJECT += $(BIN_DIR)/phase_convergence
########### SOURCE FILES ###################
SUBDIRS := $(shell ls -F | grep "\/" )

SOURCE_FILES := $(wildcard *.cxx)
CORE_OBJECTS := $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)
SOURCE_FILES += $(foreach dir, $(SUBDIRS), $(wildcard $(dir)*.cxx))

########### OBJ AND DEP FILES ###################

OBJECTS := $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

############ HEADERS ######################
FFTW3_DIR = $(shell dirname $(shell dirname $(shell which fftw-wisdom)))
HEADERS = -I$(SRCDIR) -I$(BASEDIR)/common -I$(SWFFT_DIR) -I$(BASEDIR)/initializer -I$(SWFFT_DIR)

LIBS += $(SWFFT_DIR)/build/distribution.o -lfftw3

LDFLAGS=-L$(FFTW3_DIR) -L$(SWFFT_DIR)/build

###########################################

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
 LDFLAGS +=
 LIBS += ${GIO_LIBS} 
endif
###########################################

DEPENDENCY_OPTIONS = -MM

#Create .d files

$(BUILD_DIR)/%.cxx.d : %.cxx
	$(CPP) $(DEPENDENCY_OPTIONS) $(HEADERS) $< -MT "$(BUILD_DIR)/$*.o $(BUILD_DIR)/$*.d" -MF $(BUILD_DIR)/$*.cxx.d 

all : $(DEPS) $(PROJECT)
#linking

$(BIN_DIR)/massProfile : $(BUILD_DIR)/massProfile.cxx.o 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/massProfile.cxx.o $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/phase_convergence : $(BUILD_DIR)/phase_convergence.cxx.o 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/phase_convergence.cxx.o $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

ifneq "$(strip $(DEPS))" ""
-include $(DEPS)
endif

$(BUILD_DIR)/%.cxx.o: %.cxx
	$(CPP) $(CPPFLAGS) $(HEADERS) -c $< -o $@

.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(PROJECT) $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(DEPS)

clean-all: clean depclean
