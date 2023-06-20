CPP := mpicxx
CPPFLAGS := -O3 -std=c++11 -DFFTW3 -fopenmp -MMD -MP

BUILD_BASE := $(BASEDIR)/build
BUILD_DIR := $(BUILD_BASE)/simulation
BIN_DIR := $(BASEDIR)/bin
SRCDIR := .

PROJECT = $(BIN_DIR)/powerspecICs $(BIN_DIR)/fft_test $(BIN_DIR)/constGravPot
PROJECT += $(BIN_DIR)/gaussian_blob $(BIN_DIR)/fullTest $(BIN_DIR)/sod $(BIN_DIR)/maxGrav
PROJECT += $(BIN_DIR)/noGravity $(BIN_DIR)/plane_wave
PROJECT += $(BIN_DIR)/specc 
PROJECT += $(BIN_DIR)/coherent
PROJECT += $(BIN_DIR)/wavepacket
PROJECT += $(BIN_DIR)/poisson_test
PROJECT += $(BIN_DIR)/zero_velocity
PROJECT += $(BIN_DIR)/gio_test
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
HEADERS = -I$(SRCDIR) -I$(BASEDIR)/common -I$(SWFFT_DIR) -I$(BASEDIR)/wavefunction -I$(BASEDIR)/initializer

LIBS += $(SWFFT_DIR)/build/distribution.o 
LIBS += $(BUILD_BASE)/wavefunction/Wavefunction.cxx.o 
LIBS += $(BUILD_BASE)/common/Timer.cxx.o
LIBS += $(BUILD_BASE)/initializer/Initializer.cxx.o
LIBS += $(BUILD_BASE)/initializer/CBRNG_Random.cxx.o
LIBS += $(BUILD_BASE)/initializer/initread.cxx.o
LIBS += $(BUILD_BASE)/initializer/CosmoClass.cxx.o

#LIBS += $(FFTW3_DIR)/libgsl.a 

LDFLAGS= -L$(FFTW3_DIR) -lfftw3 -lgsl -lgslcblas -L$(SWFFT_DIR)/build 

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
	@echo "[ $(CPP) ] $@"

all : $(DEPS) $(PROJECT)
gaussian: $(BIN_DIR)/gaussian_blob

#linking
$(BIN_DIR)/powerspecICs: $(BUILD_DIR)/tests/powerSpec.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/powerSpec.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/fft_test: $(BUILD_DIR)/tests/fft_test.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/fft_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/constGravPot: $(BUILD_DIR)/tests/constGravPot.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/constGravPot.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/gaussian_blob: $(BUILD_DIR)/tests/gaussian_blob.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/gaussian_blob.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/fullTest: $(BUILD_DIR)/tests/fullTest.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/fullTest.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/sod: $(BUILD_DIR)/tests/sod.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/sod.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/maxGrav: $(BUILD_DIR)/tests/maxGrav.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/maxGrav.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/noGravity: $(BUILD_DIR)/tests/noGravity.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/noGravity.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/plane_wave: $(BUILD_DIR)/tests/plane_wave.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/plane_wave.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/specc: $(BUILD_DIR)/tests/specc.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/specc.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/zero_velocity: $(BUILD_DIR)/tests/zero_velocity.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/zero_velocity.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/coherent: $(BUILD_DIR)/tests/coherent.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/coherent.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/wavepacket: $(BUILD_DIR)/tests/wavepacket.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/wavepacket.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/poisson_test: $(BUILD_DIR)/tests/poisson_test.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/poisson_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
	@echo "[ $(CPP) ] $@"

$(BIN_DIR)/gio_test: $(BUILD_DIR)/tests/gio_test.cxx.o $(CORE_OBJECTS) 
	@$(CPP) $(CPPFLAGS) $(BUILD_DIR)/tests/gio_test.cxx.o $(CORE_OBJECTS) $(LIBS) -o $@ $(LDFLAGS)
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

.PHONY: qclean
qclean:
	rm -r $(BIN_DIR)/powerspecICs $(BUILD_DIR)/tests/powerSpec.cxx.o

.PHONY: depclean
depclean:
	rm -f $(DEPS)

.PHONY: test
test:
	echo $(SWFFT_DIR)
	echo $(SOURCE_FILES)
	echo $(OBJECTS)

clean-all: clean depclean
