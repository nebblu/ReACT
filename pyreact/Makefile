COPTER_DIR = ../reactions
COPTER_LIB = $(COPTER_DIR)/lib/libcopter.a

WRAPPER_LIB = libreact_wrapper.so

CFLAGS += -fPIC
CXXFLAGS += -Wall -fPIC -std=c++11
LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial

# clang uses a different libstdc++ than GCC. Also remove OpenMP when using clang
ifneq (,$(findstring clang,$(shell '$(CXX)' -v 2>&1)))
LDFLAGS += -lc++
else
LDFLAGS += -lstdc++
CXXFLAGS += -fopenmp
endif

all: $(WRAPPER_LIB)

$(WRAPPER_LIB): react_wrapper.cpp $(COPTER_LIB)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -shared -I$(COPTER_DIR)/include

$(COPTER_LIB):
	@echo "\nBuilding COPTER\n"
	cd $(COPTER_DIR) && CXXFLAGS="$(CXXFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" ./build_copter.sh

# Clean up
.PHONY: clean
clean:
	rm -f $(WRAPPER_LIB)
	rm -rf $(WRAPPER_LIB).dSYM/
	@echo "Entering $(COPTER_DIR) and running make clean.\n"
	cd $(COPTER_DIR) && make clean && make uninstall && rm Makefile




