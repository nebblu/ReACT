COPTER_DIR = ../reactions
COPTER_LIB = $(COPTER_DIR)/lib/libcopter.a

WRAPPER_LIB = libreact_wrapper.so

CXXFLAGS += -I$(COPTER_DIR)/include -Wall

LDFLAGS += -lgsl -lgslcblas  -lstdc++ -lsundials_cvode -lsundials_nvecserial -fPIC -shared 

all: $(WRAPPER_LIB)

$(WRAPPER_LIB): react_wrapper.cpp $(COPTER_LIB)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) 

$(COPTER_LIB):
	@echo "\nBuilding COPTER\n"
	cd $(COPTER_DIR) && ./build_copter.sh

# Clean up
.PHONY: clean
clean:
	rm -f $(WRAPPER_LIB)
	rm -rf $(WRAPPER_LIB).dSYM/
	@echo "Entering $(COPTER_DIR) and running make clean.\n"
	cd $(COPTER_DIR) && make uninstall && make clean && rm Makefile




