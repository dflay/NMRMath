CXX=g++
SRC_DIR = ./src
INC_DIR = ./include
LIB_DIR = ./lib
BUILD_DIR = ./build
BIN_DIR = ./bin

# Locations and default files
PREFIX=/usr/local
LIBINSTALLDIR=$(PREFIX)/lib
INCINSTALLDIR=$(PREFIX)/include/g2field

DEPS=$(INC_DIR)/NMRMath.hh
# CFLAGS = -I$(INC_DIR) -I$(DRIVER_INC_DIR) -Wl,-rpath,$(LIB_DIR)
CFLAGS = -I$(INC_DIR) -I$(DRIVER_INC_DIR) -Wall
LIBFLAGS = -L$(LIB_DIR) 

all: $(BIN_DIR)/Test  

# $(LIB_DIR)/libNMRMath.so: $(BUILD_DIR)/NMRMath.o 
# 	$(CXX) -shared -fPIC -o $@ $^ 
# 
# $(BUILD_DIR)/NMRMath.o: $(SRC_DIR)/NMRMath.cxx $(DEPS)
# 	$(CXX) -o $@ -fPIC -c $< $(CFLAGS)

$(BUILD_DIR)/Test.o: $(SRC_DIR)/Test.cxx $(DEPS)
	$(CXX) -o $@ -fPIC -c $< $(CFLAGS) `root-config --cflags`

$(BIN_DIR)/Test: $(BUILD_DIR)/Test.o 
	$(CXX) -o $@ $< $(CFLAGS)  $(LIBFLAGS) `root-config --libs`

.PHONY: clean
.PHONY: install

clean:
	rm $(BUILD_DIR)/* $(BIN_DIR)/* $(LIB_DIR)/*

install:
	mkdir -p $(LIBINSTALLDIR)
	install $(LIB_DIR)/* $(LIBINSTALLDIR)
	install include/*.hh $(INCINSTALLDIR)
