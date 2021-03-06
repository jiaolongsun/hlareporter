CXX=		g++
CXXFLAGS=	-Wall -O3
LIBS=		-lz
UTILITIES_DIR = ../common/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/


# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/BamTools -I$(UTILITIES_DIR)/Utilities

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= BamToFastqMain.cpp BamToFastq.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=BamReader.o BGZF.o SequenceUtilities.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= bamToFastq


all: $(PROGRAM)

.PHONY: all


$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)
	
$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@echo "  * compiling dependencies" 
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/BamTools/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/Utilities/
		
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
