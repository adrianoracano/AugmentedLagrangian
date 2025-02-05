CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS += -O3 
CPPFLAGS += -w
# CPPFLAGS += -Wall -pedantic 
CPPFLAGS += -I. -I$(PACS_ROOT)/include -I$(mkEigenInc) -I$(PACS_ROOT)/src/Utilities -I$(mkBoostInc)

LDFLAGS ?= -L$(PACS_ROOT)/lib -L$(PACS_ROOT)/lib64 
LDFLAGS += -Wl,-rpath,$(PACS_ROOT)/lib # muparser is shared
LDFLAGS += -L$(mkBoostLib)

LDLIBS  ?= -lmuparser 
LDLIBS  += -lboost_iostreams -lboost_system -lboost_filesystem # to use gnuplot-iostream


LINK.o := $(LINK.cc) # Use C++ linker.

DEPEND = make.dep

EXEC ?= main

# Source files and object files
COMMON_SRCS = $(filter-out main.cpp main_test.cpp, $(wildcard *.cpp))
COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)

# Phony targets
.PHONY: all test clean distclean $(DEPEND)

# Default target
all: $(EXEC)

# Main target
main: main.o $(COMMON_OBJS)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Specialized targets (only compile necessary sources)
test: main_test

# Compile only the relevant main file
main_test: main_test.o $(COMMON_OBJS)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Object file rules
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Clean up object files and dependencies
clean:
	$(RM) $(DEPEND)
	$(RM) *.o

# Clean up everything including executables
distclean: clean
	$(RM) main main_test
	$(RM) *.csv *.out *.bak *~

# Generate dependencies
$(DEPEND): $(COMMON_SRCS)
	@$(RM) $(DEPEND)
	@for file in $(COMMON_SRCS); do \
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

# Include the dependency file
-include $(DEPEND)