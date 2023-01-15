CC = g++

CFLAGS = -std=c++17 -Wall
LFLAGS = -I. -I./include -lASImage -lMinuit
ROOTLIBS = `root-config --cflags --ldflags --glibs`
LIBSMAIN = $(CFLAGS) $(LFLAGS) $(ROOTLIBS)

SRCDIR   = src
OBJDIR   = obj

SOURCES     := $(wildcard $(SRCDIR)/*.cxx)
OBJECTS     := $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)

SOURCESEXE  := $(wildcard $(SRCDIR)/*.cc)
OBJECTSEXE  := $(SOURCESEXE:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
BINS        := $(SOURCESEXE:$(SRCDIR)/%_exe.cc=%)

all: $(BINS)
	@echo "--> Successfully created all executables."

$(BINS): $(OBJECTS) $(OBJECTSEXE)
	@echo "--> Creating executables $@"
	@$(CC) -c $(SRCDIR)/$@_exe.cc -o $(OBJDIR)/$@_exe.o $(LIBSMAIN)
	@$(CC) $(SRCDIR)/$@_exe.cc -o $@ $(OBJECTS) $(LIBSMAIN)

$(OBJECTSEXE): $(OBJDIR)/%.o : $(SRCDIR)/%.cc
	@echo "--> Creating object $@"
	@mkdir -p $(OBJDIR)
	@$(CC) -c $< -o $@ $(LIBSMAIN)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	@echo "--> Creating object $@"
	@mkdir -p $(OBJDIR)
	@$(CC) -c $< -o $@ $(LIBSMAIN)


clean:
	@rm -f $(wildcard $(OBJDIR)/*.o) $(BINS)

.PHONY: all clean
