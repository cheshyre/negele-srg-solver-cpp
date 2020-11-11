CXX			= g++
RM			= rm -rf
CXXFLAGS	= -Wall -O3 -std=c++17
LDFLAGS		= -fopenmp -lopenblas

BUILDDIR		= build
LIBSOURCEDIR	= negele
EXECSOURCEDIR	= exec
HEADERDIR		= ./
# TESTSOURCEDIR	= tests
# TESTEXECDIR		= tests/build
# TESTCATCHDIR	= tests/catch_build

# TEST3BODYSOURCEDIR		= tests3body
# TEST3BODYEXECDIR		= tests3body/build

# DATADIR			= tests/data
# DATA3BODYDIR	= tests3body/data

LIBSOURCES 		= $(shell find $(LIBSOURCEDIR) -name '*.cpp')
LIBOBJECTS		= $(addprefix $(BUILDDIR)/,$(LIBSOURCES:%.cpp=%.o))
EXECSOURCES 	= $(shell find $(EXECSOURCEDIR) -name '*.cpp')
EXECBINARIES	= $(addprefix $(BUILDDIR)/,$(EXECSOURCES:%.cpp=%.out))
# TESTSOURCES 	= $(shell find $(TESTSOURCEDIR) -name '*.cpp')
# TESTBINARIES 	= $(addprefix $(TESTEXECDIR)/,$(TESTSOURCES:%.cpp=%.out))
# TEST3BODYSOURCES 	= $(shell find $(TEST3BODYSOURCEDIR) -name '*.cpp')
# TEST3BODYBINARIES 	= $(addprefix $(TEST3BODYEXECDIR)/,$(TEST3BODYSOURCES:%.cpp=%.out))
# DATASOURCES			= $(shell find $(DATADIR) -name '*.gz')
# DATAUNPACKED		= $(DATASOURCES:%.gz=%)
# DATA3BODYSOURCES	= $(shell find $(DATA3BODYDIR) -name '*.gz')
# DATA3BODYUNPACKED	= $(DATA3BODYSOURCES:%.gz=%)

# EXTLIBCATCH		= $(TESTCATCHDIR)/main.o
# EXTLIBWIGNER	= lib/wigxjpf/libwigxjpf.a
# LIBS = $(EXTLIBCATCH) $(EXTLIBWIGNER)

all: bin

bin: $(EXECBINARIES)

# test: test_bin
# 	for test in $(TESTBINARIES); do echo $$test; $$test; done

# test3body: test3body_bin
# 	for test in $(TEST3BODYBINARIES); do echo $$test; $$test; done

# test_bin: $(TESTBINARIES) $(DATAUNPACKED)

# test3body_bin: $(TEST3BODYBINARIES) $(DATA3BODYUNPACKED)

.SECONDARY:

$(BUILDDIR)/%.out: %.cpp $(LIBOBJECTS)
	@mkdir -p $(@D)
	$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -I$(HEADERDIR) -o $@

$(BUILDDIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(HEADERDIR) -c $< -o $@

# $(TESTEXECDIR)/%.out: %.cpp $(TESTCATCHDIR)/main.o $(LIBOBJECTS)
# 	@mkdir -p $(@D)
# 	$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -I$(HEADERDIR) -o $@

# $(TEST3BODYEXECDIR)/%.out: %.cpp $(TESTCATCHDIR)/main.o $(LIBOBJECTS)
# 	@mkdir -p $(@D)
# 	$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -I$(HEADERDIR) -o $@

# $(TESTCATCHDIR)/main.o: $(TESTSOURCEDIR)/main.cc
# 	@mkdir -p $(@D)
# 	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(HEADERDIR) -c $< -o $@

# $(DATADIR)/%: $(DATADIR)/%.gz
# 	gunzip -c $< > $@

# $(DATA3BODYDIR)/%: $(DATA3BODYDIR)/%.gz
# 	gunzip -c $< > $@

clean:
	# $(RM) $(DATAUNPACKED)
	# $(RM) $(DATA3BODYUNPACKED)
	$(RM) $(LIBOBJECTS)
	# $(RM) $(TESTBINARIES)
	# $(RM) $(TEST3BODYBINARIES)
	$(RM) $(EXECBINARIES)

# clean_lib:
# 	$(RM) $(LIBS)

# lib/wigxjpf/libwigxjpf.a:
# 	make -C lib/wigxjpf/ clean
# 	make -C lib/wigxjpf/
# 	cp lib/wigxjpf/lib/libwigxjpf.a lib/wigxjpf/
