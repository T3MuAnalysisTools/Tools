include make.inc
#################################################
### Defining Compiling variables              ###
#################################################
CXX           = c++ -Wno-deprecated -Wall -std=c++0x  -std=c++1y
LD            = g++ -Wno-deprecated -Wall

INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS = -L./CommonUtils/lib -lSimpleFits
SHAREDCXXFLAGS = -I./CommonUtils
SUBDIRS = CommonUtils/ 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  

ROOTLIBS      =  -L$(ROOTSYS)/lib $(shell root-config --glibs)  -lTMVA -lMinuit2


# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Tables \
		Plots \
		DoubleEventRemoval \
		PDG_Var \
		SkimConfig \
		HistoConfig \
		Logger \
		PDGInfo	\
		PUReweight \
		Example 



CINTTARGETS = 

# now include conditional options set by flags in make.inc

ifdef USE_cherepanov 
	TARGETS += MyTest  ThreeMu  FillTMVATrees  MCBackgroundStudy MakeMVATree MakeMVACategoryII CommonSelector  VertexFits DebugFit SignalVertexSelector SimpleTauSelector ZTau3MuTauh ZTau3MuTaue ZTau3MuTaumu
	DEFS += -DUSE_cherepanov=1
endif


ifdef USE_joshi
	TARGETS += MuonPionTree
	DEFS += -DUSE_joshi=1
endif

ifdef USE_wang
	TARGETS += MyTest
	DEFS += -DUSE_wang=1
endif


ifdef USE_menendez
	TARGETS += MyTest DsToPhiPi SyncDsPhiPi
	DEFS += -DUSE_menendez=1
endif

ifdef USE_madhu
	TARGETS += MyTest MCBackgroundStudy ThreeMuonDecay CommonSelector BDTSelector MakeMVATree SignalVertexSelector DebugFit SimpleTauSelector SimpleTauMuSelector SimpleTauESelector ZTau3MuTauh_Efficiency ZTau3MuTaue_Efficiency ZTau3MuTaumu_Efficiency ZTau3MuTauh_Skimmer ZTau3MuTaue_Skimmer ZTau3MuTaumu_Skimmer ZTau3MuTauh ZTau3MuTaue ZTau3MuTaumu
	DEFS += -DUSE_madhu=1
endif



HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) libs/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) libs/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:libs

ifdef USE_cherepanov
	vpath %.cxx cherepanov
	vpath %.h cherepanov
endif


ifdef USE_joshi
	vpath %.cxx joshi
	vpath %.h joshi
endif


ifdef USE_wang
	vpath %.cxx wang
	vpath %.h wang
endif

ifdef USE_menendez
	vpath %.cxx menendez
	vpath %.h menendez
endif


ifdef USE_madhu
	vpath %.cxx madhu
	vpath %.h madhu
endif

vpath %.cxx CommonFiles
vpath %.h CommonFiles


##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o libs/$@ 

.PHONY: clean cleanall cleandf all install sharedlib 

install: Analysis.exe


clean:
	@rm libs/*.o
	@rm Analysis.exe

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@rm libs/*.o
	@rm Analysis.exe

all: sharedlib  install



sharedlib:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir echo $(CXXFLAGS) && $(MAKE) -f Makefile.standalone CXXFLAGS+='$(CXXFLAGS)' ; \
	done
	$(install)
