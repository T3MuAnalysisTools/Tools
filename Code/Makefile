include make.inc
#################################################
### Defining Compiling variables              ###
#################################################
#-Wno-parentheses
CXX           = c++ -Wno-deprecated -Wall -std=c++0x
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS = #-L./CommonUtils/lib -lTauNtuple -L./DataFormats/lib -lDataFormats  -L./CommonUtils/lib -lClassicSVfit  -L./CommonUtils/lib -lTDCILib 
SHAREDCXXFLAGS = #-I./CommonUtils -I./DataFormats
SUBDIRS = CommonUtils/ 
TAUSPINERFLAG = 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -lTreePlayer -lRIO -lMathCore -rdynamic -Wl,--rpath $(ROOTSYS)/lib

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
		Tools \
		SkimConfig \
		HistoConfig \
		Logger \
		Example 



CINTTARGETS = 



# now include conditional options set by flags in make.inc

ifdef USE_cherepanov
	TARGETS += MyTest
	DEFS += -DUSE_cherepanov=1
endif


ifdef USE_joshi
	TARGETS += MyTest
	DEFS += -DUSE_joshi=1
endif



HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux

ifdef USE_cherepanov
	vpath %.cxx cherepanov
	vpath %.h cherepanov
endif



ifdef USE_joshi
	vpath %.cxx joshi
	vpath %.h joshi
endif

vpath %.cxx CommonFiles
vpath %.h CommonFiles

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall cleandf all dataformats install sharedlib 

install: Analysis.exe


#dataformats: 
#	@cd DataFormats; gmake install; cd ../

clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

#cleandf:#=
#	@cd DataFormats; gmake clean; cd ../

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@cd DataFormats; gmake clean; cd ../
	@rm i386_linux/*.o
	@rm Analysis.exe

all: sharedlib dataformats install
