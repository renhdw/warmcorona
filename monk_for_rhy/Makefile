CXX = g++
MPICXX = mpic++

STD = -std=c++17
OPTIM = -O2
WALL = -Wall
PROFILE = -g -pg

CUDA = false

SIM5INC = .//sim5/src
SIM5OBJ = ./sim5/src/sim5lib.o
DEBUG = -UDEBUG

# compilation flags
CFLAGS = $(STD) $(OPTIM) $(WALL) -I$(SIM5INC)
LDFLAGS = $(SIM5OBJ) -lstdc++fs -lCCfits -lcfitsio

# directories
BINDIR = ./bin
OBJDIR = ./obj
TRASHDIR = ~/.trash
INSTALLDIR = ~/.local/bin

.PHONY: bins $(BINS) objs $(OBJS) clean clinstall 
# recipes for object files
OBJS = detector diskobj geoinf gridgeod kerr quadroots scatter tridgeo calhotcross
OBJSRCS = $(patsubst %, %.cpp, $(OBJS))
OBJFILES = $(patsubst %, $(OBJDIR)/%.o, $(OBJS))
objs: $(OBJS)
$(OBJS): %: %.cpp | $(OBJDIR)
	$(CXX) $(CFLAGS) $< -c -o $(OBJDIR)/$@.o

# recipes for binaries
BINS = calspec sphere
BINFILES = $(patsubst %, $(BINDIR)/%, $(BINS))
INSTALLEDBINS = $(patsubst %, $(INSTALLDIR)/%, $(BINS))
$(BINS): %: %.cpp $(OBJSRCS)
	$(CXX) $(CFLAGS) $< $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/$@
	cp $(BINDIR)/$@ $(INSTALLDIR)

3dcorona:
	$(CXX) $(CFLAGS) -fopenmp 3dcorona.cpp ../obj/scatter.o ../obj/kerr.o ../obj/gridgeod_openmp.o ../obj/scatter_openmp.o ../obj/geoinf.o ../obj/diskobj.o ../obj/quadroots.o ../obj/tridgeo.o ../obj/calhotcross.o $(LDFLAGS) -o ../bin/3dcorona
	cp ../bin/3dcorona /home/wdzhang/.local/bin

# mpi version of 3dcorona
3dcorona_h5:
	$(H5CXX) $(CFLAGS) 3dcorona_h5.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/3dcorona_h5
	cp $(BINDIR)/3dcorona_h5 $(INSTALLDIR)
3dcorona_mpi:
	$(MPICXX) $(CFLAGS) 3dcorona_mpi.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/3dcorona_mpi
	cp $(BINDIR)/3dcorona_mpi $(INSTALLDIR)
calspec_h5:
	$(H5CXX) $(CFLAGS) calspec_h5.cpp $(OBJFILES) $(LDFLAGS) -o $(BINDIR)/calspec_h5
	cp $(BINDIR)/calspec_h5 $(INSTALLDIR)
calhotcross_bin:
	$(CXX) $(CFLAGS) calhotcross.cpp ../obj/tridgeo.o ../obj/diskobj.o ../obj/scatter.o ../obj/geoinf.o ../obj/kerr.o ../obj/quadroots.o ../obj/gridgeod.o $(LDFLAGS) -o $(BINDIR)/calhotcross
	cp $(BINDIR)/calhotcross $(INSTALLDIR)

# verner95 is different from other binaries as it requires CCfits
# verner95:
#	$(CXX) -I$(SIM5INC) verner95.cpp -lstdc++fs -lCCfits -lcfitsio $(SIM5OBJ) -o $(BINDIR)/verner95
#	cp $(BINDIR)/verner95 $(INSTALLDIR)

clean:
	mv $(OBJDIR)/*.o $(TRASHDIR)
	mv $(BINDIR)/* $(TRASHDIR)
clinstall: $(INSTALLEDBINS)
	mv $^ $(TRASHDIR)
