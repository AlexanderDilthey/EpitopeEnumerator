INCS = 
LIBS = 


# an alternative line (courtesy Peter Humburg, not working for me but for him) is
# LIBS = /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_random.so /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_filesystem.so /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_system.so /home/dilthey/bamtools/bamtools/lib/libbamtools.so /home/dilthey/bamtools/bamtools/lib/libbamtools-utils.a -lz

MKDIR_P = mkdir -p

.PHONY: directories
	
# END LIBRARY SETTINGS

#
# object and binary dirs  
#

DIR_OBJ = ../obj
DIR_BIN = ../bin

CXX    = g++
COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS)
VPATH = 

OBJS = \
        $(DIR_OBJ)/Util.o \
        $(DIR_OBJ)/tests.o \
        $(DIR_OBJ)/readFiles.o \
        $(DIR_OBJ)/enumerateEpitopes_noHaplotypePairs.o \
        $(DIR_OBJ)/enumerateEpitopes_diff_pairs.o \
        $(DIR_OBJ)/findEpitopeDifferences.o \
        $(DIR_OBJ)/enumerateEpitopes_haplotypePairs.o 
        
#
# list executable file names
#
EXECS = EpitopeEnumerator

OUT_DIR = ../obj ../bin

directories: ${OUT_DIR}


#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

all: directories $(EXECS)

$(EXECS): $(OBJS)
	$(foreach EX, $(EXECS), $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, $(EXECS), $(COMPILE) $(OBJS) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

$(DIR_OBJ)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@


#
# odds and ends
#
clean:
	/bin/rm $(DIR_OBJ)/*

realclean: clean
	/bin/rm $(DIR_BIN)/*

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

