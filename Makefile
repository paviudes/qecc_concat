# Use the compiler options and link commands in file:///opt/intel/documentation_2019/en/mkl/common/mkl_link_line_advisor.htm .
# The LAPACK variables can be defined in shell by running the following command.
#	source /opt/intel/compilers_and_libraries_2019/mac/bin/compilervars.sh intel64
# The Random number generator is defined in mt19937ar.c .
# See https://gist.github.com/xuhdev/1873316.
# See also for conditional statements in Makefile: https://www.gnu.org/software/make/manual/html_node/Conditional-Syntax.html
MODE=RUN
ifdef db
	MODE=DEBUG
endif
ifeq ($(MODE), DEBUG)
	CC = gcc-10
	OPTS = -O${db} -g
	REPORT = $()
	TARGET = bmark
	LDFLAGS = $()
	LIBS_MATH = -lm
else
	CC = icc
	OPTS = -O3
	# -xavx # only works on icc
	REPORT = -qopt-report-phase=vec -qopt-report=5
	TARGET = bmark.so
	LDFLAGS = -shared
endif
### To do: the following 3 lines should not be apllicable for the debug mode.
#ifeq ($(strip $(shell icc --version)),)
#	CC = icc
#	OPTS = -O3
#endif

CFLAGS = -fPIC -Wall -Wextra -std=c11 $(OPTS) # $(REPORT)
RM = rm
SRC_DIR = src

# Detecting the OS type.
OS := $(shell uname -s)
$(info Make is being run in ${MODE} mode on the ${OS} OS.)

ifeq ($(OS), Darwin)
	CFLAGS_MKL = -I${MKLROOT}/include
	LIBS_MKL = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread $(LIBS) -ldl
else ifeq ($(OS), Linux)
	ifeq ($(MKLROOT),)
		MKLROOT="/mnt/c/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.4.311/windows/mkl/"
	endif
	MKL_DIR = ${MKLROOT}
	CFLAGS_MKL = -m64 -I${MKL_DIR}/include
	LIBS_MKL =  -L${MKL_DIR}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl
else
	# MKL_DIR = "c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.4.311/windows/mkl/lib/intel64_win"
	#Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.4.311/windows/mkl/lib/intel64_win
	MKL_DIR = "c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020/windows"
	CFLAGS_MKL =  -I${MKL_DIR}/mkl/include
	LIBS_DIR = ${MKL_DIR}/mkl/lib/intel64_win
	#LIBS_MKL = -L ${LIBS_DIR}/mkl_intel_lp64_dll.lib ${LIBS_DIR}/mkl_intel_thread_dll.lib ${LIBS_DIR}/mkl_core_dll.lib
	LIBS_MKL = -L${LIBS_DIR}/mkl_rt.lib
	# ${LIBS_DIR}/libiomp5md.lib
	#${MKL_DIR}/lib/intel64_win/mkl_lapack95_ilp64.lib
endif

$(shell mkdir -p obj/)

$(TARGET):obj/main.o obj/rand.o obj/reps.o obj/utils.o obj/sampling.o obj/constants.o obj/printfuns.o obj/mt19937ar.o obj/checks.o obj/logmetrics.o obj/memory.o obj/qecc.o obj/decode.o obj/effective.o obj/benchmark.o obj/hybrid.o obj/linalg.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) obj/main.o obj/reps.o obj/utils.o obj/rand.o obj/sampling.o obj/constants.o obj/printfuns.o obj/mt19937ar.o obj/logmetrics.o obj/checks.o obj/memory.o obj/qecc.o obj/decode.o obj/effective.o obj/benchmark.o obj/hybrid.o $(LIBS_MKL) obj/linalg.o

obj/main.o: $(SRC_DIR)/main.c Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/main.c -o obj/main.o $(LIBS_MATH)

obj/mt19937ar.o: $(SRC_DIR)/mt19937/mt19937ar.c $(SRC_DIR)/mt19937/mt19937ar.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/mt19937/mt19937ar.c -o obj/mt19937ar.o $(LIBS_MATH)

obj/rand.o: $(SRC_DIR)/rand.c $(SRC_DIR)/rand.h Makefile
	$(CC) $(CFLAGS) $(CFLAGS_MKL) $(LIBS_MKL) -c $(SRC_DIR)/rand.c -o obj/rand.o $(LIBS_MATH)

obj/linalg.o: $(SRC_DIR)/linalg.c $(SRC_DIR)/linalg.h Makefile
	$(CC) $(CFLAGS) $(CFLAGS_MKL) $(LIBS_MKL) -c $(SRC_DIR)/linalg.c -o obj/linalg.o $(LIBS_MATH)

obj/benchmark.o: $(SRC_DIR)/benchmark.c $(SRC_DIR)/benchmark.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/benchmark.c -o obj/benchmark.o $(LIBS_MATH)

obj/hybrid.o: $(SRC_DIR)/hybrid.c $(SRC_DIR)/hybrid.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/hybrid.c -o obj/hybrid.o $(LIBS_MATH)

obj/effective.o: $(SRC_DIR)/effective.c $(SRC_DIR)/effective.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/effective.c -o obj/effective.o $(LIBS_MATH)

obj/memory.o: $(SRC_DIR)/memory.c $(SRC_DIR)/memory.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/memory.c -o obj/memory.o $(LIBS_MATH)

obj/qecc.o: $(SRC_DIR)/qecc.c $(SRC_DIR)/qecc.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/qecc.c -o obj/qecc.o $(LIBS_MATH)

obj/decode.o: $(SRC_DIR)/decode.c $(SRC_DIR)/decode.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/decode.c -o obj/decode.o $(LIBS_MATH)

obj/logmetrics.o: $(SRC_DIR)/logmetrics.c $(SRC_DIR)/logmetrics.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/logmetrics.c -o obj/logmetrics.o $(LIBS_MATH)

obj/checks.o: $(SRC_DIR)/checks.c $(SRC_DIR)/checks.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/checks.c -o obj/checks.o $(LIBS_MATH)

obj/printfuns.o: $(SRC_DIR)/printfuns.c $(SRC_DIR)/printfuns.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/printfuns.c -o obj/printfuns.o $(LIBS_MATH)

obj/constants.o: $(SRC_DIR)/constants.c $(SRC_DIR)/constants.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/constants.c -o obj/constants.o $(LIBS_MATH)

obj/utils.o: $(SRC_DIR)/utils.c $(SRC_DIR)/utils.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/utils.c -o obj/utils.o $(LIBS_MATH)

obj/reps.o: $(SRC_DIR)/reps.c $(SRC_DIR)/reps.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/reps.c -o obj/reps.o $(LIBS_MATH)

obj/sampling.o: $(SRC_DIR)/sampling.c $(SRC_DIR)/sampling.h Makefile
	$(CC) $(CFLAGS) -c $(SRC_DIR)/sampling.c -o obj/sampling.o $(LIBS_MATH)

clean:
	$(RM) bmark bmark.so obj/main.o obj/rand.o obj/reps.o obj/utils.o obj/sampling.o obj/constants.o obj/printfuns.o obj/mt19937ar.o obj/checks.o obj/logmetrics.o obj/memory.o obj/qecc.o obj/decode.o obj/effective.o obj/benchmark.o obj/hybrid.o obj/linalg.o
