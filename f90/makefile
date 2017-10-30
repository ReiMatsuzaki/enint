include ../local.mk

# ==== Common ====
FF:=${FF}
MODS0=const write_title err_handle timer sys convert istream fjson math timestep 

# -- Directories --
SRC=${NAEWDY_PATH}/src
BUILD=${NAEWDY_PATH}/build
TESTS=${NAEWDY_PATH}/tests

# -- compile --
%.o: %.f
	${FC} ${FF} -c $< -o $@

%.o: %.f90
	${FC} ${FF} -c $< -o $@

# -- Clean --
clean:
	rm -f *.o
	rm -f *.mod

# ==== utest(nshel) ====
TARGET=utest_nshel
MODULES=err_handle timer const math istream fjson utest nshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/${TARGET}: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl
check_nshel: ${BUILD}/${TARGET}
	${BUILD}/${TARGET}
