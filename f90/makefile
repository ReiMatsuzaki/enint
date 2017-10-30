include ../local.mk

# ==== Common ====
FF:=${FF}
MODS0=err_handle timer const math istream fjson

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
MODULES=${MODS0} utest nshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_nshel: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl
check_nshel: ${BUILD}/utest_nshel
	${BUILD}/utest_nshel

# ==== Endvr ====
TARGET=endvr
MODULES=${MODS0} nshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/${TARGET}: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl
check_${TARGET}: ${BUILD}/${TARGET}
	${BUILD}/${TARGET} ../gms/hcp/out/nshel.json utest/rs.csv utest/mat.csv


