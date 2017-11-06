include ../local.mk

# ==== Common ====
FF:=${FF}
MODS0=err_handle timer const math istream fjson

# -- Directories --
SRC=${NAEWDY_PATH}/src
BUILD=${NAEWDY_PATH}/build
TESTS=${NAEWDY_PATH}/tests

# -- compile --
rnshel.o: temp_nshel.f90
	${FC} ${FF} -c $< -o $@ -DREALGTO

hnshel.o: temp_nshel.f90
	${FC} ${FF} -c $< -o $@ -DHERMITEGTO

rmolfunc.o: molfunc.f90
	${FC} ${FF} -c $< -o $@ -DREALFUNC

%.o: %.f90
	${FC} ${FF} -c $< -o $@

# -- Clean --
clean:
	rm -f *.o
	rm -f *.mod

# ==== utest(ArgParser) ====
TARGET=utest_argparser
MODULES=${MODS0} utest argparser ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_argparser: ${OBJS}
	${FC} ${FF} $^ -o $@ -cpp  -llapack -lblas
check_argparser: ${BUILD}/utest_argparser
	${BUILD}/utest_argparser -s s1s2 -i 523 -x 10.1 -y 10.2 

# ==== utest(nshel) ====
TARGET=utest_nshel
MODULES=${MODS0} utest rmolfunc rnshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_nshel: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl -cpp 
check_nshel: ${BUILD}/utest_nshel
	${BUILD}/utest_nshel

# ==== utest(molfunc) ====
TARGET=utest_molfunc
MODULES=${MODS0} utest rmolfunc ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_molfunc: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl -cpp 
check_molfunc: ${BUILD}/utest_molfunc
	${BUILD}/utest_molfunc

# ==== Enint ====
TARGET=enint
MODULES=${MODS0} nshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/${TARGET}: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl
check_${TARGET}: ${BUILD}/${TARGET}
	${BUILD}/${TARGET} na_rs ../gms/hcp/out/nshel.json utest/rs.csv utest/mat.csv


