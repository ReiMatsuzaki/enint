include ../local.mk

# ==== Common ====
FF:=${FF}
MODS0=err_handle timer const math strutil istream fjson

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

# ==== utest(StrUtil) ====
TARGET=utest_strutil
MODULES=${MODS0} utest strutil ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_strutil: ${OBJS}
	${FC} ${FF} $^ -o $@ -cpp  -llapack -lblas
check_strutil: ${BUILD}/utest_strutil
	${BUILD}/utest_strutil

# ==== utest(ArgParser) ====
TARGET=utest_argparser
MODULES=${MODS0} utest argparser ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_argparser: ${OBJS}
	${FC} ${FF} $^ -o $@ -cpp  -llapack -lblas
check_argparser: ${BUILD}/utest_argparser
	${BUILD}/utest_argparser -s s1s2 -i 523 -x 10.1 -y 10.2 10.3 10.4

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

# ==== utest(DRT) ====
MODULES=${MODS0} utest drt utest_drt
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_drt: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl -cpp 
check_drt: ${BUILD}/utest_drt
	${BUILD}/utest_drt

# ==== utest(Det) ====
MODULES=${MODS0} utest det utest_det
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/utest_det: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl -cpp 
check_det: ${BUILD}/utest_det
	${BUILD}/utest_det

# ==== Enint ====
TARGET=enint
MODULES=${MODS0} strutil argparser rmolfunc rnshel ${TARGET}
OBJS=$(addsuffix .o, ${MODULES})
${BUILD}/${TARGET}: ${OBJS}
	${FC} ${FF} $^ -o $@ -llapack -lblas -lfgsl
check_nars: ${BUILD}/${TARGET}
	${BUILD}/${TARGET} na_rs ../gms/hcp/out/nshel.json utest/rs.csv utest/mat.csv
check_grid: ${BUILD}/${TARGET}
	${BUILD}/${TARGET} grid -x linspace,0.0,5.0,6 -y scalar,0.0 -z scalar,0.0 -o utest/xyz.csv
check_aogrid: ${BUILD}/${TARGET} check_grid
	${BUILD}/${TARGET} ao_grid -i utest/xyz.csv -o utest/aogrid.csv -nshel ../gms/hcp/out/nshel.json -t 000


