include config.mk

OBJECTS =	random.o		\
		utility.o		\
		problem_instance.o	\
		formulas.o		\
		exact.o                 \
		pamf.o                  \
		maths.o			\
		ampl_output_reader.o	\
		instance_generator.o	\
		tabu_list.o		\
		algorithm_basic.o	\
		algorithm_unused.o	\
		algorithm.o	  	\
		algorithm_pw.o		

MyOBJ = $(OBJECTS:%=${LIB}/%)

#ALLOBJS = $(OBJECTS) $(MAINOBJECT)

MYLIB = mylib.a
MAIN3 = main
MAIN4 = makeinstance

all: k-hc-2014 makeinstance-2014

k-hc-2014: ${MYLIB} ${MAIN3}.o
	${LINKER} -o ${BIN}/$@ ${LIB}/${MAIN3}.o ${LIB}/${MYLIB} -L${LAPACKPP_LIB} -L${CPLEX_LIBS} -L${CONCERT_LIBS}

makeinstance-2014: ${MYLIB} ${MAIN4}.o
	${LINKER} -o ${BIN}/$@ ${LIB}/${MAIN4}.o ${LIB}/${MYLIB} -L${LAPACKPP_LIB} -L${CPLEX_LIBS} -L${CONCERT_LIBS}

${MYLIB}: ${MyOBJ}
	cd ${LIB}; \
	ar ruv $@ ${MyOBJ}

%.o:    %.cpp
	${COMPILER} -c -I${BOOST_INCLUDE} -I${LAPACKPP_INCLUDE} -I${CPLEX_INCLUDE} -I${CONCERT_INCLUDE} -o ${LIB}/$@ $<

${LIB}/%.o:    %.cpp
	${COMPILER} -c -I${BOOST_INCLUDE} -I${LAPACKPP_INCLUDE}  -I${CPLEX_INCLUDE} -I${CONCERT_INCLUDE} -o $@ $<

clean:
	rm -r ${LIB}/*.o ${LIB}/${MYLIB}

tar:
	tar -zcvvf hyperplane_clustering.tar.gz algorithm_apr.cpp ampl.cpp formulas.h main.cpp maths.h problem_instance.h tabu_list.h algorithm.cpp ampl.h exact.cpp inclusions.h makefile pamf.cpp random.cpp timer.h algorithm.h ampl_output_reader.cpp exact.h instance_generator.cpp makeinstance.cpp pamf.h random.h utility.cpp algorithm_pw.cpp ampl_output_reader.h formulas.cpp instance_generator.h maths.cpp problem_instance.cpp tabu_list.cpp utility.h
