#OPTFLAG         = -Wno-deprecated -Wall -O3 -m64 -funroll-all-loops -fPIC -fexceptions -DIL_STD -DILO_LINUX -DILO64 -g #-DNDEBUG                             
OPTFLAG		= -O3 #-O0 -Wall -pedantic -DNDEBUG
DBGFLAG		= -g
#LDFLAGS         = -lpthread
LDFLAGS         = -lm -lpthread -fPIC -m64 -fexceptions -DIL_STD


COMPILER        = g++ ${OPTFLAG} ${DBGFLAG}
LINKER          = g++ ${LDFLAGS} ${DBGFLAG}

# Directory for my files
SRC             = .
MYHOME          = ..
BIN             = ${MYHOME}/bin
INCLUDE         = ${MYHOME}/include 
LIB             = ${MYHOME}/lib

# Directory for Lapack++
LAPACKPP_DIR		= /usr/local
LAPACKPP_INCLUDE 	= ${LAPACKPP_DIR}/include/lapackpp
LAPACKPP_LIB	    = ${LAPACKPP_DIR}/lib -llapackpp

# Directory for Boost
BOOST_INCLUDE		= /opt/boost_1_47_0


# # Directory for Ilog Cplex solver and Ilog Concert
# CPLEX_HOME      = /opt/ibm/ILOG/CPLEX_Studio_Academic123/cplex
# CPLEX_INCLUDE   = ${CPLEX_HOME}/include
# CPLEX_LIBS      = ${CPLEX_HOME}/lib/x86-64_sles10_4.1/static_pic -lilocplex -lcplex
# CONCERT_HOME      = /opt/ibm/ILOG/CPLEX_Studio_Academic123/concert
# CONCERT_INCLUDE   = ${CONCERT_HOME}/include
# CONCERT_LIBS      = ${CONCERT_HOME}/lib/x86-64_sles10_4.1/static_pic -lconcert

# # Directory for Ilog Cplex solver and Ilog Concert
CPLEX_HOME      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
CPLEX_INCLUDE   = ${CPLEX_HOME}/include
CPLEX_LIBS      = ${CPLEX_HOME}/lib/x86-64_linux/static_pic -lilocplex -lcplex
CONCERT_HOME      = /opt/ibm/ILOG/CPLEX_Studio126/concert
CONCERT_INCLUDE   = ${CONCERT_HOME}/include
CONCERT_LIBS      = ${CONCERT_HOME}/lib/x86-64_linux/static_pic -lconcert
