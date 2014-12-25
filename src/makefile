CC = g++
CFLAGS = -Wall -g
BOOSTINCLUDE = -I /home/nii/.boost/boost_1_57_0
LDFLAGS = -lm

cell.o : cell.cpp cell.h
	${CC} ${CFLAGS} -c cell.cpp

constants.o : constants.cpp
	${CC} ${CFLAGS} -c constants.cpp

functions.o : functions.cpp
	${CC} ${CFLAGS} -c ${BOOSTINCLUDE} functions.cpp

sed_calculate.o : sed_calculate.cpp cell.cpp constants.cpp
	${CC} ${CFLAGS} -c ${BOOSTINCLUDE} sed_calculate.cpp

sed_calculate : sed_calculate.o functions.o constants.o cell.o
	${CC} ${CFLAGS} ${BOOSTINCLUDE} cell.o constants.o functions.o sed_calculate.o ${LDFLAGS} -o sed_calculate

clean :
	rm cell.o functions.o sed_calculate.o sed_calculate
