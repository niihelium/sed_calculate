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

data_reader.o : data_reader.cpp data_reader.h cell.cpp functions.cpp
	${CC} ${CFLAGS} -c ${BOOSTINCLUDE} data_reader.cpp	

sed_calculate.o : sed_calculate.cpp functions.cpp data_reader.cpp
	${CC} ${CFLAGS} -c ${BOOSTINCLUDE} sed_calculate.cpp

sed_calculate : sed_calculate.o data_reader.o functions.o constants.o cell.o
	${CC} ${CFLAGS} ${BOOSTINCLUDE} cell.o constants.o functions.o data_reader.o sed_calculate.o ${LDFLAGS} -o sed_calculate

clean :
	rm cell.o functions.o data_reader.o constants.o sed_calculate.o sed_calculate
