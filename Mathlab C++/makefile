SEProjectG21: CMatrixFiles

CMatrixFiles: Main.o CMatrix.o
	g++ Main.o CMatrix.o -o CmatrixFiles

Main.o: Main.cpp
	g++ -c Main.cpp



CMatrix.o: CMatrix.cpp
	g++ -c CMatrix.cpp

clean:
	rm -f *o CmatrixFiles 
