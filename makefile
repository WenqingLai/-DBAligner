DBAligner: ArrayList.o Configure.o DNA.o Graph.o HashTable.o Main.o PEReads.o Reference.o Sam.o SubRead.o
	g++ -o DBAligner ArrayList.o Configure.o DNA.o Graph.o HashTable.o Main.o PEReads.o Reference.o Sam.o SubRead.o
ArrayList.o: ArrayList.h Configure.h ArrayList.cpp
	g++ -c ArrayList.h Configure.h ArrayList.cpp
Configure.o: Configure.h Configure.cpp
	g++ -c Configure.h Configure.cpp
DNA.o: DNA.h DNA.cpp
	g++ -c DNA.h DNA.cpp
Graph.o: HashTable.h Reference.h Configure.h Sam.h DNA.h Graph.cpp
	g++ -c HashTable.h Reference.h Configure.h Sam.h DNA.h Graph.cpp
HashTable.o: Reference.h Configure.h HashTable.h HashTable.cpp
	g++ -c Reference.h Configure.h HashTable.h HashTable.cpp
Main.o: Reference.h DNA.h Sam.h Configure.h SubRead.h Graph.h HashTable.h PEReads.h Main.cpp
	g++ -c Reference.h DNA.h Sam.h Configure.h SubRead.h Graph.h HashTable.h PEReads.h Main.cpp
PEReads.o: Graph.h SubRead.h DNA.h PEReads.h PEReads.cpp
	g++ -c Graph.h SubRead.h DNA.h PEReads.h PEReads.cpp
Reference.o: Configure.h Sam.h Reference.h Reference.cpp
	g++ -c Configure.h Sam.h Reference.h Reference.cpp
Sam.o: Sam.h Configure.h Sam.cpp
	g++ -c Sam.h Configure.h Sam.cpp
SubRead.o: Reference.h Sam.h HashTable.h ArrayList.h Graph.h DNA.h SubRead.h SubRead.cpp
	g++ -c Reference.h Sam.h HashTable.h ArrayList.h Graph.h DNA.h SubRead.h SubRead.cpp
clean:
	rm -f DBAligner *.o
