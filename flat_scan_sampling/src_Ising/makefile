fss: Ising.o Fss_Functions.o
	g++ fss_Ising.cpp -o fss_Ising -Ofast -std=c++17 -m64 -march=native Ising.o Fss_Functions.o

Ising.o: Ising.cpp
	g++ -c Ising.cpp -o Ising.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o

clean:
	rm *.o
