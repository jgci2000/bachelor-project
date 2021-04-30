fss: Fss_Functions.o
	g++ fss_Ising_SpinS.cpp -o fss_Ising_SpinS -Ofast -std=c++17 -lstdc++fs -m64 -march=native Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o

clean_fss:
	rm fss_Ising_SpinS
	rm *.o

clean_results:
	rm -r ./data/**
