fss: Fss_Functions.o
	g++ fss_Ising_SpinS.cpp -o fss_Ising_SpinS -Ofast -std=c++17 -m64 Fss_Functions.o
	g++ fss_Ising_SpinS_new_scan.cpp -o fss_Ising_SpinS_new_scan -Ofast -std=c++17 -m64 Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o -Ofast -std=c++17 -m64

clean: clean_fss

clean_fss:
	rm fss_Ising_SpinS
	rm *.o

clean_results:
	rm -r ./data/**
