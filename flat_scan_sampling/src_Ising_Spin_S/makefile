fss: Fss_Functions.o
	clang++ fss_Ising_SpinS.cpp -o fss_Ising_SpinS -Ofast -std=c++17 -stdlib=libc++ -m64 Fss_Functions.o
	clang++ fss_Ising_SpinS_new_scan.cpp -o fss_Ising_SpinS_new_scan -Ofast -std=c++17 -stdlib=libc++ -m64 Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	clang++ -c Fss_Functions.cpp -o Fss_Functions.o

clean_fss:
	rm fss_Ising_SpinS
	rm *.o

clean_results:
	rm -r ./data/**
