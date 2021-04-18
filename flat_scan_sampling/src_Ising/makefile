fss: Fss_Functions.o
	g++ fss_Ising.cpp -o fss_Ising -Ofast -std=c++17 -m64 -march=native Fss_Functions.o

fss_mpi: Fss_Functions.o
	mpic++ fss_MPI_Ising.cpp -o fss_MPI_Ising -Ofast -std=c++17 -march=native -m64 Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o

clear_fss_mpi:
	rm *.o
	rm fss_MPI_Ising

clear_fss:
	rm *.o
	rm fss_Ising

clear_results:
	rm -r ./data/**
