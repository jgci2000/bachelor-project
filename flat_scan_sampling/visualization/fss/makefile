fss: Fss_Functions.o
	g++ fss_Ising.cpp -o fss_Ising -Ofast -std=c++17 -m64 Fss_Functions.o

fss_JDOS: Fss_Functions.o
	g++ fss_Ising_JDOS.cpp -o fss_Ising -Ofast -std=c++17 -m64 Fss_Functions.o

fss_mpi: Fss_Functions.o
	mpic++ fss_MPI_Ising.cpp -o fss_MPI_Ising -Ofast -std=c++17 -m64 Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o

clean_fss_mpi:
	rm fss_MPI_Ising
	rm *.o

clean_fss:
	rm fss_Ising
	rm *.o

clean_results:
	rm -r ./data/**
