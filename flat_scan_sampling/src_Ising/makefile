all: fss fss_mpi

fss: Fss_Functions.o
	g++ fss_Ising.cpp -o fss_Ising -Ofast -std=c++17 -m64 Fss_Functions.o
fss_mpi: Fss_Functions.o
	mpic++ fss_MPI_Ising.cpp -o fss_MPI_Ising -Ofast -std=c++17 -m64 Fss_Functions.o
	mpic++ fss_MPI_Ising_all_worker.cpp -o fss_MPI_Ising_all_worker -Ofast -std=c++17 -m64 Fss_Functions.o

Fss_Functions.o: Fss_Functions.cpp
	g++ -c Fss_Functions.cpp -o Fss_Functions.o -Ofast -std=c++17 -m64

clean: clean_fss_mpi clean_fss clean_o

clean_fss_mpi:
	rm fss_MPI_Ising
	rm fss_MPI_Ising_all_worker
clean_fss:
	rm fss_Ising
clean_o:
	rm *.o

clean_results:
	rm -r ./data/**
