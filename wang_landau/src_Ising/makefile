all: wl wl_S wl_t

wl: WL_Functions.o
	g++ wl_Ising.cpp -o wl_Ising2 -Ofast -std=c++17 -m64 WL_Functions.o

wl_S: WL_Functions.o
	g++ wl_S_Ising.cpp -o wl_S_Ising -Ofast -std=c++17 -m64 WL_Functions.o

wl_t: WL_Functions.o
	g++ wl_t_Ising.cpp -o wl_t_Ising -Ofast -std=c++17 -m64 WL_Functions.o

WL_Functions.o: WL_Functions.cpp
	g++ -c WL_Functions.cpp -o WL_Functions.o

clean_wl:
	rm *.o
	rm wl_Ising
	rm wl_S_Ising
	rm wl_t_Ising

clean_results:
	rm -r ./data/**
