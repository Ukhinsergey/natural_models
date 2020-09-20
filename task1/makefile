main: main.cpp
	mpicxx -std=c++14 -O3 -o main main.cpp 


report1: main
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result1 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 2 -t 1 -w 0:30 --stdout result2 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 4 -t 1 -w 0:30 --stdout result4 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 8 -t 1 -w 0:30 --stdout result8 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 16 -t 1 -w 0:30 --stdout result16 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 32 -t 1 -w 0:30 --stdout result32 main -- 0 100 50 0.5 4000000
	mpisubmit.pl -p 64 -t 1 -w 0:30 --stdout result64 main -- 0 100 50 0.5 4000000


report2: main
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result1 main -- 0 100 50 0.5 10
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result2 main -- 0 100 50 0.5 100
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result4 main -- 0 100 50 0.5 1000
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result8 main -- 0 100 50 0.5 10000
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result16 main -- 0 100 50 0.5 100000
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result32 main -- 0 100 50 0.5 1000000
	mpisubmit.pl -p 10 -t 1 -w 0:30 --stdout result64 main -- 0 100 50 0.5 10000000


report3: main
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result1 main -- 0 1000 500 0.5 1000
	mpisubmit.pl -p 2 -t 1 -w 0:30 --stdout result2 main -- 0 1000 500 0.5 2000
	mpisubmit.pl -p 4 -t 1 -w 0:30 --stdout result4 main -- 0 1000 500 0.5 4000
	mpisubmit.pl -p 8 -t 1 -w 0:30 --stdout result8 main -- 0 1000 500 0.5 8000
	mpisubmit.pl -p 16 -t 1 -w 0:30 --stdout result16 main -- 0 1000 500 0.5 16000
	mpisubmit.pl -p 32 -t 1 -w 0:30 --stdout result32 main -- 0 1000 500 0.5 32000

report4: main
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result1 main -- 0 100 50 0.5 10
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result2 main -- 0 100 50 0.5 100
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result4 main -- 0 100 50 0.5 1000
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result8 main -- 0 100 50 0.5 10000
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result16 main -- 0 100 50 0.5 100000
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result32 main -- 0 100 50 0.5 1000000
	mpisubmit.pl -p 1 -t 1 -w 0:30 --stdout result64 main -- 0 100 50 0.5 10000000

clean :
	rm -f *.err result* *.txt main