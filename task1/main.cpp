#include <iostream>
#include <mpi.h>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;



pair<int, long long int> dowalk(int a, int b, int x, double p) {
    long long dt = 0;
    while (a < x && x < b) {
        double r = ((double) rand()) / RAND_MAX;
        if (r < p) {
            ++x;
        } else {
            --x;
        }
        ++dt;
    }
    return pair<int , long long int>(x, dt);
}


int main(int argc, char **argv) {
    if(argc != 6) {
        cout << "input: a b x p n" << endl;
        return 1;
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int x = atoi(argv[3]);
    double p = atof(argv[4]);
    int n = atoi(argv[5]);
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int myn = n / size;
    srand(time(NULL) + rank);
    int w = 0;
    long long t = 0;
    double start_time, end_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    for(int i = 0; i < myn; ++i) {
        pair<int, long long int> temp = dowalk(a, b, x, p);
        if (temp.first == b) {
            ++w;
        }
        t += temp.second;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    double time = end_time - start_time;
    double buf[2];
    double total[2];
    total[0] = 0;
    total[1] = 0;
    buf[0] = ((double) w) / myn;
    buf[1] = ((double) t) / myn;
    MPI_Reduce(buf, total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        ofstream fout("output.txt", std::ofstream::out | std::ofstream::app);
        fout << total[0]/size << ' ' << total[1]/size << endl;
        fout << "-------------------------------------------------------------------------------------------" << endl;
        fout.close();
        fout.open("stat.txt", std::ofstream::out | std::ofstream::app);
        fout << "time = " << time << " P = " << size << endl;
        fout << a << ' ' << b << ' ' << x << ' ' << p << ' ' << n <<endl;
        fout << "-------------------------------------------------------------------------------------------" << endl;
        fout.close();
    }
    return 0;

}