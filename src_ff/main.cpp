#include <ff/ff.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

/* Stampare Matrice */
void print_M (std::vector<double> *M, uint64_t* N) {
    for (uint64_t i=0; i<*N; ++i) {
        std::cout << "\n";
        for (uint64_t j=0; j<*N; ++j){
            //std::cout << (*M)[i*(*N)+j] << "   ";
            std::cout << std::setprecision(15) << (*M)[i*(*N)+j] << "   ";
        }
    }
}

int main(int argc, char *argv[])
{
    //uint64_t N = 512; /* size of M (NxN) */
    uint64_t N = 512; /* size of M (NxN) */

    std::vector<double> M(N * N, 0);
    std::vector<double> M_test(N * N, 0); 

    for (uint64_t i=0; i<N; ++i){ /* initialize magior diagonal */
                M[i*N+i] = ((double)i+1)/(double)N; 
                M_test[i*N+i] = ((double)i+1)/(double)N; 
    }

    //int c = 0;
    //for (uint64_t i=0; i<N; ++i) {
    //    for (uint64_t j=0; j<N; ++j){
    //        M_test[i*N+j] = c++;
    //    }
    //}

    print_M(&M_test,&N);
    std::cout << "\n"; 
    std::cout << "\n"; 
    std::cout << "\n"; 
    
    for (uint64_t i=1; i<N; ++i){ 
        for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */
            
            double res = 0;
            for (uint64_t w=0; w<i; ++w){ /* dot product */

                res += M_test[j*N+j+i-w-1] /* row */ * M_test[(j+1+w)*N+j+i] /* column */; 
            }

            M_test[j*N+j+i] = std::cbrt(res); 
        }
        std::cout << "\n"; 
    }

    print_M(&M_test,&N);
    std::cout << "\n"; 
    std::cout << "\n"; 
    std::cout << "\n"; 
    
    return EXIT_SUCCESS;
}