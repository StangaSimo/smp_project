#include <ff/ff.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ff;

void print_M (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=0; i<N; ++i) {
        std::cout << "\n";
        for (uint64_t j=0; j<N; ++j){
            //std::cout << (*M)[i*(*N)+j] << "   ";
            std::cout << std::setprecision(15) << (*M)[i*N+j] << "   ";
        }
    }
}

/* initialize magior diagonal with e^0_(m,m) = (m+1)/N */
void init_M (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=0; i<N; ++i){ 
                (*M)[i*N+i] = ((double)i+1)/(double)N; 
    }
}

void simple_init (std::vector<double> *M, uint64_t N) {
    int c = 0;
    for (uint64_t i=0; i<N; ++i) {
        for (uint64_t j=0; j<N; ++j){
            (*M)[i*N+j] = c++;
        }
    }
}

/* the seq algorithm */
void compute_seq (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=1; i<N; ++i){ 
        for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */
            
            double res = 0;
            for (uint64_t w=0; w<i; ++w){ /* dot product */

                res += /* row */ (*M)[j*N+j+i-w-1] * (*M)[(j+1+w)*N+j+i] /* column */; 

            }

            (*M)[j*N+j+i] = std::cbrt(res); 
        }
    }
}

struct SourceSink: ff_monode_t<int> {
    std::vector<int> *data;
    int N_workers;
    int N_data;

    SourceSink(int N_workers,int N_data, std::vector<int> *data) {
        this->N_workers = N_workers;
        this->N_data = N_data;
        this->data = data;
    }

    int* svc(int* n) {
        if (n == nullptr) {
            int last = get_num_outchannels();
            std::cout << "allora " << last << "\n";
            //for (int i=0; i<N_data; i++) {
            //    ff_send_out_to(&(*data)[i],0);
            //    std::cout << "mando " << (*data)[i] << "\n";
            //}
            ff_send_out_to(&(*data)[1],0);
            ff_send_out_to(&(*data)[1],1);
            ff_send_out_to(&(*data)[1],2);
            ff_send_out_to(&(*data)[1],3);
			return GO_ON;
		}
        return GO_ON;
    }
};

struct Worker: ff_node_t<int> {
    int id = ff_gettid();
    int* svc(int* n) {
        std::cout << id << ": mi è arrivato :" << *n << "\n";
        return GO_ON;
    }
};

int main(int argc, char *argv[])
{
    std::cout << "*******************************\n";
    std::cout << "*           hello             *\n";
    std::cout << "*******************************\n";
    uint64_t N = 16; /* size of M (NxN) */
    std::vector<double> M(N * N, 0);
    std::vector<double> M_test(N * N, 0); 

    init_M(&M,N);
    init_M(&M_test,N);

    ffTime(START_TIME);
    compute_seq(&M_test,N);
    ffTime(STOP_TIME);
    std::printf("Time %f (ms)\n",ffTime(GET_TIME));

    /***** ff *****/

    int N_data = 10;
    std::vector<int> data(N_data, 1);
    for (int i=0; i<N_data; ++i) {
        std::cout << data[i] << "\n";
    }

    int N_workers = 4; /* size of M (NxN) */

    SourceSink sourceSink(N_workers,N_data,&data);

    ff_Farm<int> farm([&]()
                    {
						std::vector<std::unique_ptr<ff_node> > W;
						for(int i=0;i<4;++i)
						    W.push_back(make_unique<Worker>());
						return W; }()
                        , sourceSink);

    farm.remove_collector(); 
    farm.wrap_around(); 
    farm.set_scheduling_ondemand(); 

    if (farm.run_and_wait_end() < 0){
        error("running farm");
        return -1;
    }

    return EXIT_SUCCESS;
}