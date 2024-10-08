#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <ff/ff.hpp>
using namespace ff;

//#define SEQ /* comment for stopping the sequential test */
//#define ONDEMAND 

struct Task_t {
    uint64_t i;
    uint64_t j;
};

std::vector<double> *Mt;
uint64_t N;

void print_M (std::vector<double> *M) {
    for (uint64_t i=0; i<N; ++i) {
        std::cout << "\n";
        for (uint64_t j=0; j<N; ++j){
            //std::cout << (*M)[i*(*N)+j] << "   ";
            std::cout << std::setprecision(15) << (*M)[i*N+j] << "   ";
        }
    }
    std::cout << "\n";
}

/* initialize magior diagonal with e^0_(m,m) = (m+1)/N */
void init_M (std::vector<double> *M) {
    for (uint64_t i=0; i<N; ++i){ 
        (*M)[i*N+i] = ((double)i+1)/(double)N; 
    }
}

void simple_init (std::vector<double> *M) {
    int c = 0;
    for (uint64_t i=0; i<N; ++i) {
        for (uint64_t j=0; j<N; ++j){
            (*M)[i*N+j] = c++;
        }
    }
}

/* the seq algorithm */
void compute_seq (std::vector<double> *M) {
    for (uint64_t i=1; i<N; ++i){ 
        for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */
            
            double res = 0;
            for (uint64_t w=0; w<i; ++w) /* dot product */
                res += /* row */ (*M)[j*N+j+i-w-1] * (*M)[(j+1+w)*N+j+i] /* column */; 

            (*M)[j*N+j+i] = std::cbrt(res); 
        }
    }
}

/* compare Matrix */
bool compare_M (std::vector<double> *A,std::vector<double> *B) {
    bool c = true;
    for (uint64_t i=0; i<N; ++i) 
        for (uint64_t j=0; j<N; ++j)
            if ((*A)[i*N+j] != (*B)[i*N+j])
                c = false;
    return c;
}

struct SourceSink: ff_monode_t<int,int> {  /* in, out */
    std::vector<Task_t> task;
    int N_workers;
    bool stop;
    uint64_t i_c,j_c;
    uint64_t N_jobs_done;
    uint64_t N_jobs_given;
    uint64_t curr_diagonal;

    int svc_init() {
        task.resize(N_workers);
        curr_diagonal = N-1;
        i_c = 1;
        j_c = 0;

        /* for the first loop */
        N_jobs_done = 0; 
        N_jobs_given = N_workers; 

        /* for the last thread that compute the last element */
        stop = false;
        return 0;
    }

    SourceSink(int N_workers) {
        this->N_workers = N_workers;
    }

    int* svc(int* n) { /* out, in */
        if (n == nullptr) { /* initial case */
            for (int i=0; i<N_workers; i++){
                task_maker(&task,i);
                ff_send_out_to(&(task)[i],i);
            }
			return GO_ON;
		} 

        if (stop) 
            return GO_ON;

        ssize_t wid = get_channel_id();
        N_jobs_done++;

        if (N_jobs_done == curr_diagonal){ /* diagonal finished */
            N_jobs_done = 0;
            N_jobs_given = 0;
            curr_diagonal--;

            /* diagonal with less element than workers */
            if ((int)curr_diagonal < N_workers) { 

                /* final case */
                if (curr_diagonal == 1) { 
                    task_maker(&task,0);
                    ff_send_out_to(&(task)[0],0);
                    stop = true;
                    broadcast_task(EOS); 
                    return GO_ON;
                }
                /* send out the right number of tasks */
                for (int i=0; i<(int)curr_diagonal; i++) { 
                    task_maker(&task,i);
                    ff_send_out_to(&(task)[i],i);
                    N_jobs_given++;
                }
            } else { /* new diagonal */
                for (int i=0; i<N_workers; i++) {
                    task_maker(&task,(int)i);
                    ff_send_out_to(&(task)[i],i);
                    N_jobs_given++;
                } 
            }
        /* ready worker will go on with the current diagonal */
        } else { 
            /* wait if the diagonal isn't complete but the jobs has been sent */
            if (N_jobs_given < curr_diagonal) { 
                N_jobs_given++;
                task_maker(&task,(int)wid);
                ff_send_out_to(&(task)[wid],wid);
            } 
        }
        return GO_ON;
    }
    
    void task_maker(std::vector<Task_t> *task,int i) { 
        if (j_c >= N-i_c) {
            i_c++; 
            j_c = 0; 
        }
        (*task)[i].i = i_c;
        (*task)[i].j = j_c;
        j_c++;
    }
};

struct Worker: ff_node_t<Task_t, int> { /* in, out */
    int u = 1;
    int* svc(Task_t* task) {  /* out, in */
        task_compute(task);
        ff_send_out(&u);
        return GO_ON;
    }

    void task_compute(Task_t *n) {
        double res = 0;
        for (uint64_t w=0; w<(n->i); ++w)
            res += (*Mt)[(n->j)*N+(n->j)+(n->i)-w-1] * (*Mt)[((n->j)+1+w)*N+(n->j)+(n->i)]; 
        (*Mt)[n->j*N+n->j+n->i] = std::cbrt(res);
    }
};

int main(int argc, char *argv[])
{
    N =std::stol(argv[1]);
    int N_workers =std::stol(argv[2]);

    std::cout << "\nCompute a matrix " << N << " * " << N << " with " << N_workers << " threads\n" ; 

    std::vector<double> M(N * N, 0);
    init_M(&M);

#ifdef SEQ
    std::vector<double> M_test(N * N, 0); 
    init_M(&M_test);
    ffTime(START_TIME);
    compute_seq(&M_test);
    ffTime(STOP_TIME);
    std::printf("\nCompute Seq Time %f (ms)\n",ffTime(GET_TIME));
#endif

    /******************** ff ********************/

    Mt = &M;

    SourceSink sourcesink(N_workers);

    ff_Farm<int> farm([&]()
                    {
						std::vector<std::unique_ptr<ff_node> > W;
						for(int i=0;i<N_workers;++i)
						    W.push_back(make_unique<Worker>());
						return W; }()
                        , sourcesink);

    farm.remove_collector(); 
    farm.wrap_around(); 

#ifdef ONDEMAND
    farm.set_scheduling_ondemand(); 
#endif

    ffTime(START_TIME);
    if (farm.run_and_wait_end() < 0){
        error("running farm");
        return -1;
    }
    ffTime(STOP_TIME);
    std::printf("Compute Parallel Time %f (ms)\n",ffTime(GET_TIME));

#ifdef SEQ
    if (compare_M(Mt,&M_test)) 
        std::cout << "\nTest Passed " << "\n";
    else
        std::cout << "\nTest Failed " << "\n";
#endif

    return EXIT_SUCCESS;
}