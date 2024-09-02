#include <cstdio>
#include <mpi.h>
#include <unistd.h>
#include <cstdint>
#include <vector>
#include <string>
#include <cmath>

#define SEQ 

struct Task_t {
	uint64_t i;
    uint64_t j;
};

void print_M (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=0; i<N; ++i) {
        std::printf("\n");
        for (uint64_t j=0; j<N; ++j){
            std::printf("%f  ",(*M)[i*N+j]);
        }
    }
    std::printf("\n");
}

uint64_t index_compute (uint64_t n) {
	uint64_t res = 0;
    for (uint64_t i=0; i<n; ++i)
		res+=2;
	return res;
}

void init_M (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=0; i<N; ++i){ 
        (*M)[i*N+i] = ((double)i+1)/(double)N; 
    }
}

bool compare_M (std::vector<double> *A,std::vector<double> *B, uint64_t N) {
    bool c = true;
    for (uint64_t i=0; i<N; ++i) 
        for (uint64_t j=0; j<N; ++j)
            if ((*A)[i*N+j] != (*B)[i*N+j])
                c = false;
    return c;
}

void compute_seq (std::vector<double> *M, uint64_t N) {
    for (uint64_t i=1; i<N; ++i){ 
        for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */
            
            double res = 0;
            for (uint64_t w=0; w<i; ++w) /* dot product */
                res += /* row */ (*M)[j*N+j+i-w-1] * (*M)[(j+1+w)*N+j+i] /* column */; 

            (*M)[j*N+j+i] = std::cbrt(res); 
        }
    }
}

int main (int argc, char *argv[]){

	MPI_Init(&argc,&argv);

	int n_No; 
	MPI_Comm_size(MPI_COMM_WORLD,&n_No);
	uint64_t n_Nodes = (uint64_t)n_No;
	int myId;
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);

	uint64_t N = std::stol(argv[1]);




#ifdef SEQ 
	std::vector<double> M_test(N * N, 0);
	init_M(&M_test,N);
	compute_seq(&M_test,N);
#endif

	if (myId == 0) { /* task emitter */
		std::vector<double> M(N * N, 0);
		init_M(&M,N);

		//std::printf("Emitter Online\n");

		uint64_t victim = 1;
		uint64_t tasks[N*2];
		std::vector<MPI_Request> requests_send(n_Nodes-1, 0);
		std::vector<MPI_Request> requests_receive(N, 0);
		std::vector<double> result(N-1, 0);
		//double result[N-1];

		for (uint64_t i=1; i<N; ++i){ 
        	for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */ 
				uint64_t index = index_compute(j);
				tasks[index] = i;
				tasks[index+1] = j;
				MPI_Isend(&tasks[index], 2, MPI_UINT64_T, victim, 0, MPI_COMM_WORLD, &requests_send[victim-1]);
				//std::printf("HO MANDATO TASK a %ld\n",victim);

				MPI_Irecv(&result[j], 1, MPI_DOUBLE, victim, 0, MPI_COMM_WORLD, &requests_receive[j]);

				victim++;

				if (victim == n_Nodes) {
					for (uint64_t w=0; w<(n_Nodes-1); w++)
						MPI_Wait(&requests_send[w], MPI_STATUS_IGNORE);
					//std::printf("HO ASPETTATO TUTTI GLI WORKERS\n");
					victim = 1;
				}
			} /* diagonal completed */

			/* if the worker aren't even with number of element in the diagonal */
			/* TODO: test se ha senso o si puo togliere */

			if (victim != 1) {
				for (uint64_t w=0; w<victim-1; w++)
					MPI_Wait(&requests_send[w], MPI_STATUS_IGNORE);
	 			//std::printf("ASPETTATO ANCHE GLI ULTIMI \n");
				victim = 1;
			} 	

			/* TODO: test i resize se hanno senso o no */

        	for (uint64_t j=0; j<N-i; ++j) /* wait for all the diagonal result */
				MPI_Wait(&requests_receive[j], MPI_STATUS_IGNORE);
	 		//std::printf("RISULTATI arrivati \n");

        	for (uint64_t j=0; j<N-i; ++j) {/* buffer creation for the broadcast */ 
				M[j*N+j+i] = result[j];
			}

			MPI_Request requests_broadcast;
			MPI_Ibcast(&result[0], N-i, MPI_DOUBLE, myId, MPI_COMM_WORLD, &requests_broadcast);
			MPI_Wait(&requests_broadcast, MPI_STATUS_IGNORE);
			MPI_Barrier(MPI_COMM_WORLD);

			/*TODO: maybe barrirer??*/

			result.clear();
			result.resize(N-i-1);
	 		//std::printf("BROADCAST ARRIVATO \n");
		}

	 	//std::printf("EMITTER ESCE \n");
#ifdef SEQ
    if (compare_M(&M,&M_test,N)) 
        std::printf("Test Passed\n");
    else
        std::printf("Test Failed\n");
#endif

	} else { /* worker Node */
	 	std::vector<double> M(N * N, 0);
		init_M(&M,N);

		int ready_receive;
		int ready_broadcast;
		bool c = true;
		double res = 0;
		uint64_t n_broadcast = 0;
		uint64_t task[2] = {0,0};
    	//double result[N-1];
		std::vector<double> result(N-1, 0);

		while (n_broadcast < N-1) {
			MPI_Request requests_receive;
			MPI_Request requests_broadcast;
			if (c) {
				MPI_Irecv(task, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive);
				c = true;
			}
			MPI_Ibcast(&result[0], N-n_broadcast, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests_broadcast);

			while (true) {
				MPI_Test(&requests_broadcast, &ready_broadcast, MPI_STATUS_IGNORE);
				if (ready_broadcast) { 
					++n_broadcast;

					if (myId == 1)	
						std::printf("\n\n%d e il numero  %ld N = %ld \n",myId,n_broadcast,N);
					for (uint64_t j=0; j<N-n_broadcast; ++j) {
						M[j*N+j+n_broadcast] = result[j];
						if (myId == 1)	
							std::printf("\n\n%d e j = %ld e scrive %f \n",myId,j,result[j]);
					}

					MPI_Barrier(MPI_COMM_WORLD); /* mitigate the broadcast bug */
					result.clear();
					result.resize(N-n_broadcast-1);

					
					//std::printf("%d E' arrivato il broadcast %ld\n",myId, n_broadcast);
					break;
				}

				MPI_Test(&requests_receive, &ready_receive, MPI_STATUS_IGNORE);
				if (ready_receive) { /* task arrived */
					MPI_Request requests_result;
					uint64_t i = task[0];
					uint64_t j = task[1];

        			for (uint64_t w=0; w<i; ++w)
        			    res += M[j*N+j+i-w-1] * M[(j+1+w)*N+j+i]; 
        			res = std::cbrt(res);
					if (myId == 1)	
						print_M(&M,N);
					std::printf("%d mi è arrivata una task %ld %ld e il risultato %f\n",myId, task[0], task[1],res);
					MPI_Isend(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests_result);
					MPI_Wait(&requests_result, MPI_STATUS_IGNORE);
					res = 0;
					//std::printf("%d ho mandato risultato \n",myId);
					/* next possible task */
					MPI_Irecv(task, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive);
					c = false; 
				}

				
			}
		}
#ifdef SEQ
    if (compare_M(&M,&M_test,N)) 
        std::printf("Test Passed\n");
    else
        std::printf("Test Failed\n");
#endif
		//std::printf("\n\nSONO USCITO\n\n");
	}

	//print_M(&M,N);
 
	if (myId == 0)	{
		std::printf("\n\nM TEST \n\n");
		print_M(&M_test,N);
	}



	MPI_Finalize();
	return 0;
}