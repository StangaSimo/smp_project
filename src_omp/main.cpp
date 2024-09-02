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

bool compare_M (std::vector<double> *A,std::vector<double> *B) {
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
	std::vector<double> M(N * N, 0);

	init_M(&M,N);

#ifdef SEQ 
	compute_seq(&M,N);
#endif

	if (myId == 0) { /* task emitter */

		//std::printf("Emitter Online\n");

		uint64_t victim = 1;
		uint64_t tasks[N*2];
		std::vector<MPI_Request> requests_send(n_Nodes-1, 0);
		std::vector<MPI_Request> requests_receive(N, 0);
		double result[N-1];

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

        	for (uint64_t j=0; j<N-i; ++j) /* buffer creation for the broadcast */ 
				M[j*N+j+i] = result[j];

			MPI_Request requests_broadcast;
			MPI_Ibcast(result, N-i, MPI_DOUBLE, myId, MPI_COMM_WORLD, &requests_broadcast);
			MPI_Wait(&requests_broadcast, MPI_STATUS_IGNORE);
	 		//std::printf("BROADCAST ARRIVATO \n");
		}

	 	//std::printf("EMITTER ESCE \n");

	} else { /* worker Node */

		int ready_receive;
		int ready_broadcast;
		bool c = true;
		double res = 99;
		uint64_t n_broadcast = 0;
		uint64_t task[2] = {0,0};
    	double result[N-1];

		while (n_broadcast < N-1) {
			MPI_Request requests_receive;
			MPI_Request requests_broadcast;
			if (c) {
				MPI_Irecv(task, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive);
				c = true;
			}
			MPI_Ibcast(result, N-n_broadcast, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests_broadcast);

			while (true) {
				MPI_Test(&requests_receive, &ready_receive, MPI_STATUS_IGNORE);
				if (ready_receive) { /* task arrived */
					MPI_Request requests_result;

					uint64_t i = task[0];
					uint64_t j = task[1];
        			for (uint64_t w=0; w<i; ++w)
        			    res += M[j*N+j+i-w-1] * M[(j+1+w)*N+j+i]; 
        			M[j*N+j+i] = std::cbrt(res);
					res = M[j*N+j+i];

					//std::printf("%d mi è arrivata una task %ld %ld\n",myId, tasks[0], tasks[1]);
					MPI_Isend(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests_result);
					MPI_Wait(&requests_result, MPI_STATUS_IGNORE);
					//std::printf("%d ho mandato risultato \n",myId);
					/* next possible task */
					MPI_Irecv(task, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive);
					c = false; 
				}				

				MPI_Test(&requests_broadcast, &ready_broadcast, MPI_STATUS_IGNORE);
				if (ready_broadcast) { 
					n_broadcast++;

					for (uint64_t j=0; j<N-n_broadcast; j++) 
						M[j*N+j+n_broadcast] = result[j];

					//std::printf("%d E' arrivato il broadcast %ld\n",myId, n_broadcast);
					break;
				}
			}
		}
		//std::printf("\n\nSONO USCITO\n\n");
	}
	
	//std::printf("\n\nTUTTI USCITI\n\n");
	MPI_Finalize();
	return 0;
}