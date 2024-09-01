#include <cstdio>
#include <mpi.h>
#include <unistd.h>
#include <cstdint>
#include <vector>
#include <string>

// https://rookiehpc.org/mpi/docs/mpi_recv/index.html
// https://rookiehpc.org/mpi/docs/mpi_ibcast/index.html
// https://rookiehpc.org/mpi/docs/mpi_isend/index.html
// https://rookiehpc.org/mpi/docs/mpi_irecv/index.html
// https://rookiehpc.org/mpi/docs/mpi_test/index.html

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

	if (myId == 0) { /* task emitter */

		std::printf("Emitter Online\n");

		uint64_t victim = 1;
		std::vector<uint64_t> tasks(N*2, 0);
		std::vector<MPI_Request> requests_send(n_Nodes-1, 0);
		std::vector<MPI_Request> requests_receive(N, 0);
		double result[N-1];

		for (uint64_t i=1; i<N; ++i){ 
        	for (uint64_t j=0; j<N-i; ++j){ /* e^k_(i,j) */ 
				uint64_t index = index_compute(j);
				tasks[index] = i;
				tasks[index+1] = j;
				MPI_Isend(&tasks[index], 2, MPI_UINT64_T, victim, 0, MPI_COMM_WORLD, &(requests_send[victim-1]));
				std::printf("HO MANDATO TASK a %ld\n",victim);

				MPI_Irecv(&result[j], 1, MPI_DOUBLE, victim, 0, MPI_COMM_WORLD, &(requests_receive[j]));
				//MPI_Wait(&requests_receive[j], MPI_STATUS_IGNORE);
				//std::printf("RISULTATO ARRIVATO\n");

				victim++;
				if (victim == (n_Nodes)) {
					for (uint64_t w=0; w<(n_Nodes-1); w++)
						MPI_Wait(&(requests_send[w]), MPI_STATUS_IGNORE);
					std::printf("HO ASPETTATO TUTTI GLI WORKERS\n");
					victim = 1;
				}

				
				
				
			} /* diagonal completed */

			/* if the worker aren't even with number of element in the diagonal */
			/* TODO: test se ha senso o si puo togliere */

			if (victim != 1) {
				for (uint64_t w=0; w<victim-1; w++)
					MPI_Wait(&(requests_send[w]), MPI_STATUS_IGNORE);
	 			std::printf("ASPETTATO ANCHE GLI ULTIMI \n");
				victim = 1;
			} 	
			/* TODO: test i resize se hanno senso o no */

	 		std::printf("ASPETTATO i risultati \n");
        	for (uint64_t j=0; j<N-i; ++j){ /* wait for all the diagonal result */
	 			std::printf("ASPETTO LA NUMERO  \n");
				MPI_Wait(&(requests_receive[j]), MPI_STATUS_IGNORE);
			} 

	 		std::printf("risultati arrivati \n");

        	//for (uint64_t j=0; j<N-i; ++j) /* buffer creation for the broadcast */ 
			//	result[j] = M[j*N+j+i];

	  		MPI_Bcast(&result, N-i, MPI_DOUBLE, myId, MPI_COMM_WORLD);
		}

	} else { /* worker Node */

		int ready_receive;
		int ready_broadcast;
		double res = 99;
		uint64_t n_broadcast = 0;
		uint64_t tasks[2] = {0,0};
		std::vector<double> result(N-n_broadcast, 0);


		while (n_broadcast < N-1) {

			MPI_Request requests_receive;
			MPI_Request requests_broadcast;
			MPI_Irecv(&tasks, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive);
			MPI_Ibcast(&result, N-n_broadcast-1, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests_broadcast);

			while (false) {
				MPI_Test(&requests_receive, &ready_receive, MPI_STATUS_IGNORE);
				if (ready_receive) { /* task arrivata */
					MPI_Request requests_result;
					std::printf("mi è arrivata una task %ld %ld\n", tasks[0], tasks[1]);
					MPI_Isend(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests_result);
					MPI_Wait(&requests_result, MPI_STATUS_IGNORE);
					std::printf("ho mandato il risultato\n");

					/* compute task */	
					MPI_Irecv(&tasks, 2, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD, &requests_receive); /* next tasks */
				}				

				MPI_Test(&requests_broadcast, &ready_broadcast, MPI_STATUS_IGNORE);
				if (ready_broadcast) { 
					std::printf("E' arrivato il BROADCAST\n");
					return 0;
					/* scrive risultati M[j*N+j+i]*/
					n_broadcast++;
					break;
				}
			}
		}
	}
	
	MPI_Finalize();
	return 0;
}