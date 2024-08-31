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

void init_M (std::vector<double> *M, int N) {
    for (uint64_t i=0; i<N; ++i){ 
        (*M)[i*N+i] = ((double)i+1)/(double)N; 
    }
}

void buffer_broadcast_creator (){ 

}

int main (int argc, char *argv[]){

	MPI_Init(&argc,&argv);

	int n_Nodes; 
	MPI_Comm_size(MPI_COMM_WORLD,&n_Nodes);
	int myId;
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);

	int N = std::stol(argv[1]);
	std::vector<double> M(N * N, 0);
	init_M(&M,N);

	if (myId == 0) { /* task emitter */

		std::printf("Emitter Online\n");

		/* TODO: si manda task */
		/* TODO: si aspetta la ricezione delle task */
		/* TODO: si aspetta i risultati */
		/* TODO: si mandano ripete fino al compimento della diagonale */

		int victim = 0;
		int total_jobs = 0;
		std::vector<MPI_Request> requests_send(n_Nodes-1, 0);
		std::vector<MPI_Request> requests_receive(N, 0);

		for (int i=1; i<N; ++i){ 
        	for (int j=0; j<N-i; ++j){ /* e^k_(i,j) */ 

				if (victim == (n_Nodes-1)) {
					/* gestire il numero non even di worker con la diagonale */
					for (int w=0; w<n_Nodes-1; w++)
						MPI_Wait(&(requests[w]), MPI_STATUS_IGNORE);
					victim = 0;

					/* wait messagges and result */
				}

				total_jobs++;
				MPI_Isend(&buffer_sent, 2, MPI_INT, victim, 0, MPI_COMM_WORLD, &request);
				MPI_iRecv(&(requests_receive[j]),1, MPI_DOUBLE, victim, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				/* send task to worker */

			} /* diagonal completed */
        	for (int j=0; j<N-i; ++j) /*wait for all the diagonal result */
				MPI_Wait(&(requests_receive[j]), MPI_STATUS_IGNORE);

			total_jobs = 0;
			/* TODO:broacast */
			/* TODO: si manda in broadcast i risultati ordinati. */ 

		}

		

		/* il broadcast è come una barriera prima versione con broadcast alla fine, 
		poi si controllare se ha senso aspettare il numero degli worker. */
		/* e si rinizia fino a che la matrice non è finita.*/

      //buffer = 12345;
	  //MPI_Bcast(&buffer, 1, MPI_INT, myId, MPI_COMM_WORLD);
	  //std::printf("1: HO FINITO e vado avanti con la mia vita \n");

		/* TODO: si protocollo tramite task per comunicare al worker di fare task, aspettare risultati, 
		oppure uscire.*/
	} else { /* worker Node */

			


	}
	
	MPI_Finalize();
	return 0;
}