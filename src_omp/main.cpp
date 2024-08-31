#include <cstdio>
#include <mpi.h>
#include <unistd.h>
#include <cstdint>

// https://rookiehpc.org/mpi/docs/mpi_recv/index.html
// https://rookiehpc.org/mpi/docs/mpi_ibcast/index.html
// https://rookiehpc.org/mpi/docs/mpi_isend/index.html
// https://rookiehpc.org/mpi/docs/mpi_irecv/index.html
// https://rookiehpc.org/mpi/docs/mpi_test/index.html

void init_M (std::vector<double> *M) {
    for (uint64_t i=0; i<N; ++i){ 
        (*M)[i*N+i] = ((double)i+1)/(double)N; 
    }
}

int main (int argc, char *argv[]){

	MPI_Init(&argc,&argv);

	int n_Nodes; 
	MPI_Comm_size(MPI_COMM_WORLD,&n_Nodes);
	int myId;
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);

	int N = std::stol(argv[1]);
	std::vector<double> M(N * N, 0);
	init_M(&M);

	if (myId == 0) { /* task emitter */

		std::printf("Emitter Online\n");
		
		/* TODO: si manda task */
		/* TODO: si aspetta la ricezione delle task */
		/* TODO: si aspetta i risultati */
		/* TODO: si mandano ripete fino al compimento della diagonale */
		/* TODO: si manda in broadcast i risultati ordinati. */ 

		/* il broadcast è come una barriera prima versione con broadcast alla fine, 
		poi si controllare se ha senso aspettare il numero degli worker. */

		/* e si rinizia fino a che la matrice non è finita.*/

      //buffer = 12345;
	  //MPI_Bcast(&buffer, 1, MPI_INT, myId, MPI_COMM_WORLD);
	  //std::printf("1: HO FINITO e vado avanti con la mia vita \n");

		/* TODO: si protocollo tramite task per comunicare al worker di fare task, aspettare risultati, oppure uscire. */
	} else { /* worker Node */

			


	}
	
	MPI_Finalize();
	return 0;
}