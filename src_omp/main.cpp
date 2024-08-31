#include <cstdio>
#include <mpi.h>

/* idea:
 * si usa un emitter per le task, tutti fanno il broadcast dei risultati nel broadcast, poi si dice ready all'emitter */

int main (int argc, char *argv[]){
	MPI_Init(&argc,&argv);
	
	int number_of_Nodes; 
	MPI_Comm_size(MPI_COMM_WORLD,&number_of_Nodes);
	int myId;
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);
	
	std::printf("From process %d out of %d, Hello, world!\n",
		    myId, number_of_Nodes);
	
	MPI_Finalize();
	return 0;
}