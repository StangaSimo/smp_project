#include <cstdio>
#include <mpi.h>

int main (int argc, char *argv[]){
	MPI_Init(&argc,&argv);
	
	int numP; 
	MPI_Comm_size(MPI_COMM_WORLD,&numP);
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
	MPI_Get_processor_name(processor_name,&namelen);
	int myId;
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);
	
	std::printf("From process %d out of %d, running on node %s: Hello, world!\n",
		    myId, numP,processor_name);
	
	// Terminate MPI
	MPI_Finalize();
	return 0;
}