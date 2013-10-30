#include <fstream>
#include <iostream>
#include <string.h>
#include <mpi.h>

// mpic++ -o  mpiManagementExample mpiManagementExample.cpp 
// mpirun -np 4 --host localhost mpiManagementExample

#define FILE_NAME "fileExample-01.dat"
#define NUMBER 10

int main(int argc,char **argv)
{
	int mpiResult;
	int numtasks;
	int rank;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int len;

	double val[NUMBER];
	int lupe;
	//MPI_Status theStatus;

  // items to write to the file
  struct output 
  {
    double x;
    int    i;
  };
  output basicInfo;
  char   buffer[1024];

  // File related stuffs
  MPI_File mpiFileHandle;


	mpiResult = MPI_Init (&argc,&argv);
	if(mpiResult!= MPI_SUCCESS)
		{
			std::cout << "MPI not started. Terminating the process." << std::endl;
			MPI_Abort(MPI_COMM_WORLD,mpiResult);
		}

	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(hostname, &len);
	//std::cout << "Number of tasks= " <<  numtasks
	//					<< " My rank= " << rank
	//					<< " Running on " << hostname
	//					<< std::endl;

	// Initialize the buffer
	for(lupe=0;lupe<NUMBER;++lupe)
		{
			val[lupe] = (double)(lupe+rank);
		}

  // open the file
  std::cout << "opening " << FILE_NAME << std::endl;
  MPI_Status status;
  char err_buffer[MPI_MAX_ERROR_STRING];
  int resultlen;
  int ierr = MPI_File_open(MPI_COMM_WORLD,FILE_NAME, 
                           MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, 
                           MPI_INFO_NULL, 
                           &mpiFileHandle);
  std::cout << "Open: " << ierr << "," << MPI_SUCCESS << std::endl;
  std::cout << "Status: " << status.MPI_ERROR << "," 
            << status.MPI_SOURCE << "," 
            << status.MPI_TAG << std::endl;
  MPI_Error_string(ierr,err_buffer,&resultlen);
  std::cout << "Error: " << err_buffer << std::endl;


	// print out what was passed
	std::cout << "Process " << rank << " writing: ";
	for(lupe=0;lupe<NUMBER;++lupe)
		std::cout << val[lupe] << ", ";
	std::cout << std::endl;

  // Write out the basic information to the file.
  MPI_File_seek(mpiFileHandle,rank*sizeof(basicInfo),MPI_SEEK_SET);

  // Set the data, copy it to the buffer, and print the results.
  basicInfo.x = val[NUMBER-1];
  basicInfo.i = rank*2+1;
  memset(buffer,0,1024);
  memcpy(buffer,&basicInfo,sizeof(basicInfo));
  std::cout << "rank: " << rank << " moving pointer to "
           << rank*sizeof(basicInfo) 
            << " writing " << buffer << std::endl;

  // write the date at the current pointer
  MPI_File_write(mpiFileHandle,buffer,sizeof(basicInfo)/sizeof(char),
                 MPI_CHAR, &status );

  MPI_File_close(&mpiFileHandle);
	MPI_Finalize();
}
