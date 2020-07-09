#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    const int nx=1500;  // number of data per process
    float *dat;         // local data
    int i;
    float datamin, datamax, datasum, datamean;
    float globmin, globmax, globsum, globmean;
    int ierr;
    int rank, size;
    MPI_Status status;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /*
     * generate random data at each process
     */

    dat = (float *)malloc(nx * sizeof(float));
    srand(rank);
    for (i=0;i<nx;i++) {
        dat[i] = 2*((float)rand()/RAND_MAX)-1.;
    }

    /*
     * find local min/sum/max
     */ 

    datamin = 1e+19;
    datamax =-1e+19;
    datasum = 0;
    
    for (i=0;i<nx;i++) {
        if (dat[i] < datamin) datamin=dat[i];
        if (dat[i] > datamax) datamax=dat[i];
        datasum += dat[i];
    }
    datamean = datasum/nx;
    free(dat);

    printf("Proc %d min/mean/max = %f,%f,%f\n",
           rank, datamin, datamean, datamax);

    /* 
     * combine local results 
     */
    ierr = MPI_Allreduce(&datamin, &globmin, 1, MPI_FLOAT, MPI_MIN,
		    MPI_COMM_WORLD);
    /* to just send to task 0:
    ierr = MPI_Reduce(&datamin, &globmin, 1, MPI_FLOAT, MPI_MIN,
                      0, MPI_COMM_WORLD); 
     */
    ierr = MPI_Allreduce(&datasum, &globsum, 1, MPI_FLOAT, MPI_SUM,
		         MPI_COMM_WORLD);
    ierr = MPI_Allreduce(&datamax, &globmax, 1, MPI_FLOAT, MPI_MAX,
		         MPI_COMM_WORLD);
    globmean = globsum/(size*nx);

    if (rank == 0) {
       printf("Global min/mean/max = %f, %f, %f\n",
	       globmin, globmean, globmax);
    }
 
    ierr = MPI_Finalize();
    return 0;
}
