#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int merge(double *ina, int lena, double *inb, int lenb, double *out) {
    int i,j;
    int outcount=0;

    for (i=0,j=0; i<lena; i++) {
        while ((inb[j] < ina[i]) && j < lenb) {
            out[outcount++] = inb[j++];
        }
        out[outcount++] = ina[i];
    }
    while (j<lenb)
        out[outcount++] = inb[j++];

    return 0;
}

int domerge_sort(double *a, int start, int end, double *b) {
    if ((end - start) <= 1) return 0;

    int mid = (end+start)/2;
    domerge_sort(a, start, mid, b);
    domerge_sort(a, mid,   end, b);
    merge(&(a[start]), mid-start, &(a[mid]), end-mid, &(b[start]));
    for (int i=start; i<end; i++)
        a[i] = b[i];

    return 0;
}

int merge_sort(int n, double *a) {
    double b[n];
    domerge_sort(a, 0, n, b);
    return 0;
}

void printstat(int rank, int iter, char *txt, double *la, int n) {
  //  printf("[%d] %s iter %d: <", rank, txt, iter);
	printf("\nSorted array:\n");
    for (int j=0; j<n-1; j++)
        printf("%6.3lf,",la[j]);
    printf("%6.3lf\n", la[n-1]);
}

void PREZ(int localn, double *locala, int sendrank, int recvrank,
                           MPI_Comm comm) {

   
    int rank;
    double remote[localn];
    double all[2*localn];
    const int mergetag = 1;
    const int sortedtag = 2;

    MPI_Comm_rank(comm, &rank);
    if (rank == sendrank) {
        MPI_Send(locala, localn, MPI_DOUBLE, recvrank, mergetag, MPI_COMM_WORLD);
        MPI_Recv(locala, localn, MPI_DOUBLE, recvrank, sortedtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(remote, localn, MPI_DOUBLE, sendrank, mergetag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        merge(locala, localn, remote, localn, all);

        int theirstart = 0, mystart = localn;
        if (sendrank > rank) {
            theirstart = localn;
            mystart = 0;
        }
        MPI_Send(&(all[theirstart]), localn, MPI_DOUBLE, sendrank, sortedtag, MPI_COMM_WORLD);
        for (int i=mystart; i<mystart+localn; i++)
            locala[i-mystart] = all[i];
    }
}

int PCM(int n, double *a, int root, MPI_Comm comm)
{
    int rank, size, i;
    double *local_a;

// get rank and size of comm
    MPI_Comm_rank(comm, &rank); //&rank = address of rank
    MPI_Comm_size(comm, &size);

    local_a = (double *) calloc(n / size, sizeof(double));


// scatter the array a to local_a
    MPI_Scatter(a, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE,
        root, comm);
// sort local_a
    merge_sort(n / size, local_a);

//odd-even part
    for (i = 1; i <= size; i++) {

       // printstat(rank, i, "before", local_a, n/size);

        if ((i + rank) % 2 == 0) {  // means i and rank have same nature
            if (rank < size - 1) {
                PREZ(n / size, local_a, rank, rank + 1, comm);
            }
        } else if (rank > 0) {
            PREZ(n / size, local_a, rank - 1, rank, comm);
        }

    }

  //  printstat(rank, i-1, "after", local_a, n/size);

// gather local_a to a
    MPI_Gather(local_a, n / size, MPI_DOUBLE, a, n / size, MPI_DOUBLE,
           root, comm);

    if (rank == root)
        printstat(rank, i, " Sorted Array:  ", a, n);

    return MPI_SUCCESS;
}

int main(int argc, char **argv) {
int rank;
    MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int n = argc-1;
    double a[n];
    for (int i=0; i<n; i++)
        a[i] = atof(argv[i+1]);
double t1=MPI_Wtime();
    PCM(n, a, 0, MPI_COMM_WORLD);
double t2=MPI_Wtime();
if(rank==0)
printf("Time taken for the entire computation: %lf ms\n",t2-t1);
    MPI_Finalize();

    return 0;
}
