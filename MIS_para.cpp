#include <cstdlib>
#include "mpi.h"
#include <set>
#define INDEX(x,y) (uint64_t)(x * num_nodes + y)

// std::set <int> MIS(double *mat, int num_nodes){

// }

void print_arr(int *array, int len) {
    for (int i = 0; i < len; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

void assign_rand_val(int** nodes, int num_nodes, int num_proc){

    MPI_Init(NULL, NULL);

    int *recv = (int*) malloc(sizeof(int) * num_nodes);
    // Send buffer of length 1
    int *send = (int*) malloc(sizeof(int));
    int rank;
    int size;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(rank);

    // send stores a random value
    send[0] = rand();
    // printf("sending %d\n", send[0]);

    MPI_Allgather(send, 1, MPI_INT, recv, 1, MPI_INT, MPI_COMM_WORLD);

    //printf("This is root and the values are assigned \n");
    free(*nodes);
    *nodes = (int*)malloc(5 * sizeof(int));
    *nodes = recv;
    
    free(send);
    //free(recv);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

}

void add_and_record(double* mat, int* rand_vals, int num_nodes, int num_procs, std::set<int> nodes_to_delete){

    MPI_Init(NULL, NULL);

    int *recv = (int*) malloc(sizeof(int) * num_nodes);
    // Send buffer of length 1
    int *send = (int*) malloc(sizeof(int));
    int rank;
    int size;
    bool flag = 1;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Indivisible, use Scatterv
    int nodes_per_proc = num_nodes / num_procs;
    int nodes_left = num_nodes % num_procs;
    int* send_counts = (int*)malloc(sizeof(int) * (num_procs + 1));
    int* displs = (int*)malloc(sizeof(int) * (num_procs + 1));
    for (int i = 0; i < num_procs + 1; i++)
    {
        if (i < num_procs)
            send_counts[i] = nodes_per_proc;
        else
            send_counts[i] = nodes_left;
        displs[i] = i * nodes_per_proc;
    }
    MPI_Scatterv(mat, send_counts, displs, MPI_DOUBLE, recv, 
        nodes_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < nodes_per_proc; i++)
    {
        for (int j = 0; j < num_nodes; j++)
        {
            if (rand_vals[i] < rand_vals[j] && mat[INDEX(i, j)] > 0)
                flag = 0;
        }
        if (flag)
            printf("%d\n", i);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
}

int main(int argc, char** argv) {
    
    int *rand_vals = (int*) malloc(sizeof(int) * 5);

    assign_rand_val(&rand_vals, 5, 5);

    print_arr(rand_vals, 5);

    free(rand_vals);

}