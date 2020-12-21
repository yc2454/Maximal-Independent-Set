#include <cstdlib>
#include "mpi.h"
#include <set>
#include <vector>
#include <iostream>
#define INDEX(x,y) (uint64_t)(x * num_nodes + y)

// std::set <int> MIS(double *mat, int num_nodes){

// }

void print_arr(double *array, int len) {
    for (int i = 0; i < len; i++)
    {
        printf("%f ", array[i]);
    }
    printf("\n");
}

// void assign_rand_val(int** nodes, int num_nodes, int num_proc){

//     MPI_Init(NULL, NULL);

//     int *recv = (int*) malloc(sizeof(int) * num_nodes);
//     // Send buffer of length 1
//     int *send = (int*) malloc(sizeof(int));
//     int rank;
//     int size;

//     // Get rank and size of this process
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     srand(rank);

//     // send stores a random value
//     send[0] = rand();
//     // printf("sending %d\n", send[0]);

//     MPI_Allgather(send, 1, MPI_INT, recv, 1, MPI_INT, MPI_COMM_WORLD);

//     //printf("This is root and the values are assigned \n");
//     free(*nodes);
//     *nodes = (int*)malloc(5 * sizeof(int));
//     *nodes = recv;
    
//     free(send);
//     //free(recv);

//     MPI_Barrier(MPI_COMM_WORLD);
//     MPI_Finalize();

// }

void assign_rand_vals(std::vector<int> &rand_vals, int num_nodes){
    
    for (int i = 0; i < num_nodes; i++)
    {
        rand_vals.push_back(rand() % 10000);
    }
    
}

void add_and_record(double* mat, std::vector<int> rand_vals, int num_nodes, int num_procs, std::set<int> &active_nodes, std::set<int> &M){

    //MPI_Init(NULL, NULL);
    int rank;
    int size;
    bool flag = 1;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nodes_left = num_nodes % num_procs;
    int nodes_per_proc = (num_nodes - nodes_left) / num_procs;
    double *recv = (double*) malloc(sizeof(double) * num_nodes * nodes_per_proc);
    int len_grp_size = (nodes_left == 0) ? num_procs : (num_procs + 1);
    int* send_counts = (int*)malloc(sizeof(int) * len_grp_size);
    int* displs = (int*)malloc(sizeof(int) * len_grp_size);

    for (int i = 0; i < len_grp_size; i++)
    {
        if (nodes_left)
        {
            if (i < len_grp_size - 1)
                send_counts[i] = nodes_per_proc * num_nodes;
            else
                send_counts[i] = nodes_left * num_nodes;
        }
        else
            send_counts[i] = nodes_per_proc * num_nodes;

        displs[i] = i * nodes_per_proc * num_nodes;
    }

    // print_arr(send_counts, num_procs);
    // print_arr(displs, num_procs);
    // MPI_Scatterv(mat, send_counts, displs, MPI_DOUBLE, recv, 
    //     nodes_per_proc * num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(mat, nodes_per_proc * num_nodes, MPI_DOUBLE, recv, nodes_per_proc * num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::cout << "Current rank is: " << rank << std::endl;
    print_arr(recv, nodes_per_proc * num_nodes);

    // SEG FAULT here
    for (int i = 0; i < nodes_per_proc; i++)
    {
        //check whether node i has bigger random value than all its neighbors
        for (int j = 0; j < num_nodes; j++)
        {
            if (rand_vals[i] < rand_vals[j] && recv[INDEX(i, j)] > 0)
                flag = 0;
        }
        //if it does, add i to M and delete i and its neighbors from A
        if (flag){
            for (int j = 0; j < num_nodes; j++)
            {
                if(recv[INDEX(i, j)] > 0)
                    active_nodes.erase(j);
            }
            active_nodes.erase(i);
            M.insert(i);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Finalize();
    
}

std::set<int> MIS(double* mat, int num_nodes, int num_procs) {

    std::set<int> A, M, M_temp;

    for (int i = 0; i < num_nodes; i++)
        A.insert(i);
    
    std::vector<int> rand_vals;

    MPI_Init(NULL, NULL);

    // int count;

    while (!A.empty())
    {
        //rand_vals = (int*) malloc(sizeof(int) * A.size());
        //rand_vals = assign_rand_vals(5);
        assign_rand_vals(rand_vals, num_nodes);

        M_temp.clear();

        for (std::set<int>::const_iterator i = A.begin(); i != A.end(); ++i)
            add_and_record(mat, rand_vals, num_nodes, num_procs, A, M);

        rand_vals.clear();

        // count++;
        // if(count > 100){
        //     std::cout << "infinite loop :(\n";
        //     break;
        // }
    }

    MPI_Finalize();

    return M;

}

int main(int argc, char* argv[]){
    
    double *test;
    test = new double [36];

    for (int i = 0; i < 6; i++)
        test[i*6 + i] = 0;
    
    test[1] = 0;
    for (int i = 2; i < 6; i++)
        test[i] = 0;

    test[6] = 1;
    test[8] = 1;
    test[9] = 1;
    test[10] = 0;
    test[11] = 1;

    test[12] = 0;
    test[13] = 1;
    test[15] = 0;
    test[16] = 1;
    test[17] = 1;

    for (int i = 18; i < 30; i++)
        test[i] = 0;
    test[19] = 1;
    test[26] = 1;
    test[29] = 1;

    for (int i = 30; i < 35; i++)
        test[i] = 1;
    test[30] = 0;
    test[33] = 0;

    std::set<int> I = MIS(test, 6, 6);

    for(std::set<int>::const_iterator i = I.begin(); i != I.end(); i++)
        std::cout << *i << ' ';

}