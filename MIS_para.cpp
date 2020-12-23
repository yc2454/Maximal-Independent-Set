#include <cstdlib>
#include "mpi.h"
#include <set>
#include <vector>
#include <iostream>
#include <ctime>
//#define INDEX(x,y) (uint64_t)(x * num_nodes + y)

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

//     int *mat_recv = (int*) malloc(sizeof(int) * num_nodes);
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

//     MPI_Allgather(send, 1, MPI_INT, mat_recv, 1, MPI_INT, MPI_COMM_WORLD);

//     //printf("This is root and the values are assigned \n");
//     free(*nodes);
//     *nodes = (int*)malloc(5 * sizeof(int));
//     *nodes = mat_recv;
    
//     free(send);
//     //free(mat_recv);

//     MPI_Barrier(MPI_COMM_WORLD);
//     MPI_Finalize();

// }

void assign_rand_vals(std::vector<int> &rand_vals, int num_nodes){
    
    for (int i = 0; i < num_nodes; i++)
    {
        rand_vals.push_back(rand() % 10000);
    }
    
}

void add_and_record(double* mat, std::vector<int> rand_vals, 
                    int num_nodes, int num_procs, std::set<int> &neighbors, 
                    std::set<int> &M, std::set<int> &A)
{
    int rank;
    int size;
    bool flag = 1;

    // Translate A into an array, which will be the send buffer
    int* alive = (int*) malloc(sizeof(int) * A.size());
    int idx = 0;
    if (rank == 0) {
        for(std::set<int>::const_iterator a = A.begin(); a != A.end(); a++){
            alive[idx] = *a;
            idx++;
        }
    }

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // For MPI_Scatterv
    int nodes_left = num_nodes % num_procs;
    int nodes_per_proc = (num_nodes - nodes_left) / num_procs;
    int len_grp_size = (nodes_left == 0) ? num_procs : (num_procs + 1);
    int* send_counts = (int*)malloc(sizeof(int) * len_grp_size);
    int* displs = (int*)malloc(sizeof(int) * len_grp_size);

    // Buffer to hold the nodes scattered to each process
    double *mat_recv = (double*) malloc(sizeof(double) * num_nodes * nodes_per_proc);
    int *sub_alive = (int*) malloc(sizeof(int) * nodes_per_proc);
    
    // The node in consideration
    int cur_node;

    // Compute send_counts and displacements
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
    // MPI_Scatterv(mat, send_counts, displs, MPI_DOUBLE, mat_recv, 
    //     nodes_per_proc * num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(mat, nodes_per_proc * num_nodes, MPI_DOUBLE, mat_recv, nodes_per_proc * num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // print_arr(mat_recv, nodes_per_proc * num_nodes);

    // Buffer for nodes to add
    int* sub_M = (int*)malloc(sizeof(int) * nodes_per_proc);
    int* sub_nbrs = (int*)malloc(sizeof(int) * num_nodes);
    int count_M = 0;
    int count_ngb = 0;

    // Luby's Algo
    for (int i = 0; i < nodes_per_proc; i++)
    {
        cur_node = rank * nodes_per_proc + i;
        //check whether node i has bigger random value than all its neighbors
        for (int j = 0; j < num_nodes; j++)
        {
            if (rand_vals[cur_node] < rand_vals[j] && mat_recv[i * num_nodes + j] > 0)
                flag = 0;
        }
        //if it does, add i to M'
        if (flag){
            for (int j = 0; j < num_nodes; j++)
            {
                // Delete the neighbors of the current node
                if(mat_recv[i * num_nodes + j] > 0) {
                    neighbors.insert(j);
                    sub_nbrs[count_ngb] = j;
                    count_ngb++;
                }
            }
            // active_nodes.erase(cur_node);
            M.insert(cur_node);
            sub_M[count_M] = cur_node;
            count_M++;
        }
    }

    if (rank == 1) {
        for(std::vector<int>::const_iterator i = rand_vals.begin(); i != rand_vals.end(); i++)
            std::cout << *i << ' ';
        std::cout << std::endl;
        for(std::set<int>::const_iterator i = M.begin(); i != M.end(); i++)
            std::cout << *i << ' ';
        std::cout << "*" << std::endl;
        for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
            std::cout << *i << ' ';
        std::cout << "*" << std::endl;
    }

    // should gather deleting nodes here, use MPI_Gatherv
    // Buffer for the nodes to add
    int* all_M = NULL;
    if (rank == 0)
        all_M = (int*)malloc(sizeof(int) * num_nodes);
    MPI_Gather(sub_M, 1, MPI_INT, all_M, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // The neighbors of the nodes
    int* all_nbrs = NULL;
    if (rank == 0)
        all_nbrs = (int*)malloc(sizeof(int) * num_nodes);
    MPI_Gather(sub_nbrs, 1, MPI_INT, all_nbrs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Sum up the size of data
    int all_count_M;
    int all_count_ngb;
    MPI_Reduce(&count_M, &all_count_M, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&count_ngb, &all_count_ngb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    // std::cout << "size of M' this round: " << all_count_M << std::endl;

    // Add gathered data into M and its neighbors
    if(rank == 0) {
        for (int i = 0; i < all_count_M; i++)
            M.insert(all_M[i]);
        for (int i = 0; i < all_count_ngb; i++)
            neighbors.insert(all_nbrs[i]);
    }

    free(mat_recv);
    free(send_counts);
    free(displs);
    if (rank == 0){
        free(all_nbrs);
        free(all_M);
    }
    free(sub_M);
    free(sub_nbrs);
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Finalize();
    
}

std::set<int> MIS(double* mat, int num_nodes, int num_procs) {

    std::set<int> A, M, M_temp, neighbors;

    for (int i = 0; i < num_nodes; i++)
        A.insert(i);
    
    std::vector<int> rand_vals;

    //MPI_Init(NULL, NULL);

    int count;

    while (!A.empty())
    {
        //rand_vals = (int*) malloc(sizeof(int) * A.size());
        //rand_vals = assign_rand_vals(5);
        
        assign_rand_vals(rand_vals, num_nodes);
        M_temp.clear();
        neighbors.clear();
        // Add fitting nodes to M'
        add_and_record(mat, rand_vals, num_nodes, num_procs, neighbors, M_temp, A);
        // M = union(M, M')
        M.insert(M_temp.begin(), M_temp.end());
        // Delete M' and neighbors from A
        for(std::set<int>::const_iterator i = M_temp.begin(); i != M_temp.end(); i++)
            A.erase(*i);
        for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
            A.erase(*i);
        rand_vals.clear();

        // std::cout << "active nodes in this round:" << std::endl;
        // for(std::set<int>::const_iterator i = A.begin(); i != A.end(); i++)
        //     std::cout << *i << ' ';
        // std::cout << std::endl;

        count++;
        if(count > 10){
            break;
        }
    }

    //MPI_Finalize();

    return M;

}

void erase_test(std::set<int> &S) {
    S.erase(3);
}

int main(int argc, char* argv[]){

    MPI_Init(NULL, NULL);
    
    srand((unsigned) time(0));
    
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

    if (!I.empty())
        for(std::set<int>::const_iterator i = I.begin(); i != I.end(); i++)
            std::cout << *i << ' ';

    MPI_Finalize();

}