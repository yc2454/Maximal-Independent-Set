#include <cstdlib>
#include "mpi.h"
#include <set>
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>

void print_arr(int *array, int len) {
    for (int i = 0; i < len; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

void assign_rand_vals(std::vector<int> &rand_vals, int num_nodes){
    
    for (int i = 0; i < num_nodes; i++)
    {
        rand_vals.push_back(rand() % 10000);
    }
    
}

void add_and_record(int rank, double* mat, std::vector<int> rand_vals, 
                    int num_nodes, int num_procs, std::set<int> &neighbors, 
                    std::set<int> &M, int* alive, int num_alive)
{
    
    bool flag = 1;

    // Number of nodes the first (num_procs - 1) should hold
    int nodes_per_proc = (num_alive >= num_procs) ? round((1.0 * num_alive) / num_procs) : 1;
    // Number of nodes the last process should hold
    int nodes_left = num_alive - nodes_per_proc * (num_procs - 1);

    // For MPI_Scatterv
    int* send_counts = (int*)malloc(sizeof(int) * num_procs);
    int* displs = (int*)malloc(sizeof(int) * num_procs);
    
    // The node in consideration
    int cur_node;

    // Compute send_counts and displacements
    for (int i = 0; i < num_procs; i++)
    {
        if (num_alive >= num_procs)
        {
            if (i < num_procs - 1)
                send_counts[i] = nodes_per_proc;
            else
                send_counts[i] = nodes_left;
            displs[i] = i * nodes_per_proc;
        }
        else
        {
            send_counts[i] = i < num_alive ? 1 : 0;
            if (i == 0)
                displs[i] = 0;
            else
                displs[i] = i < num_alive ? (displs[i-1] + 1) : displs[i-1];
        }
    }
    
    int recving = (send_counts[rank] < 0 || rank >= num_procs) ? 0 : send_counts[rank];

    // Buffer to hold the nodes scattered to each process
    int *sub_alive = (int*) malloc(sizeof(int) * recving);
    // std::cout << "rank" << rank << ' ' << recving << " receiving " << std::endl;
    MPI_Scatterv(alive, send_counts, displs, 
        MPI_INT, sub_alive, recving, MPI_INT, 0, MPI_COMM_WORLD);
    
    // MPI_Barrier(MPI_COMM_WORLD);

    // Buffer for nodes to add
    int* sub_M = (int*) malloc(sizeof(int) * recving);
    int* sub_nbrs = (int*) malloc(sizeof(int) * num_nodes);
    int count_M = 0;
    int count_ngb = 0;

    // Luby's Algo
    if (recving > 0)
        for (int i = 0; i < send_counts[rank]; i++)
        {
            cur_node = sub_alive[i];
            //check whether node i has bigger random value than all its neighbors
            for (int j = 0; j < num_nodes; j++)
            {
                if (rand_vals[cur_node] < rand_vals[j] && mat[cur_node * num_nodes + j] > 0)
                    flag = 0;
            }
            //if it does, add i to M'
            if (flag){
                for (int j = 0; j < num_nodes; j++)
                {
                    // Delete the neighbors of the current node
                    if(mat[cur_node * num_nodes + j] > 0 
                        /*(A.find(j) != A.end())*/) {
                        if (neighbors.find(j) == neighbors.end()){
                            sub_nbrs[count_ngb] = j;
                            count_ngb++;
                        }
                        neighbors.insert(j);
                    }
                }
                // active_nodes.erase(cur_node);
                M.insert(cur_node);
                sub_M[count_M] = cur_node;
                count_M++;
            }
        }
    
    // Printing status
    for(std::vector<int>::const_iterator i = rand_vals.begin(); i != rand_vals.end(); i++)
            std::cout << *i << ' ';
    std::cout << ' ' << rank << std::endl;
    for(std::set<int>::const_iterator i = M.begin(); i != M.end(); i++)
        std::cout << *i << ' ';
    std::cout << "*" << std::endl;
    for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
        std::cout << *i << ' ';
    std::cout << "*" << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD);

    // Buffer for the nodes to add
    int* all_M = NULL;
    if (rank == 0)
        all_M = (int*)malloc(sizeof(int) * num_nodes);

    int sending = M.size();

    // Gather info for recv_counts and displs, arguments for Gatherv
    int recv_counts[num_procs];
    MPI_Gather(&sending, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, 
        MPI_COMM_WORLD);
    displs[0] = 0;
    for (int i = 1; i < num_procs; i++)
        displs[i] = displs[i-1] + recv_counts[i-1];
    
    // Send M' in each process to root
    MPI_Gatherv(sub_M, sending, MPI_INT, all_M, 
        recv_counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // The neighbors of the nodes
    int* all_nbrs = NULL;
    if (rank == 0)
        all_nbrs = (int*)malloc(sizeof(int) * num_nodes);
    
    sending = neighbors.size();
    MPI_Gather(&sending, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, 
        MPI_COMM_WORLD);
    for (int i = 1; i < num_procs; i++)
        displs[i] = displs[i-1] + recv_counts[i-1];
    
    MPI_Gatherv(sub_nbrs, sending, MPI_INT, all_nbrs, 
        recv_counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    // Sum up the size of data
    int all_count_M;
    int all_count_ngb;
    MPI_Reduce(&count_M, &all_count_M, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&count_ngb, &all_count_ngb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);

    // Add gathered data into M and its neighbors
    if(rank == 0) {
        for (int i = 0; i < all_count_M; i++)
            M.insert(all_M[i]);
        for (int i = 0; i < all_count_ngb; i++)
            neighbors.insert(all_nbrs[i]);
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    free(sub_alive);
    free(send_counts);
    free(displs);
    if (rank == 0){
        free(all_nbrs);
        free(all_M);
    }
    free(sub_nbrs);
    free(sub_M);

    MPI_Barrier(MPI_COMM_WORLD);

}

std::set<int> MIS(double* mat, int num_nodes, int num_procs) {

    std::set<int> A, M, M_temp, neighbors;

    int rank;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < num_nodes; i++)
        A.insert(i);

    // Translate the set into an array
    int* alive = (int*) malloc(sizeof(int) * A.size());
    int idx = 0;
    int num_alive = A.size();
    for(std::set<int>::const_iterator a = A.begin(); a != A.end(); a++){
        alive[idx] = *a;
        idx++;
    }
    
    std::vector<int> rand_vals;

    bool root_done = 0;
    
    while (!root_done)
    {   
        assign_rand_vals(rand_vals, num_nodes);
        M_temp.clear();
        neighbors.clear();

        // Add fitting nodes to M'
        add_and_record(rank, mat, rand_vals, num_nodes, num_procs, neighbors, M_temp, alive, num_alive);

        // M = union(M, M')
        M.insert(M_temp.begin(), M_temp.end());
        
        // Reset random values
        rand_vals.clear();

        // MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0)
        {
            // Delete M' and neighbors from A in root process
            for(std::set<int>::const_iterator i = M_temp.begin(); i != M_temp.end(); i++)
                A.erase(*i);
            for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
                A.erase(*i);
            root_done = A.empty();

            // New translate of alive nodes
            idx = 0;
            num_alive = A.size();

            // free(alive);
            alive = (int*) malloc(sizeof(int) * A.size());
            for(std::set<int>::const_iterator a = A.begin(); a != A.end(); a++){
                alive[idx] = *a;
                idx++;
            }

            // Printing
            // std::cout << "-----active nodes in this round------" << std::endl;
            // for(std::set<int>::const_iterator i = A.begin(); i != A.end(); i++)
            //     std::cout << *i << ' ';
            // std::cout << std::endl;
            // std::cout << "-------------------------------------" << std::endl;
        }
        
        // Broadcast whether A is empty in root process
        MPI_Bcast(&root_done, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

        // Broadcast the new set of active nodes to each process
        MPI_Bcast(&num_alive, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(alive, num_alive, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }


    return M;

}

int main(int argc, char* argv[]){

    MPI_Init(NULL, NULL);

    int rank, size;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    srand((unsigned) time(0));
    
    double *test;
    test = new double [49];

    for (int i = 0; i < 6; i++)
        test[i*6 + i] = 0;
    
    test[1] = 1;
    for (int i = 2; i < 6; i++)
        test[i] = 0;
    test[4] = 1;

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
    test[23] = 1;
    test[24] = 1;

    for (int i = 30; i < 36; i++)
        test[i] = 1;
    test[30] = 0;
    

    std::set<int> I = MIS(test, 6, size);

    if (!I.empty() && rank == 0)
        for(std::set<int>::const_iterator i = I.begin(); i != I.end(); i++)
            std::cout << *i << ' ';

    MPI_Finalize();

}