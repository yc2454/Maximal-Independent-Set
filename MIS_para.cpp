#include <cstdlib>
#include "mpi.h"
#include <set>
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <unistd.h>
#include <pthread.h>
#include <utility>

void print_arr(int *array, int len) {
    for (int i = 0; i < len; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

void print_arr_d(double *array, int len) {
    for (int i = 0; i < len; i++)
    {
        printf("%f ", array[i]);
    }
    printf("\n");
}

void assign_rv(int** rv, int num_nodes){

    *rv = (int*) malloc(sizeof(int) * num_nodes);
    for (int i = 0; i < num_nodes; i++)
        (*rv)[i] = rand();
    
}

// Translate the elements in a set into an array
void translate(std::set<int> A, int** arr) {

    *arr = (int*) malloc(sizeof(int) * A.size());
    int idx = 0;

    for(std::set<int>::const_iterator a = A.begin(); a != A.end(); a++){
        (*arr)[idx] = *a;
        idx++;
    }

}

void divide_and_send(int num, int num_procs, int** counts, int ** displs) {

    int per = (num >= num_procs) ? round((1.0 * num) / num_procs) : 1;
    int left = num - per * (num_procs - 1);

    (*counts) = (int*)malloc(sizeof(int) * num);
    (*displs) = (int*)malloc(sizeof(int) * num);

    for (int i = 0; i < num_procs; i++)
    {
        if (num >= num_procs) {
            if (i < num_procs - 1)
                (*counts)[i] = per;
            else
                (*counts)[i] = left;
            (*displs)[i] = i * per;
        }
        else {
            (*counts)[i] = i < num ? 1 : 0;
            if (i == 0)
                (*displs)[i] = 0;
            else
                (*displs)[i] = i < num ? ((*displs)[i-1] + 1) : (*displs)[i-1];
        }
        
    }

}

void add_and_record(int rank, double* mat, int num_nodes, int num_procs, 
                    std::set<int> &neighbors, std::set<int> &M, 
                    int* alive, int num_alive, int nodes_per_proc)
{
    // This variable will be used in Luby's Algo
    bool flag = 1;

    // For MPI_Scatterv
    int* send_counts;
    int* displs;
    divide_and_send(num_alive, num_procs, &send_counts, &displs);

    // For MPI_Allgatherv
    int* r_c_rv;
    int* displs_rv;
    divide_and_send(num_nodes, num_procs, &r_c_rv, &displs_rv);

    // The node in consideration
    int cur_node;
    // The rank of the cur_node
    int cur_node_rank;
    // The index of cur_node in mat
    int cur_node_idx;
    
    int recving = (send_counts[rank] < 0 || rank >= num_procs) ? 0 : send_counts[rank];

    // Buffer to hold the nodes scattered to each process
    // int *sub_alive = (int*) malloc(sizeof(int) * recving);
    // std::cout << "rank" << rank << ' ' << recving << " receiving " << std::endl;
    // MPI_Scatterv(alive, send_counts, displs, 
    //     MPI_INT, sub_alive, recving, MPI_INT, 0, MPI_COMM_WORLD);

    int* rand_vals = (int*) malloc(sizeof(int) * num_nodes);
    int* sub_rv = (int*) malloc(sizeof(int) * recving);
    
    assign_rv(&sub_rv, r_c_rv[rank]);

    MPI_Allgatherv(sub_rv, r_c_rv[rank], MPI_INT, rand_vals, 
        r_c_rv, displs_rv, MPI_INT, MPI_COMM_WORLD);

    // Buffer for nodes to add
    int* sub_M = (int*) malloc(sizeof(int) * recving);
    int* sub_nbrs = (int*) malloc(sizeof(int) * num_nodes);
    int count_M = 0;
    int count_ngb = 0;

    // Luby's Algo
    // if (recving > 0)
    for (int i = 0; i < num_alive; i++)
    {
        cur_node = alive[i];
        cur_node_rank = cur_node / nodes_per_proc;
        cur_node_idx = cur_node - (cur_node_rank * nodes_per_proc);

        // When the mat in the current rank doesn't have cur_node,
        // we can do nothing
        if (cur_node_rank != rank)
            continue;
        else {
            //check whether node i has bigger random value than all its neighbors
            for (int j = 0; j < num_nodes; j++)
            {
                if (rand_vals[cur_node] < rand_vals[j] && mat[cur_node_idx * num_nodes + j] > 0)
                    flag = 0;
            }
            //if it does, add i to M'
            if (flag){
                for (int j = 0; j < num_nodes; j++)
                {
                    // Delete the neighbors of the current node
                    if(mat[cur_node_idx * num_nodes + j] > 0 
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
    }
    
    //Printing status
    // std::cout << "rank " << rank << std::endl;
    // print_arr(rand_vals, num_nodes);
    // for(std::set<int>::const_iterator i = M.begin(); i != M.end(); i++)
    //     std::cout << *i << ' ';
    // std::cout << "*" << std::endl;
    // for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
    //     std::cout << *i << ' ';
    // std::cout << "*" << std::endl;
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

    // Add gathered data into M and its neighbors
    if(rank == 0) {
        for (int i = 0; i < all_count_M; i++)
            M.insert(all_M[i]);
        for (int i = 0; i < all_count_ngb; i++)
            neighbors.insert(all_nbrs[i]);
    }

    free(send_counts);
    free(displs);
    if (rank == 0){
        free(all_nbrs);
        free(all_M);
    }
    free(sub_nbrs);
    free(sub_M);
    free(sub_rv);
    free(rand_vals);

    MPI_Barrier(MPI_COMM_WORLD);

}

struct MIS_para
{
    double* mat;
    int num_nodes;
    int num_procs;
    int nodes_per_proc;
};

void *MIS(void *args) {
    
    struct MIS_para *pstru;
    pstru = (struct MIS_para *) args;
    int num_nodes = pstru->num_nodes;
    int num_procs = pstru->num_procs;
    double* mat = pstru->mat;
    int nodes_per_proc = pstru->nodes_per_proc;

    std::set<int> nodes_remain, A, M, M_temp, neighbors;

    int rank;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < num_nodes; i++)
        nodes_remain.insert(i);
    
    // Translate the sets into arrays
    int* rem;
    int num_rem = nodes_remain.size();
    translate(nodes_remain, &rem);

    // Keeping check of alive nodes at each round
    int* alive;
    int num_alive;

    // This variable indicates whether we have found ONE MIS
    bool root_done;

    // This variable indicates whether we have found ALL MIS
    bool all_done = 0;

    // message buffer to send back to parent thread
    int* message;
    int msg_size;
    int msg_count = 0;
    // bool msg_sent = true;
    
    while (!all_done) { 

        num_alive = num_rem;
        alive = (int*) malloc(sizeof(int) * num_rem);
        alive = rem;
        A.clear();
        for (int i = 0; i < num_alive; i++)
            A.insert(alive[i]);
        
        root_done = 0;

        while (!root_done)
        {   
            M_temp.clear();
            neighbors.clear();

            // Add fitting nodes to M'
            add_and_record(rank, mat, num_nodes, num_procs,
                neighbors, M_temp, alive, num_alive, nodes_per_proc);

            MPI_Barrier(MPI_COMM_WORLD);

            // M = union(M, M')
            M.insert(M_temp.begin(), M_temp.end());

            if (rank == 0)
            {
                // Delete M' and neighbors from A in root process
                for(std::set<int>::const_iterator i = M_temp.begin(); i != M_temp.end(); i++)
                    A.erase(*i);
                for(std::set<int>::const_iterator i = neighbors.begin(); i != neighbors.end(); i++)
                    A.erase(*i);
                root_done = A.empty();

                // New translate of alive nodes
                num_alive = A.size();
                translate(A, &alive);
                // printf("The alive array at this round: \n");
                // print_arr(alive, num_alive);

            }
            
            // Broadcast whether A is empty in root process
            MPI_Bcast(&root_done, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

            // Broadcast the new set of active nodes to each process
            MPI_Bcast(&num_alive, 1, MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Bcast(alive, num_alive, MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);
        }

        if (rank == 0) {

            for (std::set<int>::const_iterator i = M.begin(); i != M.end(); i++)
                nodes_remain.erase(*i);

            std::cout << "-----generated set in this round-----" << std::endl;
            for(std::set<int>::const_iterator i = M.begin(); i != M.end(); i++)
                std::cout << *i << ' ';
            std::cout << std::endl;
            std::cout << "-------------------------------------" << std::endl;
            msg_count++;

            translate(M, &message);
            msg_size = M.size();
            for (int i = 0; i < num_procs; i++) {
                // MPI_Send(&msg_sent, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
                MPI_Send(message, msg_size, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&msg_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            
            all_done = nodes_remain.empty();
            num_rem = nodes_remain.size();
            translate(nodes_remain, &rem);
            
        }
        
        if (rank == 0) {
        // msg_sent = false;
        for (int i = 0; i < num_procs; i++) 
            MPI_Send(&msg_count, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        }

        M.clear();

        // Broadcast whether there is any nodes left in root process
        MPI_Bcast(&all_done, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        // Broadcast the new set of remaining nodes to each process
        MPI_Bcast(&num_rem, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(rem, num_rem, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(alive);
    // return M;
    return((void *)0);

}

int main(int argc, char* argv[]){

    // MPI_Init(NULL, NULL);
    int provided;

    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Warning MPI did not provide MPI_THREAD_MULTIPLE\n");
    }

    int rank, size;

    // Get rank and size of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank);
    
    double *test;
    test = new double [36];

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
    
    double* test1 = new double[12]; 
    double* test2 = new double[12];
    double* test3 = new double[12];
    for (int i = 0; i < 36; i++)
    {
        if (i < 12) 
            test1[i] = test[i];
        else if (i < 24)
            test2[i - 12] = test[i];
        else
            test3[i - 24] = test[i];
    }

    std::vector<double*> tests;
    tests.push_back(test1);
    tests.push_back(test2);
    tests.push_back(test3);
    
    int npp = round((1.0 * 6) / size);

    std::vector<MIS_para> args;
    MIS_para arg;
    for (int i = 0; i < 3; i++)
    {
        arg.mat = tests[i];
        arg.num_nodes = 6;
        arg.num_procs = size;
        arg.nodes_per_proc = npp;
        args.push_back(arg);
    }
    
    int err;
    pthread_t ntid;
    err = pthread_create(&ntid, NULL, &MIS, &(args[rank]));
    if (err != 0)
        printf("can't create thread: %s\n", strerror(err));

    // pthread_join(ntid, NULL);
    int num_nodes = 6;
    int* msg = (int*) malloc(sizeof(int) * num_nodes);
    int msg_size;
    // bool msg_sent = true;
    int msg_count;

    MPI_Recv(&msg_count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("We will send message %d times\n", msg_count);

    for (int i = 0; i < msg_count; i++){
        MPI_Recv(msg, num_nodes, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&msg_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank == 0)
            printf("MIS received in main thread:\n");
        print_arr(msg, msg_size);
        // MPI_Recv(&msg_sent, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    pthread_join(ntid, NULL);

    MPI_Finalize();

}