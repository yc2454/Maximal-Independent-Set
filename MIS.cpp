#include <cstdlib>
#include <set>
#include <vector>
#include <iostream>
#include <tuple>

std::set< std::set<int> > MIS(double **mat, int num_nodes){

    std::set< std::set<int> > mis;
    std::set<int> I;
    //int cur_len = 0;
    int num_remain = num_nodes;
    std::set<int> nodes_remain;
    for (int i = 0; i < num_nodes; i++)
    {
        nodes_remain.insert(i);
    }
    int cur_node = 0;
    std::set<int> neighbors;

    //neighbors.insert(0);
    I.insert(cur_node);

    //int count = 0;

    while (!nodes_remain.empty())
    {
        while (num_remain > 0)
        {
            for (int i = 0; i < num_nodes; i++)
            {
                // removing current node's neighbors
                if ((mat[cur_node][i] > 0) && (nodes_remain.find(i) != nodes_remain.end())){
                    //std::cout << i << std::endl;
                    num_remain--;
                    neighbors.insert(i);
                    nodes_remain.erase(i);
                }
                
            }
            // std::cout << "node " << cur_node << " done" << std::endl;
            // remove current node
            nodes_remain.erase(cur_node);
            num_remain--;
            // choose a new current node, insert
            if(!nodes_remain.empty()){
                cur_node = *(nodes_remain.begin());
                I.insert(cur_node);
            }
        }

        mis.insert(I);
        // Code for testing, commented out:
        // std::cout << "What neighbors look like:";
        // for (std::set<int>::const_iterator j = neighbors.begin(); j != neighbors.end(); ++j)
        //     std::cout << *j << ' ';
        // std::cout << std::endl;
        // std::cout << "What this INDEPEN look like:";
        // for (std::set<int>::const_iterator j = I.begin(); j != I.end(); ++j)
        //     std::cout << *j << ' ';
        // std::cout << std::endl;
        I.clear();
        // std::cout << "One set DONE" << std::endl;
        // std::cout <<std::endl;

        num_remain = neighbors.size();
        nodes_remain = neighbors;
        if(!nodes_remain.empty()){
            cur_node = *(nodes_remain.begin());
            I.insert(cur_node);
        }
        neighbors.clear();

        // count++;
        // if (count > 5)
        //     break;

    }
    
    return mis;
}

// main, for testing purpose
int main(int argc, char* argv[]){
    
    double **test;
    test = new double *[6];
    double *val0 = (double*)malloc(6 * sizeof(double));
    *val0 = 0;
    *(val0 + 1) = 1;
    *(val0 + 2) = 0;
    *(val0 + 3) = 0;
    *(val0 + 4) = 0;
    *(val0 + 5) = 0;

    double *val1 = (double*)malloc(6 * sizeof(double));
    *val1 = 1;
    *(val1 + 1) = 0;
    *(val1 + 2) = 1;
    *(val1 + 3) = 1;
    *(val1 + 4) = 0;
    *(val1 + 5) = 1;

    double *val2 = (double*)malloc(6 * sizeof(double));
    *val2 = 0;
    *(val2 + 1) = 1;
    *(val2 + 2) = 0;
    *(val2 + 3) = 0;
    *(val2 + 4) = 1;
    *(val2 + 5) = 1;

    double *val3 = (double*)malloc(6 * sizeof(double));
    *val3 = 0;
    *(val3 + 1) = 1;
    *(val3 + 2) = 0;
    *(val3 + 3) = 0;
    *(val3 + 4) = 0;
    *(val3 + 5) = 0;
    
    double *val4 = (double*)malloc(6 * sizeof(double));
    *val4 = 0;
    *(val4 + 1) = 0;
    *(val4 + 2) = 1;
    *(val4 + 3) = 0;
    *(val4 + 4) = 0;
    *(val4 + 5) = 1;

    double *val5 = (double*)malloc(6 * sizeof(double));
    *val5 = 0;
    *(val5 + 1) = 1;
    *(val5 + 2) = 1;
    *(val5 + 3) = 0;
    *(val5 + 4) = 1;
    *(val5 + 5) = 0;

    test[0] = val0;
    test[1] = val1;
    test[2] = val2;
    test[3] = val3;
    test[4] = val4;
    test[5] = val5;
    
    std::set< std::set<int> > I = MIS(test, 6);

    // for (int i = 0; i < 5; i++)
    // {
    //     for (int j = 0; j < 5; j++)
    //     {
    //         std::cout << test[i][j] << ' ';
    //     }
    //     std::cout << '\n';
    // }
    
    for (std::set< std::set<int> >::const_iterator i = I.begin(); i != I.end(); ++i){
        for (std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j)
            std::cout << *j << ' ';
        std::cout << "one set done!" << '\n';
    }

    // for (std::set<int>::const_iterator j = I2.begin(); j != I2.end(); ++j)
    //     std::cout << *j << ' ';
        
}