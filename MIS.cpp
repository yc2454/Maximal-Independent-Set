#include <cstdlib>
#include <set>
#include <iostream>

std::set< std::set<int> > MIS(double **mat, int num_nodes){

    std::set< std::set<int> > mis;
    std::set<int> I;
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
        I.clear();
        // Code for testing, commented out:
        // std::cout << "What neighbors look like:";
        // for (std::set<int>::const_iterator j = neighbors.begin(); j != neighbors.end(); ++j)
        //     std::cout << *j << ' ';
        // std::cout << std::endl;
        // std::cout << "What this INDEPEN look like:";
        // for (std::set<int>::const_iterator j = I.begin(); j != I.end(); ++j)
        //     std::cout << *j << ' ';
        // std::cout << std::endl;
        // std::cout << "One set DONE" << std::endl;
        // std::cout <<std::endl;
        num_remain = neighbors.size();
        nodes_remain = neighbors;
        if(!nodes_remain.empty()){
            cur_node = *(nodes_remain.begin());
            I.insert(cur_node);
        }
        neighbors.clear();

    }
    
    return mis;
}
