#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <climits>
#include <cfloat>
#include <stack>
#include <set>
#include <algorithm>
#include <stdlib.h>


Graph::Graph(const char* const & edgelist_csv_fn) {
    // TODO
    //char buffer[256];
    string buffer;
    ifstream input_file(edgelist_csv_fn, ios::in);

    if(!input_file.is_open()){
        exit(1);
    }

    while(!input_file.eof()){
        // store a line to buffer
        //input_file.getline(buffer, 100);
        getline(input_file, buffer);
        if(buffer.size()==0){
            break;
        }
        buffer = buffer + '\n';

        // create istringstream of current line to extract substrings
        istringstream ss(buffer);
        string first, second, third;
        getline(ss, first, ',');
        getline(ss, second, ',');
        getline(ss, third, '\n');
        double w = stod(third);

        // store elements of an edge to vector named edge
        //vector<string> edge;
        //edge.push_back(first);
        //edge.push_back(second);
        //edge.push_back(third);

        // put complete edge {u:{v: weight}} and {v:{u: weight}} into map
        // tuple: (cummulative weight, weight, curr node, prev node)
        if(graph_map.find(first) == graph_map.end()){
            unordered_map<string, double> inner_map;
            inner_map[second] = w; //make_tuple(DBL_MAX, w, second, "empty");
            graph_map[first] = make_tuple(DBL_MAX, "empty", false, inner_map);
            edge_count++;
        }else{
            get<3>(graph_map.at(first))[second] = w;
            edge_count++;
        }

        if(graph_map.find(second) == graph_map.end()){
            unordered_map<string, double> inner_map;
            inner_map[first] = w; //make_tuple(DBL_MAX, w, first,"empty");
            graph_map[second] = make_tuple(DBL_MAX, "empty", false, inner_map);
        }else{
            get<3>(graph_map.at(second))[first] = w;
        }

    }
    input_file.close();
}

unsigned int Graph::num_nodes() {
    // TODO
    return graph_map.size();
}

vector<string> Graph::nodes() {
    // TODO
    vector<string> nodes;
    for(const auto &pair : graph_map){
        nodes.push_back(pair.first);
    }
    return nodes;
}

unsigned int Graph::num_edges() {
    // TODO
    return edge_count;
}

unsigned int Graph::num_neighbors(string const & node_label) {
    // TODO
    unsigned int num = get<3>(graph_map.at(node_label)).size();
    return num;
}

double Graph::edge_weight(string const & u_label, string const & v_label) {
    // TODO
    if(graph_map.find(u_label) == graph_map.end()){
        return -1;
    }else if(get<3>(graph_map.at(u_label)).find(v_label) == get<3>(graph_map.at(u_label)).end()){
        return -1;
    }else if(u_label == v_label){
        return -1;
    }
    double weight = get<3>(graph_map.at(u_label)).at(v_label);
    return weight;
}

vector<string> Graph::neighbors(string const & node_label) {
    // TODO
    vector<string> neighbors;
    for(const auto &pair : get<3>(graph_map.at(node_label))){
        neighbors.push_back(pair.first);
    }
    return neighbors;
}

vector<string> Graph::shortest_path_unweighted(string const & start_label, string const & end_label) {
    // TODO
    vector<string> mynodes = nodes();
    vector<string> path;
    // this part should check if start_label and end_label are in mynodes
    if(find(mynodes.begin(), mynodes.end(), start_label) == mynodes.end()){
        return path;
    }else if(find(mynodes.begin(), mynodes.end(), end_label) == mynodes.end()){
        return path;
    }

    // return starting node if start and end are equal
    if(start_label == end_label){
        path.push_back(start_label);
        return path;
    }

    // Dijkstra
    priority_queue<tuple<double, string, string>, 
                        vector<tuple<double, string, string>>, 
                        greater<tuple<double, string, string>>> myqueue;
    //queue<tuple<string, string>> myqueue;
    tuple<double, string, string> start_tuple = make_tuple(0, start_label, "empty");
    unordered_map<string, string> prev_map;    // this is a map that store previous node (curr:prev)
    string curr;
    string prev;
    myqueue.push(start_tuple);

    // set initial node's distance to 0
    get<0>(graph_map.at(start_label)) = 0;

    while(!myqueue.empty()){
        // update curr node and prev node
        tuple<double, string, string> curr_tuple = myqueue.top();
        myqueue.pop();
        curr = get<1>(curr_tuple);
        prev = get<2>(curr_tuple);
        bool visited = get<2>(graph_map.at(curr));
        double curr_distance = get<0>(graph_map.at(curr));
        unordered_map<string, double> edges_map = get<3>(graph_map.at(curr));
        if(!visited){
            // if not visited, mark it as visited.
            get<2>(graph_map.at(curr)) = true;
            // loop through all edges
            for(const auto &pair : edges_map){
                if(pair.first != prev && curr_distance + 1 < get<0>(graph_map.at(pair.first))){ // && get<2>(graph_map.at(pair.first)) == false
                    get<0>(graph_map.at(pair.first)) = curr_distance + 1;
                    myqueue.push(make_tuple(curr_distance + 1, pair.first, curr));
                    prev_map[pair.first] = curr;
                }
            }
        }
    }

    // restore path
    string temp = end_label;
    path.push_back(temp);
    // if cant find temp in prev_map, then this is temp node is not connected with starting node
    while(prev_map.find(temp) != prev_map.end() && temp != start_label){
        path.push_back(prev_map[temp]);
        temp = prev_map[temp];
    }
    if(path.size() == 1){
        path.pop_back();
        for(unsigned int i=0; i<num_nodes(); i++){
            get<0>(graph_map.at(mynodes[i])) = DBL_MAX;
            get<1>(graph_map.at(mynodes[i])) = "empty";
            get<2>(graph_map.at(mynodes[i])) = false;
        }
        return path;
    }
    reverse(path.begin(), path.end());



    // reset graph's tuples
    for(unsigned int i=0; i<num_nodes(); i++){
        get<0>(graph_map.at(mynodes[i])) = DBL_MAX;
        get<1>(graph_map.at(mynodes[i])) = "empty";
        get<2>(graph_map.at(mynodes[i])) = false;
    }

    return path;
}




vector<tuple<string,string,double>> Graph::shortest_path_weighted(string const & start_label, string const & end_label) {
    // TODO
    vector<string> mynodes = nodes();
    vector<string> path;
    // create a vector of tuples (`from_label`, `to_label`, `edge_weight`)
    vector<tuple<string, string, double>> toreturn;
    // this part should check if start_label and end_label are in mynodes
    if(find(mynodes.begin(), mynodes.end(), start_label) == mynodes.end()){
        return toreturn;
    }else if(find(mynodes.begin(), mynodes.end(), end_label) == mynodes.end()){
        return toreturn;
    }

    // return starting node if start and end are equal
    if(start_label == end_label){
        toreturn.push_back(make_tuple(start_label, start_label, -1));
        return toreturn;
    }

    // BFS
    priority_queue<tuple<double, string, string>, 
                        vector<tuple<double, string, string>>, 
                        greater<tuple<double, string, string>>> myqueue;
    //queue<tuple<string, string>> myqueue;
    tuple<double, string, string> start_tuple = make_tuple(0, start_label, "empty");
    unordered_map<string, string> prev_map;    // this is a map that store previous node (curr:prev)
    string curr;
    string prev;
    myqueue.push(start_tuple);

    // set initial node's distance to 0
    get<0>(graph_map.at(start_label)) = 0;

    while(!myqueue.empty()){
        // update curr node and prev node
        tuple<double, string, string> curr_tuple = myqueue.top();
        myqueue.pop();
        curr = get<1>(curr_tuple);
        prev = get<2>(curr_tuple);
        bool visited = get<2>(graph_map.at(curr));
        double curr_distance = get<0>(graph_map.at(curr));
        unordered_map<string, double> edges_map = get<3>(graph_map.at(curr));
        if(!visited){
            // if not visited, mark it as visited.
            get<2>(graph_map.at(curr)) = true;
            // loop through all edges
            for(const auto &pair : edges_map){
                // store edge weight
                double weight = pair.second;
                if(pair.first != prev && curr_distance + weight < get<0>(graph_map.at(pair.first))){ // && get<2>(graph_map.at(pair.first)) == false
                    get<0>(graph_map.at(pair.first)) = curr_distance + weight;
                    myqueue.push(make_tuple(curr_distance + weight, pair.first, curr));
                    prev_map[pair.first] = curr;
                }
            }
        }
    }

    // restore path
    string temp = end_label;
    path.push_back(temp);
    // if cant find temp in prev_map, then this is temp node is not connected with starting node
    while(prev_map.find(temp) != prev_map.end() && temp != start_label){
        path.push_back(prev_map[temp]);
        temp = prev_map[temp];
    }
    // if such path dose not exist, return empty vector;
    if(path.size() == 1){
        path.pop_back();
        for(unsigned int i=0; i<num_nodes(); i++){
            get<0>(graph_map.at(mynodes[i])) = DBL_MAX;
            get<1>(graph_map.at(mynodes[i])) = "empty";
            get<2>(graph_map.at(mynodes[i])) = false;
        }
        return toreturn;
    }
    reverse(path.begin(), path.end());

    // construct a vector of tuples to return from path
    for(unsigned int i=0; i<path.size()-1; i++){
        string curr = path[i];
        string next = path[i+1];
        unordered_map<string, double> curr_map = get<3>(graph_map.at(curr));
        double edge_weight = curr_map[next];
        toreturn.push_back(make_tuple(path[i], path[i+1], edge_weight)); //get<3>(graph_map.at(path[i]))[path[i+1]]
    }
    

    for(unsigned int i=0; i<num_nodes(); i++){
        get<0>(graph_map.at(mynodes[i])) = DBL_MAX;
        get<1>(graph_map.at(mynodes[i])) = "empty";
        get<2>(graph_map.at(mynodes[i])) = false;
    }

    return toreturn;
}

vector<vector<string>> Graph::connected_components(double const & threshold) {
    // TODO
    vector<vector<string>> toreturn;

    // this for-loop remove edges greater than threshold in graph
    for(auto &node : graph_map){
        unordered_map<string, double> temp_map = get<3>(node.second);
        for(auto &pair : temp_map){
            if(get<3>(node.second).find(pair.first) != get<3>(node.second).end()){
                if(pair.second > threshold){
                // node is the pair, node.second is the inner tuple
                // get<3>(node.second) is the inner map
                get<3>(node.second).erase(pair.first);
                } 
            }
        }
    }

    
    // BFS to traverse graph
    vector<string> mynodes = nodes();
    int num = num_nodes();
    vector<string> path;
    unordered_map<string, string> prev_map;    // this is a map that store previous node (curr:prev)
    
    queue<tuple<string, string>> myqueue;
    string curr;
    string prev;
    int count = 0;

    for(const auto &node : graph_map){
        tuple<string, string> start_tuple = make_tuple(node.first, "empty");
        set<string> inner;
        //inner.insert(node.first);
        myqueue.push(start_tuple);

        // set initial node's distance to 0
        get<0>(graph_map.at(node.first)) = 0;

        // traverse graph, bfs
        while(!myqueue.empty()){
            tuple<string, string> curr_tuple = myqueue.front();
            myqueue.pop();
            curr = get<0>(curr_tuple);
            prev = get<1>(curr_tuple);
            bool visited = get<2>(graph_map.at(curr));
            //double curr_distance = get<0>(graph_map.at(curr));
            unordered_map<string, double> edges_map = get<3>(graph_map.at(curr));
            if(!visited){
                get<2>(graph_map.at(curr)) = true;
                // where should this line be
                inner.insert(curr);
                for(const auto &pair : edges_map){
                    if(pair.first != prev && get<2>(graph_map.at(pair.first)) == false){
                        myqueue.push(make_tuple(pair.first, curr));
                        // that line was here
                    }
                }
            }
        }

        vector<string> myvect;
        for(string s : inner){
            myvect.push_back(s);
            count++;
        }
        if(myvect.size() != 0){
            toreturn.push_back(myvect);
        }
        if(count == num){
            break;
        }
    }

    for(unsigned int i=0; i<num_nodes(); i++){
        get<0>(graph_map.at(mynodes[i])) = DBL_MAX;
        get<1>(graph_map.at(mynodes[i])) = "empty";
        get<2>(graph_map.at(mynodes[i])) = false;
    }

    return toreturn;
}

double Graph::smallest_connecting_threshold(string const & start_label, string const & end_label) {
    // TODO
    double val = 0;
    return val;
}
