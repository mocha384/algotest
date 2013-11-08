//This program calculates the average shorted path in a graph

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility> //pairing

using namespace std;

//create namespace
namespace algo {

//#define VERBOSE_DEBUG_INFO

int const VERTEX_UNDEFINED = numeric_limits<int>::max();
double const DISTANCE_INFINITY = numeric_limits<double>::infinity();

//checks the infinity
template <typename T>
inline bool is_any_infinity(T value) {
    return !( numeric_limits<T>::min() <= value && 
              value <= numeric_limits<T>::max() );
}
///////////////////////////////////////////////////////
//////////////////// class Graph //////////////////////
//////////////////////////////////////////////////////
class Graph {
public:
//////////////////// Internal class Edge ///////////////////////
    class Edge {
    public:
        int mConnectedVertex; //the edge has connected vertex
        double mDistance; //weight on the edge
        Edge(double dist, int connected_vertex) {
            mDistance = dist;
            mConnectedVertex = connected_vertex;
        }
        inline int get_connected_vertex() {
            return mConnectedVertex;
        }
        inline double get_edge_value() {
            return mDistance;
        }
        bool operator<(const Edge& e) const {
            return mDistance < e.mDistance;
        }
    }; //class Edge

//////////////////// Internal class Vertex /////////////////////
    class Vertex {
    public:
        vector <class Edge> mEdgeList; //list of connected vertices with their distances
        //returns the connected vertices to a given vertex
        inline int get_degrees() {
            return mEdgeList.size(); 
        }
        //get list of neighbours
        vector <class Edge>& get_neighbors() {
            return mEdgeList;
         }
        //find a connection exists
        bool adjacent(int connected_vertex) {
            if (0 < get_degrees()) {
                vector<class Edge>::iterator it; 
                for ( it = mEdgeList.begin(); it != mEdgeList.end(); ++it){
                    if (connected_vertex == it->mConnectedVertex)
                        return true;
                }
            }
            return false;
        }
        //get an edge value 
        double get_edge_value(int connected_vertex) {
             if (0 < get_degrees()) {
                vector<class Edge>::iterator it; 
                for ( it = mEdgeList.begin(); it != mEdgeList.end(); ++it){
                    if (connected_vertex == it->mConnectedVertex)
                        return it->get_edge_value();
                }
            }
            return 0;
        }
        //print connected vertices
        void print() {
             if (0 < get_degrees()) {
                vector<class Edge>::iterator it; 
                for ( it = mEdgeList.begin(); it != mEdgeList.end(); ++it){
                    cout << "(" << it->get_connected_vertex() << "," << it->get_edge_value() << ")" << " ";
                }
            }
        }
        //add an edge
        void add_edge(Edge& edge) {
            mEdgeList.push_back(edge);
        }
       //sort all the neighbours in their ascending order of distance
        void sort_neighbors() { 
            std::sort (mEdgeList.begin(), mEdgeList.end());
        }
    }; //class Vertex

//////////////////////////////////////////////////////
    //list of constructors
    Graph(int vertex) {
        mVertexSize = vertex;
        for(int i=0; i<mVertexSize; i++) {
            Vertex v;
            mVertexList.push_back(v);
        }
    }
    //V(G) : returns numbers of vertices in a graph
    inline int get_total_vertex() {
        return mVertexSize;
    }
    //E(G) : returns number of edges in a graph
    inline int total_edges() {
        return mEdgeSize;
    }
    //Test whether x has an edge to y
    bool adjacent(int x, int y) {
        return mVertexList[x].adjacent(y);
    }
    //return the value associated with edge
    int get_edge_value(int x, int y) {
        return mVertexList[x].get_edge_value(y);
    }
    //sets the value to an edge from vertex x to vertex y
    void add_edge(int x, int y, double dist) {
        Edge edge(dist,y);
        mVertexList[x].add_edge(edge);
    }
    //randomly generates a distance on a given minimum and maximum distance
    double get_distance(double dist_min, double dist_max) {
        double random_num = (double)rand() / (double)RAND_MAX;
        return ((dist_max-dist_min)*random_num + dist_min); //fit the given generated random number into the range
    }
    //get the indexed vertex
    inline Vertex& get_vertex(int index) {
        return mVertexList[index];
    }
    //pick the random vertex
    int get_random_vertex() {
        return rand() % mVertexSize;
    }
    
    //generate random edges for a graph with a given distance ranges
    void generate_edges(double density, int dist_min, int dist_max) {
        int max_possible_edges = mVertexSize * (mVertexSize-1) / 2; //complete graph has n*(n-1)/2 edges
        mDensity = density;
        mDistMin = dist_min;
        mDistMax = dist_max;
        mEdgeSize = max_possible_edges * mDensity;
        for (int i=0; i<mEdgeSize; i++) {
            int x,y;
            while (1) {
                x = get_random_vertex(); // get a vertex
                y = get_random_vertex(); // get another vertex to make a connection
                if (x == y) //both are same, find another pair
                    continue;
                if (adjacent(x,y) == true) //connection is already made, find another pair
                    continue;
                break;
            }
            double dist = get_distance(dist_min,dist_max);
            //cout << "Making pair " << x << "," << y << " : " << dist << endl;
            add_edge(x,y,dist);
            add_edge(y,x,dist); //undirected graph has edge from y to x too
        }   
    }
     void sort() {
        vector <class Vertex>::iterator it;
        for (it = mVertexList.begin(); it != mVertexList.end(); ++it) {
            it->sort_neighbors();
        }
    }   
    void print() {
        cout << "Total Edges : " << mEdgeSize << endl;
        cout << "Vertex size : " << mVertexSize << " Density " << mDensity << endl;
        cout << "Distance range is (" << mDistMin << "," << mDistMax << ")" << endl;
#ifdef VERBOSE_DEBUG_INFO
        vector <class Vertex>::iterator it;
        int i;
        for (it = mVertexList.begin(),i=0; it != mVertexList.end(); ++it,++i) {
            cout << endl;
            cout << "V[" << i << " D(" << it->get_degrees() << ")] :";
            it->print();
        }
        cout << endl;
#endif
    }			
private:
    int mVertexSize; //total available vertices for the graph
    int mEdgeSize;
    int mSource;
    double mDensity;
    double mDistMin, mDistMax;	
    vector <class Vertex> mVertexList;
}; //class Graph
///////////////////////////////////////////////////////
//////////////////// class ShortedPathAlgo ////////////
//////////////////////////////////////////////////////
 
class ShortestPathAlgo {
public:
    //constructor
    ShortestPathAlgo(Graph& G):mGraph(G) {
        int total_vertex = G.get_total_vertex();
        //initialize variables to the default state
        for (int i = 0; i < total_vertex; ++i) {
            mVisited.push_back(false);
            mParentList.push_back(VERTEX_UNDEFINED);
            mDistanceList.push_back(DISTANCE_INFINITY);
        }
    }
    //it prepares the dijkstra structures which can be used later for finding shortest path
    void dijkstra(int source) {
        mSource = source;
        mDistanceList[source] = 0;

        //create a set which holds the pair of (distance, vertex)
        std::set< pair<double,int> > pqueue; 
        pqueue.insert(make_pair(mDistanceList[source],source)); //start with the source vertex first

        while( !pqueue.empty() ) {
            double current_distance = pqueue.begin()->first;
            int current_vertex = pqueue.begin()->second;
            pqueue.erase(pqueue.begin());

            mVisited[current_vertex] = true; //we have visited this vertex now
            
	    //find all the connected vertices to the current_vertex   
            Graph::Vertex v = mGraph.get_vertex(current_vertex); //get the Vertex from the index
            vector <class Graph::Edge> elist = v.get_neighbors(); //get the list of neighbours of this Vertex
            
	    vector <class Graph::Edge>::iterator it;
            for (it = elist.begin();it != elist.end(); ++it) {
                double alt_dist = 0.0;
                int connected_vertex = it->get_connected_vertex(); //connected vertex
                double connected_distance = it->get_edge_value(); //distance from current to connected vertex

                alt_dist = current_distance + connected_distance; //total distance to reach connected_vertex
                if (alt_dist < mDistanceList[connected_vertex]) { //do we have better distance?
                    pqueue.erase(make_pair(mDistanceList[connected_vertex],connected_vertex)); //remove the old pair

                    mDistanceList[connected_vertex] = alt_dist;
                    mParentList[connected_vertex] = current_vertex;

                    pqueue.insert(make_pair(mDistanceList[connected_vertex],connected_vertex)); //add new pair with new distance
                }
            }
        }
    }
    //gets the list of vertex in the path form for a shortest distance
    vector <int> get_shortest_path(int target) {
        int u = target;
        vector <int> path;
        while (mParentList[u] != VERTEX_UNDEFINED) {
            path.push_back(u);
            u = mParentList[u];
        }
        return path;
    }
    //gets the shorted distance from source to destination vertex
    //0 in case of no direct path
    inline double get_shortest_distance(int target) { 
       return mDistanceList[target];
    }
    void print() {
        cout << "Source node : " << mSource << endl;
#ifdef VERBOSE_DEBUG_INFO
        cout << "Parents :" << endl;
        vector<int>::iterator it;
        for (it = mParentList.begin(); it != mParentList.end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
        cout << "Distances :" << endl;
        vector<double>::iterator itd;
        for (itd = mDistanceList.begin(); itd != mDistanceList.end(); ++itd) {
            cout << *itd << " ";
        }
        cout << endl;
        cout << "Visted :" << endl;
        vector<bool>::iterator itb;
        for (itb = mVisited.begin(); itb != mVisited.end(); ++itb) {
            cout << *itb << " ";
        }
        cout << endl;
#endif
    }
private:
    Graph& mGraph;
    int mSource;
    vector <double> mDistanceList;
    vector <bool> mVisited;
    vector <int> mParentList;
}; // class ShortestPathAlgo

class MonteCarlo {
public:
	MonteCarlo(int max_vertex, double density, double dist_min, double dist_max):
           mMaxVertex(max_vertex),mDensity(density),mDistMin(dist_min),mDistMax(dist_max) {
	}

	//run simulator
	void run_simulation() {
	    Graph graph(mMaxVertex);
	    graph.generate_edges(mDensity,mDistMin,mDistMax);
	    graph.print();
	    graph.sort();

	    ShortestPathAlgo spa(graph);
	    int source_vertex = 0;
	    Graph::Vertex v = graph.get_vertex(source_vertex);
	    while (!v.get_degrees()) {
	        ++source_vertex;
	        v = graph.get_vertex(source_vertex);
	    }
	    spa.dijkstra(source_vertex);
	    spa.print();
	    int avg_counter = 0;
	    double total_dist = 0.0;
	    for (int i = source_vertex+1; i < mMaxVertex; i++) {
	        double dist = spa.get_shortest_distance(i); 
	        if (!is_any_infinity(dist)) { 
	            ++avg_counter;
	            total_dist += dist;
	        }
	    }
	    cout << "Average shortest distance of a graph is " << total_dist / avg_counter << endl;
	}
private:
	int mMaxVertex;
    double mDensity;
	double mDistMin;
	double mDistMax;
};

}//namespace algo

using namespace algo;

int main () {
    //initialize the seed for random number generation
    srand((unsigned)time(0));
    const int max_vertex = 50;
    const double dist_min = 1.0;
    const double dist_max = 10.0;
    double density = 0.2;

    //density with 0.2
    MonteCarlo mc_1(max_vertex,density,dist_min,dist_max);
    mc_1.run_simulation();    

    density = 0.4;
    //desnsity with 04.
    MonteCarlo mc_2(max_vertex,density,dist_min,dist_max);
    mc_2.run_simulation();

    return 0;
}
