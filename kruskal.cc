//This program calculates the average shorted path in a graph

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility> //pairing

//readind from file
#include <fstream>
#include <sstream>
#include <string>

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
    	init(vertex);
    }
    //read from a file
    Graph(char* filename) {
    	ifstream file(filename);
    	string str;
    	if (getline(file,str)) {
        	istringstream iss(str);
        	int v = 0;
        	iss >> v;
        	init(v);
        	cout << "Found vertex : " << v << endl;
    	}
    	while (getline(file,str)) {
        	istringstream iss(str);

        	int x=0, y=0, dist=0;
        	iss >> x;
        	iss >> y;
        	iss >> dist;
        	cout << "Found x : " << x <<  " y : "  <<  y  << " dist : " << dist << endl;
        	add_edge(x,y,dist);
    	}
    }
    inline void init(int vertex) {
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
//////////////////// class UnionFind ////////////
//////////////////////////////////////////////////////


//Implement union find data structure which quickly tells whether
//two vertices are connected or not
class UnionFind {
public:
	vector<int> mIdList;
	vector<bool> mExistList;
	int mCount; //component counter
	UnionFind(int n) {
		mCount = n;
		for(int i = 0; i < n; i++) {
			mIdList.push_back(i);
			mExistList.push_back(false);
		}
	}
	//checks whether x and y are connected
	bool connected( int x, int y) {
		return (mIdList[x] == mIdList[y]);
	}
	bool union_add(int x, int y) {
		if (mIdList[x] == mIdList[y])
			return false;
		//check whether x and y are of the same component
		if (!mExistList[x] && !mExistList[y]) {
			mIdList[x] = mCount;
			mIdList[y] = mCount;
			mExistList[x] = true;
			mExistList[y] = true;
			++mCount;
		} else if (mExistList[x] && mExistList[y]) { //different component but not connected yet
			int replace_value, replace_with;
			//consider the lower value
			if (mIdList[x] < mIdList[y]) {
				replace_with = mIdList[x];
				replace_value = mIdList[y];
			} else {
				replace_with = mIdList[y];
				replace_value = mIdList[x];
			}
			//find all the component values with previously found one
			for (unsigned int i = 0; i < mIdList.size(); i++) {
				if (mIdList[i] == replace_value)
					mIdList[i] = replace_with;
			}
		} else if (mExistList[x]){
			mExistList[y] = true;
			mIdList[y] = mIdList[x]; //assign component of x to y
		} else {
			mExistList[x] = true;
			mIdList[x] = mIdList[y]; //assign component of y to x
		}
		return true;
	}
	void print() {
		for(unsigned int i = 0; i < mIdList.size(); i++) {
			cout << mIdList[i] << " ";
		}
		cout << endl;
	}
};
///////////////////////////////////////////////////////
//////////////////// class KruskalMST ////////////
//////////////////////////////////////////////////////

class KruskalMST {
public:
	//Internal class to hold the full edge from x to y with its distance
	class FullEdge {
	public:
		int mX, mY;
		double mDistance;
		FullEdge() {
			mX = 0;
			mY = 0;
			mDistance = 0.0;
		}
		FullEdge(int x,int y,double dist) {
			mX = x;
			mY = y;
			mDistance = dist;
		}
		//assignmnet overloading
		FullEdge& operator=(const FullEdge &rhs) {
			if (this != &rhs) {
				mX = rhs.mX;
				mY = rhs.mY;
				mDistance = rhs.mDistance;
			}
			return *this;
		}
		//< operator overloading for sorting
        bool operator<(const FullEdge& fe) const {
            return mDistance < fe.mDistance;
        }
	};
	KruskalMST(Graph& graph) {
		mVertexSize = graph.get_total_vertex() ;
		mTotalCost = 0;

		//prepare the x,y,dist tuple for processing MST
		for (unsigned int i = 0; i < mVertexSize; i++) {
			Graph::Vertex v = graph.get_vertex(i); //get the Vertex from the index
			vector <class Graph::Edge> elist = v.get_neighbors(); //get the list of neighbours of this Vertex

			vector <class Graph::Edge>::iterator it;
			for (it = elist.begin();it != elist.end(); ++it) {
				int connected_vertex = it->get_connected_vertex(); //connected vertex
				double connected_distance = it->get_edge_value(); //distance from current to connected vertex
				FullEdge fe(i,connected_vertex,connected_distance);
				mFullEdgeList.push_back(fe);
			}
		}
		//sort tuples according to the dist
		std::sort(mFullEdgeList.begin(),mFullEdgeList.end());
	}
	//implement the algorithm
	void runMST() {
		UnionFind uf(mVertexSize); //union find data structure

		mVistedMST.clear();
		for (unsigned int e = 0; e < mFullEdgeList.size() &&
					    mVistedMST.size() < mVertexSize - 1; ++e) {
			FullEdge fe = mFullEdgeList[e];
			cout << setw(2) << fe.mX << " to " << fe.mY << " dist = " << fe.mDistance << " ### ";
			if (!uf.connected(fe.mX, fe.mY)) {
				if (uf.union_add(fe.mX, fe.mY)) {
					mVistedMST.push_back(fe);
#ifdef VERBOSE_DEBUG_INFO
					cout << "Adding (" << fe.mX << "," << fe.mY << ") ";
#endif
				}
			}
			uf.print();
		}
	}
	void print() {
		mTotalCost = 0;
		for (unsigned int e = 0; e < mVistedMST.size(); ++e) {
			FullEdge fe = mVistedMST[e];
			mTotalCost += fe.mDistance;
			cout << fe.mX << " to " << fe.mY << " dist = " << fe.mDistance << endl;
		}
		cout << "Total MST cost is " << mTotalCost << endl;
	}
private:
	unsigned int mVertexSize; //total vertices
	double mTotalCost; //cost of MST
	vector<class FullEdge> mFullEdgeList; //maintain the list of edge tuples
	vector<class FullEdge> mVistedMST; //visited nodes of MST
};

}//namespace algo

using namespace algo;

int main (int argc, char* argv[]) {

    if (argc != 2) {
    	cout << "Error ## Usage : program input_file" << endl;
    	return 1;
    }

    Graph graph(argv[1]);
#ifdef VERBOSE_DEBUG_INFO
    graph.print();
#endif

    KruskalMST km(graph);
    km.runMST();
    km.print();

    return 0;
}
