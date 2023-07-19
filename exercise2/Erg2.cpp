#include <iostream>
#include <time.h>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/shortest_path.h>

using namespace std;
using namespace boost;
using namespace leda;

// This function is used to generate a random directed graph 
typedef property<edge_weight_t, int> EdgeWeightProperty;

// Define the directed graph as an adjacency list
typedef adjacency_list<vecS, vecS, directedS, no_property, EdgeWeightProperty, no_property, listS> DirectedGraph;

// Define the vertex class as vertex_descriptor
typedef graph_traits<DirectedGraph>::vertex_descriptor Vertex;

// Define the edge class as edge_descriptor
typedef graph_traits<DirectedGraph>::edge_descriptor Edge;

// Define the edge iterator as edge_iterator
typedef graph_traits<DirectedGraph>::edge_iterator EdgeIterator;

// Define the edge iterator as edge_iterator
typedef graph_traits<DirectedGraph>::vertex_iterator VertexIterator;

// Define the edge weight map as property_map
typedef property_map<DirectedGraph, edge_weight_t>::type EdgeWeightMap;


// This function is used to generate a random directed graph
void CopyLedaGraphToBoostGraph(DirectedGraph& BoostDirectedGraph, leda::graph& LedaGraph, edge_array<int>& LedaEdgeWeightMap)
{
	// Create a boost directed graph with the same number of nodes as the leda directed graph
	DirectedGraph boostGraph(LedaGraph.number_of_nodes());

	// Create a leda edge
	leda::edge tempEdge;

	
	forall_edges(tempEdge, LedaGraph)
	{
		// Get the source node of the edge
		node source = LedaGraph.source(tempEdge);

		// Get the target node of the edge
		node target = LedaGraph.target(tempEdge);

		// Get the weight of the edge
		int edgeWeight = LedaEdgeWeightMap[tempEdge];

		// Add the edge in the boost directed graph
		add_edge(LedaGraph.index(source), LedaGraph.index(target), edgeWeight, boostGraph);
	}

	// Set the boost directed graph as the input parameter
	BoostDirectedGraph = boostGraph;
}


 // Applies the Bellman Ford algorithm to the  directedGraph

bool BellmanFord(DirectedGraph& directedGraph, Vertex startingVertex)
{
	// Create a property map that contains the nodes's costs
	std::map<Vertex, int> nodeCostMap;

	// Create a property map that contains the edges's weights
	EdgeWeightMap boostEdgeWeightMap = get(edge_weight, directedGraph);

	// 	Initialise the boost vertex iterators
	VertexIterator vertexIteratorBegin, vertexIteratorEnd;

	
	for(tie(vertexIteratorBegin, vertexIteratorEnd) = vertices(directedGraph); vertexIteratorBegin != vertexIteratorEnd; vertexIteratorBegin++)
	{
		//
		nodeCostMap.insert(pair<Vertex,int>(*vertexIteratorBegin, INT_MAX));
	}

	// Set the starting vertex cost to 0
	nodeCostMap[startingVertex] = 0;

	// Initialise the boost edge iterators
	EdgeIterator edgeIteratorBegin, edgeIteratorEnd;

	
	for(tie(vertexIteratorBegin, vertexIteratorEnd) = vertices(directedGraph); vertexIteratorBegin != vertexIteratorEnd; vertexIteratorBegin++)
	{
		
		for(tie(edgeIteratorBegin, edgeIteratorEnd) = edges(directedGraph); edgeIteratorBegin != edgeIteratorEnd; edgeIteratorBegin++)
		{
			// Get the current edge's source node
			int sourceNodeCost = nodeCostMap[source(*edgeIteratorBegin, directedGraph)];

			// Get the current edge's target node
			int targetNodeCost = nodeCostMap[target(*edgeIteratorBegin, directedGraph)];

			// Get the current edge's weight
			int edgeWeight = boostEdgeWeightMap[*edgeIteratorBegin];
			
			// If the current edge's source node cost is not INT_MAX and the current edge's source node cost + the current edge's weight is less than the current edge's target node cost
			if(sourceNodeCost != INT_MAX && (sourceNodeCost + edgeWeight < targetNodeCost))
			{
				nodeCostMap[target(*edgeIteratorBegin, directedGraph)] = sourceNodeCost + edgeWeight;
			}
	 	}
	}

	
	for(tie(edgeIteratorBegin, edgeIteratorEnd) = edges(directedGraph); edgeIteratorBegin != edgeIteratorEnd; edgeIteratorBegin++)
	{
		// Get the current edge's source node
		int sourceNodeCost = nodeCostMap[source(*edgeIteratorBegin, directedGraph)];

		// Get the current edge's target node
		int targetNodeCost = nodeCostMap[target(*edgeIteratorBegin, directedGraph)];

		// Get the current edge's weight
		int edgeWeight = boostEdgeWeightMap[*edgeIteratorBegin];

		// If the current edge's source node cost is not INT_MAX and the current edge's source node cost + the current edge's weight is less than the current edge's target node cost
		if(sourceNodeCost != INT_MAX && (sourceNodeCost + edgeWeight < targetNodeCost))
		{
			// Return false since the directed graph contains a negative weight cycle
			return false;
		}
	}

	// Return true since the directed graph doesn't contain a negative weight cycle
	return true;
}

// Main function
int main()
{
	#pragma region Initialization

	// Create an empty boost directed graph
	DirectedGraph boostDirectedGraph;

	// Create an empty leda directed graph
	leda::graph ledaDirectedGraph;

	//Create an empty edge array
	edge_array<int> ledaEdgeWeightArray;
	
	// User graph option
	std::string graphOption;

	// Number of nodes
	int numberOfNodes;

	cout << "Choose the testing graph between grid or random." << endl; 

	// Read the graph type
	cin >> graphOption;

	cout << "Enter the number of nodes." << endl;

	// Read the number of nodes
	cin >> numberOfNodes;

	
	if(graphOption == "grid")
	{
		// Create a grid graph
		grid_graph(ledaDirectedGraph, numberOfNodes);

		// Create an edge array
		edge_array<int> edgeWeightArray(ledaDirectedGraph);

		// Copy the edge array
		ledaEdgeWeightArray = edgeWeightArray;

		// Initialise a random seed
		srand(time(NULL));

		// Create a leda edge
		leda::edge tempEdge;

		// Iterate through all the edges of the leda directed graph
		forall_edges(tempEdge, ledaDirectedGraph)
		{
			// Set the edge weight to a random number between 1 and 10
			int tempEdgeSourceNodeIndex = ledaDirectedGraph.index(ledaDirectedGraph.source(tempEdge));

			// Get the point representation of the edge's target node index
			int tempEdgeTargetNodeIndex = ledaDirectedGraph.index(ledaDirectedGraph.target(tempEdge));

			// Get the point representation of the edge's source node index
			div_t tempEdgeSourceNodeIndexDivResult = div(tempEdgeSourceNodeIndex, numberOfNodes);

			// Get the point representation of the edge's target node index
			div_t tempEdgeTargetNodeIndexDivResult = div(tempEdgeTargetNodeIndex, numberOfNodes);

		    // Check if the edge is a vertical edge that belong in the third quarter
			bool verticalEdgeThirdQuarterPresence = (tempEdgeSourceNodeIndexDivResult.quot >= (numberOfNodes/2)) && (tempEdgeSourceNodeIndexDivResult.rem <= (numberOfNodes/2)) && (tempEdgeTargetNodeIndexDivResult.quot > (numberOfNodes/2)) && (tempEdgeTargetNodeIndexDivResult.rem < (numberOfNodes/2));

			// Check if the edge is a horizontal edge that belong in the third quarter
			bool horizontalEdgeThirdQuarterPresence = (tempEdgeSourceNodeIndexDivResult.quot > (numberOfNodes/2)) && (tempEdgeSourceNodeIndexDivResult.rem < (numberOfNodes/2)) && (tempEdgeTargetNodeIndexDivResult.quot >= (numberOfNodes/2)) && (tempEdgeTargetNodeIndexDivResult.rem <= (numberOfNodes/2));

			// Check if the edge is a vertical edge that belong in the third quarter or a horizontal edge that belong in the third quarter
			if(verticalEdgeThirdQuarterPresence || horizontalEdgeThirdQuarterPresence)
			{
				// Check if the edge is the special horizontal edge
				if((rand() % 2) == 0)
				{
					// Reverse the current edge
					ledaDirectedGraph.rev_edge(tempEdge);
				}
				
				// Check if the edge is the special horizontal edge
				bool specialThirdQuarterHorizontalEdge = (tempEdgeSourceNodeIndexDivResult.quot == (numberOfNodes/2 + 1)) && (tempEdgeSourceNodeIndexDivResult.rem == (numberOfNodes/2 - 1)) && (tempEdgeTargetNodeIndexDivResult.quot == (numberOfNodes/2 + 1)) && (tempEdgeTargetNodeIndexDivResult.rem == (numberOfNodes/2));

				// Check if the edge is the special vertical edge
				bool specialThirdQuarterVerticalEdge = (tempEdgeSourceNodeIndexDivResult.quot == (numberOfNodes/2)) && (tempEdgeSourceNodeIndexDivResult.rem == (numberOfNodes/2 - 1)) && (tempEdgeTargetNodeIndexDivResult.quot == (numberOfNodes/2 + 1)) && (tempEdgeTargetNodeIndexDivResult.rem == (numberOfNodes/2 - 1));

				// Check if the edge is the special horizontal edge or the special vertical edge
				if(specialThirdQuarterVerticalEdge || specialThirdQuarterHorizontalEdge)
				{
					// Assign a negative cost to the special edge
					ledaEdgeWeightArray[tempEdge] = -100000;
				}
			}
			else
			{
				// Assign random integer values as costs between 0 and 10000
				ledaEdgeWeightArray[tempEdge] = (rand() % 10000);
			}
		}
		// Get the grid graph number of nodes
		numberOfNodes = ledaDirectedGraph.number_of_nodes();

		// Copy the leda directed graph to the boost directed graph
		CopyLedaGraphToBoostGraph(boostDirectedGraph, ledaDirectedGraph, ledaEdgeWeightArray);

		// Initialise a property map that contain the boost graph edges weights
		EdgeWeightMap boostEdgeWeightMap = get(edge_weight, boostDirectedGraph);
	}
	else
	{
		
		if(graphOption == "random")
		{
			// Calculate the number of edges
			int numberOfEdges = ceil(20 * numberOfNodes * log2(numberOfNodes));

			// Generate a random directed graph
			random_graph(ledaDirectedGraph, numberOfNodes, numberOfEdges, false, true, true);

			// Make the graph cohesive
			Make_Connected(ledaDirectedGraph);

			// Intialise an edge array that will contain the leda graph edges weights
			edge_array<int> edgeWeightArray(ledaDirectedGraph);

			// Copy the edge array
			ledaEdgeWeightArray = edgeWeightArray;

			// Initialise a random seed
			srand(time(NULL));

			// Edge that will be used for the iteration
			leda::edge tempEdge;

			
			forall_edges(tempEdge, ledaDirectedGraph)
			{
				// Assign random integer values as costs between 10 and 10000
				ledaEdgeWeightArray[tempEdge] = (rand() % 10100) - 100;
			}

			// Copy the leda directed graph to the boost directed graph
			CopyLedaGraphToBoostGraph(boostDirectedGraph, ledaDirectedGraph, ledaEdgeWeightArray);

			// Initialise a property map that contain the boost graph edges weights
			EdgeWeightMap boostEdgeWeightMap = get(edge_weight, boostDirectedGraph);
		}
		else
		{
			cout << "Choose between grid or random." << endl;

			// Exit the program
			exit(0);
		}
	}

	// Initialise a node array the will contain the last edge on a shortest path from the starting node to a node
	node randomLedaNode = ledaDirectedGraph.choose_node();

	// Get the point representation of the random node
	Vertex randomBoostVertex = vertex(ledaDirectedGraph.index(randomLedaNode), boostDirectedGraph);
	
	// Initialise a node array the will contain the last edge on a shortest path from the starting node to a node
	node_array<leda::edge> ledaPredecessorNodeArray(ledaDirectedGraph);

	// Initialise a node array the will contain the distance from the starting node to a node
	node_array<int> ledaDistanceNodeArray(ledaDirectedGraph);

	// Initialize a property map that contain the edges's weights
	EdgeWeightMap boostEdgeWeightMap = get(edge_weight, boostDirectedGraph);

	// Initialize a property map that contain the vertices's distances
	std::vector<int> boostDistanceVector(numberOfNodes,INT_MAX);

	// Set the distance of the chosen vertex to 0
	boostDistanceVector[randomBoostVertex] = 0;

	// Initialize a vector that will contain the predeccesor of each node
	std::vector<std::size_t> boostPredeccesorVector(numberOfNodes);

	for(int index = 0; index < numberOfNodes; index++)
	{
		// Set the current vertex predeccesor to itself
		boostPredeccesorVector[index] = index;
	}
	
	#pragma endregion Initialization

	#pragma region Simulation

	// Initialise the starting CPU time
	float CPUTime = used_time();

	// Execute the user defined Bellman Ford algorithm for the boost directed graph using the random node
	bool negativeWeightCircleNotFound = BellmanFord(boostDirectedGraph, randomBoostVertex);

	//	Print the user defined Bellman Ford function execution time
	cout << "User defined Bellman Ford function execution time: " << used_time(CPUTime) << " seconds."<< endl;

	// Execute the Leda Bellman Ford algorithm for the Leda directed graph using the random node and the defined arrays
	negativeWeightCircleNotFound = BELLMAN_FORD(ledaDirectedGraph, randomLedaNode, ledaEdgeWeightArray, ledaDistanceNodeArray, ledaPredecessorNodeArray);

	// Print the Leda Bellman Ford function execution time
	cout << "Leda Bellman Ford function execution time: " << used_time(CPUTime) << " seconds."<< endl;

	// If the directed graph doesn't contain a negative weight circle...
	if(negativeWeightCircleNotFound)
	{
		// Execute the Boost Bellman Ford algorithm for the boost directed graph using the defined maps
		bellman_ford_shortest_paths(boostDirectedGraph, numberOfNodes, weight_map(boostEdgeWeightMap).distance_map(&boostDistanceVector[0]).distance_map(&boostPredeccesorVector[0]));

		// Print the Boost Bellman Ford function execution time
		cout << "Boost Bellman Ford function execution time: " << used_time(CPUTime) << " seconds."<< endl;

		cout << "The directed graph doesn't contain a negative weight cycle." << endl; 
	}
	else
	{
		cout << "The directed graph contains a negative weight cycle." << endl;
	}

	#pragma endregion Simulation

	// Return 0
	return 0;
}