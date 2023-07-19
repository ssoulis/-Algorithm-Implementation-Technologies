#include <iostream>
#include <time.h>
#include <LEDA/graph/graph_misc.h>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/graph_alg.h>
#include <LEDA/system/basic.h>
#include <LEDA/core/dynamic_trees.h>
#include <algorithm>
#include <cstdio>
#include <LEDA/graph/graph_gen.h>


using namespace std;
using namespace leda;

// Define the graph type

class TreeCohesiveComponent
{
	#pragma region Properties

	public:

		list<node> NodeList;

	private:

		int Size;

	
		node firstnodepointer;

		node lastnodepointer;

	#pragma endregion Properties

	#pragma region Public Methods

	public:

		// Default constructor

		TreeCohesiveComponent()
		{
		
			Size = 0;

			firstnodepointer = NULL;

			lastnodepointer = NULL;
		}


		// Constructor for a single node

		TreeCohesiveComponent(node Node)
		{
			
			Size = 1;

			firstnodepointer = Node;

			lastnodepointer = Node;

			NodeList.append(Node);
		}

		// Constructor for a list of nodes
        
		int GetSize()
		{
			return Size;
		}

       

		node GetFirstNode()
		{
			return firstnodepointer;
		}

		

		node GetLastNode()
		{
			return lastnodepointer;
		}

		//  Node size list getter

		void SetSize(int ListSize)
		{
			Size = ListSize;
		}

		// First node pointer setter

		void Setfirstnodepointer(node NodePointer)
		{
			firstnodepointer = NodePointer;
		}

		// Last node pointer setter

		void Setlastnodepointer(node NodePointer)
		{
			lastnodepointer = NodePointer;
		}


        // Overload the == operator

		friend bool operator == (TreeCohesiveComponent const &LeftOperand, TreeCohesiveComponent const &RightOperand) 
		{
			// Check if the two operands are the same

         	return (LeftOperand.firstnodepointer == RightOperand.firstnodepointer) && (LeftOperand.lastnodepointer == RightOperand.lastnodepointer) && (LeftOperand.Size == RightOperand.Size); 
    	}
		
	#pragma endregion Public Methods
};


// Compare the cost of two edges

int EdgeCostCompare(const edge &a, const edge &b)
{
	
	GRAPH<int, int> testGraph;


	if(testGraph.inf(a) < testGraph.inf(b))
	{
		// Set a before b

		return -1;
	}	
	else
	{
		// Set b before a
		
		return 0;
	}
}



// Search the tree cohesive component list for a node

TreeCohesiveComponent SearchNode(list<TreeCohesiveComponent> &ComponentList, node Node)
{
	
	TreeCohesiveComponent tempTreeCohesiveComponent;

	node tempNode;

	forall(tempTreeCohesiveComponent, ComponentList)
	{
		forall(tempNode, tempTreeCohesiveComponent.NodeList)
		{
			if(tempNode == Node)
			{
				return tempTreeCohesiveComponent;
			}
		}
	}

	return tempTreeCohesiveComponent;
}

// Merge two tree cohesive components

TreeCohesiveComponent MergeTreeCohesiveComponents(TreeCohesiveComponent FirstTreeCohesiveComponent, TreeCohesiveComponent SecondTreeCohesiveComponent)
{

	node tempNode;

	// Check which tree cohesive component is smaller

	if(FirstTreeCohesiveComponent.GetSize() >= SecondTreeCohesiveComponent.GetSize())
	{
		
		forall(tempNode, SecondTreeCohesiveComponent.NodeList)
		{
		
			FirstTreeCohesiveComponent.NodeList.append(tempNode);
		}

		// Set the new node list size

		FirstTreeCohesiveComponent.SetSize(FirstTreeCohesiveComponent.NodeList.size());

		// Set the new last node pointer

		FirstTreeCohesiveComponent.Setlastnodepointer(FirstTreeCohesiveComponent.NodeList.tail());

		// Return the merged tree cohesive component

		return FirstTreeCohesiveComponent;
	}
	else
	{
		// Merge the second tree cohesive component into the first one

		forall(tempNode, FirstTreeCohesiveComponent.NodeList)
		{
			// Append the node to the first tree cohesive component

			SecondTreeCohesiveComponent.NodeList.append(tempNode);
		}

		// Set the new node list size

		SecondTreeCohesiveComponent.SetSize(SecondTreeCohesiveComponent.NodeList.size());

		// Set the new last node pointer

		SecondTreeCohesiveComponent.Setlastnodepointer(SecondTreeCohesiveComponent.NodeList.tail());

		
		return SecondTreeCohesiveComponent;
	}
}

// Add an edge to the tree cohesive component list

bool TreeCohesiveComponentAddEdge(list<TreeCohesiveComponent>& TreeCohesiveComponentList, node EdgeSourceNode, node EdgeTargetNode)
{
	// Search the tree cohesive component list for the edge source node

	TreeCohesiveComponent sourceNodeTreeCohesiveComponent = SearchNode(TreeCohesiveComponentList, EdgeSourceNode);

    // Search the tree cohesive component list for the edge target node

	TreeCohesiveComponent targetNodeTreeCohesiveComponent = SearchNode(TreeCohesiveComponentList, EdgeTargetNode);

	//  Check if the edge source and target nodes are in the same tree cohesive component

	if(sourceNodeTreeCohesiveComponent.GetFirstNode() == targetNodeTreeCohesiveComponent.GetFirstNode())
	{
		// Return true because a circle was created

		return true;
	}

	else
	{
		// Merge the two tree cohesive components

		TreeCohesiveComponentList.append(MergeTreeCohesiveComponents(sourceNodeTreeCohesiveComponent, targetNodeTreeCohesiveComponent));

		//  Remove the old tree cohesive component that contained the edge's source node
		
		TreeCohesiveComponentList.remove(sourceNodeTreeCohesiveComponent);

		//  Remove the old tree cohesive component that contained the edge's target node

		TreeCohesiveComponentList.remove(targetNodeTreeCohesiveComponent);
	}
	
	//  Return false because no circle was created

	return false;
}

// Kruscal's algorithm

GRAPH<int, int> Kruscal(GRAPH<int, int>& UndirectedGraph, list<edge>& CircleEdgeList)
{
	// Minimum spanning tree

	GRAPH<int, int> MinSpanningTree(UndirectedGraph);

	// List that will contain the edges of the minimum spanning tree
    
	list<edge> edgeList = MinSpanningTree.all_edges();

	// Delete all edges from the minimum spanning tree

	MinSpanningTree.del_all_edges();

	// Sort the edge list by cost

	edgeList.sort(EdgeCostCompare);

	// List that will contain the tree cohesive components

	list<TreeCohesiveComponent> TreeCohesiveComponentList;

	// Node that will be used for the iteration

	node tempNode;

	
	forall_nodes(tempNode, MinSpanningTree)
	{
		
		TreeCohesiveComponent treeNode(tempNode);

		
		TreeCohesiveComponentList.append(treeNode);
	}

	
	edge tempEdge;

	// List that will contain the edges of the minimum spanning tree
    
	list<edge> MinSpanningTreeEdgesList;
    

	forall(tempEdge, edgeList)
	{
		// Get the edge's source node

		node tempEdgeSourceNode = UndirectedGraph.source(tempEdge);

		//  Get the edge's target node

		node tempEdgeTargetNode = UndirectedGraph.target(tempEdge);

		// Get the edge's cost

		int tempEdgeCost = UndirectedGraph.inf(tempEdge);

		// Check if the edge's source and target nodes are in the same tree cohesive component

		bool graphCircleFound = TreeCohesiveComponentAddEdge(TreeCohesiveComponentList, tempEdgeSourceNode, tempEdgeTargetNode);

		// Check if a circle was created

		if(!graphCircleFound)
		{
			
			MinSpanningTreeEdgesList.append(tempEdge);
		}
		else
		{
			
			CircleEdgeList.append(tempEdge);
		}
	}


	// Add the edges of the minimum spanning tree to the minimum spanning tree

	forall(tempEdge,MinSpanningTreeEdgesList)
	{
		
		MinSpanningTree.new_edge(MinSpanningTree.source(tempEdge), MinSpanningTree.target(tempEdge), MinSpanningTree.inf(tempEdge));
	}

	return MinSpanningTree;
}

// Minimum spanning tree validator

bool MinSpanningTreeValidator(GRAPH<int, int>& UndirectedGraph, GRAPH<int, int>& MinSpanningTree, list<edge>& CircleEdgeList)
{
	// List that will contain the edges of the minimum spanning tree

	list<edge> MinSpanningTreeEdgeList = MinSpanningTree.all_edges();

	// Delete all edges from the minimum spanning tree

	leda_dynamic_trees dynamicTree;

	//  Edge that will be used for the iteration

	edge tempEdge;

	// Node that will be used for the iteration

	node tempNode;

	// Map that will associate a vertex for each node in the minimum spanning tree

	node_map<vertex> nodeVertexMap(MinSpanningTree, NULL);



	forall_nodes(tempNode, MinSpanningTree)
	{


		nodeVertexMap[tempNode] = dynamicTree.make();
	}



	forall(tempEdge, MinSpanningTreeEdgeList)
	{
		
		dynamicTree.evert(nodeVertexMap[MinSpanningTree.source(tempEdge)]);

	
		dynamicTree.evert(nodeVertexMap[MinSpanningTree.target(tempEdge)]);

		
		dynamicTree.link(nodeVertexMap[MinSpanningTree.source(tempEdge)], nodeVertexMap[MinSpanningTree.target(tempEdge)], MinSpanningTree.inf(tempEdge));
	}
	
	
	forall(tempEdge, CircleEdgeList)
	{

		vertex edgeSourceNodeVertex = nodeVertexMap[MinSpanningTree.source(tempEdge)];


		vertex edgeTargetNodeVertex = nodeVertexMap[MinSpanningTree.target(tempEdge)];


		vertex leastCommonAncestorVertex = dynamicTree.lca(edgeSourceNodeVertex,edgeTargetNodeVertex);
		
		
        // Search the path between the edge source node vertex and the least common ancestor vertex

		while(edgeSourceNodeVertex != leastCommonAncestorVertex)
		{
			// 
			if(MinSpanningTree.inf(tempEdge) >= dynamicTree.cost(edgeSourceNodeVertex))
			{
				
				edgeSourceNodeVertex = dynamicTree.parent(edgeSourceNodeVertex);
			}
			else
			{
				
				return false;
			}
		}

		// Search the path between the edge target node vertex and the least common ancestor vertex

		while(edgeTargetNodeVertex != leastCommonAncestorVertex)
		{
			
			if(MinSpanningTree.inf(tempEdge) >= dynamicTree.cost(edgeTargetNodeVertex))
			{
				
				edgeTargetNodeVertex = dynamicTree.parent(edgeTargetNodeVertex);
			}
			else
			{
				
				return false;
			}
		}
	}
	
	
	return true;
}




// Main function

int main()
{
	#pragma region Initialization

	// Undirected graph

	GRAPH<int, int> UndirectedGraphNew;
    

	// Minimum spanning tree

	list<edge> circleEdgeList;

    #pragma endregion

	std::string GraphOptions;
            

	int NumNodes;

	cout << "Choose the graph you want to test (options : grid , rand , synth )" << endl; 

	cin >> GraphOptions;

	cout << "Enter the number of nodes" << endl;

    // Read the number of nodes

	cin >> NumNodes;

	edge tempEdge;



    // If the grid graph is selected

	if(GraphOptions == "grid")
	{

		grid_graph(UndirectedGraphNew, NumNodes);
		
	
		srand(time(NULL));

	
		forall_edges(tempEdge, UndirectedGraphNew)
		{
            // random number between 10 and 1000

			UndirectedGraphNew.assign(tempEdge, (rand() % 990) + 10);
		}
	}
	else if (GraphOptions == "rand")	// If rand graph is selected
	{
		
				GRAPH InitGraph;
                        
				// 2nlog2(n) edges

				int NumEdges = ceil(2 * NumNodes * log2(NumNodes));

				// Create a rand undirected graph 
				rand_simple_undirected_graph(InitGraph, NumNodes, NumEdges);
				
				// Make the graph cohesive
				Make_Biconnected(InitGraph);

				// Copy the GRAPH to the GRAPH<int, int>
				CopyGraph(UndirectedGraphNew, InitGraph);

				// Initialize a rand seed
				srand(time(NULL));

				// For every edge in the undirected graph...
				forall_edges(tempEdge, UndirectedGraphNew)
				{
					// Assign rand integer values as costs between 10 and 10000
					UndirectedGraphNew.assign(tempEdge, (rand() % 990) + 10);
				}
			}
	else if (GraphOptions == "synth") 			// If the synth graph is selected
	{
			

			
			complete_GRAPH(UndirectedGraphNew, NumNodes);

		
			forall_edges(tempEdge,UndirectedGraphNew)
			{
				// worst case is P-1

				UndirectedGraphNew.assign(tempEdge, UndirectedGraphNew.index(UndirectedGraphNew.target(tempEdge)) - 1);
			}
	}
	else
	{
				cout << " Try again, choose between : (grid, synth , rand) " << endl;
		
				exit(0);
	}
	
	
	// Initialize the edge cost array

	edge_array<int> edgeCostArray(UndirectedGraphNew);

	#pragma endregion Initialization


	// minimun spanning tree	

	GRAPH<int, int> MinSpanningTree = Kruscal(UndirectedGraphNew, circleEdgeList);



	// Initialize the execution time

	float MyTime = used_time();



	// kruscal execution time

	cout << "User defined Kruscal function execution time: " << used_time(MyTime) << endl;

	// minimum spanning tree 

	list<edge> min_tree_edges = MIN_SPANNING_TREE(UndirectedGraphNew,edgeCostArray);

	// Print Kruscal execution time

	cout << "Leda minimum spanning tree function execution time: " << used_time(MyTime) << endl;
	
	// Print the minimum spanning tree
	
	bool MinSpanningTreeValidation = MinSpanningTreeValidator(UndirectedGraphNew, MinSpanningTree, circleEdgeList);

	// if the minimum spanning tree is valid

	if(MinSpanningTreeValidation)
	{
		cout << "The minimum spanning tree was validated." << endl; 
	}
	else
	{
		cout << "The minimum spanning tree is invalid." << endl;
	}
	
	
	return 0;
}