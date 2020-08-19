#include "ShortestPath.h"

using namespace std;

void TestShortestPath()
{
  // Create a graph (taken from website below)
  // http://ciips.ee.uwa.edu.au/~morris/Year2/PLDS210/dij-op.html
  unsigned int AI[] = {0,2,4,5,7,10};
  unsigned int A[] =  {1,  4,  2,  4,  3,  2,  0,  1,  2,  3};
  unsigned int w[] =  {10, 5,  1,  2,  4,  6,  7,  3,  9,  2};
  
  // Create shortest pather
  DijkstraShortestPath<unsigned int> sp(5, AI, A, w);

  // Compute shortest path from starting point
  sp.ComputePathsFromSource(0);

  // Reconstruct the shortest path to each vertex
  for(unsigned int i=1; i<5; i++)
    {
    std::cout << "Path to 0 from " << i << " : ";
    unsigned int j = i;
    while(j != 0)
      {
      std::cout << j << " - ";
      j = sp.GetPredecessorArray()[j];
      }
    std::cout << "0" << std::endl;
    }
}

/* 
int main(int argc, char *argv[])
{
  TestShortestPath();
}
*/ 


