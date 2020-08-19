#ifndef __BinaryHeap_h_
#define __BinaryHeap_h_

#include <iostream>
#include <cassert>

/**
 * This class defines the binary heap data structure. It is an array
 * data structure that has the following property at all times:
 *
 * A[i] <= A[2(i+1)-1] and A[i] <= A[2(i+1)]
 *
 * where i is a 0-based index (this is more elegant for 1-based indices,
 * but this is a C program, so be it).
 *
 * The heap is allocated at the beginning to hold N elements and the 
 * weights of the N elements are supplied in the constructor. This is
 * different from the STL approach to priority queues. The nice thing
 * is that there is no memory allocation outside of the constructor.
 *
 * Another big difference between this heap and STL priority queue is that
 * this heap allows you to change the weight of an element already in the
 * heap, as well as to check whether an element belongs to the heap.
 *
 * This heap is intended to be used with the implementation of
 * Dijkstra's shortest path algorithm
 */
template<class TWeight>
class BinaryHeap {
public:
  /** 
   * Allocate the memory for the heap, passing in the array
   * of weights. The weights can be changed later, but that
   * should be done using the UpdateWeight method, not directly,
   * as the heap property would be violated.
   */
  BinaryHeap(unsigned int nWeights, TWeight *inWeightArray)
    {
    m_WeightArray = inWeightArray;
    m_ReserveSize = nWeights;
    m_HeapIndex = new int[nWeights];
    m_Heap = new unsigned int[nWeights];
    m_HeapSize = 0;
    }

  ~BinaryHeap() 
    {
    delete[] m_Heap;
    delete[] m_HeapIndex;
    }

  /** 
   * Reinitialize the heap to full size of the weights array. This is used in 
   * conjunction with Dijkstra's algorithm. The weights are all set to the 
   * specified value in order to maintain the heap property.
   *
   * This operation is O(n)
   */
  void InsertAllElementsWithEqualWeights(TWeight weight)
    { 
    m_HeapSize = m_ReserveSize;
    for(int i=0;i<m_ReserveSize;i++)
      {
      m_WeightArray[i] = weight;
      Put(i,i);
      }
    }

  /** 
   * Insert an element into the heap. 
   * 
   * This operation is O(log n) 
   */
  void InsertElement(unsigned int iElement)
    { SiftUp(iElement, m_HeapSize++); }

  /** 
   * Extract the smallest element from the heap, removing it 
   * 
   * This operation is O(log n) 
   */
  unsigned int PopMinimum()
    {
    assert(m_HeapSize > 0);
    unsigned int rtn = m_Heap[0];
    Put(0, m_Heap[m_HeapSize-1]);
    m_HeapSize--;
    SiftDown(0,false);

    // Update the position of the element in the heap
    m_HeapIndex[rtn] = m_ReserveSize;

    return rtn;
    }

  /** 
   * Lower the weight of an element and reorder the heap accordingly.
   * The parameter to the call is the new weight of the element, not
   * the difference. It's important that the new weight is less than
   * the old weight! It's also important that the element is in the 
   * heap (use ContainsElement to check)
   *
   * This operation is O(log n) with small constant term
   */
  void DecreaseElementWeight(unsigned int iElement, TWeight xNewWeight)
    { 
    // Check the predicates
    assert(m_HeapIndex[iElement] < m_HeapSize
      && xNewWeight <= m_WeightArray[iElement]);

    // Change the weight
    m_WeightArray[iElement] = xNewWeight;

    // Find the place for the element upstream
    SiftUp(iElement, m_HeapIndex[iElement]);
    }

  /** 
   * Increase the weight of an element. 
   * Same usage as DecreaseElementWeight
   */
  void IncreaseElementWeight(unsigned int iElement, TWeight xNewWeight)
    { 
    // Check the predicates
    assert(m_HeapIndex[iElement] < m_HeapSize
      && xNewWeight >= m_WeightArray[iElement]);

    // Change the weight
    m_WeightArray[iElement] = xNewWeight;

    // Find the place for the element down stream
    SiftDown(m_HeapIndex[iElement]);
    }

  /** 
   * Change the weight of an element and reorder the heap accordingly.
   * This is a bit less efficient than calling DecreaseElementWeight or 
   * IncreaseElementWeight directly if you know apriori if the weight
   * is increased or decreassed.
   * 
   * This operation is O(log n) however 
   */
  void UpdateElementWeight(unsigned int iElement, TWeight weight)
    {
    // Check the direction
    if(weight < m_WeightArray[iElement])
      DecreaseElementWeight(iElement, weight);
    else
      IncreaseElementWeight(iElement, weight);

    // Push the lement down to the leaf
    // int iPos = SiftDown(m_HeapIndex[iElement], true);

    // Update the weight
    // m_WeightArray[iElement] = weight;

    // Perform insertion on the element from that position
    // SiftUp(iElement, iPos);
    }

  // Print the heap (prints on one line, not useful for big heaps)
  void PrintHeap(std::ostream &out)
    {
    for(unsigned int j=0;j<m_ReserveSize;j++)
      std::cout << "[" << j << "," << m_WeightArray[j] << "] ";
    std::cout << std::endl;
    for(unsigned int i=0;i<m_HeapSize;i++)
      std::cout << "[" << i << "," << m_Heap[i] << "," 
        << m_WeightArray[m_Heap[i]] << "] ";
    std::cout << std::endl;
    }

  /** 
   * Get the number of elements currenly in the heap. This is
   * not the same as the capacity of the heap, i.e., number of
   * weights
   */
  unsigned int GetSize()
    { return m_HeapSize; }

  /**
   * Check if an element is in the heap or not.
   *
   * This operation is O(1)
   */
  bool ContainsElement(unsigned int iPos)
    { return m_HeapIndex[iPos] < m_HeapSize; }

private:
  // The number of elements allocated
  int m_ReserveSize;

  // The number of elements in the heap
  int m_HeapSize;

  // The weights associated with the heap
  TWeight *m_WeightArray;

  // The position in the heap of each element
  int *m_HeapIndex;

  // The heap array
  unsigned int *m_Heap;

  // Get the parent of an element
  inline int Parent(int x) 
    { return ((x + 1) >> 1) - 1; }

  // Get the element to the right
  inline int Right(int x)
    { return ((x + 1) << 1); }

  // Get the element to the left
  inline int Left(int x)
    { return ((x + 1) << 1) - 1; }

  // Put an element in the heap at position x
  inline void Put(unsigned int iPos, unsigned int iElt)
    {
    m_Heap[iPos] = iElt;
    m_HeapIndex[iElt] = iPos;
    }

  // Get the weight at heap index
  TWeight WeightAt(unsigned int iPos)
    { return m_WeightArray[m_Heap[iPos]]; }
  
  /** Push an element from position iPos up until it fits */
  void SiftUp(unsigned int iElement, int iPos)
    {
    while(iPos > 0 && WeightAt(Parent(iPos)) > m_WeightArray[iElement])
      {
      Put(iPos, m_Heap[Parent(iPos)]);
      iPos = Parent(iPos);
      }
    
    Put(iPos, iElement);
    }

  /** Push the element iPos down until it fits. If infWeight is set
   * the element comes down to the bottom. Returns position where the
   * element settled.*/
  int SiftDown(int iPos, bool infWeight)
    {
    int l = Left(iPos);
    int r = Right(iPos);
    int largest = iPos;
    
    if(l < m_HeapSize && (infWeight || WeightAt(l) < WeightAt(iPos)))
      {
      largest = l;
      }
    
    if(r < m_HeapSize && (infWeight || WeightAt(r) < WeightAt(largest)))
      {
      largest = r;
      }
    
    if(largest != iPos)
      {
      // Swap iPos, largest
      unsigned int tmp = m_Heap[iPos];
      Put(iPos, m_Heap[largest]);
      Put(largest, tmp);

      // Recurse
      return SiftDown(largest,infWeight);
      }
    else return iPos;
    }

};

#endif
