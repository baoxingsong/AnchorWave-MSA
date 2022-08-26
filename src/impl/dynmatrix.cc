#include "dynmatrix.h"
#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <map>
#include <set>
#include <regex>

double MIN_DISTANCE = 0.0;
double Max_DISTANCE =std::numeric_limits<double>::max();
void addCluster(ClusterNode *&head, ClusterNode *&tail, const std::string& name)
// adds a cluster (at the tail) and the corresponding row and column to data structure
// distance of all added DistanceNodes should be initialized to 0.0
// at the end, tail should point to the newly added ClusterNode
{
  // Initialize ClusterNode with name and the number of clusters inside
  ClusterNode *newNode = new ClusterNode;
  newNode->name = name;
  newNode->numClusters = 1; // Used to keep track of how many clusters have been combined into one ClusterNode

  if (head == NULL) { // If we are adding a cluster to an empty DynMatrix
    DistanceNode *firstDistance = new DistanceNode;
    firstDistance->distance = Max_DISTANCE;

    head = tail = newNode;
    head->column = head->row = firstDistance;
  }

  else {
    // Reattach pointers of the added ClusterNode and set tail to newly created ClusterNode
    newNode->prev = tail;
    tail->next = newNode;
    tail = newNode;
    // Ensure tail->next is set to NULL
    tail->next = NULL;

    // Initialize first DistanceNodes
    DistanceNode *firstCol = new DistanceNode;
    DistanceNode *firstRow = new DistanceNode;
    firstCol->distance = Max_DISTANCE;
    firstRow->distance = Max_DISTANCE;

    // Ensure pointers are set to NULL to prevent pointing to wrong locations
    firstCol->nextInColumn = NULL;
    firstCol->nextInRow = NULL;

    firstRow->nextInColumn = NULL;
    firstRow->nextInRow = NULL;

    tail->column = firstCol;
    tail->row = firstRow;

    // Attach first DistanceNodes to first column and row of newly added ClusterNode
    DistanceNode *currCol = tail->column;
    DistanceNode *currRow = tail->row;

    DistanceNode *prevCol = tail->prev->column;
    DistanceNode *prevRow = tail->prev->row;

    prevCol->nextInRow = currCol;
    prevRow->nextInColumn = currRow;

    // Go through each DistanceNode and attach pointers like a singly-linked list
    while (prevCol != prevRow) {
      // Initialize DistanceNodes for use
      DistanceNode *newCol = new DistanceNode;
      DistanceNode *newRow = new DistanceNode;
      newCol->distance = Max_DISTANCE;
      newRow->distance = Max_DISTANCE;

      // Ensure pointers are set to NULL to prevent pointing to wrong locations
      newCol->nextInColumn = NULL;
      newCol->nextInRow = NULL;

      newRow->nextInColumn = NULL;
      newRow->nextInRow = NULL;

      currCol->nextInColumn = newCol;
      currRow->nextInRow = newRow;

      currCol = currCol->nextInColumn;
      currRow = currRow->nextInRow;

      prevCol = prevCol->nextInColumn;
      prevRow = prevRow->nextInRow;

      prevCol->nextInRow = currCol;
      prevRow->nextInColumn = currRow;

    }
    // Initialize and attach the last DistanceNode in newly added ClusterNode
    DistanceNode *edgeNode = new DistanceNode;
    edgeNode->distance = Max_DISTANCE;

    currCol->nextInColumn = currRow->nextInRow = edgeNode;
  }
}

void removeCluster(ClusterNode *&head,ClusterNode *&tail,ClusterNode *toBeRemoved)
// removes a cluster pointed to by toBeRemoved and the corresponding row and column
{

  // If there is only one ClusterNode left in the list, delete the only node and set pointers to NULL for addCluster to function correctly
  if (toBeRemoved == head && toBeRemoved == tail) {
    ClusterNode *tempNode = head;

    tempNode->next = NULL;
    tempNode->prev = NULL;

    head = NULL;
    tail = NULL;
    delete tempNode;
    return;
  }
  if (toBeRemoved == head) {
    // Steps:
    // reassign head to head->next
    // go down ->nextInColumn until NULL, then go down ->nextInRow until NULL and remove pointers

    ClusterNode *tempHead = head;
    ClusterNode *tempHeadRow;
    DistanceNode *tempVerticalDistance;
    DistanceNode *tempHeadCol;
    head = head->next;
    head->prev->next = NULL;
    head->prev = NULL;

    for (tempVerticalDistance = tempHead->column; tempVerticalDistance != NULL; tempVerticalDistance = tempVerticalDistance->nextInColumn) {
      tempVerticalDistance->nextInRow = NULL;
    }

    tempHeadCol = head->column->nextInColumn;

    for (tempHeadRow = head; tempHeadRow != NULL; tempHeadRow = tempHeadRow->next) {
      tempHeadRow->column = tempHeadRow->column->nextInColumn;
      tempHeadRow->row = tempHeadCol;

      tempHeadCol = tempHeadCol->nextInColumn;
    }
    tempHead = NULL;
    return;
  }
  else if (toBeRemoved == tail) {
    // Steps:
    // go down ->column->nextInColumn until NULL and remove pointers; also set pointers pointing to tail list to NULL
    // go down ->row->nextInRow until NULL and remove pointers; also set pointers pointing to tail list NULL;

    tail = tail->prev;
    tail->next = NULL;
    DistanceNode *tempDistanceRow;
    DistanceNode *tempDistanceCol;

    for (tempDistanceRow = tail->row, tempDistanceCol = tail->column; tempDistanceRow != tempDistanceCol; tempDistanceRow = tempDistanceRow->nextInRow, tempDistanceCol = tempDistanceCol->nextInColumn) {
      tempDistanceCol->nextInRow = NULL;
      tempDistanceRow->nextInColumn = NULL;

    }
    tempDistanceCol->nextInColumn = tempDistanceCol->nextInRow = NULL;


  }
  else {
    // Steps:
    // reattach pointers to cluster node
    // Go down nextInRow and nextInColumn until pointers are equal and remove + reattach pointers
    // go to ->prev and ->next cluster node and traverse down list
    ClusterNode *prevNode = toBeRemoved->prev;
    ClusterNode *nextNode = toBeRemoved->next;

    toBeRemoved->prev->next = toBeRemoved->next;
    toBeRemoved->next->prev = toBeRemoved->prev;

    DistanceNode *tempColPrev = prevNode->column;
    DistanceNode *tempColNext = nextNode->column;

    DistanceNode *tempRowPrev = prevNode->row;
    DistanceNode *tempRowNext = nextNode->row;

    for (;tempColPrev != NULL; tempColPrev = tempColPrev->nextInColumn, tempColNext = tempColNext->nextInColumn) {
      tempColPrev->nextInRow = tempColNext;
    }

    for (;tempRowPrev != NULL; tempRowPrev = tempRowPrev->nextInRow, tempRowNext = tempRowNext->nextInRow) {
      tempRowPrev->nextInColumn = tempRowNext;
    }
  }

}

/*
void findMinimum(ClusterNode *head,ClusterNode *&C,ClusterNode *&D)
// finds the minimum distance (between two different clusters) in the data structure 
// and returns the two clusters via C and D
{
  // Create temporary Cluster and Distance Nodes
  ClusterNode *curr = NULL;
  DistanceNode *tempStart = NULL;

  ClusterNode *secondCurr = NULL;

  ClusterNode *minClusterNode = NULL;
  DistanceNode *minDistanceNode = NULL;
  // Set minDistance to max Double value
  double minDistance = std::numeric_limits<double>::max();
  // Go through each DistanceNode and check for minimum value
  for (curr = head; curr != NULL; curr = curr->next) {
      int i=0;
    for (tempStart = curr->column; tempStart != NULL; tempStart = tempStart->nextInColumn) {
        std::cout<<"209"<<std::endl;
        std::cout<<"tempStart->distance   "<<tempStart->distance<<std::endl;
//      if (tempStart->distance < minDistance && tempStart->distance != 0) {
      if (tempStart->distance < minDistance) {
        minDistance = tempStart->distance;
        std::cout<<"minDistance   "<<minDistance<<std::endl;
        minClusterNode = curr;
        minDistanceNode = tempStart;
      }
      std::cout<<"216"<<std::endl;
    }
      std::cout<<"218"<<std::endl;
  }
  // set corresponding ClusterNode of minimum DistanceNode to given ClusterNode C
  C = minClusterNode;
    std::cout<<"c  "<<C<<std::endl;
  // Go through each DistanceNode and find matching minimum value
  for (secondCurr = minClusterNode->next; secondCurr != NULL; secondCurr = secondCurr->next) {
    for(tempStart = secondCurr->column; tempStart != NULL; tempStart = tempStart->nextInColumn) {
      if (tempStart->distance == minDistanceNode->distance) {
        // Set minimum distance value to global variable
        MIN_DISTANCE = minDistanceNode->distance;
        D = secondCurr;
        return;
      }
    }
  }
}
*/

void findMinimum(ClusterNode *head,ClusterNode *&C,ClusterNode *&D)
// finds the minimum distance (between two different clusters) in the data structure
// and returns the two clusters via C and D
{
    ClusterNode *tempC = head;
    ClusterNode *tempD = head;
    DistanceNode *temp;
    double globalMin = std::numeric_limits<double>::max();
    while(tempC){
        temp = tempC->row;
        while(tempD) {
            if(tempC==tempD){
                temp->distance=std::numeric_limits<double>::max();
            }
//            if(tempC != tempD && temp->distance != 0.0 && temp->distance < globalMin){
            if(tempC != tempD  && temp->distance < globalMin){
                globalMin = temp->distance;
                C = tempC;
                D = tempD;
            }
            temp = temp->nextInRow;
            tempD = tempD->next;
        }
        tempC = tempC->next;
        tempD = head;
    }
}



void UPGMA(ClusterNode *&head, ClusterNode *&tail,std::map<std::string, std::string> sequences,  std::vector<std::string> &Rnames,std::vector<std::string> &Qnames,    std::vector<std::string> &dna_refs,std::vector<std::string> &dna_queries)
// Implement the UPGMA Algorithm
{
    // Create temporary Cluster and Distance Nodes
    ClusterNode *tempHead = head;

    ClusterNode *clusterOne = NULL;
    ClusterNode *clusterTwo = NULL;

    DistanceNode *distanceOne;
    DistanceNode *distanceTwo;
    int j=0;
    // Create vectors to hold distances of DistanceNodes
    std::vector<double> clusterOneDistances;
    std::vector<double> clusterTwoDistances;
    std::vector<double> resultantDistances;

    // Call findMinimum function to locate two ClusterNodes containing the minimum distance in DynMatrix
    // clusterOne and clusterTwo are the ClusterNodes containing the minimum distance
    findMinimum(tempHead, clusterOne, clusterTwo);
    std::cout << clusterOne->name << "\tline 247\t" << clusterTwo->name << std::endl;


    // Store the distances in the two ClusterNodes' clusterOne and clusterTwo in vectors
    for (distanceOne = clusterOne->column; distanceOne != NULL; distanceOne = distanceOne->nextInColumn) {
        clusterOneDistances.push_back(distanceOne->distance);
    }
    for (distanceTwo = clusterTwo->column; distanceTwo != NULL,distanceTwo->distance!=Max_DISTANCE; distanceTwo = distanceTwo->nextInColumn) {
        clusterTwoDistances.push_back(distanceTwo->distance);
    }
    // Use useFormula function to compute average distances and store inside resultant vector
    resultantDistances = useFormula(clusterOne, clusterTwo, clusterOneDistances, clusterTwoDistances);
    std::smatch result;
    std::regex  pattern("[^(,)]+");

    std::string::const_iterator iterStartR = clusterOne->name.begin();
    std::string::const_iterator iterEndR = clusterOne->name.end();
    std::string::const_iterator iterStartQ = clusterTwo->name.begin();
    std::string::const_iterator iterEndQ = clusterTwo->name.end();

   while (regex_search(iterStartR, iterEndR, result, pattern)){
       Rnames.push_back(result[0]);
       iterStartR = result[0].second;
   }
   while (regex_search(iterStartQ, iterEndQ, result, pattern)){
       Qnames.push_back(result[0]);
       iterStartQ = result[0].second;
   }
    // Call combineCluster function to remove the two ClusterNodes and add the new ClusterNode into the DynMatrix
    combineCluster(head, tail, clusterOne, clusterTwo, resultantDistances);
    for (int j = 0; j < Rnames.size(); ++j) {
        for (std::map<std::string, std::string>::iterator it1 = sequences.begin(); it1 != sequences.end(); ++it1) {
            if (sequences.find(Rnames[j]) == it1) {
                dna_refs.push_back(it1->second);
            }
        }
    }
    for (int k = 0; k < Qnames.size(); ++k) {
        for (std::map<std::string, std::string>::iterator it2 = sequences.begin(); it2 != sequences.end(); ++it2) {
            if (sequences.find(Qnames[k]) == it2) {
                dna_queries.push_back(it2->second);
            }
        }
    }
    for (int i=0;i<dna_refs.size();++i){
        std::cout<<"cin "<<"dna_refs "<<dna_refs[i]<<std::endl;
    }
    std::cout<<std::endl;
    for (int i = 0; i < dna_queries.size(); ++i) {
        std::cout << "cin  "<<"dna_queries "<<dna_queries[i] << std::endl;
    }

}
void combineCluster(ClusterNode *&head,ClusterNode *&tail,ClusterNode *&C,ClusterNode *&D, std::vector<double> values)
// Adds a cluster using a combination of C->name and D->name, then fills the DistanceNodes of the resulting ClusterNode using the values of the given vector parameter
{
  // Set variables to make new ClusterNode name
  std::string leftparen = "(";
  std::string comma = ",";
  std::string rightparen = ")";
  // Assign values of vector parameter to new vector
  std::vector<double> upgmaValues = values;
  // Construct new ClusterNode name
  std::string name = leftparen + C->name +comma + D->name + rightparen;
  std::cout<<"name    "<<name<<std::endl;
  DistanceNode *firstDNCol;
  DistanceNode *firstDNRow;
  ClusterNode *tempAddedClusterNode = NULL;
  int i = 0; // initialize vector index
  // Remove the two ClusterNodes containing the minimum distance
  removeCluster(head, tail, C);
  removeCluster(head, tail, D);
  // Add a new cluster to the matrix, with 0.0 distance values for rows and columns.
  addCluster(head, tail, name);
  // Update numClusters of newly added ClusterNode to keep track of how many ClusterNodes are in the newly added one
  tempAddedClusterNode = tail;
  tempAddedClusterNode->numClusters = C->numClusters + D->numClusters;
  // Fill distance values of DistanceNodes in newly added ClusterNode with values from given vector
  for (firstDNCol = tempAddedClusterNode->column, firstDNRow = tempAddedClusterNode->row; firstDNCol != firstDNRow; firstDNCol = firstDNCol->nextInColumn, firstDNRow = firstDNRow->nextInRow) {
    firstDNCol->distance = upgmaValues[i];
    firstDNRow->distance = upgmaValues[i];
    i++;
  }
  // Ensure pointers are NULL to prevent accessing wrongly linked locations
  firstDNCol->nextInColumn = firstDNRow->nextInRow = NULL;
  firstDNCol->nextInRow = firstDNRow->nextInColumn = NULL;

}

std::vector<double> useFormula(ClusterNode *clusterOne, ClusterNode *clusterTwo, std::vector<double> firstValues, std::vector<double> secondValues)
// Uses given formula to calculate values to be put into each DistanceNode of added ClusterNode
{
  // Initialize and set variables
  std::vector<double> result;
  DistanceNode *tempOne;
  DistanceNode *tempTwo;
  int numClusterOne = clusterOne->numClusters;
  int numClusterTwo = clusterTwo->numClusters;
  double numerator = 0.0;
  double denominator = 0.0;

  // Loop through each DistanceNode
  for (tempOne = clusterOne->column, tempTwo = clusterTwo->column; tempOne != NULL; tempOne = tempOne->nextInColumn, tempTwo = tempTwo->nextInColumn) {
    // Only compute values that are not the minimum distance or zero
    // if ((tempOne->distance != MIN_DISTANCE && tempOne->distance != 0) && (tempTwo->distance != MIN_DISTANCE && tempTwo->distance != 0) ) {
    if ( (tempOne->distance!= Max_DISTANCE)&&(tempTwo->distance!=Max_DISTANCE) ) {
       // Compute numerator and denominator of formula
       numerator = (numClusterOne * tempOne->distance) + (numClusterTwo * tempTwo->distance);
       denominator = numClusterOne + numClusterTwo;
       // Push computed value into vector to be returned
       result.push_back(numerator/denominator);
    }
  }

  return result;
}


void printRowByRow(ClusterNode *head)
// Print DynMatrix row by row for debugging to ensure correctly linked Cluster and Distance Nodes
{
  while(head) {
    std::cout << std::setw(10) << head->name <<":\t";
    DistanceNode *curD = head->row;
    while(curD) {
      std::cout << curD->distance << "\t";
      curD = curD->nextInRow;
    }
    std::cout <<std::endl;
    head = head->next;

  }
}

void printColumnByColumn(ClusterNode *head)
// Print DynMatrix column by column for debugging to ensure correctly linked Cluster and Distance Nodes
{
  while(head) {
    std::cout << std::setw(10) << head->name <<":\t";
    DistanceNode *curD = head->column;
    while(curD) {
      std::cout << curD->distance << "\t";
      curD = curD->nextInColumn;
    }
    std::cout <<std::endl;
    head = head->next;
  }
}



void sequencesToDifferDistanceMatrix( std::map<std::string, std::string>& sequences,std::vector<std::string> & seqNames,float **differ_distanceMatrix) {
    size_t i = 0;
    std::map<std::string/*sequence name*/, std::map<std::string/*sequence name*/, int16_t> > differMatrix;
    std::map<std::string, std::string>::iterator  it;
    std::map<std::string, std::string>::iterator  it1;
    for(std::string seqName:seqNames) {
        it = sequences.find(seqName);
        size_t j = 0;
        for(std::string seqName:seqNames) {
            it1 = sequences.find(seqName);
            for (int k = 0; k < it->second.size(); ++k) {
                if (it->second[k] != it1->second[k]) {
                    differMatrix[it->first][it1->first]++;
                }
            }
            differ_distanceMatrix[i][j] = differMatrix[it->first][it1->first];
            j++;
        }
    i++;
    }
}


void updataSequences(std::map<std::string, std::string> &sequences, std::vector<std::string> &_alignment_ds,std::vector<std::string> &_alignment_qs,std::vector<std::stack<char>> &SRs,std::vector<std::stack<char>> &SQs, std::vector<std::string> &Rnames,std::vector<std::string> &Qnames,std::vector<std::string> &dna_refs,std::vector<std::string> &dna_queries){

    while (!SQs[0].empty()) {
        for( int k=0; k<dna_refs.size(); ++k ){
            _alignment_ds[k] += SRs[k].top();
            SRs[k].pop();
        }
        for( int l=0; l<dna_queries.size(); ++l ){
            _alignment_qs[l] += SQs[l].top();
            SQs[l].pop();
        }
    }

    for(int m=0 ,z=0;m<Rnames.size();++m,++z){
        for (std::map<std::string, std::string>::iterator it3 = sequences.begin();it3 != sequences.end(); ++it3) {
            if (it3->first == Rnames[m]){
                sequences[Rnames[m]] = _alignment_ds[z];
            }
        }
    }
    for(int n=0,z=0;n<Qnames.size();++n,++z){
        for (std::map<std::string, std::string>::iterator it4 = sequences.begin();it4 != sequences.end(); ++it4) {
            if (it4->first == Qnames[n]) {
                sequences[Qnames[n]] = _alignment_qs[z];
            }
        }
    }
    for (std::map<std::string,std::string>::iterator it6=sequences.begin();it6!=sequences.end();++it6) {
        std::cout << it6->first <<"  "<< it6->second << std::endl;
    }
}

