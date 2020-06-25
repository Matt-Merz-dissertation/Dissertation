/*
** This is a program written by Matthias Merzenich to search for spherical
** pictures over certain relative group presentations.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define P_MAXVERTS 0
#define P_MAXFACES 1
#define P_MAXNEWVERTS 2
#define P_MAXHALFEDGES 3
#define P_GROUPORDER 4
#define P_DEGREE 5
#define P_SYMMETRYTYPE 6
#define P_NUMPICS 7
#define P_OUTFREQ 8
#define P_OUTLARGEST 9

#define NUM_PARAMS 10

#define MAXDEGREE 10

//#define TABLE_DEBUG

const char *wordType;
const char *startWord;

//int letterValue[10] = {-3,1,0,0,-1,3};

int letterValue[MAXDEGREE];

int groupOrder;

int symmetryType;

int degree;

int params[NUM_PARAMS];

int* extraFace;         //This tracks whether an extra face was added at a particular step.
int* onlyOneEdge;       //This tracks whether only one edge was seen at this step.

int totalNodes;
int totalFreeEdges;

int numPics = 0;

int** canFollow;

/* ========================== */
/* Data structure definitions */
/* ========================== */

// Vertices are added to the picture as part of a face.  The orientation of a
// vertex represents the corner of the vertex that is included in the new face.
// 
// seenCorner[] gives the first step at which a corner became part of a face.
// seenCorner[0] represents the first unused corner when going around clockwise.
// 
// seenEdge[] gives the first step at which an edge became part of a face.
// seenEdge[0] represents the first unused edge when going around clockwise.
// 
// Note that one corner and two edges are already used when a new face is added.

typedef struct vertex {
   int* seenCorner;
   int* seenEdge;
   uint8_t orient;
} vertex;

// This can be used to store the location of any single edge in the picture.
// (note that we do not keep a list of all edges)
// The root of an edge is the first vertex added that includes that edge.
// The edgeID identifies which edge of the root we are looking at.
typedef struct edge {
   vertex* root;
   uint8_t edgeID;
} edge;

// This can be used to store the location of any single corner in the picture.
// (note that we do not keep a list of all corners)
// The root of a corner is the vertex associated to the corner.
// The cornerID identifies which corner of the root we are looking at.
// (consider seenCorner[cornerID])
typedef struct corner {
   vertex* root;
   uint8_t cornerID;
} corner;

vertex* verts;          //Stack of vertices in the current picture
vertex* top;            //top of vertex stack
vertex* bottom;         //bottom of vertex stack

int* corners;
int* topCorner;

int* edges;
int* topEdge;

vertex* getVertexFromEdge(int* theEdge){
   return bottom + (int)((theEdge - edges) / (degree - 2));
}

vertex* getVertexFromCorner(int* theCorner){
   return bottom + (int)((theCorner - corners) / (degree - 1));
}

/* ======================= */
/* Vertex stack operations */
/* ======================= */

unsigned int totalVerts;      //Total number of vertices in current picture

void push(uint8_t orientation){
   int i;
   ++top;
   ++totalVerts;
   totalFreeEdges += degree - 2;
   top->orient = orientation;
   top->seenCorner = topCorner + 1;
   top->seenEdge = topEdge + 1;
   topCorner += degree - 1;
   topEdge += degree - 2;
   for(i = 0; i < degree - 2; ++i){
      top->seenCorner[i] = 0;
      top->seenEdge[i] = 0;
   }
   top->seenCorner[degree - 2] = 0;    //need to initialize one more corner than edge
}

vertex* pop(){
   --top;
   --totalNodes;
   totalFreeEdges -= degree - 2;
   topCorner -= degree - 1;
   topEdge -= degree - 2;
   return (top + 1);
}

// Orientation is simply a choice of corner.  The possible orientations are
// numbers from 0 to degree - 1.  The first degree orientations correspond to
// ..., d, c, b, a.  The remaining orientations go in opposite order
// A, B, C, D, ...

uint8_t orientFromLetter(char theLetter){
   if(theLetter < 'a') return (degree + theLetter - 'A');
   else return (uint8_t)(degree + 'a' - theLetter - 1);
}

char letterFromOrient(uint8_t theOrient){
   if(theOrient >= degree) return (theOrient + 'A' - degree);
   else return (char)(degree + 'a' - theOrient - 1);
}

char getCornerLetter(vertex* theVertex, uint8_t theCorner){
   int i = (theVertex->orient + theCorner + 1) % degree;
   if(theVertex->orient >= degree) i += degree;
   return letterFromOrient(i);
}

int getCornerOrient(vertex* theVertex, uint8_t theCorner){
   int i = (theVertex->orient + theCorner + 1) % degree;
   if(theVertex->orient >= degree) i += degree;
   return i;
}


/* ===================================== */
/* Data-printing functions (for testing) */
/* ===================================== */


// Print entire vertex list
void printVerts(){
   vertex* i;
   int j;
   long long int c = 0;
   printf(" Vertex ID | orient | letter | seenCorner? | seenEdge? \n");
   printf("-----------+--------+--------+-------------+-----------\n");
   for(i = bottom; i <= top; ++i){
      ++c;
      printf(" %9lli |     %2d |    %c   | ",c,i->orient,letterFromOrient(i->orient));
      for(j = 0; j < degree - 1; ++j){
         //only tells if corner has been seen, not when corner was seen:
         putchar((i->seenCorner[j] ? '1' : '0'));
      }
      printf("%*s | ",12 - degree," ");
      for(j = 0; j < degree - 2; ++j){
         //only tells if edge has been seen, not when edge was seen:
         putchar((i->seenEdge[j] ? '1' : '0'));
      }
      printf("\n");
   }
}

int getVertexNumber(vertex* theVertex){
   int c = 0;
   vertex* i = bottom;
   while(i < theVertex){
      ++i;
      ++c;
   }
   return c;
}

int theWordType;

int* coeffInit;
int* coeffTerm;

int** nextEdgeValue;

void buildFollow(){
   int i,j;
   canFollow  = (int **)malloc(sizeof(int *) * degree * 2);
   canFollow[0] = (int *)malloc(sizeof(int) * degree * degree * 4);
   for(i = 0; i < 2*degree; ++i) canFollow[i] = (*canFollow + 2 * degree * i);
   
   nextEdgeValue  = (int **)malloc(sizeof(int *) * degree * 2);
   nextEdgeValue[0] = (int *)malloc(sizeof(int) * degree * degree * 4);
   for(i = 0; i < 2*degree; ++i) nextEdgeValue[i] = (*nextEdgeValue + 2 * degree * i);
   
   coeffInit = malloc(2*degree*sizeof(*coeffInit));
   coeffTerm = malloc(2*degree*sizeof(*coeffTerm));
   
   // The following is for words of the form xaxbxcxd...
   for(i = 1; i < degree; ++i){
      coeffTerm[i] = (((theWordType >> (i - 1)) & 1) ? 0 : 1);
   }
   for(i = 0; i < degree; ++i){
      coeffInit[i] = (((theWordType >> i) & 1) ? 1 : 0);
   }
   coeffTerm[0] = (((theWordType >> (degree - 1)) & 1) ? 0 : 1);
   
   for(i = 0; i < degree - 1; ++i){
      coeffInit[i + degree] = ((theWordType >> ((degree - 2 - i)) & 1) ? 0 : 1);
   }
   for(i = 0; i < degree; ++i){
      coeffTerm[i + degree] = ((theWordType >> ((degree - 1 - i)) & 1) ? 1 : 0);
   }
   coeffInit[2*degree - 1] = (((theWordType >> (degree - 1)) & 1) ? 0 : 1);
   
   for(i = 0; i < 2*degree; ++i){
      for(j = 0; j < 2*degree; ++j){
         if(coeffTerm[i] == coeffInit[j] && i != 2*degree - 1 - j) canFollow[i][j] = 1;
         else canFollow[i][j] = 0;
         
#ifdef TABLE_DEBUG
         printf("%c %c %d   ",letterFromOrient(i),letterFromOrient(j),canFollow[i][j]);
#endif
         
      }
#ifdef TABLE_DEBUG
      printf("\n");
#endif
   }
}

typedef struct starEdge starEdge;

struct starEdge {
   starEdge** nextEdge;
   uint8_t orientation;
};

starEdge* starGraph;
starEdge** starEdges;

void buildStarGraph(){
   
   int i,j,c;
   
   starGraph = malloc(2 * degree * sizeof(*starGraph));
   starEdges = malloc(2 * degree * (degree - 1) * sizeof(*starEdges));
   
   for(i = 0; i < 2*degree; ++i){
      starGraph[i].orientation = (uint8_t)i;
      starGraph[i].nextEdge = starEdges + i * (degree - 1);
   }
   
#ifdef TABLE_DEBUG
   printf("\n\n");
#endif
   
   for(i = 0; i < 2*degree; ++i){
      for(j = 0; j < 2*degree; ++j){
         
#ifdef TABLE_DEBUG
         printf("%c %c %d   ",letterFromOrient(i),letterFromOrient(j),canFollow[i][j]);
#endif
         
      }
#ifdef TABLE_DEBUG
      printf("\n");
#endif
   }
   for(i = 0; i < 2*degree; ++i){
      c = 0;
      for(j = 0; j < 2*degree; ++j){
         nextEdgeValue[i][j] = -1;
         if(canFollow[i][j]){
            starGraph[i].nextEdge[c] = starGraph + j;
            nextEdgeValue[i][j] = c;
            
#ifdef TABLE_DEBUG
            printf("%d  %d  %d\n",i,j,c);
#endif
            ++c;
         }
      }
   }
   
#ifdef TABLE_DEBUG
   for(i = 0; i < 2*degree; ++i){
      printf("%c:  ",letterFromOrient(starGraph[i].orientation));
      for(j = 0; j < degree - 1; ++j){
         printf("%c ",letterFromOrient((starGraph[i].nextEdge[j])->orientation));
         fflush(stdout);
      }
      printf("\n");
      fflush(stdout);
   }
#endif
   
}

long long int myPow(int base, int exponent){
   int i;
   long long int result = 1;
   for(i = 0; i < exponent; ++i) result *= base;
   return result;
}

uint32_t* newPaths;
uint8_t* newPathLengths;
uint32_t* newPathInd;

// Test if a path in the star graph is a valid closed loop representing a
// freely-reduced word representing the identity in G.
// initCornerSum does NOT include firstLetter or lastLetter
int goodPath(int firstLetter, int lastLetter, int initCornerSum, uint32_t thePath, int theLength){
   int i,j,c;
   int sum = initCornerSum;
   
   int totalZeros = 0;
   int totalNonZero = 0;
   starEdge* y = starGraph + lastLetter;    //y is our starting point in the star graph
   
   starEdge* x = y;
   
   for(i = 0; i < theLength; ++i){
      c = ((long long int)thePath / myPow(degree - 1,i)) % (degree - 1);
      x = x->nextEdge[c];
      
      if(letterValue[x->orientation] == 0) ++totalZeros;
      else ++totalNonZero;
      
      sum += letterValue[x->orientation];
   }
   
   if(totalNonZero == 0 && totalZeros > 1) return 0;
   if(canFollow[x->orientation][firstLetter] && sum % groupOrder == 0){
      x = y;
#ifdef TABLE_DEBUG
      printf("%d  %c ",initCornerSum,letterFromOrient(x->orientation));
#endif
      for(i = 0; i < theLength; ++i){
         c = (thePath / myPow(degree - 1,i)) % (degree - 1);
         x = x->nextEdge[c];
#ifdef TABLE_DEBUG
         printf("%c",letterFromOrient(x->orientation));
#endif
      }
#ifdef TABLE_DEBUG
      printf(" %c",letterFromOrient(firstLetter));
      printf("\n");
#endif
      
      return 1;
   }
   else return 0;
}

//paths in newPaths are stored "backwards".  That is, the least significant "bit" represents the first step in the star graph

void buildTables2(){
   long long int c = 0;
   long long int i,j;
   int cornerSum;
   int length,path;
   int startSymb,endSymb;
   int maxInputs = (2*degree) * (2*degree) * groupOrder;
   long long int maxPaths = myPow(degree - 1, params[P_MAXNEWVERTS]) * maxInputs;
   
   uint32_t* newPathsTemp;
   uint8_t* newPathLengthsTemp;
   
   newPathsTemp = malloc(maxPaths * sizeof(*newPathsTemp));
   newPathLengthsTemp = malloc(maxPaths * sizeof(*newPathLengthsTemp));
   
   uint32_t theIndex;
   
   uint32_t* pathCount;  //pathCount[i] == number of valid paths for input i;
   
   pathCount = malloc(maxInputs * sizeof(*pathCount));
   newPathInd = malloc((maxInputs + 1) * sizeof(*newPathInd));
   for(i = 0; i < maxInputs; ++i){
      pathCount[i] = 0;
      newPathInd[i] = 0;
   }
   
   newPathInd[maxInputs] = 0;
   
   i = -1;
   
#ifdef TABLE_DEBUG
   printf("corner values:\n");
   for(j = 0; j < 2*degree; ++j){
      printf("%c: %d\n",letterFromOrient(j), letterValue[j]);
   }
#endif
   
   for(cornerSum = 0; cornerSum < groupOrder; ++cornerSum){
      for(startSymb = 0; startSymb < 2*degree; ++startSymb){
         for(endSymb = 0; endSymb < 2*degree; ++endSymb){
#ifdef TABLE_DEBUG
            printf("%d  %c %c %d\n",theIndex,letterFromOrient(startSymb),letterFromOrient(endSymb),cornerSum);
#endif
            theIndex = cornerSum*(2*degree)*(2*degree) + startSymb*(2*degree) + endSymb;
            for(length = 0; length <= params[P_MAXNEWVERTS]; ++length){
               for(path = 0; path < myPow(degree - 1, length); ++path){
                  if(goodPath(startSymb,endSymb,cornerSum,path,length)){
                     ++i;
                     ++pathCount[theIndex];
                     newPathsTemp[i] = path;
                     newPathLengthsTemp[i] = length;
                  }
               }
            }
         }
      }
   }
   
   newPaths = malloc(i * sizeof(*newPaths));
   newPathLengths = malloc(i * sizeof(*newPathLengths));
   
   for(j = 0; j <= i; ++j){
      newPaths[j] = newPathsTemp[j];
      newPathLengths[j] = newPathLengthsTemp[j];
   }
   newPathInd[0] = 0;
   for(i = 1; i <= maxInputs; ++i){
      newPathInd[i] = newPathInd[i - 1] + pathCount[i - 1];
   }
   free(pathCount);
   free(newPathsTemp);
   free(newPathLengthsTemp);
}

//firstLetter does not count as part of theLength
void printWordFromPath(int firstLetter, uint32_t thePath, int theLength){
   int i,c;
   printf("%c",letterFromOrient(firstLetter));
   starEdge* y = starGraph + firstLetter;    //y is our starting point in the star graph
   
   for(i = 0; i < theLength; ++i){
      c = ((long long int)thePath / myPow(degree - 1,i)) % (degree - 1);
      y = y->nextEdge[c];
      printf("%c",letterFromOrient(y->orientation));
   }
   printf("\n");
   
}

void printWord(uint32_t startPath, uint32_t endPath, int startLength, int endLength){
   
   int i;
   
   int firstLetter = startPath / myPow(degree - 1,startLength - 1);
   
   printf("%c",letterFromOrient(firstLetter));
   
   starEdge* y = starGraph + firstLetter;
   
   for(i = 0; i < startLength - 1; ++i){
      y = y->nextEdge[(startPath / myPow(degree - 1,startLength - 2 - i)) % (degree - 1)];
      printf("%c",letterFromOrient(y->orientation));
   }
   
   printf(" ");
   
   for(i = 0; i < endLength; ++i){
      y = y->nextEdge[(endPath / myPow(degree - 1,i)) % (degree - 1)];
      printf("%c",letterFromOrient(y->orientation));
   }
   
   printf("\n");
   
}

typedef struct picFace picFace;

struct picFace {
   uint32_t path;
   uint8_t startLetter;
   uint8_t length;
};

picFace* currentPic;

/* === getFirstUnseenEdge() & getLastUnseenEdge() ===
// These functions find the first and last edges that have not yet been used.
// It uses the global list of vertices and edges to determine what's been seen.
// 
// Input: None (works on global variables)
// 
// Output: returns pointer to element in edge list
//
// Does not change global variable values
*/

int* getFirstUnseenEdge(){
   int* theEdge = edges;
   while(*theEdge != 0 && theEdge < topEdge) ++theEdge;
   return theEdge;
}

int* getLastUnseenEdge(){
   int* theEdge = topEdge;
   while(*theEdge != 0 && theEdge > edges) --theEdge;
   return theEdge;
}

int* getNextUnseenEdge(int* theEdge){
   theEdge++;
   while(*theEdge != 0 && theEdge < topEdge) ++theEdge;
   return theEdge;
}

int* getPreviousUnseenEdge(int* theEdge){
   theEdge--;
   while(*theEdge != 0 && theEdge > edges) --theEdge;
   return theEdge;
}

int* getFirstUnseenCorner(){
   int* theCorner = corners;
   while(*theCorner != 0 && theCorner < topCorner) ++theCorner;
   return theCorner;
}

// The following function is for a sanity check.
// It counts the current number of unseen (free) edges.
int getNumberOfUnseenEdges(){
   int theCount = 0;
   int* theEdge = edges;
   while(theEdge <= topEdge){
      if(*theEdge == 0) ++theCount;
      ++theEdge;
   }
   return theCount;
}

int getOrientFromCorner(int* theCorner){
   int c = (theCorner - corners) % (degree - 1);
   return getCornerOrient(getVertexFromCorner(theCorner), c);
}

//want to "return" the length along with the path
int currentStartPathLength;

int step;
int* numAddedVerts;

/* === getCurrentStartPath() ===
// This function finds the partial word to be used in the next added face.
// It uses the global list of vertices and edges to determine what's been seen.
// 
// Input: None (works on global variables)
// 
// Output: Integer representing the first and last corners and the corner sum
//
// Changes global variable values
//   Edges and corners in edge and corner lists are set to the current step.
//   When the program backs up at a particular step, it will set these
//   values back to 0.
*/

uint32_t getCurrentStartPath(){
   
   int cornerSum = 0;
   
   int* firstEdge = getFirstUnseenEdge();
   int* lastEdge = getLastUnseenEdge();
   
   if(firstEdge == lastEdge) onlyOneEdge[step] = 1;   //This is needed.
   else onlyOneEdge[step] = 0;                        //This might be redundant.
   
   int* firstCorner;
   int* lastCorner;
   
   *firstEdge = step;
   *lastEdge = step;
   
   vertex* firstVertex = getVertexFromEdge(firstEdge);
   // The location mod (degree - 2) of firstEdge in the edge list tells us which
   // edge of the vertex we have
   int firstEdgeNum = (firstEdge - edges) % (degree - 2);
   firstCorner = (firstVertex->seenCorner) + firstEdgeNum;
   
   vertex* lastVertex = getVertexFromEdge(lastEdge);
   int lastEdgeNum = (lastEdge - edges) % (degree - 2);
   lastCorner = (lastVertex->seenCorner) + lastEdgeNum + 1;
   int* x = firstCorner;
   int* firstUnseenCorner;
   int* lastUnseenCorner;
   int foundFirstYet = 0;
   int i = 0;
   
   while(x >= corners){
      if(*x == 0){
         if(!foundFirstYet){
            foundFirstYet = 1;
            firstUnseenCorner = x;
         }
         lastUnseenCorner = x;
         cornerSum += letterValue[getOrientFromCorner(x)];
         *x = step;
         ++i;
      }
      --x;
   }
   x = topCorner;
   
   while(x >= lastCorner){
      if(*x == 0){
         if(!foundFirstYet){
            foundFirstYet = 1;
            firstUnseenCorner = x;    //This line might be unnecessary
         }
         lastUnseenCorner = x;
         
         cornerSum += letterValue[getOrientFromCorner(x)];
         *x = step;
         ++i;
      }
      --x;
   }
   cornerSum = (cornerSum % groupOrder);
   if(cornerSum < 0) cornerSum += groupOrder;
   return cornerSum*(2*degree)*(2*degree) + getOrientFromCorner(firstUnseenCorner)*(2*degree) + getOrientFromCorner(lastUnseenCorner);
}


int canComplete(){
   int cornerSum = 0;
   
   // Again, first edge is later in list, which is confusing:
   int* lastEdge = getFirstUnseenEdge();
   
   // but this matches the process in getCurrentStartPath():
   int* firstEdge = getNextUnseenEdge(lastEdge);  
   int* firstCorner;
   int* lastCorner;
   
   *firstEdge = step;
   *lastEdge = step;
   
   vertex* firstVertex = getVertexFromEdge(firstEdge);
   // The location mod (degree - 2) of firstEdge in the edge list tells us which
   // edge of the vertex we have
   int firstEdgeNum = (firstEdge - edges) % (degree - 2);
   firstCorner = (firstVertex->seenCorner) + firstEdgeNum;
   
   vertex* lastVertex = getVertexFromEdge(lastEdge);
   int lastEdgeNum = (lastEdge - edges) % (degree - 2);
   lastCorner = (lastVertex->seenCorner) + lastEdgeNum + 1;
   int* x = firstCorner;
   int* firstUnseenCorner;
   int* lastUnseenCorner;
   int foundFirstYet = 0;
   int i = 0;
   
   while(x >= lastCorner){
      if(*x == 0){
         if(!foundFirstYet){
            foundFirstYet = 1;
            firstUnseenCorner = x;
         }
         lastUnseenCorner = x;
         
         cornerSum += letterValue[getOrientFromCorner(x)];
         
         *x = step;
      }
      x--;
   }
   
   cornerSum = (cornerSum % groupOrder);
   if(cornerSum < 0) cornerSum += groupOrder;
   
   uint32_t theStartPath = cornerSum*(2*degree)*(2*degree) + getOrientFromCorner(firstUnseenCorner)*(2*degree) + getOrientFromCorner(lastUnseenCorner);
   
   uint32_t theStartIndex = newPathInd[theStartPath + 1];
   int numberOfPaths = newPathInd[theStartPath + 1] - newPathInd[theStartPath];
   
   int theNewLength = newPathLengths[theStartIndex - numberOfPaths];
   
   if(cornerSum == 0 && numberOfPaths > 0 && theNewLength == 0) return 1;
   else{
      
      *firstEdge = 0;         //"back up" over the improperly added face
      *lastEdge = 0;
      
      int* x = corners;
      
      while(x <= topCorner){
         if(*x == step){
            *x = 0;
         }
         ++x;
      }
      
      return 0;
   }
      
   
   
}

int addNewVerts(uint32_t startPath, uint32_t thePath, int theLength){
   
   int i;
   
   totalFreeEdges -= 2;
   
   int firstLetter = startPath % (2*degree);
   int lastLetter = (startPath / (2*degree)) % (2*degree);
   
   starEdge *x = starGraph + firstLetter;
   starEdge *y = starGraph + firstLetter;
   
   if(onlyOneEdge[step] == 1){
      if(theLength < 2) return 0;
      
      x = x->nextEdge[(thePath / myPow(degree - 1,i)) % (degree - 1)];
      int requiredFinalOrient = (x->orientation) + 1;
      if(requiredFinalOrient%degree == 0) requiredFinalOrient -= degree;
      
      for(i = 1; i < theLength; i++){
         x = x->nextEdge[(thePath / myPow(degree - 1,i)) % (degree - 1)];
      }
      
      if(x->orientation == requiredFinalOrient) theLength--;
      else return 0;
      
   }
   for(i = 0; i < theLength; i++){
      y = y->nextEdge[(thePath / myPow(degree - 1,i)) % (degree - 1)];
      push(y->orientation);
   }
   
   if(onlyOneEdge[step] == 1){
      *getFirstUnseenCorner() = -1;
      *getFirstUnseenEdge() = -1;
   }
   
   numAddedVerts[step] = theLength;
   return 1;
   
}

void unseeVertsAndEdges(){
   totalFreeEdges += 2;
   int* x = corners;
   
   while(x <= topCorner){
      if(*x == step){
         *x = 0;
      }
      ++x;
   }
   
   x = edges;
   
   while(x <= topEdge){
      if(*x == step){
         *x = 0;
      }
      ++x;
   }
   
}

void backUp(){
   unseeVertsAndEdges(step);
   int i;
   for(i = 0; i < numAddedVerts[step - 1]; ++i) pop();
}

int printPictureFace(int theStep){
   
   printf("%04d  ",theStep+1);
   
   int* firstEdge = edges;
   while(*firstEdge != theStep && firstEdge < topEdge) ++firstEdge;
   
   int* lastEdge = topEdge;
   while(*lastEdge != theStep && lastEdge > edges) --lastEdge;
   
   int* firstCorner;
   int* lastCorner;
   
   int wordCounter = 0;
   vertex* firstVertex = getVertexFromEdge(firstEdge);
   int firstEdgeNum = (firstEdge - edges) % (degree - 2);
   firstCorner = (firstVertex->seenCorner) + firstEdgeNum;
   
   vertex* lastVertex = getVertexFromEdge(lastEdge);
   int lastEdgeNum = (lastEdge - edges) % (degree - 2);
   lastCorner = (lastVertex->seenCorner) + lastEdgeNum + 1;
   
   int* x = firstCorner;
   int lastCornerOrient;
   
   while(x >= corners){
      if(*x == theStep){
         lastCornerOrient = getOrientFromCorner(x);
         printf("%c",letterFromOrient(getOrientFromCorner(x)));
         wordCounter++;
      }
      --x;
   }
   x = topCorner;
   
   while(x >= lastCorner){
      if(*x == theStep){
         lastCornerOrient = getOrientFromCorner(x);
         printf("%c",letterFromOrient(getOrientFromCorner(x)));
         wordCounter++;
      }
      --x;
   }
   
   printf(" ");
   
   starEdge* y = starGraph + lastCornerOrient;
   int i;
   
   for(i = 0; i < currentPic[theStep].length; ++i){
      y = y->nextEdge[(currentPic[theStep].path / myPow(degree - 1,i)) % (degree - 1)];
      printf("%c",letterFromOrient(y->orientation));
      wordCounter++;
   }
   
   printf("\n");
   fflush(stdout);
   
   return wordCounter;
}


int printPictureFaceHop(int theStep){
   
   printf("%04d  ",theStep+1);
   
   int wordCounter = 0;
   
   int* x = topCorner;
   
   while(x >= corners){
      if(*x == theStep){
         printf("%c",letterFromOrient(getOrientFromCorner(x)));
         wordCounter++;
      }
      x--;
   }
   printf(" (hop)\n");
   fflush(stdout);
   
   return wordCounter;
}

printPictureFaceHop2(int theStep){
   
   printf("%04d  ",theStep+1);
   
   theStep = 0;
   
   int* x = topCorner;
   
   while(x >= corners){
      if(*x == theStep){
         printf("%c",letterFromOrient(getOrientFromCorner(x)));
      }
      x--;
   }
   fflush(stdout);
}




void printGraph(int firstFaceSize){
   
   int i,j;
   
   printf("\n\nEL <- matrix(c(");
   for(i = 0; i < firstFaceSize - 1; ++i){
      printf("%d,%d, ",i+1,i+2);
   }
   printf("%d,1, ",i+1);
   int *x;
   int *y;
   int b = firstFaceSize;
   
   for(i = 1; i < step; ++i){
      
      x = edges;
      while(*x != i) ++x;
      y = x + 1;
      while(*y != i) ++y;
      
      if(numAddedVerts[i] == 0){
         printf("%d,%d, ",getVertexFromEdge(x) - bottom + 1,getVertexFromEdge(y) - bottom + 1);
         continue;
      }
      
      printf("%d,%d, ",getVertexFromEdge(y) - bottom + 1,b + 1);
      for(j = 1; j < numAddedVerts[i]; ++j){
         printf("%d,%d, ",b + 1,b + 2);
         ++b;
      }
      printf("%d,%d, ",b + 1,getVertexFromEdge(x) - bottom + 1);
      
      ++b;
      
   }
   
   printf("),ncol=2,byrow=TRUE)\n\n");
   
   printf("c(");
   for(i = 1; i < step; ++i) printf("%d,",numAddedVerts[i]);
   printf(")\n\n");
   
   fflush(stdout);
}


//The stack structure for the main search
uint32_t* startPath;
uint32_t* startIndex;
int* startPathLength;
int* pathsRemaining;

void printPicture(int initialNumberOfVerts,int symmetryType,int isCompletePicture){
   
   int i;
   
   int maxWordLength = 0;
   int wordLength;
   
   vertex* x = bottom;
   
   printf("%d-fold symmetry\n\n",symmetryType);
   printf("0001  ");
   
   
   for(i = 0; i < initialNumberOfVerts; i++){
      printf("%c",letterFromOrient(x->orient));
      x++;
   }
   printf(" (x%d)",symmetryType);
   printf("\n");
   
   for(i = 1; i < step; i++){
      if(extraFace[i] == 1) wordLength = printPictureFaceHop(i);
      else wordLength = printPictureFace(i);
      if(wordLength > maxWordLength) maxWordLength = wordLength;
   }
   if(isCompletePicture){
      printPictureFaceHop2(step);
      printf(" (x%d)\n",symmetryType);
      printf("\n");
      printf("Max Word Length: %d",maxWordLength);
   }
}

void printSymmGraph(int firstFaceSize, int n){
   
   int i,j,k;
   
   printf("\n\nEL <- matrix(c(");
   for(k = 0; k < n; k++){
      for(i = 0; i < firstFaceSize - 1; i++){
         printf("%d,%d, ", n*i + k, n*(i+1) + k);
      }
      printf("%d,%d, ", n*i + k, (k+1) % n);
   }
   int *x;
   int *y;
   int b = firstFaceSize;
   
   
   for(i = 1; i < step; i++){
      
      x = edges;
      while(*x != i) x++;
      if(onlyOneEdge[i] == 1) y = x;
      else{
         y = x + 1;
         while(*y != i) y++;
      }
      
      if(extraFace[i] == 1){
         for(k = 0; k < n; k++){
            printf("%d,%d, ",n*(getVertexFromEdge(x) - bottom) + k, n*(getVertexFromEdge(y) - bottom) + k);
         }
         continue;
      }
      if(numAddedVerts[i] == 0){
         for(k = 0; k < n; k++){
            printf("%d,%d, ", n*(getVertexFromEdge(y) - bottom) + k,n*(getVertexFromEdge(x) - bottom) + ((k+1)%n));
         }
         continue;
      }
      
      for(k = 0; k < n; k++){
         printf("%d,%d, ", n*(getVertexFromEdge(y) - bottom) + k , n*b + k);
      }
      
      int temp = b;
      
      for(j = 1; j < numAddedVerts[i]; ++j){
         for(k = 0; k < n; k++){
            printf("%d,%d, ",n*b + k, n*(b+1) + k);
         }
         ++b;
      }
      
      if(onlyOneEdge[i] == 1){
         for(k = 0; k < n; k++){
            printf("%d,%d, ",n*b + k, n*temp + ((k+1)%n));
         }
      }
      else{
         for(k = 0; k < n; k++){
            printf("%d,%d, ",n*b + k, n*(getVertexFromEdge(x) - bottom) + ((k+1)%n));
         }
      }
      
      ++b;
   }
   
   printf(")+1,ncol=2,byrow=TRUE)\n\n");
   
   printf("c(");
   for(i = 1; i < step; ++i) printf("%d,",n*numAddedVerts[i]);
   printf(")\n\n");
   
   fflush(stdout);
}

void echoParams(){
   
   int i;
   
   printf("Max vertices: %d\n", params[P_MAXVERTS]);
   printf("Max faces: %d\n", params[P_MAXFACES]);
   printf("Max new vertice: %d\n", params[P_MAXNEWVERTS]);
   printf("Max half edges: %d\n", params[P_MAXHALFEDGES]);
   printf("Num pictures: %d\n", params[P_NUMPICS]);
   printf("Output frequency: %d\n", params[P_OUTFREQ]);
   if(params[P_OUTLARGEST]) printf("Print largest partial picture\n");
   printf("Coefficient group order: %d\n", groupOrder);
   printf("Degree: %d\n", degree);
   
   printf("theWordType: ");
   for(i = degree - 1; i >= 0; i--){
      putchar((theWordType & (1 << i)) ? 'x' : 'X');
   }
   printf("\n");
   
   
   printf("Exponents: ");
   for(i = 0; i < degree; i++){
      printf("%d ",letterValue[i]);
   }
   printf("\n");
   
   printf("Start Word: %s\n",startWord);
   printf("Symmetry Type: %d\n", symmetryType);
   
   printf("\n");
}

int gcd(int a, int b) {
   if (a > b) return gcd(b,a);
   else if (a == 0) return b;
   else return gcd(b-a,a);
}

int main(int argc, char *argv[]){
   
   long long int i;
   
   for (i = 0; i < argc; i++)
      printf(" %s", argv[i]) ;
   printf("\n\n");
   
   params[P_MAXVERTS] = 1 << 20;
   params[P_MAXFACES] = 1 << 20;
   params[P_MAXNEWVERTS] = 0;
   params[P_MAXHALFEDGES] = 0;
   params[P_NUMPICS] = 1;
   params[P_OUTFREQ] = 26;
   params[P_OUTLARGEST] = 0;
   
   totalFreeEdges = 0;
   
   degree = 3;
   theWordType = 0b110;
   
   
   int s;
   
   int j;
   
   
   
   
   if (argc > 1) while(--argc > 0) {
      if ((*++argv)[0] == '-'){
        switch ((*argv)[1]) {
           case 'v': case 'V':
              --argc;
              sscanf(*++argv, "%d", &params[P_MAXVERTS]);
              break;
           case 'n': case 'N':
              --argc;
              sscanf(*++argv, "%d", &params[P_MAXNEWVERTS]);
              break;
           case 'h': case 'H':
              --argc;
              sscanf(*++argv, "%d", &params[P_MAXHALFEDGES]);
              break;
           case 'f': case 'F':
              --argc;
              sscanf(*++argv, "%d", &params[P_MAXFACES]);
              break;
           case 'p': case 'P':
              --argc;
              sscanf(*++argv, "%d", &params[P_NUMPICS]);
              break;
           case 'q': case 'Q':
              --argc;
              sscanf(*++argv, "%d", &params[P_OUTFREQ]);
              break;
           case 'o': case 'O':
              --argc;
              sscanf(*++argv, "%d", &params[P_GROUPORDER]);
              groupOrder = params[P_GROUPORDER];
              break;
           case 'w': case 'W':
              --argc;
              startWord = *++argv;
              break;
           case 'r': case 'R':
              --argc;
              wordType = *++argv;
              params[P_DEGREE] = strlen(wordType);
              degree = params[P_DEGREE];
              for(i = 0; i < degree; i++){
                 theWordType |= ((wordType[i] == 'x') ? 1 : 0) << (degree - i - 1);
              }
              break;
           case 'e': case 'E':
              for(i = 0; i < degree; i++){
                 --argc;
                 sscanf(*++argv, "%d", &letterValue[degree - i - 1]);
              }
              for(i = 0; i < degree; i++){
                 letterValue[degree + i] = -letterValue[degree - i - 1];
              }
              break;
           case 'l': case 'L':
              params[P_OUTLARGEST] = 1;
              break;
           default:
              printf("Unrecognized option %s\n", *argv);
              exit(1);
        }
      }
   }
   else{
      printf("This is a program for finding reduced spherical pictures over certain\n");
      printf("relative presentations. Input options are as follows:\n");
      printf("\n");
      printf("-r RR  Searches for pictures over a relative presentation with relator\n");
      printf("       shape RR specified as a string consisting of \'x\' and \'X\', where \n");
      printf("       \'X\' represents x^-1. For example, the relator xaxbx^(-1)c is\n");
      printf("       specified by xxX.\n");
      printf("-o NN  specifies the order of the coefficient group to NN\n");
      printf("-e LL  LL is a space-separated list specifying the exponent values of\n");
      printf("       the coefficient elements in the relator\n");
      printf("\n");
      printf("-v NN  Searches for pictures with at most NN vertices (default: 2^20)\n");
      printf("-f NN  Searches for pictures with at most NN faces (default: 2^20)\n");
      printf("\n");
      printf("-n NN  Add no more than NN vertices when creating a new face\n");
      printf("-h NN  Limit the number of half edges at any step to NN\n");
      printf("\n");
      printf("-p NN  Output up to NN complete spherical pictures (default: 1)\n");
      printf("-q NN  Output the current search state every 2^NN steps (default: 26)\n");
      printf("-l     Output the current search state when a new largest picture occurs\n");
      exit(0);
   }

   
   symmetryType = 0;
   for(i = 0; i < strlen(startWord); i++){
      symmetryType += letterValue[orientFromLetter(startWord[i])];
   }
   symmetryType = symmetryType % groupOrder;
   if(symmetryType < 0) symmetryType += groupOrder;
   symmetryType = groupOrder / gcd(groupOrder,symmetryType);
   params[P_SYMMETRYTYPE] = symmetryType;
   
   echoParams();
   //return 0;
   
   verts = malloc(params[P_MAXVERTS]*sizeof(*verts));
   
   corners = malloc(params[P_MAXVERTS]*(degree - 1)*sizeof(*corners));
   edges = malloc(params[P_MAXVERTS]*(degree - 2)*sizeof(*edges));
   
   
   topCorner = corners - 1;
   topEdge = edges - 1;
   
   top = verts - 1;
   bottom = verts;
   
   numAddedVerts = malloc((params[P_MAXFACES] + 5)*sizeof(*numAddedVerts));
   
   buildFollow();
   buildStarGraph();
   
   for(i = 0; i < strlen(startWord) - 1; i++){
      if(!canFollow[orientFromLetter(startWord[i])][orientFromLetter(startWord[i+1])]){
         printf("Error: start word is invalid.\n");
         exit(0);
      }
   }
   
   if(!canFollow[orientFromLetter(startWord[strlen(startWord) - 1])][orientFromLetter(startWord[0])]){
      printf("Error: start word is invalid.\n");
      exit(0);
   }
   
   
   buildTables2();
   
   
   currentPic = malloc((params[P_MAXFACES] + 5)*sizeof(*currentPic));
   
   startPath = malloc((params[P_MAXFACES] + 5)*sizeof(*startPath));
   startPathLength = malloc((params[P_MAXFACES] + 5)*sizeof(*startPathLength));
   startIndex = malloc((params[P_MAXFACES] + 5)*sizeof(*startIndex));
   pathsRemaining = malloc((params[P_MAXFACES] + 5)*sizeof(*pathsRemaining));
   
   extraFace = calloc(params[P_MAXFACES] + 5,sizeof(*extraFace));
   onlyOneEdge = calloc(params[P_MAXFACES] + 5,sizeof(*onlyOneEdge));
   
   search();
   
   return 0;
   
}

void getTheStartPath(int step){
   startPath[step] = getCurrentStartPath();
   startIndex[step] = newPathInd[startPath[step] + 1];
   pathsRemaining[step] = newPathInd[startPath[step]+1] - newPathInd[startPath[step]] + 1;
   
   currentPic[step].startLetter = startIndex[step] % (2*degree);
}


int search(){
   
   int j;
   long long int i, c;
   
   for(i = 0; i < strlen(startWord); i++){
      push(orientFromLetter(startWord[i]));
   }
   
   int initialNumberOfVerts = top - bottom + 1;

   uint32_t newPath;
   int newLength;
   
   step = 1;
   int biggestPicture = 0;
   long long int calcs;
   
   startPath[step] = getCurrentStartPath();
//   printf("startPath: %d\n",startPath[step]);
   
   startIndex[step] = newPathInd[startPath[step] + 1];
   pathsRemaining[step] = startIndex[step] - newPathInd[startPath[step]] + 1;
   
   currentPic[step].startLetter = startIndex[step] % (2*degree);
   
   calcs = 0;
   
   long long int freeEdgeError = 0;
   long long int stepError = 0;
   
   printf("Beginning search\n");
   
   for(;;){
      ++calcs;
      
      if(!(calcs & ((1 << params[P_OUTFREQ]) - 1)) || (step > biggestPicture && params[P_OUTLARGEST])){
//         printVerts();
         printf("max free edges exceeded: %lli\n",freeEdgeError);
         printf("max steps exceeded: %lli\n",stepError);
         printf("current step: %d\n",step);
         printPicture(initialNumberOfVerts,symmetryType,0);
         printf("Number of half edges: %d\n",getNumberOfUnseenEdges());
         printf("Largest Partial: %d faces\n",biggestPicture);
         fflush(stdout);
         
//         if(!(calcs & (((long long int)1 << 33) - 1))){
//            printPicture(initialNumberOfVerts,symmetryType,0);
//            printSymmGraph(strlen(startWord) + 1,symmetryType);
//            fflush(stdout);
//         }
         
         if(step > biggestPicture) biggestPicture = step;
         
      }
      
      --pathsRemaining[step];
      
      if(pathsRemaining[step] < 1){
         
         if(step == 1){
            printf("None found!\n");
            printf("Largest Partial: %d faces\n",biggestPicture);
            return 0;
         }
         
         if(extraFace[step] == 1){
            
            unseeVertsAndEdges(step);
            extraFace[step] = 0;
            
            getTheStartPath(step);
            
            

            continue;
         }
         else{
            backUp();
            onlyOneEdge[step] = 0;     //This might be redundant.
            --step;
            continue;
         }
      }
      
      newLength = newPathLengths[startIndex[step] - pathsRemaining[step]];
      newPath = newPaths[startIndex[step] - pathsRemaining[step]];
      
      
      if(getNumberOfUnseenEdges() + newLength*(degree - 2) > params[P_MAXHALFEDGES]){
         ++freeEdgeError;
         continue;
      }
      
      if(step + (newLength*(degree - 2) + getNumberOfUnseenEdges())/ 2 - 2 > params[P_MAXFACES]){
         ++stepError;
         continue;
      }
      
      currentPic[step].path = newPath;
      currentPic[step].length = newLength;
      
      if(!addNewVerts(startPath[step],newPath,newLength)) continue;
      
      ++step;
      
      if(getNumberOfUnseenEdges() == 0){
//         printVerts();
         printf("found!\n\n");
         
         
         printPicture(initialNumberOfVerts,symmetryType,1);
//         printSymmGraph(initialNumberOfVerts,symmetryType);
         fflush(stdout);
         if(params[P_NUMPICS] == 1){
            exit(0);
         }
         backUp();
         step--;
         params[P_NUMPICS]--;
         continue;
      }
      while(getNumberOfUnseenEdges() >= 4 && canComplete()){
         numAddedVerts[step] = 0;
         pathsRemaining[step] = 1;
         extraFace[step] = 1;
         ++step;
      }
      getTheStartPath(step);
   }
   return 0;
}
