//Just the prototype please look down to the function and its documentation
void commP2P(int *nextPart, int *recvCount, int *sendCount, int *part, int *recvIdx, int *sendIdx, int *offsetSnd, int *offsetRcv, int numProcs, int procRank, int startIndex, int nbToProcess);

/** Print the cutsize and the imbalance of a graph partition
 * @param N the number of vertices in the graph
 * @param adjBeg beginning pointer of the adjacency list of each node
 * @param adj the list of adjacent vertex indices
 * @param part the array indicating the part index of each vertex (0 or 1)
 */
void printCutSizeImbal(
    int N,
    int *adjBeg,
    int *adj,
    int *part);

/** Computing the cutsize of a graph partition
 * @param N the number of vertices in the graph
 * @param adjBeg beginning pointer of the adjacency list of each node
 * @param adj the list of adjacent vertex indices
 * @param part the array indicating the part index of each vertex (0 or 1)
 * @param maxIters the number of iterations for which the algorithm must run
 * @param epsilon the load imbalance parameter >= 1.0
 * @param algoName the name of the algorithm to be executed
 */
void bisectGraph(
    int N,
    int *adjBeg,
    int *adj,
    int *part,
    int maxIters,
    double epsilon,
    char *algoName)
{
  int procRank, numProcs;
  printCutSizeImbal(N, adjBeg, adj, part);
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (strcmp(algoName, "bisect-seq") == 0) { // Sequential graph bisection
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE
      int otherPart=0;
      int count[2]={0,0}; //count the number of vertices in each part
      int nextPart[N];
      for (int j=0; j<N; ++j) {
          part[j]==0 ? ++count[0]: ++count[1];
      }
      
      //Starting Iterations
      for(int i = 1; i<=maxIters;++i){
          for (int j=0; j<N; ++j) {
              otherPart=1-part[j];
              if (count[otherPart]>((N/2)*epsilon)) {
                  nextPart[j]=part[j];
              }
              else{
                  int countNeig[2]={0,0};
                  for (int l=adjBeg[j]; l<adjBeg[j+1]; ++l) {
                     part[adj[l]] ==0 ? ++countNeig[0]: ++countNeig[1];
                  }
                  if(countNeig[otherPart]>countNeig[1-otherPart]){
                      nextPart[j]=otherPart;
                      ++count[otherPart];
                      --count[1-otherPart];
                      
                  }
                  else{
                      nextPart[j]=part[j];
                  }
              }
          }
          for (int j=0; j<N; ++j) {
              part[j]=nextPart[j];
          }
          printCutSizeImbal( N, adjBeg, adj, part);
      }
  } else if (strcmp(algoName, "bisect-a2a") == 0) { // Graph bisection with all-to-all communication
    // BEGIN IMPLEMENTATION HERE
      int otherPart=0;
      int count[2]={0,0};
      int tabPartition=(N/numProcs); //value  N/P that will be round up if neccessary
      if (N%numProcs !=0){
          ++tabPartition;
      }
      int nbToProcess=tabPartition; //Number of elements handle by this processor
      if(procRank==numProcs-1){
          nbToProcess=N-((numProcs-1)*nbToProcess); //Last processor takes what remains
      }
      
      int startIndex=procRank*tabPartition;
      int nextPart[nbToProcess];
      int nbPerProc[numProcs]; //Array used by Allgatherv to know how many elements each processor is sending
      int offsets[numProcs]; //offset array for Allgatherv
      
      //initalisation of nbPerProc and offsets arrays
      for (int i=0; i<numProcs; ++i) {
          if(i==numProcs-1){
              nbPerProc[i]=N-((numProcs-1)*tabPartition);
          }else{
              nbPerProc[i]=tabPartition;
          }
          offsets[i]= i*tabPartition;
      }
      
      //Starting Iterations
      for(int i = 1; i<=maxIters;++i){
          count[0]=0;
          count[1]=0;
          for (int j=0; j<N; ++j) {
              part[j]==0 ? ++count[0]: ++count[1];
          }
          for (int j=startIndex; j<startIndex+nbToProcess; ++j) {
              otherPart=1-part[j];
              if (count[otherPart]>((N/2)*epsilon)) {
                  nextPart[j-startIndex]=part[j];
              }
              else{
                  int countNeig[2]={0,0};
                  for (int l=adjBeg[j]; l<adjBeg[j+1]; ++l) {
                     part[adj[l]] ==0 ? ++countNeig[0]: ++countNeig[1];
                  }
                  if(countNeig[otherPart]>countNeig[1-otherPart]){
                      nextPart[j-startIndex]=otherPart;
                      ++count[otherPart];
                      --count[1-otherPart];
                  }
                  else{
                      nextPart[j-startIndex]=part[j];
                  }
              }
          }
        //Communication part
        MPI_Allgatherv(
            nextPart,
            nbToProcess,
            MPI_INT,
            part,
            nbPerProc,
            offsets,
            MPI_INT,
            MPI_COMM_WORLD
          );
          printCutSizeImbal( N, adjBeg, adj, part);
      }

  } else if (strcmp(algoName, "bisect-p2p") == 0) { // Graph bisection with point-to-point communication
    // BEGIN IMPLEMENTATION HERE
      
      int tabPartition=(N/numProcs); //value  N/P that will be round up if neccessary
      if (N%numProcs !=0){ ++tabPartition;}
      
      int nbToProcess=tabPartition; //Number of elements handled by this processor
      if(procRank==numProcs-1){ nbToProcess=N-((numProcs-1)*nbToProcess);} //Le dernier proc prend ce qu'il reste
      
      int startIndex=procRank*tabPartition;
      
      int *recvCount=calloc(numProcs,sizeof(int)); // Use of calloc to have values initialized at 0 in the array
      int *visited=calloc(N,sizeof(int)); // Use of calloc to have values initialized at 0 in the array
      
      for(int i=startIndex;i<startIndex+nbToProcess;++i){
          for (int l=adjBeg[i]; l<adjBeg[i+1]; ++l) {
              int neig = adj[l];
              if (visited[neig]!=1){
                  ++recvCount[neig/tabPartition];
                  visited[neig]=1;
              }
          }
      }
      
      free(visited);

      int totalRecv=0;
      int *idx=calloc(numProcs,sizeof(int)); //idx[i] is the number of elements already added to the list recvIdx[i]
      for (int i=0;i<numProcs;++i){
          totalRecv = totalRecv+recvCount[i];
          if (i>0){
              idx[i]=idx[i-1]+recvCount[i-1];
          }
      }
      visited=calloc(N,sizeof(int)); // Use of calloc to have values initialized at 0 in the array
      
      int recvIdx[totalRecv];
      
      for(int i=startIndex;i<startIndex+nbToProcess;++i){
          for (int l=adjBeg[i]; l<adjBeg[i+1]; ++l) {
              int neig = adj[l];
              if (visited[neig]!=1){
                  recvIdx[idx[neig/tabPartition]]=neig;
                  ++idx[neig/tabPartition];
                  visited[neig]=1;
              }
          }
      }
      
      int sendCount[numProcs];
      
      MPI_Alltoall(
        recvCount,
        1,
        MPI_INT,
        sendCount,
        1,
        MPI_INT,
        MPI_COMM_WORLD
      );
      
      int offsetRcv[numProcs];
      int offsetSnd[numProcs];
      
      //initialization of the different arrays fo communication
      offsetRcv[0]=0;
      offsetSnd[0]=0;
      int totalSend=sendCount[0];
      int otherPart=0;
      int count[2]={0,0};
      for (int i=1;i<numProcs;++i){
          offsetSnd[i]=offsetSnd[i-1]+sendCount[i-1];
          offsetRcv[i]=offsetRcv[i-1]+recvCount[i-1];
          totalSend=totalSend+sendCount[i];
      }
      
      int sendIdx[totalSend];
      
      MPI_Alltoallv(
        recvIdx,
        recvCount,
        offsetRcv,
        MPI_INT,
        sendIdx,
        sendCount,
        offsetSnd,
        MPI_INT,
        MPI_COMM_WORLD
      );
      
      free(idx);
      
      int nextPart[nbToProcess];
      
      //Start iterations
      for(int i = 1; i<=maxIters;++i){
          count[0]=0;
          count[1]=0;
          for (int j=0; j<N; ++j) {
              part[j]==0 ? ++count[0]: ++count[1];
          }
          for (int j=startIndex; j<startIndex+nbToProcess; ++j) {
              otherPart=1-part[j];
              if (count[otherPart]>((N/2)*epsilon)) {
                  nextPart[j-startIndex]=part[j];
              }
              else{
                  int countNeig[2]={0,0};
                  for (int l=adjBeg[j]; l<adjBeg[j+1]; ++l) {
                     part[adj[l]] ==0 ? ++countNeig[0]: ++countNeig[1];
                  }
                  if(countNeig[otherPart]>countNeig[1-otherPart]){
                      nextPart[j-startIndex]=otherPart;
                      ++count[otherPart];
                      --count[1-otherPart];
                  }
                  else{
                      nextPart[j-startIndex]=part[j];
                  }
              }
          }
          
          //Communication part
          commP2P(nextPart, recvCount, sendCount, part, recvIdx,sendIdx, offsetSnd, offsetRcv, numProcs, procRank, startIndex, nbToProcess);
          
          printCutSizeImbal( N, adjBeg, adj, part);
      }
      
      
      free(visited);
      
      free(recvCount);
  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}


/** Do all the communication needed by the p2p algorithm
 * @param nextPart the updated part of the array "part"
 * @param recvCount array containing the number of elements to receive from each processors
 * @param sendCount array containing the number of elements to send to each processors
 * @param part the array indicating the part index of each vertex (0 or 1)
 * @param recvIdx array containing the index of elements to receive from each processors
 * @param sendIdx array containing the index of elements to send to each processors
 * @param offsetSnd array of offset to know where each sequence of elements start for which processeur in sendIdx
 * @param offsetRcv array of offset to know where each sequence of elements start for which processeur in recvIdx
 * @param numProcs number of processsors
 * @param procRank rank of the processor
 * @param startIndex index where this processor start to update
 * @param nbToProcess number of elements to process
 */
void commP2P(int *nextPart, int *recvCount, int *sendCount, int *part, int *recvIdx, int *sendIdx, int *offsetSnd, int *offsetRcv, int numProcs, int procRank, int startIndex, int nbToProcess){
    
    int *recvPart[numProcs];
    
    int *sendPart[numProcs];
    for(int i=0; i<numProcs;++i){
        recvPart[i]=calloc(recvCount[i],sizeof(int));
        sendPart[i]=calloc(sendCount[i],sizeof(int));
        for(int j=0;j<sendCount[i];++j){
            sendPart[i][j]=nextPart[sendIdx[offsetSnd[i]+j]-startIndex];
        }
    }
    
    MPI_Request sendRequest[numProcs];
    MPI_Status stat;
    for(int i=0; i<numProcs;++i){
        MPI_Isend(sendPart[i], sendCount[i], MPI_INT, i, i,
                  MPI_COMM_WORLD, &(sendRequest[i]));
    }

    for(int i=0; i<numProcs;++i){
        MPI_Recv(recvPart[i], recvCount[i], MPI_INT, i, procRank,
                 MPI_COMM_WORLD, &stat);
    }
    MPI_Waitall(numProcs, sendRequest,MPI_STATUSES_IGNORE);
    
    for (int j=startIndex; j<startIndex+nbToProcess; ++j) {
        part[j]=nextPart[j-startIndex];
    }
    for (int i=0; i<numProcs; ++i) {
        for (int j=0; j<recvCount[i]; ++j) {
            part[recvIdx[offsetRcv[i]+j]]=recvPart[i][j];
        }
    }
    for(int i=0; i<numProcs;++i){
        free(sendPart[i]);
        free(recvPart[i]);
    }
    
}
