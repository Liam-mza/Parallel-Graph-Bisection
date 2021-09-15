#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "graph-bisect-solution.c"

void readGraph(
    char *fileName,
    int *pnVtx,
    int *pnEdge,
    int **padjBeg,
    int **padj)
{
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("ERROR: Unable to open the file %s.\n", fileName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  int nVtx, nEdge;
  int *adjBeg, *adj;
  int chckScan = fscanf(file, " %d %d", pnVtx, pnEdge);
  nVtx = *pnVtx; nEdge = *pnEdge;
  int *edgeLeft, *edgeRight;
  edgeLeft = (int *) malloc(nEdge * sizeof(edgeLeft[0]));
  edgeRight = (int *) malloc(nEdge * sizeof(edgeRight[0]));
  *padjBeg = adjBeg = (int *) malloc((nVtx + 1) * sizeof(adjBeg[0]));
  memset(adjBeg, 0, (nVtx + 1) * sizeof(adjBeg[0]));
  *padj = adj = (int *) malloc(2 * nEdge * sizeof(adj[0]));
  
  int i;
  for (i = 0; i < nEdge; i++) {
    int checkScan = fscanf(file, " %d %d", edgeLeft + i, edgeRight + i);
    adjBeg[edgeLeft[i]]++;
    adjBeg[edgeRight[i]]++;
  }
  for (i = 1; i <= nVtx; i++) { adjBeg[i] += adjBeg[i - 1]; }
  for (i = 0; i < nEdge; i++) {
    adj[--adjBeg[edgeRight[i]]] = edgeLeft[i];
    adj[--adjBeg[edgeLeft[i]]] = edgeRight[i];
  }

  free(edgeLeft);
  free(edgeRight);
  fclose(file);
}

void initPart(
    int **ppart,
    int nVtx,
    int seed)
{
  int *part;
  *ppart = part = (int *) malloc(nVtx * sizeof(part[0]));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) { // Initialize the part array in the master process
    srand(seed);
    int i;
    for (i = 0; i < nVtx; i++) { part[i] = rand() % 2; }
  }
  MPI_Bcast(part, nVtx, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the part array to all other processes
}

void printUsage() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    printf("Usage: mpirun -np [num-procs] ./graph-bisect [graph-file-name] [bisect-algo-name] [max-iters] [epsilon] [random-seed]\n"
        "Arguments:\n"
        "\t[num-procs]: The number of MPI ranks to be used. For the sequential algorithm it has to be 1.\n"
        "\t[graph-file-name]: Name of the graph file.\n"
        "\t[bisect-algo-name]: Name of the graph bisection algorithm. There are three possibilities:\n"
        "\t\tbisect-seq: Sequential bisection algorithm.\n"
        "\t\tbisect-a2a: Parallel bisection algorithm using all-to-all communication routines.\n"
        "\t\tbisect-p2p: Parallel bisection algorithm using point-to-point communication routines.\n"
        "\t[max-iters]: The number of iterations for the label propagation to run.\n"
        "\t[epsilon]: Maximum allowed imbalance ratio for the label propagation.\n"
        "\t[seed]: Random seed used to initialize the part array. This parameter is optional, and you can use it to test the correctness of your algorithm. If not specified, it is set to time(0).\n"
        );
  }
}

void printCutSizeImbal(
    int N,
    int *adjBeg,
    int *adj,
    int *part)
{
  int cutSize = 0;
  int numProcs, procRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  int vtxPerProc = (N + numProcs - 1) / numProcs;
  int vtxBeg = procRank * vtxPerProc;
  int vtxEnd = (procRank == numProcs - 1) ? N : (procRank + 1) * vtxPerProc;
  int i;
  int partZeroCount = 0;
  for (i = vtxBeg; i < vtxEnd; i++) { // Make a pass over the local vertices, compute the local cutsize
    int j;
    if (part[i] == 0) { partZeroCount++; }
    for (j = adjBeg[i]; j < adjBeg[i + 1]; j++) {
      if (part[adj[j]] != part[i]) { // A cut edge is found
        cutSize++;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cutSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // Sum all cut edges
  MPI_Allreduce(MPI_IN_PLACE, &partZeroCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // Sum all part zero counts
  cutSize /= 2; // Note that every cut edge is counted exactly twice.
  double imbal;
  if (2 * partZeroCount >= N) { // There are more 0s
    imbal = (2.0 * partZeroCount) / N;
  } else { // There are more 1s
    imbal = (2.0 * (N - partZeroCount)) / N;
  }
  if (procRank == 0) {
    printf("cutsize: %d imbal: %.3e\n", cutSize, imbal);
  }
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 5) {
    printUsage();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int procRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  char *algoName = argv[2];
  int maxIters; maxIters = atoi(argv[3]);
  int nVtx, nEdge;
  int *adjBeg, *adj;
  int *part;
  double epsilon; epsilon = atof(argv[4]);
  int seed;
  if (argc == 6) { seed = atoi(argv[5]); }
  else { seed = time(0); }

  readGraph(argv[1], &nVtx, &nEdge, &adjBeg, &adj);
  initPart(&part, nVtx, seed);
  double startTime = MPI_Wtime();
  bisectGraph(nVtx, adjBeg, adj, part, maxIters, epsilon, algoName);
  if (procRank == 0) { printf("bisectGraph took %e seconds.\n", MPI_Wtime() - startTime); }
  free(part);

  MPI_Finalize();
  return 0;
}
