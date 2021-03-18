#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#include "mmio.h"

double** G;
int M, N;
//this is a test commit
void dijkstraBetweenness(int, int, int, int);
void readMatrix(int, char**);
void freeMatrix();

int main(int argc, char * argv[])
{
    double begin, end, time_spent;
    int start, B;

    readMatrix(argc, argv);

    printf("Enter the starting node: \n");
    scanf("%d", &start);
    printf("Enter the node: \n");
    scanf("%d", &B);
    
    begin = omp_get_wtime();
    dijkstraBetweenness(M, N, start, B);
    end = omp_get_wtime();
    
    time_spent = (end - begin);
    printf("%lf\n", time_spent);
    
    freeMatrix();
    return 0;
}

void dijkstraBetweenness(int m, int n,int startnode, int B)
{
	int visited[m], pred[m], count, pathCount, nextnode, i, j, btw;
    double **cost, distance[m], mindistance, betweenness;
    
    cost = (double**)malloc(sizeof(double*)*m);
    for (i = 0; i < m; ++i) {
        cost[i] = (double*)malloc(sizeof(double)*n);
    }

    omp_set_num_threads(4); 

	for(i=0;i<n;i++)
        #pragma omp parallel for private(cost)
		{
            for(j=0;j<n;j++)
            {
                if(G[i][j]==0)
                {
                    cost[i][j]=INFINITY;
                }
                else
                {
                    cost[i][j]=G[i][j];
                }
            }
        }
	#pragma omp parallel for private(distance)
    {
        for(i=0;i<n;i++)
        {
            distance[i]=cost[startnode][i];
            pred[i]=startnode;
            visited[i]=0;
        }
    }
	
	
	distance[startnode]=0;
	visited[startnode]=1;
	count=1;

	while(count<n-1)
	{
		mindistance=INFINITY;
		for(i=0;i<n;++i)
			if(distance[i]<mindistance&&!visited[i])
			{
				mindistance=distance[i];
				nextnode=i;
			}
            
        visited[nextnode]=1;
        #pragma omp parallel for private(distance, pred, visited)
        {
            for(i=0;i<n;i++)
            {   
                if(!visited[i])
                {
                    if(mindistance+cost[nextnode][i]<distance[i])
                    {
                        
                        distance[i]=mindistance+cost[nextnode][i];
                        pred[i]=nextnode;
                    }
                }
            }
                
        }
		count++;
	}
    pathCount = 0;
    btw = 0;

	for(i=0;i<n;i++)
    {
		if(i!=startnode)
		{
			printf("\nPath=%d",i);
            pathCount++;
            if(i == B)
            {
                btw++;
            }
			j=i;
			do
			{
				j=pred[j];
				printf("<-%d",j);
                if(j == B)
                {
                    btw++;
                }
			}while(j!=startnode);
	    }
    }
    betweenness = (double)btw / pathCount;
    printf("\nPath count: %d", pathCount);
    printf("\nBTW: %d", btw);
    printf("\nBetweennes centrality of %d: %f", B, betweenness);
    printf("\n");
    for (i = 0; i < m; ++i) {
        free(cost[i]);
    }
    free(cost);
}

void readMatrix(int argc, char** argv)
{
    int n, i, j, ret_code, *I, *J, nz, weighted = 0;
    MM_typecode matcode;
    FILE *f; 
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
        if (argc >= 3) {
            if (!strcmp(argv[2], "weighted"))
                weighted = 1;
        }
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));

    G = malloc( M*sizeof(double*) );
    for(i=0; i<M; i++)
    {
        G[i] = malloc( N*sizeof(double));
    }


    for(i=0; i<M; i++)
        for(int j=0; j<N;j++)
            G[i][j]=0;

    if (weighted) {
        for (i=0; i<nz; i++)
        {
            fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
            I[i]--; 
            J[i]--;
            printf("%d: %d %d\n",i, I[i], J[i]);
        }
    }
    else {
        for (i=0; i<nz; i++)
        {
            fscanf(f, "%d %d\n", &I[i], &J[i]);
            I[i]--; 
            J[i]--;
            val[i] = 1;
            printf("%d: %d %d\n",i, I[i], J[i]);
        }
    }

    if (f !=stdin) fclose(f);

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i=0; i<nz; i++)
    {     
        G[I[i]][J[i]] = fabs(val[i]);   
    }
    
    for (i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%.3f\t", G[i][j]);
        }
        printf("\n");
    }
    free(I);
    free(J);
}

void freeMatrix()
{
    int i;
    for(i=0; i<M; i++)
    {
        free(G[i]);
    }
    free(G);
}
