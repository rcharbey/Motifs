/*
Info:
Feel free to use these lines as you wish.
This program computes a lower bound and an upper bound on the number of k-motifs that can be induced by a node and k-1 of its neighbors. I call such a motif a k-star for simplicity as there is a star on k nodes in it.

The program runs in linear time and memory.

The lower and upper bounds are empiricaly good ((upper-lower)/lower<0.1) for k>4 on some real-world graphs.

To compile:
"gcc kstarLBUB.c -O3 -o kstarLBUB".

To execute:
"./kstarLBUB k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print an lower-bound and upper-bound on the number of k-stars.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define NLINKS 100000000 //maximum number of edges in the input graph (used for memory allocation), will automatically increase if needed!

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg ;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	unsigned md;//maximum degree

	unsigned *map;//map[i]=original label of node i

	unsigned *d;//d[i]=degree of node i
	unsigned *cd;//cumulative degree: (start with 0) length=n+1
	unsigned *adj;//list of neighbors with lower degree
} sparse;

void freesparse(sparse *g){
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

sparse* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {//Add one edge
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		g->e++;
		if (g->e==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//for degree ordering
int compare_nodedeg(void const *a, void const *b){
	nodedeg const *pa = a;
	nodedeg const *pb = b;
	return (pa->deg < pb->deg) ? 1 : -1;/////////
}

//sorting nodes in non-increasing order of degree and relabeling
void degord(sparse *g) {
	unsigned i;
	unsigned *newlabel;
	nodedeg *nodedeglist=malloc(g->n*sizeof(nodedeg));
	for (i=0;i<g->n;i++) {
		nodedeglist[i].node=i;
		nodedeglist[i].deg=0;
	}
	for (i=0;i<g->e;i++) {
		nodedeglist[g->edges[i].s].deg++;
		nodedeglist[g->edges[i].t].deg++;
	}
	qsort(nodedeglist,g->n,sizeof(nodedeg),compare_nodedeg);
	newlabel=malloc(g->n*sizeof(unsigned));
	g->map=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++) {
			newlabel[nodedeglist[i].node]=i;
			g->map[i]=nodedeglist[i].node;
	}
	free(nodedeglist);

	for (i=0;i<g->e;i++) {
		g->edges[i].s=newlabel[g->edges[i].s];
		g->edges[i].t=newlabel[g->edges[i].t];
	}

}

//for future use in qsort
int cmpfunc (const void * a, const void * b){
	if (*(unsigned*)a>*(unsigned*)b){
		return 1;
	}
	return -1;
}

//Building the sparse graph structure
void mksparse(sparse *g){
	unsigned i,j;

	g->d=calloc(g->n,sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->d[g->edges[i].s]++;
		g->d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	g->md=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+g->d[i-1];
		g->md=(g->md>g->d[i-1])?g->md:g->d[i-1];
		g->d[i-1]=0;
	}

	g->adj=malloc(2*g->e*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->adj[ g->cd[g->edges[i].s] + g->d[ g->edges[i].s ]++ ]=g->edges[i].t;
		g->adj[ g->cd[g->edges[i].t] + g->d[ g->edges[i].t ]++ ]=g->edges[i].s;
	}

	for (i=0;i<g->n;i++) {//sorting neighbors
		qsort(&g->adj[g->cd[i]],g->d[i],sizeof(unsigned),cmpfunc);
		
		/*
		printf("%u\n",d[i]);
		for (j=g->cd[i];j<g->cd[i+1];j++){
			printf("%u ",d[g->adj[j]]);
		}
		printf("\n");
		exit(1);
		*/
	}

	free(g->edges);
}

inline unsigned long long choose(unsigned n, unsigned char k){
    if (k == 0)
	 return 1;
    return (n * choose(n - 1, k - 1)) / k;
}


//Computing a lower bound on number of kstars
void kstar(sparse *g, unsigned char k, unsigned long long *ub, unsigned long long *lb) {
	unsigned u,v,i,d;
	bool *tab=calloc(g->n,sizeof(bool));

	*lb=0;
	*ub=0;

	for (u=0;u<g->n;u++){

		d=g->d[u];
		if (d>=k-1){
			(*ub)+=choose(d,k-1);
		}
		if (d>=k-1){

			for (i=g->cd[u];i<g->cd[u+1];i++) {
				tab[g->adj[i]]=1;
			}

			for (i=g->cd[u];i<g->cd[u+1];i++) {
				v=g->adj[i];
				if (v<=u){
					d--;
					if (d<k-1){
						break;
					}
				}
			}

			if (d>=k-1){
				(*lb)+=choose(d,k-1);
			}

			for (i=g->cd[u];i<g->cd[u+1];i++) {
				tab[g->adj[i]]=0;
			}

		}
	}

	free(tab);

}

int main(int argc,char** argv){
	sparse* g;
	unsigned char k;
	unsigned long long ub,lb;

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[2]);

	g=readedgelist(argv[2]);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);

	printf("Sorting nodes in non-increasing order of degree\n");
	degord(g);

	printf("Building the graph structure\n");

	mksparse(g);

	printf("Maximum degree = %u\n",g->md);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing lower-bound on the number of %s-stars\n",argv[1]);
	k=atoi(argv[1]);

	kstar(g, k,&ub,&lb);

	printf("Lower-bound on the number of %u-stars: %llu\n",k,lb);
	printf("Upper-bound on the number of %u-stars: %llu\n",k,ub);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	freesparse(g);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}

