/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-motifs.
It is inspired by the 1985 algorithm of Chiba And Nishizeki to list all k-cliques detailed in "Arboricity and subgraph listing".

To compile:
"gcc kmotif.c -O3 -o kmotif".

To execute:
"./kmotif edgelist.txt k".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the number of k-motifs.
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

	unsigned *cd;//cumulative degree: (start with 0) length=n+1
	unsigned *adj;//list of neighbors with lower degree
} sparse;

void freesparse(sparse *g){
	free(g->cd);
	free(g->adj);
	free(g);
}

//special structure to buid the k-motifs from the graph
typedef struct {
	unsigned char k;//size of k-motifs to list
	unsigned *p;//nodes[p[l]..q[l]]=list of nodes to consider for the l-iest node of the current k-motif
	unsigned *q;
	unsigned *list;//list of candidates
	bool *r;//r[i]==1 iff i should not be considered
	unsigned *motif;//list of nodes in the k-motif
} special;

special *allocspecial(unsigned n,unsigned char k){
	special *s=malloc(sizeof(special));
	s->k=k;
	s->p=malloc(k*sizeof(unsigned));
	s->q=malloc(k*sizeof(unsigned));
	s->list=malloc(n*sizeof(unsigned));
	s->r=calloc(n,sizeof(bool));
	s->motif=malloc(k*sizeof(unsigned));
	return s;
}

void freespecial(special *s){
	free(s->p);
	free(s->list);
	free(s->r);
	free(s->motif);
	free(s);
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
	unsigned *d=calloc(g->n,sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->cd[0]=0;
	g->md=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		g->md=(g->md>d[i-1])?g->md:d[i-1];
		d[i-1]=0;
	}

	g->adj=malloc(2*g->e*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->adj[ g->cd[g->edges[i].s] + d[ g->edges[i].s ]++ ]=g->edges[i].t;
		g->adj[ g->cd[g->edges[i].t] + d[ g->edges[i].t ]++ ]=g->edges[i].s;
	}

	for (i=0;i<g->n;i++) {//sorting neighbors
		qsort(&g->adj[g->cd[i]],d[i],sizeof(unsigned),cmpfunc);
		
		/*
		printf("%u\n",d[i]);
		for (j=g->cd[i];j<g->cd[i+1];j++){
			printf("%u ",d[g->adj[j]]);
		}
		printf("\n");
		exit(1);
		*/
	}

	free(d);
	free(g->edges);
}

//The heart of the function kmotif <3
void kmotif_rec(sparse *g,unsigned char k, unsigned char l, special *s, unsigned long long *n) {
	unsigned i,j,u,v;

	if(l==k-1){

		(*n)+=s->q[l-1]-s->p[l-1];
		/*
		//to enumerate (slower)
		for(i=s->p[l-1]; i<s->q[l-1]; i++){//list all nodes
			u=s->list[i];
			s->motif[l]=u;//here motif[0..l] contains the nodes in the k-motif.
			(*n)++;
			
			//to print all k-motifs
			//printf("%u",g->map[s->motif[0]]);
			//for (j=1;j<k;j++){
			//	printf(" %u",g->map[s->motif[j]]);
			//}
			//printf("\n");
		}
		//*/
	}

	else{
		for(i=s->p[l-1]; i<s->q[l-1]; i++){
			u=s->list[i];
			s->motif[l]=u;
			s->p[l]=i+1;
			s->q[l]=s->q[l-1];
			for (j=g->cd[u];j<g->cd[u+1];j++) {
				v=g->adj[j];
				if (s->r[v]==0){
					s->r[v]=1;
					s->list[s->q[l]++]=v;
				}
			}
	
			kmotif_rec(g,k,l+1,s,n);
	
			for (j=s->q[l-1];j<s->q[l];j++){
				s->r[s->list[j]]=0;
			}
		}
	}
}

//initializing kmotif
unsigned long long kmotif(sparse *g, unsigned char k) {
	unsigned u,v,i;
	unsigned long long n=0;
	special* s=allocspecial(g->n,k);

	for (u=0;u<g->n;u++){
		s->r[u]=1;
		s->motif[0]=u;

		s->p[0]=0;
		s->q[0]=0;
		for (i=g->cd[u];i<g->cd[u+1];i++) {
			v=g->adj[i];
			if (s->r[v]==0){
				s->r[v]=1;
				s->list[s->q[0]++]=v;
			}
		}

		kmotif_rec(g,k,1,s,&n);

		for (i=0;i<s->q[0];i++){
			s->r[s->list[i]]=0;
		}

	}

	freespecial(s);
	return n;
}

int main(int argc,char** argv){
	sparse* g;
	unsigned char k;
	unsigned long long n;

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

	printf("Iterates over all %s-motifs\n",argv[1]);
	k=atoi(argv[1]);

	n=kmotif(g, k);

	printf("Number of %u-motifs: %llu\n",k,n);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	freesparse(g);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}

