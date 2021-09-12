/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 25, 2010 
 * Usage: This is part of the bicluster package. Use, redistribute, modify
 *        without limitations.
 *
 * Produces two graphs sequentially, derived from microarray data.
 *  
 * The first graph is generated from the raw data where an edge is defined
 * as two genes having common components from same condition, the score on
 * the edge being the number of the same components. The edges in the first
 * graph, are used as vertices in the second graph, where edges are defined
 * as the common columns between two co-vertex edges in the first graph,
 * with scores defined as the number of same columns for all three genes.
 *
 */

#include "make_graph.h"
/*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/
static const int HEAP_SIZE = 20000000;
/*extern double uber_ratio;*/
/**************************************************************************/

/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge in the first graph*/
static int str_intersect_r (const discrete *s1, const discrete *s2)
{
	int common_cnt = 0;
	/* s1 and s2 of equal length, so we check s1 only */
	int i;
	for (i=0;i<cols;i++)
	{
		if (*s1==*s2 && (*s1!=0)) 
			common_cnt++;
		s1++; 
		s2++;
	}
	return common_cnt;
}

/**************************************************************************/

/*seed_deduct is not used in the whole program while it may be useful someday*/
void seed_deduct (const discrete *s)
/* remove a row from the profile */
{
	int i;
	discrete ss;
	for (i=0; i<cols; i++)
	{
		ss = s[i];
		profile[i][ss]--;
	}
}

void seed_update (const discrete *s)
{
	int i;
	for (i = 0; i < cols; i++) 
		profile[i][s[i]]++;
}

/*****************************************************************************/

/* Fibonacci heap related subroutines */
static int edge_cmpr(void *a, void *b)
{
	int score_a, score_b;
	score_a = ((Edge *)a)->score;
	score_b = ((Edge *)b)->score;

	if (score_a < score_b) return -1;
	if (score_a == score_b) return 0;
	return 1;
}

/* Maintain a fixed size heap */
static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min)
{
	if (a->fh_n < HEAP_SIZE) 
	{
		fh_insert(a, (void *)i);
	}
	else
	{
		if (edge_cmpr(cur_min, i) < 0)
		{
			/* Remove least value and renew */
			fh_extractmin(a);
			fh_insert(a, (void *)i);
			/* Keep a memory of the current min */
			*cur_min = (Edge *)fh_min(a);
		}
	}
}

/*sort the edges in decrease order so that e1 is the largest edge*/
static void fh_dump(struct fibheap *a, Edge **res)
{
	int i;
	int n = a->fh_n;
	for (i=n-1; i>=0; i--)
		res[i] = (Edge *) fh_extractmin(a);
}

/**************************************************************************/

void make_graph (const char* fn)
{
	FILE *fw = mustOpen(fn, "w");
	int i, j,  k, cnt,block_size, block_border;
	int block_num = 10; 
	int rec_num = 0;
	if (po->COL_WIDTH == 2) 
		po->COL_WIDTH = MAX(cols/20, 4);
	
	/* finding too many zeros genes 
	int too_many_zero_array[rows];
	find_too_many_zero_genes(too_many_zero_array, arr_c);*/
	
	/* edge_ptr describe edges */
	AllocArray(edge_list, HEAP_SIZE);

	/* Allocating heap structure */
	struct fibheap *heap;
	heap = fh_makeheap();
	fh_setcmp(heap, edge_cmpr);

	/* Generating seed list and push into heap */
	progress("Generating seed list (minimum weight %d)", po->COL_WIDTH);	

	Edge __cur_min = {0, 0, po->COL_WIDTH};
	Edge *_cur_min = &__cur_min;
	Edge **cur_min = & _cur_min;
	/* iterate over all genes to retrieve all edges */
	block_size = rows/block_num;
	for (k = 0; k < block_num; k++)
	{
		block_border = (k+1) * block_size;
		if (k == block_num - 1)
			block_border = rows;		
		for (i = k* block_size; i < block_border; i++)
		{
			/* if(too_many_zero_array[i] == 1) continue;*/
			for (j = i+1; j < block_border; j++)
			{
			    cnt = str_intersect_r(arr_c[i], arr_c[j]);
		    	cnt += uber_score(i, j) * uber_ratio * cnt;     /*uber_ratio*/
			    if (cnt < (*cur_min)->score) continue;
		
		    	AllocVar(edge_ptr);
		    	edge_ptr -> gene_one = i;
			    edge_ptr -> gene_two = j;
			    edge_ptr -> score = cnt;	
			    fh_insert_fixed(heap, edge_ptr, cur_min);
			}
		}
	}
	rec_num = heap->fh_n;
	if (rec_num == 0)
		errAbort("Not enough overlap between genes");

	/* sort the seeds */
	uglyTime("%d seeds generated", rec_num);
	ReAllocArray(edge_list, rec_num);
	fh_dump(heap, edge_list);

	/* bi-clustering */
        int n_blocks = 0;
	progress("Modules generated");
	n_blocks = cluster(fw, edge_list, rec_num);
	uglyTime("%d modules are written to %s", n_blocks, fn);

	/* clean up */
	for (i=0; i<rec_num; i++)
		free(edge_list[i]);
	free(edge_list);
}

/***************************************************************************/
