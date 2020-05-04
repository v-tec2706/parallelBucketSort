#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "list_utils.h"

#define SIZE 10000000
#define INT_MAX 2147483647
#define DEBUG 0

void randNumbers(int* numbers){
  int i;

  #pragma omp parallel
  {
  unsigned int myseed = (unsigned int) omp_get_wtime() + omp_get_thread_num();
  #pragma omp for schedule(auto) private(i)
   for(i=0; i<SIZE; i++)
    {
      numbers[i] = (int) rand_r(&myseed)%200000;
    }
  }
}

void validateIfSorted(int* numbers){
  int i;
  int res = 1;
   #pragma omp parallel for shared(res) private(i)
   for(i=0; i<SIZE-1; i++)
    {
      if (numbers[i+1] - numbers[i] < 0){
        res *= 0;
      }
    }
  if(res == 0){
    printf("Number are not sorted properly! :\()");
  }
}

void splitIntoBuckets(int startBucket, int endBucket, int* numbers, Node** buckets, float bucketSize, int min, omp_lock_t* locks){
  int i;
  for(i=startBucket; i<endBucket; i++){
    int val = numbers[i];
    int relativeVal = val - min;
    int bucketNumber = (relativeVal / bucketSize) ;
    if(bucketNumber >= SIZE){
      bucketNumber = SIZE-1;
    }
    omp_set_lock(&locks[bucketNumber]);
    push( &buckets[bucketNumber%SIZE], val);
    omp_unset_lock(&locks[bucketNumber]);
  }
}

void sortBuckets(int startBucket, int endBucket, Node** buckets){
  int i;
  for(i=startBucket; i<endBucket; i++){
      quickSort(&buckets[i]);
  }
}

void findNumberOfItems(int startBucket, int endBucket, Node** buckets, int* bucketsSizeInBunch){
  int totalSize = 0;
  int i;
  for(i=startBucket; i<endBucket; i++){
    struct Node *node = buckets[i];
    while(node != NULL){
      node = node->next;
      totalSize++;
    }
  }
  bucketsSizeInBunch[omp_get_thread_num()]=totalSize;
}

void setWritingOffset(int n_threads, int* bucketsSizeInBunch){
  #pragma omp single
  {
      int i;
      int total = 0;
      for(i=0;i<n_threads;i++){
        int tmp = bucketsSizeInBunch[i];
        bucketsSizeInBunch[i]=total;
        total += tmp;
        if (DEBUG) printf("Thread no: %d has %d elems to perform starting at %d\n", i, tmp, bucketsSizeInBunch[i]);
      }
      if(total != SIZE){
        printf("Not all elements from buckets are covered in thread's subtasks!\n");
      }
  }
}

int findMax(int* A, int i, int j)
{
    int idx;
    int max_val = 0;

    #pragma omp parallel for reduction(max:max_val)
    for (idx = i; idx < j; idx++)
       max_val = max_val > A[idx] ? max_val : A[idx];

    return max_val;
}


int findMin(int* A, int i, int j)
{
    int idx;
    int min_val = INT_MAX;

    #pragma omp parallel for reduction(min:min_val)
    for (idx = i; idx < j; idx++)
       min_val = min_val < A[idx] ? min_val : A[idx];

    return min_val;
}

void writeBucketsIntoFinalTable(int startBucket, int endBucket,Node** buckets,int* numbers, int* bucketsSizeInBunch){
  int offset = bucketsSizeInBunch[omp_get_thread_num()];
  int i;
  for(i=startBucket; i<endBucket; i++){
    struct Node *node = buckets[i];
    while (node != NULL)
    {
        numbers[offset] = node->data;
        offset++;
        node = node->next;
    }
  }
}

int main(int argc, char *argv[])
{
    long i;
    double start, end, rand_start, rand_end, sort_start, sort_end, split_start, split_end, rewrite_start, rewrite_end;
    float bucketSize;

    int* numbers = (int*) calloc(SIZE,sizeof(int));
    Node ** buckets = (Node**) calloc(SIZE, sizeof(struct Node*));
    omp_lock_t* locks = (omp_lock_t*) calloc(SIZE, sizeof(omp_lock_t));
    int * bucketsSizeInBunch = (int*) calloc(omp_get_num_threads(), sizeof(int));

    for (i = 0; i < SIZE; i++){
        omp_init_lock(&locks[i]);
    }

    start = omp_get_wtime();

    rand_start = omp_get_wtime();
    randNumbers(numbers);
    rand_end = omp_get_wtime();

    float min=(float) findMin(numbers, 0, SIZE);
    float max=(float) findMax(numbers, 0, SIZE);

    int n_threads, id;
    long ppt, sp;

    #pragma omp parallel private(id, ppt, sp) \
        shared(n_threads, buckets, numbers, max, min, bucketsSizeInBunch)
    {

        #pragma omp single
        {
            bucketSize = (max-min) / SIZE;
            if (DEBUG) printf("\nmin: %.6f, max: %.6f, bucketSize: %.6f \n", min, max, bucketSize);
            n_threads = omp_get_num_threads();
        }

        id = omp_get_thread_num();

        ppt = (long)((SIZE + n_threads - 1) / n_threads);
        sp = id * ppt;

        if (id == n_threads - 1)
           ppt = SIZE - sp;

        if (DEBUG) printf("Thread no: %d / %d, start point: %ld, end point: %ld\n", id, n_threads - 1, sp, sp+ppt-1);

        split_start = omp_get_wtime();
        splitIntoBuckets(sp, sp+ppt, numbers, buckets, bucketSize, min, locks);
        split_end = omp_get_wtime();

        #pragma omp barrier

        sort_start = omp_get_wtime();

        sortBuckets(sp, sp + ppt, buckets);

        sort_end = omp_get_wtime();

        findNumberOfItems(sp, sp+ppt, buckets, bucketsSizeInBunch);

        #pragma omp barrier
        setWritingOffset(n_threads, bucketsSizeInBunch);

        rewrite_start = omp_get_wtime();
        writeBucketsIntoFinalTable(sp, sp+ppt, buckets, numbers, bucketsSizeInBunch);
        rewrite_end = omp_get_wtime();

        validateIfSorted(numbers);

    }

  end = omp_get_wtime();

  if(DEBUG){
    printf("\n");
    for(i=0; i<SIZE; i++)
     {
       printf("%d ", numbers[i]);
     }
  }

  printf("\nRand time is: %fs\n", rand_end - rand_start);
  printf("\nSplit time is: %fs\n", split_end - split_start);
  printf("\nSort time is: %fs\n", sort_end - sort_start);
  printf("\nWrite time is: %fs\n", rewrite_end - rewrite_start);
  printf("\nTotal time is: %fs\n", end - start);

  free(numbers);
  free(buckets);
  free(locks);
  free(bucketsSizeInBunch);

  return 0;
}
