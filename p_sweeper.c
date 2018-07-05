#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

/*
 * compile: 
 *
 *    $ pgcc -o p_sweeper -mp p_sweeper.c
 *    $ gcc -o p_sweeper -std=c99 -fopenmp p_sweeper.c
 *
 * execution:
 *    (where 32 is number of threads)
 *
 *    $ /usr/bin/time p_sweeper 32
 *
 * p_sweeper.c is an OpenMP memory sweeper that will attempt to
 * allocate NBUFFERS gigabytes of RAM in an attempt to clear
 * Linux file/disk cache prior to model execution.
 * If it can allocate NBUFFERS gigabytes of RAM, it will terminate
 * normally, if it attempts to allocate NBUFFERS gigabytes and
 * it is not available, 'kswapd' will terminate execution, but
 * the net result will still be cleared buffer cache.
 */

//#define NBUFFERS 254   /* 1GB buffers */
#define NBUFFERS 121   /* 1GB buffers */
#define GB 1024*1024*1024

int main(int argc, char** argv) {
    char* buffers[NBUFFERS];
    int thread_count = strtol(argv[1], NULL, 10);

#pragma omp parallel for num_threads(thread_count)
    for(int i = 0; i < NBUFFERS; i++) {
      buffers[i] = malloc(GB);
      memset(buffers[i], 0, GB);
    }

    return 0;
}
