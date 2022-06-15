#include <stdio.h>
#include <unistd.h>
#include <pthread.h>

int value = 20;
pthread_t tidA, tidB, tidC;

void* fC(void*ptr) {
  printf("%d\n", *(int*)ptr );
  pthread_exit(NULL);
}

void* fB(void*ptr) {
  value++;
  pthread_create(&tidC, 0, fC, ptr);
  return NULL;
}

void* fA(void*ptr) {
   value++;
   pthread_create(&tidB, 0, fB, ptr);
   pthread_exit(NULL);
}

int c_test() {
  value++;
  pthread_create(&tidA, 0, fA, &value);
  pthread_exit(NULL);
  return 0;
}