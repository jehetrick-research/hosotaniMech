#include <stdio.h>

int main() {
   int i;
   float zap[2];

   for(i=0; i<10; i++) {
      zap[i] = i*i;
      printf("%d %f\n", i, zap[i]);
   }
}
