#include "tools.hh"
#include "TROOT.h"
#include <cstdio>
#include <cstdlib>
int main(int argc, char* argv[]) {
  int in = 0;
  if (argc > 1) in = atoi(argv[1]);
  int testint = testfunc(in);
  printf("bonjour et %i\n", testint);
  return 0;
}
