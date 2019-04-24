#include <stdio.h>

#include "lib/third_party/samtools/sam.h"
#include "lib/third_party/samtools/faidx.h"


int main(int argc, char** argv) {

  fai_build(argv[1]);
  
  return 0;
}
