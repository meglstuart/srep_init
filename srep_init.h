#ifndef __SREP_INIT_H
#define __SREP_INIT_H

#include <string>

class srep_init
{
public:
  double dt = 0.1f;
  double smoothAmount = 0.1f;
  int maxIter = 10;
  std::string inputMesh = "../test_data/bunny.off";
  srep_init(double d, double smooth, int max);

  int process_inputMesh();
  int load_mesh();
  int step_forwardflow();
};

#endif
