#ifndef DD_MAIN_H
#define DD_MAIN_H

#include <gmpxx.h>
#include <memory>
#include <string>
#include <vector>

#include "Counter.h"

namespace gpmc_wrapper {

class GPMCWrapper {
public:
  GPMCWrapper(){};
  // Initializes and reads input.
  void main_dd(std::string input_file_path);
  void PrintDeltaCNF(std::string out_file_path);

  void InitSolver(GPMC::Configuration config) {
    S_ = std::make_unique<GPMC::Counter<mpz_class>>(config);
  }

  void InitProbDD(int heuristic, std::string log_file_path,
                             std::string counter_name);
  void UpdateProbDD(bool test_result);
  void SelectProbDDTest(std::string tested_cnf_output_path);
  void PrintReducedCNF(std::string reduced_cnf_output_path);
  bool IsProbDDFinished();

  std::unique_ptr<GPMC::Counter<mpz_class>> S_;
};

} // namespace gpmc_wrapper

#endif
