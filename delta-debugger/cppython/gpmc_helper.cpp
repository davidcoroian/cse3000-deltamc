/*
* ProbDD implementation -- Copyright (c) 2025, David Coroian

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "gpmc_helper.h"

#include <string>

#include "core/dd_main.h"

namespace gpmc_helper {

using gpmc_wrapper::GPMCWrapper;

void GPMCHelper::InitGPMCWrapper() {
  gpmc_wrapper_ = std::make_unique<GPMCWrapper>();
}

void GPMCHelper::InitProbDD(int heuristic, std::string log_file_path,
                             std::string counter_name) {
  gpmc_wrapper_->InitProbDD(heuristic, log_file_path, counter_name);
}

void GPMCHelper::UpdateProbDD(bool test_result) {
  gpmc_wrapper_->UpdateProbDD(test_result);
}

void GPMCHelper::PrintFinalCNF(std::string file_path) {
  gpmc_wrapper_->PrintReducedCNF(file_path);
}

void GPMCHelper::PrintTestCNF(std::string file_path) {
  gpmc_wrapper_->SelectProbDDTest(file_path);
}

bool GPMCHelper::Finished() { return gpmc_wrapper_->IsProbDDFinished(); }

void GPMCHelper::ReadCNF(std::string file_path) {
  // Initializes GPMC and reads from file_path.
  gpmc_wrapper_->main_dd(file_path);
}

void GPMCHelper::PrintCNF(std::string file_path) {
  gpmc_wrapper_->PrintDeltaCNF(file_path);
}

} // namespace gpmc_helper
