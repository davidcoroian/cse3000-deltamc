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

#ifndef DELTA_DEBUGGER_H
#define DELTA_DEBUGGER_H

#include <memory>
#include <string>

#include "gpmc_helper.h"

namespace delta_debugger {

class DeltaDebugger {
public:
  DeltaDebugger(){};
  DeltaDebugger(int size, int prob) : size_(size), prob_(prob){};

  // dummy cython test
  int GetSize();
  void PrintSize();
  void UpdateSize(int scale);

  // DD
  void InitDD();
  void InitProbDD(int heuristic, std::string log_file_path,
                             std::string counter_name);

  void UpdateProbDD(bool test_result);

  bool Finished();

  void ReadCNF(std::string file_path);
  void PrintCNF(std::string out_file_path);
  void PrintTestCNF(std::string out_file_path);
  void PrintFinalCNF(std::string out_file_path);

private:
  int size_, prob_;
  std::unique_ptr<gpmc_helper::GPMCHelper> gpmc_helper_;
};

} // namespace delta_debugger

#endif