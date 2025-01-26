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

#include "delta_debugger.h"

#include <iostream>
#include <memory>

#include "gpmc_helper.h"

namespace {

// hidden helpers

} // namespace

namespace delta_debugger {

using gpmc_helper::GPMCHelper;

int DeltaDebugger::GetSize() { return size_ * prob_; }

void DeltaDebugger::UpdateSize(int scale) { size_ *= scale; }

void DeltaDebugger::PrintSize() { std::cout << "Size is " << size_ << "\n"; }

void DeltaDebugger::InitProbDD(int heuristic, std::string log_file_path,
                             std::string counter_name) {
  gpmc_helper_->InitProbDD(heuristic, log_file_path, counter_name);
}

void DeltaDebugger::UpdateProbDD(bool test_result) {
  gpmc_helper_->UpdateProbDD(test_result);
}

void DeltaDebugger::PrintTestCNF(std::string out_file_path) {
  gpmc_helper_->PrintTestCNF(out_file_path);
}

void DeltaDebugger::PrintFinalCNF(std::string out_file_path) {
  gpmc_helper_->PrintFinalCNF(out_file_path);
}

bool DeltaDebugger::Finished() { return gpmc_helper_->Finished(); }

void DeltaDebugger::InitDD() {
  gpmc_helper_ = std::make_unique<GPMCHelper>();
  gpmc_helper_->InitGPMCWrapper();
}

void DeltaDebugger::ReadCNF(std::string file_path) {
  gpmc_helper_->ReadCNF(file_path);
}

void DeltaDebugger::PrintCNF(std::string out_file_path) {
  gpmc_helper_->PrintCNF(out_file_path);
}

} // namespace delta_debugger
