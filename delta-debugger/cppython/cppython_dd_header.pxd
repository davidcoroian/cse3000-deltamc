"""

ProbDD implementation -- Copyright (c) 2025, David Coroian

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

"""

# distutils: language = c++

from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "delta_debugger.h" namespace "delta_debugger":
  cdef cppclass DeltaDebugger:
        DeltaDebugger() except +
        DeltaDebugger(int, int) except +
        int GetSize()
        void PrintSize()
        void UpdateSize(int scale)
        void InitDD()
        void InitProbDD(int heuristic, string log_file_path, string counter_name)
        void UpdateProbDD(bool test_result)
        void ReadCNF(string file_path)
        void PrintCNF(string out_file_path)
        void PrintTestCNF(string out_file_path)
        void PrintFinalCNF(string out_file_path)
        bool Finished()
        int size_, prob_

cdef extern from "delta_debugger.cpp":
    pass
