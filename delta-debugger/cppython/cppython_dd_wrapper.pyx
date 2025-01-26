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

from cppython_dd_header cimport DeltaDebugger
from libcpp cimport bool

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class PyDeltaDebugger:
    # Hold a C++ instance which we're wrapping
    cdef DeltaDebugger c_dd  

    def __init__(self):
        self.c_dd = DeltaDebugger()

    def get_size(self):
        return self.c_dd.GetSize()

    def print_size(self):
        self.c_dd.PrintSize()

    def update_size(self, scale):
        self.c_dd.UpdateSize(scale)

    def read_cnf(self, file_path):
        self.c_dd.ReadCNF(file_path)

    def print_cnf(self, out_file_path):
        self.c_dd.PrintCNF(out_file_path)

    def init_dd(self):
        self.c_dd.InitDD()

    def init_probdd(self, heuristic, log_file_path, counter_name):
        self.c_dd.InitProbDD(heuristic, log_file_path, counter_name)

    def update_prob_dd(self, test_result):
        self.c_dd.UpdateProbDD(test_result)

    def finished(self):
        return self.c_dd.Finished()

    def print_test_cnf(self, out_file_path):
        self.c_dd.PrintTestCNF(out_file_path)

    def print_final_cnf(self, out_file_path):
        self.c_dd.PrintFinalCNF(out_file_path)
