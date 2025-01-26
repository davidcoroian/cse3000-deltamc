from setuptools import Extension, setup

from Cython.Build import cythonize

delta_debugger_path = "path/to/delta-debugger/"
gpmc_path = delta_debugger_path + "gpmc-main/"

srcfiles = ["cppython_dd_wrapper.pyx",
            "gpmc_helper.cpp",
            gpmc_path + "preprocessor/Preprocessor.cc",
            gpmc_path + "preprocessor/lib_sharpsat_td/subsumer.cpp",
            gpmc_path + "preprocessor/IFlowCutter.cc",
            gpmc_path + "preprocessor/TreeDecomposition.cc",
            gpmc_path + "core/ComponentCache.cc",
            gpmc_path + "core/ComponentManager.cc",
            gpmc_path + "core/Config.cc",
            gpmc_path + "core/Counter.cc",
            gpmc_path + "core/dd_main.cc",
            gpmc_path + "core/Instance.cc",
            gpmc_path + "core/Solver.cc",
            gpmc_path + "ddnnf/DecisionTree.cc",
            gpmc_path + "flow-cutter-pace17/src/cell.cpp",
            gpmc_path + "flow-cutter-pace17/src/greedy_order.cpp",
            gpmc_path + "flow-cutter-pace17/src/list_graph.cpp",
            gpmc_path + "flow-cutter-pace17/src/tree_decomposition.cpp",
            gpmc_path + "utils/Options.cc",
            gpmc_path + "utils/System.cc",
            gpmc_path + "preprocessor/TestSolver.cc",]

incldirs = [delta_debugger_path + "gpmc-main", 
            delta_debugger_path + "cppython",
            ]

extensions = [
  Extension("cppython_dd_wrapper", 
            sources=srcfiles,
            include_dirs=incldirs,
            libraries=["gmpxx", "gmp", "mpfr", "z"],
            ),
]

setup(ext_modules=cythonize(extensions))
