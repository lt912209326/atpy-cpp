# atpy-cpp
This package is for accelerator design and optimization, working on python as a library. The core code is written in C++14, wrapped by cython.
## Install
This package depends on cmake, gcc/msvc(support c++14), cython, pymoo 
1. Install Visual Studio Code, then, install the extensions of cmake and cmake tools, and build the c++ library with cmake extension.
2. run `python setup.py bdist_wheel` in the directory, the the python package will be compiled.
3. install this package as a python library `pip install ./dist/atpy.xxx.xxx.whl`, where xxx.xxx is some information about your platform.
