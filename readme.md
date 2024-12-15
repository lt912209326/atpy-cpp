# atpy
This package is for accelerator design and optimization, working on python as a library. The core code is written in C++14, wrapped by cython.
## Install
This package depends on cmake, gcc/msvc(support c++14, linux better), cython, pymoo, matplotlib, numpy and Eigen 
Pre-work: install msvc by Visual Studio 2022 生成工具(select ‘使用C++的桌面开发’ ), install
1. Install Visual Studio Code, then, install the extensions of cmake and cmake tools, and build the c++ library with cmake extension. (select your compiler, select the release mode, build the c++ project)
2. run `python setup.py bdist_wheel` in the directory, the the python package will be compiled.
3. install this package as a python library `pip install ./dist/atpy.xxx.xxx.whl`, where xxx.xxx is some information about your platform.
