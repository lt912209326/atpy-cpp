# from distutils.core import setup, Extension
from setuptools import setup, Extension,find_packages
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Compiler import Options
import platform
import numpy as np
import os
Options.generate_cleanup_code=True

# os.environ["CC"] = "gcc-12.1.0"
# os.environ["CXX"] = "gcc-12.1.0"

link_args = ['-static-libgcc',
             '-static-libstdc++',
            #  '-Wl,-Bstatic,--whole-archive',
            #  '-lwinpthread',
            #  '-lpthread',
            #  '-Wl,--no-whole-archive'
             ]

# link_args=None

#Extension args
language ='c++'
include_dirs=["atpy",
            "atpy/core",
            "atpy/core/interface",
            "atpy/core/physics",
            "atpy/core/physics/elements",
            "atpy/core/physics/beamline",
            "atpy/core/physics/utils",
            np.get_include()
            ]



if platform.system() == "Linux":
    extra_compile_args = [ "-O2","-Wall", "-std=c++14","-fopenmp"]
    extra_link_args=["-fopenmp"]
    # extra_link_args+=link_args
elif platform.system() == "Windows":
    extra_compile_args = ["/O2", "/w", "/std:c++14","/openmp"]
    extra_link_args=["/NODEFAULTLIB:libcmt.lib"]

sources=[
          {'name':"atpy.core.interface.constants"  , 'sources':['atpy/core/interface/constants.pyx'], "extra_compile_args":extra_compile_args   },
          {'name':"atpy.core.parser.lexer"  , 'sources':['atpy/core/parser/lexer.pyx'], "extra_compile_args":extra_compile_args},
          {'name':"atpy.core.utils"  , 'sources':['atpy/core/utils.pyx'], "extra_compile_args":extra_compile_args},
          {'name':"atpy.core.parser.parser"  , 'sources':['atpy/core/parser/parser.pyx'],"libraries":["BeamLine"],
          "extra_compile_args":extra_compile_args},
          {'name':"atpy.core.elements"  , 'sources':['atpy/core/elements.pyx'],"libraries":["Elements"], 
          "extra_compile_args":extra_compile_args },
          {'name':"atpy.core.beamline"  , 'sources':['atpy/core/beamline.pyx'],"libraries":["BeamLine","Elements"],
           "extra_compile_args":extra_compile_args,
           "extra_link_args":extra_link_args
          },
        ]
def set_extension():
    extension=[]
    for pkg in sources:
        extension.append( Extension(
                                    **pkg,
                                    language=language, 
                                    include_dirs=include_dirs, 
                                    library_dirs=[ "atpy/core","atpy/core/Release"]
                                    )

                        )
    return extension

# 需在conda环境中进行
# cmake -G “MinGW Makefiles”
# mingw32-make

# sys.system()

 
#Setup config
setup(name='atpy',
      version='2.0.1',
      author='LiuTao',
      packages=find_packages(),

      cmdclass = {'build_ext': build_ext},
      ext_modules=cythonize(set_extension(),
                            annotate=not True,
                            force=True,
                            verbose=not True,
                            compiler_directives={'profile':not True, 'language_level':3,
                                                'boundscheck':not True, 'cdivision':True, 
                                                'initializedcheck':False,'linetrace':not True},
                            # include_path = [".", np.get_include() ]
                            # include_path=[".","./atpy/core","./atpy/core/interface","./atpy/core/parser"] #,"physics/utils"
                            ),
      package_dir={'':'.',},
      package_data = {
            "atpy":["*.py"],
            "atpy.graphics":["*.py"],
            "atpy.core.parser":["*.py"],
            "atpy.core.Release":["*.a","*.lib"],
            "atpy.core.interface":["*.pyd"],
            "atpy.core":["*.pyd","*.a","*.lib"],
            },
        exclude_package_dir = {
            "atpy":[".pyx"],
            "atpy.graphics":[".pyx"],
            "atpy.core.parser":[".pyx"],
            "atpy.core.interface":[".pyx"],
            "atpy.core":[".pyx"],
            },
      zip_safe = False,
      
      install_requires=[
            'Cython>=3.0.0a10',
            'cymem>=2.0.2',
            'matplotlib>=3.1.0',
            'numpy>=1.19.1',
            "pymoo==0.6.1"
        ],
     )
