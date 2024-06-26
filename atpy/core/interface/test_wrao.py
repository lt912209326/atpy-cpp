
import os
import glob

import autowrap
import autowrap.Code
import autowrap.CodeGenerator
import autowrap.DeclResolver
import autowrap.Main
import autowrap.PXDParser
import autowrap.Utils


def test_full_lib():
    """
    Example with multi-file library and multi-file result.
    This shows a full run through of a case where multiple class files (A, B,
    C, D) with multiple classes in them (Aklass, A_second, etc.) need to be
    wrapped, a total of 10 different entities over 8 header files (4 hpp and 4
    pxd files). Autowrap will generate a .pxd and a .pyx file for each module.
    We decided to wrap the library into three modules, A, B and CD to show the
    capability of autowrap to do that. Note that we have perform multiple steps:
    - Step 1: parse all header files *together* - all pxd files need to be
              parsed together so that declarations are properly resolved.
    - Step 2: Map the different parsed entities to the pxd files and the
              desired modules, we use a master dict here that can be consumed
              by autowrap and specifies which pxd files and which declarations
              make up which module.
    - Step 3: Generate Cython code for each module
    - Step 4: Generate C++ code for each module (note that Step 3 has to be
              completed first before we can start to generate C++ code)
    - Step 5: Compile (run setup.py)
    Note that autowrap gives you full control how many modules you want to
    produce and which classes go into which modules. It automatically generates
    correct cimport statements in each so that dependencies between the modules
    are not an issue.
    """

    curdir = os.getcwd()
    workdir = curdir #tmpdir.strpath + "/package"
    os.makedirs(workdir)
    os.chdir(workdir)
    open("__init__.py", "a").close()

    try:

        mnames = ["constants","elements", "beamline"]

        # Step 1: parse all header files
        PY_NUM_THREADS = 1
        pxd_files = ["cppconstants.pxd", "cppstructures.pxd", "cppelements.pxd","cppcomponent.pxd", "cppparser.pxd", "cppoptimization.pxd" ,"cppbeamline.pxd"]
        full_pxd_files = [os.path.join(test_files, f) for f in pxd_files]
        decls, instance_map = autowrap.parse(
            full_pxd_files, ".", num_processes=int(PY_NUM_THREADS)
        )

        assert len(decls) == 13, len(decls)

        # Step 2: Perform mapping
        pxd_decl_mapping = {}
        for de in decls:
            tmp = pxd_decl_mapping.get(de.cpp_decl.pxd_path, [])
            tmp.append(de)
            pxd_decl_mapping[de.cpp_decl.pxd_path] = tmp

        masterDict = {}
        masterDict[mnames[0]] = {
            "decls": pxd_decl_mapping[full_pxd_files[0]],
            "addons": [],
            "files": [full_pxd_files[0]],
        }
        masterDict[mnames[1]] = {
            "decls": pxd_decl_mapping[full_pxd_files[0]]+pxd_decl_mapping[full_pxd_files[1]]+pxd_decl_mapping[full_pxd_files[2]],
            "addons": [],
            "files": [full_pxd_files[0]]+[full_pxd_files[1]]+[full_pxd_files[2]] ,
        }
        masterDict[mnames[2]] = {
            "decls":   [ pxd_decl_mapping[full_pxd_files[i]] for i in range(7)],
            "addons": [],
            "files": [full_pxd_files[i] for i in range(7) ] ,
        }

        # Step 3: Generate Cython code
        converters = []
        for modname in mnames:
            m_filename = "%s.pyx" % modname
            cimports, manual_code = autowrap.Main.collect_manual_code(
                masterDict[modname]["addons"]
            )
            autowrap.Main.register_converters(converters)
            autowrap_include_dirs = autowrap.generate_code(
                masterDict[modname]["decls"],
                instance_map,
                target=m_filename,
                debug=False,
                manual_code=manual_code,
                extra_cimports=cimports,
                include_boost=True,
                allDecl=masterDict,
                add_relative=True
            )
            masterDict[modname]["inc_dirs"] = autowrap_include_dirs
        os.chdir("..")
    finally:
        os.chdir(curdir)

test_full_lib()