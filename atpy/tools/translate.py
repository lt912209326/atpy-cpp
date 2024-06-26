
# __all__=["translate"]

def translate(lat):
    import re
    LATTICE=re.sub(r"^[ \t]+",r"",lat,flags=re.M|re.X)
    ELEM_TYPES=["Marker", "Drift","ExactDrift", "Dipole", "Quadrupole", "Sextupole", "Octupole","Multipole","Solenoid",
        "RFCavity","Corrector","Wiggler", "Line","BeamLine","Tuning"]
    elem_token='|'.join([ rf"{elem}[ \t]*\(" for elem in ELEM_TYPES ])
    # print(elem_names)
    lattice=re.sub(rf"([A-Za-z_][A-Za-z0-9_]*)[ \t]*=[ \t]*({elem_token})",r"\1=\2'\1',",LATTICE)
    
    import os
    import sys
    import traceback
    import tempfile
    #将修改后得元件定义字符串写入临时文件，因为exec() 直接执行字符串不能 traceback错误具体位置
    path=tempfile.gettempdir()
    if path not in sys.path:
        sys.path.append(path)
    fn,filename=tempfile.mkstemp(suffix=".py")
    outsock = os.fdopen(fn,'w')
    outsock.write(lattice)
    outsock.close()
    file=os.path.basename(filename)
    try:
        exec(f"from {file[:-3]} import*",{},globals() )
    except Exception:
        traceback.print_exc()
    finally:
        os.remove(filename)
        return lattice 