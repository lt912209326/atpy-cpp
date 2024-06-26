
__all__=[ "kwd_index","tws_index","loc_index","glb_index","elem_index","rcParams" ]


cdef dict tmp_index_kwd=KEYWORDS_DICT
INDEX_KWD = {key:value.decode("utf-8") for key,value in tmp_index_kwd.items()}
KWD_INDEX = {value.decode("utf-8"):key for key,value in tmp_index_kwd.items()}


cdef dict tmp_index_tws = TWISS_DICT
INDEX_TWS = {key:value.decode("utf-8") for key,value in tmp_index_tws.items()}
TWS_INDEX = {value.decode("utf-8"):key for key,value in tmp_index_tws.items()}


cdef dict tmp_index_loc = LOCAL_DICT
INDEX_LOC = {key:value.decode("utf-8") for key,value in tmp_index_loc.items()}
LOC_INDEX = {value.decode("utf-8"):key for key,value in tmp_index_loc.items()}


cdef dict tmp_index_glb = GLOBALS_DICT
INDEX_GLB = {key:value.decode("utf-8") for key,value in tmp_index_glb.items()}
GLB_INDEX = {value.decode("utf-8"):key for key,value in tmp_index_glb.items()}


ALL_INDEX={**KWD_INDEX, **LOC_INDEX, **TWS_INDEX }

ELEM_INDEX ={"Marker":MARKER, "Drift":DRIFT ,"ExactDrift":EXACTDRIFT , "Dipole":DIPOLE ,  "Quadrupole":QUADRUPOLE ,  
            "Sextupole":SEXTUPOLE , "Octupole":OCTUPOLE ,  "Girder":GIRDER ,  "RFCavity":RFCAVITY ,  
            "Wiggler":WIGGLER ,  "Crab":CRAB, "Tuning":TUNING   }

INDEX_ELEM={value:key for key,value in ELEM_INDEX.items()}

POSITION_DEPEND_NAMES= list( KWD_INDEX.keys() ) + list( LOC_INDEX.keys() ) + list( TWS_INDEX.keys() )

DEFAULT_ELEM_KARGS={"l":0, "angle":0, "k1":0, "k2":0, "k3":0, "e1":0.5, "e2":0.5, "nslice":1, 
    "dnux":0 , "dnuy":0 , "betax1":1 , "alphax1":0 , "betay1":1 , "alphay1":0,"etax1":0,"etapx1":0,
    "betax2":1 , "alphax2":0 , "betay2":1 , "alphay2":0,"etax2":0,"etapx2":0 }


kwd_index=KWD_INDEX

tws_index=TWS_INDEX

loc_index=LOC_INDEX

glb_index=GLB_INDEX

elem_index=ELEM_INDEX

rcParams={"precision":1e-10 }


