
from math import pi,log10,cos
import re
from .interface.constants import rcParams
# from constants import  KWD_INDEX, INDEX_KWD, TWS_INDEX, LOC_INDEX, GLB_INDEX, ELEM_INDEX, DEFAULT_ELEM_KARGS

# export_hash ={
#     "elegant":{L:"L", ANGLE:"ANGLE", K1:"K1", K2:"K2",K3:"K3", E1:"e1", E2:"e2",NSLICE:"N_SLICES" },

# }


def generate_tuners(ring, names):
    index,optics=ring[names,["s","betax"]  ]
    index=index if isinstance(index,list) else [index]
    index = [ i-1 if i>0 else i for i in index]
    return [generate_tuner(ring,i) for i in index]
    
def generate_tuner(ring, position):
    args = ["betax","alphax","betay","alphay","etax","etapx"]
    index,optics=ring[position,args]
    args1=','.join( [f"{arg}1={optics[ii]}" for ii,arg in enumerate(args)] )
    args2=','.join( [f"{arg}2={optics[ii]}" for ii,arg in enumerate(args)] )
    return f"Tuning(dnux=0, dnuy=0, {args1 },\n\t\t{args2} )"


#("L", "ANGLE" ,"K1" ,"K2"         ,"K3"         ,"E1"         ,"E2"         ,"NSLICE"     ,"TILT"       )

# print(KWD_INDEX )
# print(INDEX_KWD )
# print(TWS_INDEX )
# print(LOC_INDEX )
#print(Line)
#print(Element)

cdef class Line:
#    cdef:
#        readonly str name
#        readonly list line
#        readonly list expand
#        readonly list reverse
#        readonly list unit_cell
#        readonly dict elems
#        readonly dict lines
#        readonly list ordered_lines
#        readonly 
    
    def __init__(self,str name,*args, bint reverse=False ):
        self.name=name
        self.line=[]
        self.expand=[]
        self.reverse=[]
        self.elems={}
        # self.unit_cell=[]
        self.ordered_lines=[]

        # lines=args[::-1] if reverse else args
        for index,arg in enumerate(args):
            if isinstance(arg,Element):
                if arg.name in ("START","END") and (index != 0 and index != len(args)-1): 
                    raise ValueError(f"START and END are reserved for the first and last Marker element, but {index}-th used the {arg.name}!")
                # self.unit_cell.append(arg)
                self.line.append(arg)
                self.elems[arg.name]=arg

                self.expand.append(arg)
                self.reverse.append(False)
            elif isinstance(arg,Line):
                # when a line is used by BeamLine, START and END will be inserted. to reuse the line, delete these two Elements
                arg.renew()
                self.elems={**self.elems,**arg.elems}
                self.ordered_lines+=arg.ordered_lines
                self.line.append(arg )
                self.expand+=<Line>arg.expand
                self.reverse+=<Line>arg.reverse
            elif arg is None:
                raise RuntimeWarning(f"{index}-th component in {self.name} is None!")
            else:
                raise ValueError(f"{name}: {index}-th element is invalid!")
        tmp_lines=self.ordered_lines
        self.ordered_lines = list(set(self.ordered_lines) )
        self.ordered_lines.sort(key=tmp_lines.index )
        self.ordered_lines.append(self)
        

    def renew(self):
        """
        remove the START and END Marker inserted by BeamLine
        """
        if self.expand[0].name == "START": 
            self.expand.pop(0)
            self.reverse.pop(0)
        if self.expand[-1].name == "END": 
            self.expand.pop(-1)
            self.reverse.pop(-1)

    
    def __rmul__(self,int other):
    
        cdef Line newline=Line.__new__(self.__class__ )
        cdef list reverse,expand
        newline.name=f"{other}*{self.name}"
        self.renew()
        if other<0:
            reverse=-other*[not value for value in self.reverse[::-1] ]
            expand=-other*self.expand[::-1]
        elif other>0:
            reverse=other*self.reverse
            expand=other*self.expand
        else:
            raise RuntimeWarning("0*{self.name} is not suggested!")
           # return None
        newline.reverse=reverse
        newline.expand =expand
        newline.elems=self.elems 

        newline.ordered_lines=self.ordered_lines
        if self.name[0] in "+-0123456789":
            newline.line=self.line
        else:
            newline.line=[self]
        return newline



        
    def __neg__(self):
        cdef Line newline=Line.__new__(self.__class__ )
        newline.name=f"-{self.name}"
        self.renew()
        newline.reverse=[not value for value in self.reverse[::-1] ]
        newline.expand =self.expand[::-1]
        newline.ordered_lines=self.ordered_lines
        newline.elems=self.elems 
        if self.name[0] in "+-0123456789":
            newline.line=self.line
        else:
            newline.line=[self]
        return newline

    
    def __repr__(self):
        return f"<class '{self.__class__.__name__}': {self.name}>"

    def __str__(self):
        return f'{self.name}=Line("{self.name}",{", ".join([line.name for line in self.line ])})'


    
    def export(self,str filename,str filetype="atpy" ):
        # elems,lines=self.get()
        num_component=len(self.expand)
        ordered_elems = sorted(self.elems.items() ,key=lambda item: item[1].kind_code+self.expand.index(item[1])/num_component )
        with open(filename, "a+") as fn:
            for elem_name,elem in ordered_elems:
                fn.write(elem.str(filetype) )
            for line in self.ordered_lines:
                fn.write(line.str(filetype) )
        
    def str(self,str filetype):
        elems=[value.name for value in self.line ]
        width=4*(len(self.name)//4+1)
        
        if filetype=="atpy":
            long_str = f"{self.name:{width}} = Line({', '.join(elems)})\n"
        elif filetype=="sad":
            long_str = f"LINE  {self.name:{width}}=({' '.join(elems)});\n"
        elif filetype=="opa":
            long_str = f"{self.name:{width}} : {', '.join(elems)};\r"
        elif filetype=="madx":
            long_str = f"{self.name:{width}} : line=({', '.join(elems)});\n"
        elif filetype=="elegant":
            long_str = f"{self.name:{width}} : line=({', '.join(elems)})\n"
        else:
            raise ValueError(f"Unexpected Type {filetype}!")
        if filetype=="opa":
            return "\r    ".join(re.findall(r".{131}\w*[,;) \t\n]+|.{1,131};", long_str ) )+"\r\n"
        elif filetype=="elegant":
            return "&\n\t\t".join(re.findall(r".{131}\w*[,); \t\n]+|.{1,131}[;)]", long_str ) )+"\n"
        else:
            return "\n\t\t".join(re.findall(r".{131}\w*[,); \t\n]+|.{1,131}[;)]", long_str ) )+"\n"


    

# print(KWD_INDEX )
# print(TWS_INDEX )
# print(INDEX_KWD )
# print(Line)
# print(Element)




cdef class Element:
#    cdef:
#        CppElement* elem 
#        string name
#        bint       owner
#        str        kind
#        list       eids
    def __cinit__(self,str name0,*, **kargs):
        if name0 is not None:
            self.name = name0
        else:
            raise ValueError("Element name is needed!")
        cdef string name=<string>name0.encode("utf8")
        cdef str kind=self.__class__.__name__
        self.kind=kind
        cdef size_t num_kargs= len(kargs)
        cdef dict _kargs={**DEFAULT_ELEM_KARGS,**kargs}
        _float_nslice = _kargs["nslice"]%1

        if _kargs["nslice"]==0 or _float_nslice !=0 :
            raise ValueError("nslice should be int no less than 1!")
        if kind == "Element":
            raise  TypeError("Can't instantiate abstract class Element")
        elif kind=="Marker":
            self.elem = new CppMarker(name)
        elif kind=="Drift":
            self.elem = new CppDrift(name,_kargs["l"],_kargs["nslice"])
        elif kind=="ExactDrift":
            self.elem = new CppExactDrift(name,_kargs["l"],_kargs["nslice"])

        elif kind=="Dipole":
            self.elem = new CppDipole(name,_kargs["l"], _kargs["angle"], _kargs["k1"], _kargs["e1"], _kargs["e2"], _kargs["nslice"])
        elif kind=="Quadrupole":
            self.elem = new CppQuadrupole(name,_kargs["l"], _kargs["k1"], _kargs["nslice"])
        elif kind=="Sextupole":
            self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        elif kind=="Octupole":
            self.elem = new CppOctupole(name,_kargs["l"], _kargs["k3"], _kargs["nslice"])
            #self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        elif kind=="Tuning":
            self.elem = new CppTuning(name, _kargs["dnux"] , _kargs["dnuy"] , 
                    _kargs["betax1"] , _kargs["betay1"], _kargs["alphax1"]  , _kargs["alphay1"], _kargs["etax1"], _kargs["etapx1"],
                    _kargs["betax2"] , _kargs["betay2"], _kargs["alphax2"]  , _kargs["alphay2"], _kargs["etax2"], _kargs["etapx2"] )
        elif kind=="Multipole":
            pass
            #self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        elif kind=="Solenoid":
            pass
            #self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        elif kind=="RFCavity":
            pass
            #self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        elif kind=="Corrector":
            pass
            #self.elem = new CppSextupole(name,_kargs["l"], _kargs["k2"], _kargs["nslice"])
        
        self.keywords=self.elem.keywords
        self.kind_code=self.elem.kind

        for key in kargs.keys():
            if KWD_INDEX[key] not in self.keywords:
                raise ValueError(f"{key} is not a valid argument for {kind}!")

            
        
    def __init__(self,*arg, **kargs):
        pass
        # self.owner = True
    
    
    def __dealloc__(self):
        del self.elem 
        self.elem=NULL

    
    def __rmul__(self,int other):
        cdef Line newline = Line.__new__(Line )
        newline.name=f"{other}*{self.name}"
        newline.elems={self.name:self}
        newline.line=[self]
        newline.ordered_lines=[]

        if other<0:
            newline.expand = [self for i in range(-other)]
            newline.reverse= [True for i in range(-other)]
        elif other >0:
            newline.expand = [self for i in range(other)]
            newline.reverse= [False for i in range(other)]
        else:
            raise RuntimeWarning("0*{self.name} is not suggested!")
           # return None
        return newline 

    def __neg__(self):
        cdef Line newline = Line.__new__(Line )
        newline.name=f"-{self.name}"
        newline.elems = {self.name : self}
        #newline.elems[self.name]=self
        newline.line=[self]
        newline.ordered_lines=[]
        newline.expand=[self]
        newline.reverse=[True]
        return newline

    
    def __repr__(self):
        return f"<class '{self.__class__.__name__}': {self.name}>"

    def __str__(self):
        return self.str("atpy")
        
    def __getitem__(self,index not None):
        cdef dict kwds=self.elem.values
        if isinstance(index,int):
            if 0<=index<len(self.eids):
                return self.eids[index]
            elif -len(self.eids) <=index<0:
                return self.eids[index]
            else:
                raise ValueError(f'{self.name}: Index is out of Range!')
        elif index in KWD_INDEX.keys():
            return kwds[ KWD_INDEX[index] ]
        elif index == 'name':
            return self.name
        elif index == 'kind':
            return self.kind
        elif index == 'keywords':
            return [kwds[i] for i in self.keywords]
        elif index == "all":
            return [self.name,self.kind ] + [kwds[i] for i in self.keywords]
        else:
            raise ValueError(f'{self.name}: Input arg is invalid!')

    def str(self, str filetype):
        """
        根据不同的文件类型，将元素信息转换为对应的字符串格式。

        参数:
        - self: 对象自身引用。
        - filetype: 字符串，指定要转换的文件类型，可以是"atpy"、"sad"、"opa"、"madx"或"elegant"。

        返回值:
        - 根据指定的文件类型，返回相应格式的字符串，描述了元素的类型、名称和参数。
        """

        # 提取关键字字典
        cdef dict kwds=self.elem.values
        cdef dict kind2others
        cdef dict kwds2others

        if filetype=="atpy":
            # 为ATPY格式构建参数字符串
            params=[f"{INDEX_KWD[index]}={kwds[index]:.{int(-log10(rcParams['precision']))}f}" for index in self.keywords if index != NSLICE or (index == NSLICE and kwds[index] !=1 ) ]
            return f"{self.name:{4*(len(self.name)//4+1)}} = {self.kind}({', '.join(params)})\n"

        elif filetype=="sad":
            # 为SAD格式构建参数字符串
            kind2others ={"Marker":"MARK", "Drift":"DRIFT","ExactDrift":"DRIFT", "Dipole":"BEND", 
                        "Quadrupole":"QUAD", "Sextupole":"SEXT", "Octupole":"OCT", "Multipole":"MULT", 
                        "RFCavity":"CAVI", "Solenoid":"SOL", "Tuning":"MAP"  }
            kwds2others ={L:"L", ANGLE:"ANGLE", K1:"K1", K2:"K2",K3:"K3", E1:"E1", E2:"E2" }
            params=[]
            if self.kind in ("Multipole","Solenoid","RFCavity","Corrector"):
                raise TypeError(f"{self.kind} can not be translated to {filetype}!")
            for index in self.keywords:
                # 根据关键字构建参数字符串
                if index in (K1,K2):
                    param=f"{kwds2others[index]}={kwds[index]*kwds[L]}"
                elif index!=NSLICE:
                    param=f"{kwds2others[index]}={kwds[index]}"
                else:
                    continue
                params.append(param)
            if self.kind not in ("Marker","Map","Drift","ExactDrift"):
                params.append(" DISFRIN=1 ")
            return f"{kind2others[self.kind]:8}  {self.name:<{4*(len(self.name)//4+1)}}=({'  '.join(params)});\n"

        elif filetype=="opa":
            # 为OPA格式构建参数字符串
            kind2others={"Marker":"OpticsMarker", "Drift":"Drift","ExactDrift":"Drift", 
                        "Dipole":"Bending", "Quadrupole":"Quadrupole", "Sextupole":"Sextupole", 
                        "Octupole":"Multipole", "Multipole":"Multipole", "Tuning":"OpticsMarker" }
            kwds2others ={L:"L", ANGLE:"T", K1:"K", K2:"K" , K3:"K" , E1:"T1", E2:"T2" ,NSLICE:"N"}
            params=[]
            if self.kind in ("Multipole","Solenoid","RFCavity","Corrector"):
                raise TypeError(f"{self.kind} can not be translated to {filetype}!")
            for index in self.keywords:
                # 根据关键字构建参数字符串
                if index == K2:
                    param=f"{kwds2others[index]}={0.5*kwds[index]:.6f}"
                elif index == ANGLE:
                    param=f"{kwds2others[index]}={kwds[index]*180/pi:.6f}"
                elif index in (E1,E2):
                    param=f"{kwds2others[index]}={kwds[index]*kwds[ANGLE]*180/pi:.6f}"
                elif index == K3:
                    param=f"{kwds2others[index]}={0.16666666667*kwds[index]:.6f}"
                elif index==NSLICE and kwds[index]==1:
                    continue
                else:
                    param=f"{kwds2others[index]}={kwds[index]:.6f}"
                params.append(param)
            return f"{self.name:{ 4*(len(self.name)//4+1)}} : {', '.join( [kind2others[self.kind] ] + params)};\r"

        elif filetype=="madx":
            # 为MADX格式构建参数字符串
            kind2others={"Marker":"Marker", "Drift":"Drift","ExactDrift":"Drift", "Dipole":"sbend", "Quadrupole":"Quadrupole", "Sextupole":"Sextupole", "Octupole":"Octupole", "Multipole":"Multipole" }
            kwds2others ={L:"L", ANGLE:"ANGLE", K1:"K1", K2:"K2",K3:"K3", E1:"e1", E2:"e2",NSLICE:"N_SLICES" }
            params=[]
            if self.kind in ("Multipole","Solenoid","RFCavity","Corrector"):
                raise TypeError(f"{self.kind} can not be translated to {filetype}!")
            for index in self.keywords:
                # 根据关键字构建参数字符串
                if index in (E1,E2):
                    param=f"{kwds2others[index]}={kwds[index]*kwds[ANGLE] }"
                elif index == NSLICE and kwds[index]==1:
                    continue
                else:
                    param=f"{kwds2others[index]}={kwds[index]}"
                params.append(param)
            return f"{self.name:{4*(len(self.name)//4+1)}} : {', '.join( [kind2others[self.kind] ] + params)};\n"

        elif filetype=="elegant":
            # 为ELEGANT格式构建参数字符串
            kind2others={"Marker":"Marker", "Drift":"Drift","ExactDrift":"EDrift", "Dipole":"csbend", 
                        "Quadrupole":"KQUAD", "Sextupole":"KSEXT", "Octupole":"KOCT", "Multipole":"Multipole","Tuning":"Marker" }
            kwds2others ={L:"L", ANGLE:"ANGLE", K1:"K1", K2:"K2",K3:"K3", E1:"e1", E2:"e2",NSLICE:"N_SLICES" }
            params=[]
            if self.kind in ("Multipole","Solenoid","RFCavity","Corrector"):
                raise TypeError(f"{self.kind} can not be translated to {filetype}!")
            for index in self.keywords:
                # 根据关键字构建参数字符串
                if index in (E1,E2):
                    param=f"{kwds2others[index]}={kwds[index]*kwds[ANGLE] }"
                elif index == NSLICE and kwds[index]==1:
                    continue
                else:
                    param=f"{kwds2others[index]}={kwds[index]}"
                params.append(param)
            if self.kind in ("Marker","Drift","ExactDrift"):
                return f"{self.name:{4*(len(self.name)//4+1)}} : {', '.join( [kind2others[self.kind] ] + params)}\n"
            else:
                return f"{self.name:{4*(len(self.name)//4+1)}} : {', '.join( [kind2others[self.kind] ] + params)},SYNCH_RAD=1,ISR=1,N_SLICES=10\n"

cdef class Marker(Element):
    pass


cdef class Drift(Element):
    pass

cdef class ExactDrift(Element):
    pass

cdef class Dipole(Element):
    pass


cdef class Quadrupole(Element):
    pass


cdef class Sextupole(Element):
    pass


cdef class Octupole(Element):
    pass



cdef class Multipole(Element):
    pass


cdef class Girder(Element):
    pass


cdef class RFCavity(Element):
    pass


cdef class Wiggler(Line):
    def __init__(self,name,l,*,lambda0=0.2, Bmax=0, nslices=20, brho=20.0/3.0):
        m=l/lambda0
        l_slice=lambda0/nslices
        Bis=[nslices*Bmax*(cos(2*pi*(i-1)/nslices ) - cos(2*pi*i/nslices ))/(2*pi) for i in range(nslices) ]
        Ais =[Bi*l_slice/brho for Bi in Bis]
        x0p = 0
        A1s=[]
        A2s=[]
        A1=x0p
        A2=0
        for i in range(nslices):
            A1s.append( A1 )
            A2=Ais[i]-A1
            A2s.append( A2 )
            A1=-A2
        e1s = [ A1s[i]/Ais[i]  for i in range(nslices)]
        e2s = [ A2s[i]/Ais[i]  for i in range(nslices)]
        print(f"number of period: ",m)
        print(f"Max By : ",Bmax)
        print("angles: ",*Ais)
        seq= int(m)*[Dipole(f"{name}{i:0>2}",l=l_slice,angle=Ais[i],e1=e1s[i],e2=e2s[i]) for i in range(nslices) ] 
        super(Wiggler, self).__init__(name,*seq )

    pass
