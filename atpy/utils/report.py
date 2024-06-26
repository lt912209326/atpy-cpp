from ..graphics.opticsplot import OpticsPlot
import numpy as np
from datetime import datetime
import os

def summary(ring, main=True):
    concerned_param = [ 'circumference','energy', 'emitx', 'U0', 'alphac','e_spread',  'Qx', 'Qy', 'taux', 'tauy', 'tauz', 'dQx', 'dQy' ]
    if main==False:
        concerned_param =concerned_param+['d2Qx', 'd2Qy','nature_chromx', 'nature_chromy', 'RI1', 'RI2', 'RI3',  'RI4', 'RI5', 'spin', 'damp_factor', 'Jx', 'Jy', 'Jz'] 
    params=ring.parse(f"{';'.join(concerned_param)}") 
    output_= "\n".join([f"{value:<20}:{params[index]:>15.6}" for index, value in enumerate(concerned_param) ])
    print(output_)
    return { key:value for key,value in zip(concerned_param,params) }

def display(ring, name=".*", terms=None):
    if terms==None:
        terms=["s","l","betax","alphax","betay","alphay", "etax","etapx","nux","nuy" ]
    index,optics=ring[name,terms ]
    names= [ elem["name"] for elem in ring[index] ]
    print(f"{'No.  '} {'Name':<13}",":",*[f"{vi:>13}" for vi in terms])
    for i,value in enumerate(optics):
        print(f"{index[i]:<5} {names[i]:<13}",":",*[f"{vi:>13.7}" for vi in value]  ) 

def nonlinear_parameters(ring, RDTs=True, chrom=True ):
    names=[]
    if chrom==True:
        names+=["END[0].chromx","END[0].chromy","dQx","dQy","d2Qx","d2Qy","d3Qx","d3Qy","d4Qx","d4Qy", ]
    if RDTs==True:
        names+=[ 'H30000', 'H10110', 'H10020', 'H10200', 'H20001', 'H00201',"dnux_dJx","dnux_dJy","dnuy_dJy", 'H21000',
                  'H22000', 'H11110', 'H00220', 'H31000', 'H40000', 'H20110', 'H11200', 'H20020', 'H20200', 'H00310', 'H00400']
    values = ring.parse(";".join(names))
    for name,value in zip(names,values):
        if abs(value)>1e6:
            print(f"{name:<20}:{value:20.6e}")
        else:
            print(f"{name:<20}:{value:20.6f}")

def save_lattice(ring, prefix="STCF",save_type=["atpy","elegant","sad","opa"] ):
    # rcParams["precision"]=1e-8
    opti = OpticsPlot(ring)
    now = datetime.now()
    month= now.strftime('%Y_%m')
    date_time = now.strftime('%m_%d_%H_%M')
    path=f'./history/{month}/{now.strftime("%Y_%m_%d")}'
    if not os.path.exists(path):
        os.makedirs(path)
    blinename=ring.pyline.name
    
    filename=f'{path}/{prefix}_{blinename}_{date_time}'
    ring.save()

    concerned_param = [ 'circumference',"emitx","Qx","Qy","dQx","dQy","d2Qx","d2Qy","dnux_dJx","dnux_dJy","dnuy_dJy","taux","tauz","alphac","RI1","H22000*2/3.1415926","H11110/3.1415926","H00220*2/3.1415926" ]
    extra=["energy", 'U0',"nature_chromx","nature_chromy", 'alphac', 'damp_factor', 'e_spread']
    
    terms=["s","betax","alphax","betay","alphay", "etax","etapx" ,"nux","nuy"]
        
    concerned_param=extra + concerned_param
    index,optics=ring[:,terms ]
    params=ring.parse(f"{';'.join(concerned_param)}") 
    output_= "\n".join([f"{value:<20}:{params[index]:>15.6}" for index, value in enumerate(concerned_param) ])
    if "atpy" in save_type:
        ring.export(f"{filename}.py","atpy")
    if "opa" in save_type:
        ring.export(f"{filename}.opa","opa")
    if "sad" in save_type:
        ring.export(f"{filename}.sad","sad")
    if "elegant" in save_type:
        ring.export(f"{filename}.lte","elegant")
    if "madx" in save_type:
        ring.pyline.export(f"{filename}.madx","madx")
    opti.draw(["betax","betay","etax","H0"],dpi=100 )
    opti.save(f"{filename}.png",dpi=300)
    with open(f"{filename}.txt","w+") as fn:
        fn.write(output_ )
    np.savetxt(f"{filename}_twiss.csv", optics, delimiter=",",header=", ".join(terms))

def export_vars(ring,file_type="atpy"):
    elem_names = set(ring.parser.vars_elem_name )
    name_elems = [(name,ring[name]) for name in elem_names]

    ordered_elems = sorted(name_elems ,key=lambda item: item[1].kind_code )

    vars_elem = [elem.str(file_type) for name,elem in ordered_elems ]
    for elem in vars_elem:print(elem,end="")
    return vars_elem