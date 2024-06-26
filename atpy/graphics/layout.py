import numpy as np
from math import cos,sin,pi
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# __all__=["Layout, layout_datas"]

def layout_datas(ring):
    """
        return:list of [elem_type,x,y,theta,angle,length]
    """
    ring.calc()
    terms=["Gx","Gy","thetax","angle","l"]
    kinds= [ elem["kind"] for elem in ring[1:] ]
    index,optics=ring[1:, terms]
    optics=optics.tolist()
    return [ [kinds[index],*value]  for index, value in enumerate(optics) ]


class Layout:
    def __init__(self,project=""):
        self.datas=[]
        self.fig,self.ax=plt.subplots()
        self.xlim=[0,1]
        self.ylim=[0,1]
        self.ax.set_title(f"Layout of {project}" )
        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]",rotation=90)
        # self.ax.set_xticks(self.ax.get_xticks(minor=True),minor=True)



    def add_beamline(self,datas,drift_color="black",scale=1):
        """
        input:
            datas:list of [elem_type,x,y,theta,angle,length]
        """
        self.datas.append(datas)
        for index,value in enumerate(self.datas[-1]):
            elem_type,x,y,theta,angle,length=value
            if self.xlim[0]>x:
                self.xlim[0]=x
            if self.xlim[1]<x:
                self.xlim[1]=x
            if self.ylim[0]>y:
                self.ylim[0]=y
            if self.ylim[1]<y:
                self.ylim[1]=y
            if elem_type=="Marker" or elem_type=="Tuning":
                elem_color="black"
                length=0.0001
                height=0.6*scale
            elif elem_type=="Drift" or elem_type=="ExactDrift":
                elem_color=drift_color
                height=0.05*scale
            elif elem_type=="Dipole":
                elem_color="yellow"
                if (index+1)<len(self.datas[-1]):
                    theta=theta-0.5*angle
                    # 0.5*theta+0.5*self.datas[-1][index+1][3]
                height=0.5*scale
            elif elem_type=="Quadrupole":
                elem_color="blue"
                height=0.6*scale
            elif elem_type=="Sextupole":
                elem_color="lawngreen"
                height=0.6*scale
            elif elem_type=="Octupole":
                elem_color="purple"
                height=0.6*scale
            elif elem_type=="RFCavity":
                elem_color="orange"
                height=0.6*scale
            elif elem_type=="Solenoid":
                elem_color="green"
                height=0.6*scale
            else:
                raise ValueError(f'Unknown Element type {elem_type}')
            x_start,y_start=x+0.5*height*sin(theta)-length*cos(theta), y-0.5*height*cos(theta)-length*sin(theta)
            theta=theta/pi*180
            self.ax.add_patch(patches.Rectangle((x_start, y_start), length, height,angle=theta,facecolor=elem_color,linewidth=1.0,antialiased=True,alpha=1))
        self.xlim[0]-=4
        self.xlim[1]+=4
        self.ylim[0]-=4
        self.ylim[1]+=4
        self.ax.set_xlim(*self.xlim)
        self.ax.set_ylim(*self.ylim)
        #self.fig.show()
        #self.fig.savefig("stcf_lattice.png",dpi=1920)

