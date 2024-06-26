from math import cos,sin,pi
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib as mpl
from ..core.beamline import BeamLine
from ..core.elements import Element
from ..core.interface.constants import tws_index

class OpticsPlot:
    """
    Class for plot of optics functions
    """
    def __init__(self,line, sqrt_beta_ulimit=169):
        """
        input:
            lat_datas : lsit of [no.,name,kind,position,length,value ]
            optics_datas:list of element information [coord,alphax,betax,nux,etax,etapx,alphay,betay,nuy]
        """
        self.line=line
        self.sqrt_beta_ulimit=sqrt_beta_ulimit
        # self.input_optics_list=["coord",'alphax','betax','nux', 'etax','etapx', 'alphay','betay','nuy']
        # self.input_local_list=list(self.line["LOC"])
        self.input_lat_list = [ "no.","name","kind","position","length","value" ]
        self.functions= list(tws_index.keys())+["sigmax","sigmay"]
        
        # self.optics_datas=optics_datas.transpose()
        # self.lat_datas = list(map(list,zip(*lat_datas)))

        # optics_index={value:index for index,value in enumerate(self.input_optics_list)}
        # lat_index = {value:index for index,value in enumerate(self.input_lat_list) }
        # self.func_index= {**optics_index,**lat_index }
        self.subplots_num=1
        self.plot_functions={}
        self.fig=None
        self.ax=None
    
    def init_datas(self, start=0, end=None):
        self.optics_datas=self.line.optics(start=start,stop=end, key=self.input_optics_list ).transpose()
        # ["no.","name","kind","position","length","value"]
        index,position=self.line[start:end,"s" ]
        position=position.tolist()
        elem_datas=[elem["all"] for elem in self.line[start:end] ]
        lat_datas=[]
        #  [ "no.","name","kind","position","length","value" ]
        for i,value in enumerate(elem_datas):
            if value[1]=="Marker" or value[1]=="Tuning":
                tmp_data=[ index[i], value[0], value[1], position[i], 0, 0   ]
            elif value[1]=="Drift" or value[1]=="ExactDrift" :
                tmp_data=[ index[i], value[0], value[1], position[i], value[2], 0   ]
            elif value[1]=="Dipole":
                tmp_data=[ index[i], value[0], value[1], position[i], value[2], 0   ]
            elif value[1]=="Quadrupole":
                tmp_data=[ index[i], value[0], value[1], position[i], value[2], value[3]   ]
            elif value[1]=="Sextupole":
                tmp_data=[ index[i], value[0], value[1], position[i], value[2], value[3]   ]
            elif value[1]=="Octupole":
                tmp_data=[ index[i], value[0], value[1], position[i], value[2], value[3]   ]
            else:
                raise NotImplementedError(f"{ value[1]} is still not implented!")
            lat_datas.append(tmp_data )
        self.lat_datas =list(map(list,zip(*lat_datas)))
        

    def save(self,filename,*, dpi=300):
        self.fig.savefig(filename,dpi=dpi)


    def draw(self, plot_list,*, start=0, end=None,reverse=False, dpi=300):

        plt.rcParams["figure.dpi"]=dpi
        if len(plot_list)==0:
            raise ValueError("No function is inputted!")
        self.plot_functions={value:0 for value in plot_list}
        
        optics_list=[]
        for value in plot_list:
            if value == "sigmax":
                optics_list.append("betax")
            elif value == "sigmay":
                optics_list.append("betay")
            elif value  in self.functions:
                optics_list.append(value)
            else:
                raise ValueError(f"{value} is out of range of {self.functions} !")

        self.input_optics_list=["coord"]+list(set(optics_list))
        optics_index={value:index for index,value in enumerate(self.input_optics_list)}
        lat_index = {value:index for index,value in enumerate(self.input_lat_list) }
        self.func_index= {**optics_index,**lat_index }

        self.init_datas(start=start,end=end )

        subplots_num=1
        if "betax" in plot_list or "betay" in plot_list:
            subplots_num+=1
            self.plot_functions["betax"]=1 if "betax" in plot_list else 0
            self.plot_functions["betay"]=1 if "betay" in plot_list else 0
        if "sigmax" in plot_list or "sigmay" in plot_list:
            subplots_num+=1
            self.plot_functions["sigmax"]=1 if "sigmax" in plot_list else 0
            self.plot_functions["sigmay"]=1 if "sigmay" in plot_list else 0
        if "alphax" in plot_list or "alphay" in plot_list or "etapx" in plot_list:
            subplots_num+=1
            self.plot_functions["alphax"]=1 if "alphax" in plot_list else 0
            self.plot_functions["alphay"]=1 if "alphay" in plot_list else 0
            self.plot_functions["etapx"]=1 if "etapx" in plot_list else 0
        if "nux" in plot_list or "nuy" in plot_list:
            subplots_num+=1
            self.plot_functions["nux"]=1 if "nux" in plot_list else 0
            self.plot_functions["nuy"]=1 if "nuy" in plot_list else 0
        if "etax" in plot_list or "H0" in plot_list:
            subplots_num+=1
            self.plot_functions["etax"]=1 if "etax" in plot_list else 0
            self.plot_functions["H0"]=1 if "H0" in plot_list else 0
        if any([True if name[:2]=="lh" else False for name in plot_list ] ):
            subplots_num+=1
            for name in plot_list:
                if name[:2]=="lh":
                    self.plot_functions[name]=1

        
        if subplots_num!=self.subplots_num:
            self.subplots_num=subplots_num
            # if self.fig: self.fig.clf()
            self.fig, self.ax = plt.subplots(subplots_num, 1, 
                        gridspec_kw={
                            # 'width_ratios': [2, 1],
                            'height_ratios': (subplots_num-1)*[8]+[subplots_num-1] } )
        else:
            # if self.fig: self.fig.clf()
            self.fig, self.ax = plt.subplots(subplots_num, 1, 
                        gridspec_kw={
                            # 'width_ratios': [2, 1],
                            'height_ratios': (subplots_num-1)*[8]+[subplots_num-1] } )
        work_cell=0
        
        self.plot_lattice(-1)
        if self.plot_functions.get("betax") or self.plot_functions.get("betay") :
            self.plot_beta(work_cell)
            work_cell+=1
        if self.plot_functions.get("sigmax") or self.plot_functions.get("sigmay") :
            self.plot_beamsize(work_cell)
            work_cell+=1
        if self.plot_functions.get("etax") or self.plot_functions.get("H0")  :
            self.plot_dispersion(work_cell)
            work_cell+=1
        if self.plot_functions.get("alphax") or self.plot_functions.get("alphay") or self.plot_functions.get("etapx"):
            self.plot_alpha(work_cell)
            work_cell+=1
        if self.plot_functions.get("nux") or self.plot_functions.get("nuy") :
            self.plot_phaseadvance(work_cell)
            work_cell+=1
        if any([True if name[:2]=="lh" else False for name in plot_list ] ):
            self.plot_local_RDTs(work_cell)
            work_cell+=1
        if reverse:
            for axes in self.fig.get_axes():
                axes.invert_xaxis()
        plt.subplots_adjust(hspace=0.0)

        
        # def on_move(event):
        #     if event.inaxes in self.ax:
        #         xlim=event.inaxes.get_xlim()
        #         for ax in self.ax:
        #             if event.inaxes != ax:
        #                 ax.set_xlim(xlim)
        #     else:
        #         return
        #     self.fig.canvas.draw_idle()
        # self.fig.canvas.mpl_connect('motion_notify_event', on_move)


    def plot_lattice(self,cell_id):
        kind_id=self.func_index["kind"]
        position_id=self.func_index["position"]
        length_id=self.func_index["length"]
        value_id=self.func_index["value"]
        positions=self.lat_datas[position_id]
        lengthes=self.lat_datas[length_id]
        values=self.lat_datas[value_id]
        drift_height=0.2
        self.ax[cell_id].plot(positions,len(positions)*[0.5*drift_height],"k",alpha=0.1,linewidth=0.1)
        # self.ax[cell_id].set_facecolor("none")
        self.ax[cell_id].set_frame_on(True)
        self.ax[cell_id].patch.set_alpha(1)

        lat_elem_dict={"Marker":["black",1.0 ],"Tuning":["black",1.0 ], "Drift":[ "white",drift_height ],"ExactDrift":[ "white",drift_height ], "Dipole":["yellow",0.5 ], 
                    "Quadrupole": ["blue",0.6 ], "Sextupole":["lawngreen",0.8], "Octupole":["purple",0.8],
                    "RFCavity":["orange", 0.4], "Solenoid":["green",0.4]  }
        line_width=0.5
        self.ax[cell_id].add_patch(patches.Rectangle((positions[0], -0.5*drift_height), positions[-1]-positions[0], drift_height,edgecolor="black",linewidth=line_width ,facecolor="white") )
        for index,kind in enumerate(self.lat_datas[kind_id]):
            length=lengthes[index]
            # edgecolor0="none"
            if kind not in lat_elem_dict.keys():
                raise ValueError(f'Unknown Element type {kind}')
            elem_color,height=lat_elem_dict[kind]
            edge_color=elem_color
            x_start=positions[index]-lengthes[index]
            if kind=="Marker" or kind=="Tuning":
                # length=0.0001
                self.ax[cell_id].plot([x_start,x_start],[-0.5*height,0.5*height],color=elem_color, linewidth=line_width ,linestyle="-" )
                continue
            elif kind=="Drift" or kind=="ExactDrift" :
                continue
            elif kind=="Dipole":
                edge_color="black"
            if values[index]>0:
                y_start =-0.5*drift_height
            elif values[index]<0:
                y_start = 0.5*drift_height-height
            else:
                y_start=-0.5*height
            # y_start=-0.5*drift_height if values[index]>=0 else 0.5*drift_height-height
            self.ax[cell_id].add_patch(patches.Rectangle((x_start, y_start), length, height,edgecolor=edge_color,linewidth=line_width ,facecolor=elem_color))
        # xlim=self.ax[0].get_xlim()
        ylim=[-1,1]
        # self.ax[cell_id].set_xlim(*xlim)
        self.ax[cell_id].set_ylim(*ylim)
        self.ax[cell_id].set_xlabel(f"s [m]")
        self.ax[cell_id].minorticks_on()
        # self.ax[cell_id].xaxis.set_visible(False)
        # self.ax[cell_id].yaxis.set_visible(False)
        # self.ax[cell_id].spines['top'].set_visible(False)
        self.ax[cell_id].spines['right'].set_visible(False)
        self.ax[cell_id].spines['left'].set_visible(False)
        self.ax[cell_id].yaxis.set_ticklabels([])
        self.ax[cell_id].yaxis.set_ticks([])
        # self.ax[cell_id].axis("off")

    def plot_beta(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        max_beta=10
        if self.plot_functions["betax"]:
            max_betax=max(self.optics_datas[self.func_index["betax"]] )
            max_beta=max(max_betax,max_beta )
        if self.plot_functions["betay"]:
            max_betay=max(self.optics_datas[self.func_index["betay"]] )
            max_beta=max(max_betay,max_beta )
        if max_beta>self.sqrt_beta_ulimit:
            unit_str="[m^{1/2}]"
        else:
            unit_str="[m]"

        if self.plot_functions["betax"]:
            index=self.func_index["betax"]
            if max_beta> self.sqrt_beta_ulimit:
                self.ax[cell_id].plot(positions,np.sqrt(self.optics_datas[index]),"r",label=r"$\sqrt{\beta_x}$")
                ylabel.append(r"\sqrt{\beta_x}")
            else:
                self.ax[cell_id].plot(positions, self.optics_datas[index],"r",label=r"$\beta_x$")
                ylabel.append("\\beta_x")
        if self.plot_functions["betay"]:
            index=self.func_index["betay"]
            if max_beta> self.sqrt_beta_ulimit:
                self.ax[cell_id].plot(positions,np.sqrt(self.optics_datas[index]),"b",label=r"$\sqrt{\beta_y}$")
                ylabel.append(r"\sqrt{\beta_y}")
            else:
                self.ax[cell_id].plot(positions, self.optics_datas[index],"b",label=r"$\beta_y$")
                ylabel.append("\\beta_y")
        legend=self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel(f"${','.join(ylabel)}\,\,{unit_str}$")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        legend.get_frame().set_facecolor("none")
        # if cell_id < self.subplots_num-1:
        #     self.ax[cell_id].xaxis.set_ticklabels([])
    

    def plot_beamsize(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        max_beta=10
        if self.plot_functions["sigmax"]:
            max_betax=max(self.optics_datas[self.func_index["betax"]] )
            max_beta=max(max_betax,max_beta )
        if self.plot_functions["sigmay"]:
            max_betay=max(self.optics_datas[self.func_index["betay"]] )
            max_beta=max(max_betay,max_beta )
        unit_str="[m]"

        if self.plot_functions["sigmax"]:
            index=self.func_index["betax"]
            self.ax[cell_id].plot(positions, np.sqrt(self.line["emitx"]*self.optics_datas[index]),"r",label=r"$\sigma_x$")
            ylabel.append("\\sigma_x")
        if self.plot_functions["sigmay"]:
            index=self.func_index["betay"] 
            self.ax[cell_id].plot(positions, np.sqrt(self.line["STAT"]["mincouple"]*self.line["emitx"]*self.optics_datas[index]),"b",label=r"$\sigma_y$")
            ylabel.append("\\sigma_y")
        legend=self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel(f"${','.join(ylabel)}\,\,{unit_str}$")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        legend.get_frame().set_facecolor("none")


    def plot_dispersion(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        unit_str="[m]"
        if self.plot_functions["etax"]:
            index=self.func_index["etax"]
            self.ax[cell_id].plot(positions,self.optics_datas[index],"c",label=r"$\eta_x$")
            ylabel.append("\\eta_x")
        if self.plot_functions["H0"]:
            index=self.func_index["H0"]
            self.ax[cell_id].plot(positions,10.0*self.optics_datas[index],"c--",label=r"$10 \mathcal{H}$")
            ylabel.append("10 \\mathcal{H}")
        self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel(f"${','.join(ylabel)}\,\,{unit_str}$")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        legend=self.ax[cell_id].legend(loc="upper right")
        legend.get_frame().set_facecolor("none")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        # if cell_id < self.subplots_num-1:
        #     self.ax[cell_id].xaxis.set_ticklabels([])

    def plot_alpha(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        if self.plot_functions["alphax"]:
            index=self.func_index["alphax"]
            
            self.ax[cell_id].plot(positions,self.optics_datas[index],"r",label=r"$\alpha_x$")
            ylabel.append("\\alpha_x")
        if self.plot_functions["alphay"]:
            index=self.func_index["alphay"]
            self.ax[cell_id].plot(positions,self.optics_datas[index],"b",label=r"$\alpha_y$")
            ylabel.append("\\alpha_y")
        if self.plot_functions["etapx"]:
            index=self.func_index["etapx"]
            self.ax[cell_id].plot(positions,self.optics_datas[index],"c--",label=r"$\eta_x^{\prime}$")
            ylabel.append("\\eta_x^{\prime}")
        self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel(f"${','.join(ylabel)}$")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        legend=self.ax[cell_id].legend(loc="upper right")
        legend.get_frame().set_facecolor("none")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        # if cell_id < self.subplots_num-1:
        #     self.ax[cell_id].xaxis.set_ticklabels([])

    def plot_phaseadvance(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        if self.plot_functions["nux"]:
            index=self.func_index["nux"]
            self.ax[cell_id].plot(positions,self.optics_datas[index],"r",label=r"$\nu_x$")
            ylabel.append("\\nu_x")
        if self.plot_functions["nuy"]:
            index=self.func_index["nuy"]
            self.ax[cell_id].plot(positions,self.optics_datas[index],"b",label=r"$\nu_y$")
            ylabel.append("\\nu_y")
        self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel(f"${','.join(ylabel)}$")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        legend=self.ax[cell_id].legend(loc="upper right")
        legend.get_frame().set_facecolor("none")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        # if cell_id < self.subplots_num-1:
        #     self.ax[cell_id].xaxis.set_ticklabels([])

    def plot_local_RDTs(self,cell_id):
        ylabel=[]
        positions=self.optics_datas[self.func_index["coord"]]
        color_map=[ "r", "b", "y", "g", "c","r--", "b--", "y--", "g--", "c--"]
        color_index=0
        for index,name in enumerate(self.input_optics_list):
            if name[:2]!="lh":
                continue
            index=self.func_index[name ]
            self.ax[cell_id].plot(positions,self.optics_datas[index],color_map[color_index],label=name )
            ylabel.append(name )
            color_index+=1
        self.ax[cell_id].legend(loc="upper right")
        self.ax[cell_id].set_ylabel("Local RDTs")
        self.ax[cell_id].tick_params(axis="x",which="both",direction="in")
        legend=self.ax[cell_id].legend(loc="upper right")
        legend.get_frame().set_facecolor("none")
        self.ax[cell_id].grid(visible=True)
        self.ax[cell_id].sharex(self.ax[-1] )
        # if cell_id < self.subplots_num-1:
        #     self.ax[cell_id].xaxis.set_ticklabels([])