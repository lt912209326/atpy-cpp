
#[ "misaligment", "radiation", "fluctuation", "fringe", "edge", "lossmap", "fma", "slice", "combineddipole", "computedrivingterms", "leaderordertermonly",
#     "nonlineartermonly", "linear", "period", "expand", 
#     "npara",
#      "dp0", "dp"
#]

__all__=["Status", "default_status" ]

PARTICLES_MAP={"electron":ELECTRON,"positron":POSITRON, "proton":PROTON}

DEFAULT_STATUS=STAT0

STAT_INDEPEND_FIELDS=list(DEFAULT_STATUS.keys() )

default_status=DEFAULT_STATUS

cdef class Status:
    def __init__(self, **kargs):
        for key,value in kargs.items():
            if key not in STAT_INDEPEND_FIELDS:
                raise ValueError(f"Unexpected key {key}!")
        self._dict= kargs 
        kargs_keys = list(kargs.keys() )

        self.stat=STAT0

        self.bool_type_fields=["misaligment","radiation", "fluctuation", "fringe","edge" , "lossmap" ,"fma" , "slice" , "combineddipole", "computedrivingterms",  
                                "leaderordertermonly","fast_2nd_order_RDTs", "nonlineartermonly", "linear" , "period" , "expand", "second_order_chrom" , "third_order_chrom", "printout",
                                "transfermatrix", "lazy_compute","rdt_fluctuation","local_twiss", "off_momentum_rdts"]
        self.int_type_fields= ["npara"  , "particle" , "nperiods","track_lines", "track_turns","NP","max_da_range","chrom_refpt" ]

        self.float_type_fields= [ "dp0", "dp" ,"rf_dp", "energy", "monitor_dp", "mincouple","larger_monitor_dp","max_betax", "max_etax","off_rdts_observer" ]

        self.readonly_fields = ["Trx" , "Try" ,"totalslice" , "multipoleslice"]

        for key,value in kargs.items():
            self[key]=value
        
    def __setitem__(self,keys,values):

        if isinstance(keys,str):
            if keys  in self.bool_type_fields:
                if not isinstance(values,bool):raise ValueError(f"Value for {keys} should be bool type, not {values}!")
            elif keys  in self.int_type_fields:
                if not isinstance(values,int):raise ValueError(f"Value for {keys} should be int type, not {values}!")
            elif keys  in self.float_type_fields:
                if not isinstance(values,float):raise ValueError(f"Value for {keys} should be float type, not {values}!")
            else:
                raise ValueError(f"Unexpected key {keys}!")
            if "slice"==keys:
                self.stat.slice=values
            elif "lossmap"==keys:
                self.stat.lossmap=values
            elif "fma"==keys:
                self.stat.fma=values
            elif "edge"==keys:
                self.stat.edge=values
            elif "fringe"==keys:
                self.stat.fringe=values
            elif "fluctuation"==keys:
                self.stat.fluctuation=values
            elif "radiation"==keys:
                self.stat.radiation=values
            elif "misaligment"==keys:
                self.stat.misaligment=values
                if values ==True:
                    self.stat.nonlineartermonly=False
                    self.stat.expand=True
            elif "combineddipole"==keys:
                self.stat.combineddipole=values
            elif "computedrivingterms"==keys:
                # raise ValueError("Changing of computedrivingterms is not implement,which need reinitialization of srdt cache!")
                self.stat.computedrivingterms=values
                if self.stat.period !=True:
                    print("warning: period might be True?")
            elif "leaderordertermonly"==keys:
                self.stat.leaderordertermonly=values
            elif "fast_2nd_order_RDTs"==keys:
                self.stat.fast_2nd_order_RDTs=values
            elif "lazy_compute"==keys:
                self.stat.lazy_compute=values
            elif "nonlineartermonly"==keys:
                self.stat.nonlineartermonly=values
                if self.stat.period !=True:
                    print("warning: period might be True?")
            elif "linear"==keys:
                self.stat.linear=values
            elif "period"==keys:
                self.stat.period=values
                if values ==False:
                    print("warning: nonlineartermonly and computedrivingterms might be True?")
                    #self.stat.nonlineartermonly=False
                    #self.stat.computedrivingterms=False
            elif "expand"==keys:
                self.stat.expand=values
            elif "second_order_chrom"==keys:
                self.stat.second_order_chrom=values
                if values==False:
                    self.stat.third_order_chrom=values
            elif "third_order_chrom"==keys:
                self.stat.third_order_chrom=values
                if values==True:
                    self.stat.second_order_chrom=True
            elif "printout"==keys:
                self.stat.printout=values
            elif "rdt_fluctuation"==keys:
                self.stat.rdt_fluctuation=values
            elif "local_twiss" ==keys:
                self.stat.local_twiss = values
            elif "off_momentum_rdts" ==keys:
                self.stat.off_momentum_rdts = values
            elif "off_rdts_observer" ==keys:
                self.stat.off_rdts_observer = values
            elif "npara"==keys:
                self.stat.npara=values
            elif "NP"==keys:
                self.stat.NP=values
            elif "max_da_range"==keys:
                self.stat.max_da_range=values
            elif "chrom_refpt"==keys:
                self.stat.chrom_refpt=values                
            elif "particle"==keys:
                self.stat.particle=values
            elif "nperiods"==keys:
                self.stat.nperiods=values
            elif "track_lines"==keys:
                self.stat.track_lines=values
            elif "track_turns"==keys:
                self.stat.track_turns=values
            elif "dp0"==keys:
                self.stat.dp0=values
            elif "dp"==keys:
                self.stat.dp=values
            elif "rf_dp"==keys:
                self.stat.rf_dp=values
            elif "monitor_dp"==keys:
                self.stat.monitor_dp=values
            elif "larger_monitor_dp"==keys:
                self.stat.larger_monitor_dp=values
            elif "max_betax"==keys:
                self.stat.max_betax=values
            elif "max_etax"==keys:
                self.stat.max_etax=values
            elif "energy"==keys:
                self.stat.energy=values
            elif "mincouple"== keys:
                self.stat.mincouple=values
            elif "transfermatrix"==keys:
                self.stat.transfermatrix=values
            else:
                raise ValueError(f"Unexpected key {keys}!")
        elif isinstance(keys,list):
            if not isinstance(values,list) or len(keys) != len(values):
                raise ValueError(f"Values should be correspond to keys!")
            for key,value in zip(keys,values):
                if not isinstance(key,str):
                    raise ValueError(f"Unexpected key {key}, keys should be str type or list of str type!")
                self[key]=value
        else:
            raise ValueError(f"Unexpected keys {keys}, keys should be str type or list of str type!")
        self._dict=self.stat

    def export(self, str filename, str filetype="atpy"):
        with open(filename,"a+") as fn:
            fn.write(self.str(filetype) )

    def str(self,str filetype="atpy"):
        if filetype=="atpy":
            kargs_str= [f"{key}={value}" for key,value in (self._dict.items()-DEFAULT_STATUS.items()  ) if key not in self.readonly_fields ]
            return f"stat=Status( {', '.join(kargs_str)} )"
        elif filetype=="sad":
            pass
        elif filetype=="opa":
            pass
        elif filetype=="madx":
            pass

        
    def fields(self):
        return f"{', '.join([ f'{type(value).__name__}:{key}' for key,value in default_status.items() ])}"
