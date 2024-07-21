
import re
import numpy as np
import math
from cython.operator cimport dereference as deref
cimport cython
from cython.parallel cimport prange
cimport openmp
from libc.string cimport memcpy

##@cython.no_gc_clear
cdef class BeamLine:
    """
    methods:
    eval(self, str expr)
    parse(self,str code)
    save(self)              :
    evolution(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV )
    calc(self)
    findclosedorbit(self, double dp)
    highorderchromaticity(self,double dp0=0.0001 )
    compute_large_off_momentum_tunes(self)
    correctchrom(self, dQx=None, dQy =None)
    compute_off_momentum_twiss(self, list dp_range, double dp_step=1e-4, local_twiss=True)
    compute_off_momentum_RDTs(self)
    track(self, double[:,::1] beam0, int start_pos=0, int end_pos=-1, int nturn0=1,int nturn1=1)
    display(self,str token,bint detail=False)
    export(self,str filename, str filetype="atpy")
    str(self,str filetype="atpy")
    """
#    cdef:
#        CppBeamLine lat
#        readonly Parse parser
#        readonly Line pyline
#        readonly Status stat
#        readonly dict elems_index
    
    def __cinit__(self,str name0, Status stat0, arg,**kargs):
        self.name=name0
        self.kind=self.__class__.__name__
        cdef string name=<string>(self.name.encode("utf8") )
        self.stat=stat0
        self.elems_index={}
        # print("python start __cinit__:")
        self.lat=NULL
        self.parser=None
        self.lat_pool.resize(0)
        self.parser_pool = []
        self.current_worker = -1
        self.initial_worker = -1
        # print("python start __cinit__:")
        self.nkernel=1
        # cdef size_t num_args=len(args)
        # cdef Elements elem
        cdef size_t length0
        # print("python start __cinit__:")
        if isinstance(arg,str):
            pass
        elif isinstance(arg,Line):
            self.pyline=arg
            if self.pyline.expand[-1].name !="END":
                self.pyline.expand.append(Marker("END") )
                self.pyline.reverse.append( False)
            if self.pyline.expand[0].name !="START":
                self.pyline.expand.insert(0,Marker("START") )
                self.pyline.reverse.insert(0, False)
            length0=len(self.pyline.expand)-1
            self.lat=new CppBeamLine(name, self.stat.stat.particle, self.stat.stat.energy,length0,self.stat.stat)
            # self.elems_index["START"]=[0]
            # print("__cinit__:",f"{<size_t>self.lat:x}")
        else:
            raise ValueError(f"Third arg is invalid!")
        # print("python start __cinit__:")
    

    def __init__(self, str name, Status, arg,**kargs):
        cdef Element value
        cdef size_t index,position
        # print("before add element: ")
        # print(self.pyline.expand)
        #print()
        for index,value in enumerate(self.pyline.expand ):
            # 在CppBeamLine中有默认起始Marker元素 "START"
            position=index
            if value.name not in self.elems_index.keys():
                (<Element?>value).eids=[position]
                self.elems_index[value.name] = [position ]
            else:
                # print(value.name,": ",index,",",value)
                (<Element?>value).eids.append(position)
                self.elems_index[value.name].append(position)
        # first component "START" is already in lat
        for key,value0 in kargs.items():
            if key in TWS_INDEX.keys():
                # print(key,":",TWS_INDEX[key] )
                (&(self.lat.line[0].tws[0]) )[TWS_INDEX[key] ]=<double?>value0
            elif key in ("Gx","Gy","Gz","thetax", "thetay","thetaz"):
                self.lat.line[0].local[ LOC_INDEX[key] ]=<double?>value0
            else:
                raise ValueError("Unexpected initial optics or geometry parameter {key}!")
        self.parser=self.start_cppbeamline(self.lat )

    cdef Parser start_cppbeamline(self, CppBeamLine* lat ):
        cdef Element value
        cdef size_t index,position
        cdef Parser parser=Parser(self.elems_index )
        for index,value in enumerate(self.pyline.expand ):
            # 在CppBeamLine中有默认起始Marker元素 "START"
            if index==0: continue
            lat.append(value.elem ,self.pyline.reverse[index ] )
        parser._set_database(lat )
        # print("after add element:4 ")
        # self.cache_parser._set_database(lat )
        # print("after add element:5 ")
        if lat.stat.period==True:
            self._parse(parser,"CONSTRAINT,EXPR:=DIM(ABS(END[0].R11+END[0].R22),2)+DIM(ABS(END[0].R33+END[0].R44),2 );")
        self.length=lat.length+1
        lat.calculate()

        #register off-momentum RDTs to id_table, so that they're accessible for parser.
        dp0=lat.stat.dp0
        lat.register_RDTs(dp0,True)
        # lat.stat.dp0=-lat.stat.off_rdts_observer
        lat.register_RDTs(dp0-lat.stat.off_rdts_observer,False)
        # lat.stat.dp0=lat.stat.off_rdts_observer
        lat.register_RDTs(dp0+lat.stat.off_rdts_observer,False)
        # lat.stat.dp0=dp0
        # print("after add element:6 ")
        self.stat.stat=lat.stat
        lat.recover_twiss()
        return parser

    def set_parallel(self, int nkernel=-1):
        """
        set_parallel(self, Py_ssize_t nkernel=None)
            set threads of parallel
            params: 
                nkernel: (default: -1, max_threads of machine )
        """
        cdef Py_ssize_t i,num_thread = openmp.omp_get_max_threads()
        
        cdef string name=<string>(self.name.encode("utf8") )
        cdef CppBeamLine *lat=NULL 
        cdef Parser parser

        if nkernel == -1:
            self.nkernel = num_thread
        elif nkernel>num_thread:
            self.nkernel = num_thread
        elif nkernel ==0 or nkernel <0 :
            raise ValueError("nkernel should not be 0 or minus number(other than -1)!")
        else:
            self.nkernel = nkernel
        cdef Py_ssize_t __SIZEOF_DOUBLE__=8
        if len(self.parser_pool)<self.nkernel:
            for i in range(len(self.parser_pool), self.nkernel):
                lat = new CppBeamLine(name, self.stat.stat.particle, self.stat.stat.energy,self.length-1,self.stat.stat)
                memcpy( &lat.line[0].local, &self.lat.line[0].local, LOC_NUM*__SIZEOF_DOUBLE__ )
                memcpy(&(lat.line[0].tws[0]), &(self.lat.line[0].tws[0]), TWS_NUM*__SIZEOF_DOUBLE__)
                self.lat_pool.emplace_back( lat )
                parser = self.start_cppbeamline(lat)
                self.parser_pool.append( parser )

        for i in range(self.nkernel):
            self.update_parser(self.parser_pool[i], self.parser)

    def update_parser(self,Parser parser1, Parser parser0):
        """
        update_parser(self,Parser parser1, Parser parser0)
        update the parser1 from parser0, usually called internally.
        """
        len0 = len(parser0.tokens)
        len1 = len(parser1.tokens)
        if len1>len0:
            raise ValueError("first parser should bee updated by second parser!")
        for i in range(len1,len0):
            parser1.lexer.tokens = parser0.tokens[i]
            parser1.lexer.count = 0
            parser1.lexer.token_num = len(parser1.lexer.tokens)
            parser1.tokens.append(parser0.tokens[i] )
            parser1.parse()


    def set_worker(self,Py_ssize_t index=0):
        """
        set_worker(self,Py_ssize_t index=0)
        change the default worker, index <= nkernel ,then the old default worker will stored to the index worker
        """
        cdef:
            CppBeamLine* lat0 = self.lat 
            Parser parser0= self.parser
        if index > len(self.parser_pool) or index<0:
            raise ValueError(f"index should less than size of thread pool {len(self.parser_pool)} and no less than 0, while {index} is given!")

        self.lat = self.lat_pool.at(index)
        self.parser=self.parser_pool[index]
        self.update_parser(self.parser, parser0)
        self.lat_pool[index] = lat0
        self.parser_pool[index] = parser0
        
        if self.initial_worker == -1:
            # -1 mean current worker is the initial worker, so the 
            self.initial_worker = index
            self.current_worker = index
        elif self.initial_worker == index:
            # index mean intial worker is stored in pool, so the intial worker is set as current worker again,current worker is -1
            self.initial_worker = -1
            self.current_worker = -1
        else:
            # intial worker is stored in pool, so the intial worker is still stored in pool, current worker changed to index

            self.current_worker = index

    
    
    def eval(self, str expr):
        pass
    


    cdef object _parse(self,Parser parser,str code):
        parser.lexer.tokenize(code )
        parser.tokens.append(parser.lexer.tokens )
        # self.parser.codes.append(self.parser.lexer.codes )
        return parser.parse()

    def parse(self,str code):
        """
        VAR,NAME=$TTQ[1-8]$[0].k1, LOWER=-3, UPPER=3,STEP=1e-6;
        VAR,TTQ8[0].l:=2/TTQ7[0].k2;
        CONSTRAINT,EXPR:= DIM(ABS(END[0].betax-4.7428944),1e-2)/4+DIM(ABS(END[0].betay-3.0217133),1e-2)/4;
        OPTIMIZE,EXPR:=2e-0*(DIM(dnux_dJx,1e4 )+DIM(-5000,dnux_dJx ) +DIM(dnux_dJy,1e4 )+DIM(-5000,dnux_dJy )   );
        CHROM,AIM_DQX=1.0,KNOB=CCXS1;
        CHROM,AIM_DQY=1.65,KNOB=CCYS1;
        # above are the rules of set optimization parameters
        h11001;h00111;dnux_dJx;dnuy_dJy;START[0].betax #will return a list of 5 values before
        """
        # print("BeamLine.parse")
        # self.parser.lexer.tokenize(code )
        # self.parser.tokens.append(self.parser.lexer.tokens )
        # self.parser.codes.append(self.parser.lexer.codes )
        return self._parse(self.parser, code)
        #pass
    
    
    def _save(self):
        """
        internal funcction
        """"
        cdef:
            CppBeamLine* lat=NULL
            bytes name
        lat=self.lat
        for index in lat.elems:
            name=index.first
            if name.decode("utf8") in ("START","END") :continue
            lat.save((<Element?>self.pyline.elems[name.decode("utf8")]).elem )
    
    def save(self):
        """
        save the value of internal to the python objects
        """
        self._save()


    


    cdef void _update_variables(self,CppBeamLine* lat,  double* variables)nogil:
        cdef: 
            int i
            string name
            Variables* pvar=&lat.vars
        i=0
        # print("_update_variables")
        for name in pvar.ordered_independ_var_names:
            (<Var*>pvar.independent_vars[name]).left.data[0]=variables[i]
            #print((<Var*>pvar.independent_vars[name]).left.data[0])
            i+=1
        for name in pvar.ordered_depend_var_names:
            pvar.dependent_vars[name].calc()
        pvar=NULL
            #pass
        # for i in range(self.nseq):
        #     self.elems[i].update(&self.kwd_properties[i][0])
    
    
    cdef double _update_constraints(self,CppBeamLine* lat,  double* CV)nogil:
        cdef:
            size_t i
            string name
            double value=0,summary=0,cell=0
            Constraints* pconstraint=&lat.constraints
            double *tmp_ptr = NULL 
        # print("_update_constraints")
        computeRDTs = lat.stat.computedrivingterms
        lat.stat.computedrivingterms=False
        lat.calculate()
        lat.stat.computedrivingterms=computeRDTs
        if lat.stat.period:
            cell=pconstraint.values[0].calc()
        if cell>0:
            fill(CV,CV+pconstraint.num_constraint,1e10 )
            CV[0]=cell
            summary=cell*cell+1e20*(pconstraint.num_constraint-1)
        else:
            if lat.chrom_corrector.iscorr1 or lat.chrom_corrector.iscorr2:
                #computeRDTs = lat.stat.computedrivingterms
                #lat.stat.computedrivingterms=False
                #self.lat.calculate()
                lat.correctChrom(lat.chrom_corrector.aim_dQx, lat.chrom_corrector.aim_dQy )
                # calling linearoptics, integrate radiation secondly
            if lat.constraints.time_consuming_terms[DRIVING_TERMS ] and computeRDTs:
                lat.computeRDTs()
                if lat.stat.off_momentum_rdts:
                    lat.compute_off_momentum_RDTs()
                    lat.recover_twiss()
            if lat.constraints.time_consuming_terms[CHROMATIC_TERMS ] and (lat.stat.second_order_chrom or lat.stat.third_order_chrom):
                lat.computeSecondOrderChromaticities(lat.stat.dp)
            if lat.constraints.time_consuming_terms[MONITOR_OFF_MOMENTUM_TERMS ] :
                lat.compute_large_off_momentum_tunes(lat.stat.monitor_dp)
            if lat.constraints.time_consuming_terms[OFF_MOMENTUM_SUM_TERMS ] :
                lat.compute_off_momentum_sum_square(lat.stat.monitor_dp)
            if lat.constraints.time_consuming_terms[DA_TRACKING_TERMS ] :
                lat.get_DA_area( tmp_ptr )
            #更新ID表和约束
            lat.TwissPropagate()
            for name in lat.id_table.id_table:
                lat.id_table.id_dict[name].calc()
            for i in range(pconstraint.num_constraint):
                value=pconstraint.values[i].calc()
                CV[i]=value
                # print("CV: ",value)
                summary+=value*value
        pconstraint=NULL
        return summary



    cdef double _update_optima(self,CppBeamLine* lat,  double* optima, bint is_feasible)nogil:
        cdef:
            size_t i
            double value=0,summary=0
            Optima* poptima=&lat.optima
            double *tmp_ptr = NULL 
        # print("_update_optima")

        computeRDTs = lat.stat.computedrivingterms
        cdef bint lazy_compute=lat.stat.lazy_compute
        if is_feasible or not lazy_compute:
            if poptima.time_consuming_terms[DRIVING_TERMS ] and not lat.constraints.time_consuming_terms[DRIVING_TERMS ] and computeRDTs:
                lat.computeRDTs()
                if lat.stat.off_momentum_rdts:

                    lat.compute_off_momentum_RDTs()
                    lat.recover_twiss()
            if poptima.time_consuming_terms[CHROMATIC_TERMS ] and not lat.constraints.time_consuming_terms[CHROMATIC_TERMS ] and (lat.stat.second_order_chrom or lat.stat.third_order_chrom):
                lat.computeSecondOrderChromaticities(lat.stat.dp)
            if poptima.time_consuming_terms[MONITOR_OFF_MOMENTUM_TERMS ] and not lat.constraints.time_consuming_terms[MONITOR_OFF_MOMENTUM_TERMS ] :
                lat.compute_large_off_momentum_tunes(lat.stat.monitor_dp)
            if lat.stat.larger_monitor_dp >1.0:
                lat.compute_large_off_momentum_tunes(lat.stat.larger_monitor_dp*lat.stat.monitor_dp)
            if poptima.time_consuming_terms[DA_TRACKING_TERMS ] and not lat.constraints.time_consuming_terms[DA_TRACKING_TERMS ] :
                lat.compute_large_off_momentum_tunes(lat.stat.monitor_dp)
            if poptima.time_consuming_terms[OFF_MOMENTUM_SUM_TERMS ] and not lat.constraints.time_consuming_terms[OFF_MOMENTUM_SUM_TERMS ] :
                lat.compute_off_momentum_sum_square(lat.stat.monitor_dp)
            if poptima.time_consuming_terms[DA_TRACKING_TERMS ] and not lat.constraints.time_consuming_terms[DA_TRACKING_TERMS ] :
                lat.get_DA_area( tmp_ptr )


        for i in range(poptima.num_optima):
            value = poptima.minormax[i]*poptima.values[i].calc()
            optima[i]=value
            summary+=value*value
        poptima=NULL
        return summary


    cdef void _evolution(self,CppBeamLine* lat, double[:] variables, double[:] objectives , double[:] CV  )nogil:
        cdef:
            # CppBeamLine* lat=NULL
            CppStatus stat0
            bint computeRDTs
            double cv_summary
        self._update_variables(lat, &variables[0])
        # print("_evolution: ",np.asarray(CV.base) )
        cv_summary=self._update_constraints(lat, &CV[0] )
        if 0.0 == cv_summary:
            self._update_optima(lat, &objectives[0], True)
            # print(np.asarray(objectives[:]))
        else:
            self._update_optima(lat, &objectives[0], False)


    def evolution(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV ):
        """
            evolution(self,double[:,:] variables, double[:,:] objectives, double[:,:] CV )
                variables: MxN array corresponing to the VAR of parser 
                objectives: MxL array corresponing to the CONSTRAINT of parser 
                CV: MxP array corresponing to the OPTIMIZE of parser 
            returns:
                CV, objectives
        """
        cdef:
            int i,j,i_cv_best,i_objv_best, num_variables= variables.shape[1], num_optima = objectives.shape[1], num_constraints=CV.shape[1], num_pop = variables.shape[0]
            double cv_tmp,value,objv_value
            CppBeamLine* lat=NULL
        lat=self.lat

        if num_variables != lat.vars.num_independent_vars:
            raise ValueError(f"Input numpy array:{num_variables} doesn't match the variables:{lat.vars.num_independent_vars}!")
        if num_optima != lat.optima.num_optima:
            raise ValueError(f"Input numpy array:{num_optima} doesn't match the optima:{lat.optima.num_optima}!")
        if num_constraints != lat.constraints.num_constraint:
            raise ValueError(f"Input numpy array: {num_constraints} doesn't match the constraints:{lat.constraints.num_constraint}!")
        cdef Py_ssize_t job_num_per_woker,  num_full_job, nkernel=self.nkernel

        # if num_pop< nkernel, don't use all kernels
        if num_pop < nkernel: nkernel = num_pop
        if self.nkernel>1 and num_pop > 1:
            for i in prange(num_pop,nogil=True,num_threads=nkernel):
                self._evolution(self.lat_pool[openmp.omp_get_thread_num() ], variables[i,:], objectives[i,:], CV[i,:]   )
        else:
            for i in range(num_pop):
                self._evolution(lat, variables[i,:], objectives[i,:], CV[i,:]   )
        # cdef list a=_job_start_index, b=_job_end_index
        # print(a,"\n",b)
        # print(np.asarray(CV) )
    

    def __call__(self,double[:] X):
        """
        X: 1-d array of N elements corresponing to the VAR of parser 
        return:
            CV, objectives
            CV: 1-d array of P elements corresponing to the OPTIMIZE of parser 
            objectives: 1-d array L elements corresponing to the CONSTRAINT of parser 
        """
        cdef:
            int i,i_cv_best,i_objv_best, num_variables= X.shape[1]
            double cv_tmp,value,objv_value
            CppBeamLine* lat=NULL
        lat=self.lat
        
        if num_variables != lat.vars.num_independent_vars:
            raise ValueError(f"Input numpy array:{num_variables} doesn't match the variables:{lat.vars.num_independent_vars}!") 

        cdef double[:] CV=np.zeros(lat.constraints.num_constraint,dtype=float)
        
        cdef double[:] Obj=np.zeros(lat.optima.num_optima,dtype=float)
     
        self._evolution(self.lat, X, Obj, CV   )
        
        return (np.asarray(CV),np.asarray(Obj) )


    def calc(self):
        cdef:
            double value=0,summary=0,cell=0
            CppBeamLine*  lat=self.lat
            Constraints* pconstraint=&self.lat.constraints
            Optima* poptima=&self.lat.optima
            string name
        lat.calculate()
        if lat.stat.second_order_chrom or lat.stat.third_order_chrom:
            lat.computeSecondOrderChromaticities(lat.stat.dp)
        for name in lat.id_table.id_table:
            lat.id_table.id_dict[name].calc()
        for i in range(pconstraint.num_constraint):
            pconstraint.values[i].calc()
        for i in range(poptima.num_optima):
            value = poptima.minormax[i]*poptima.values[i].calc()
        lat.TwissPropagate()
    
    def findclosedorbit(self, double dp):
        """
        dp: off-momentum Delta_p/p_0
        find the closed orbit with 4-D track
        """
        cdef double[::1] r=np.array([0,0,0,0,0,dp],dtype="float")
        self.lat.findClosedOrbit( &r[0] )
    
    def highorderchromaticity(self,double dp0=0.0001 ):
        self.lat.computeSecondOrderChromaticities(dp0)

    def compute_large_off_momentum_tunes(self):
        self.lat.compute_large_off_momentum_tunes(self.lat.stat.monitor_dp)

    def compute_off_momentum_RDTs(self):
        self.lat.compute_off_momentum_RDTs()

    def correctchrom(self, dQx=None, dQy=None):
        """
        correctchrom(self, dQx=None, dQy=None)
        correct 1st order chromaticity to (dQx,dQy) use the knob set by CHROM in parser
        if dQx, dQy is None, the chromaticity is corrected to the value set by AIM_DQX, AIM_DQY
        """
        cdef:
            CppBeamLine* lat=self.lat
        if lat.chrom_corrector.iscorr1 or lat.chrom_corrector.iscorr2:
            #computeRDTs = lat.stat.computedrivingterms
            #lat.stat.computedrivingterms=False
            #self.lat.calculate()
            dqx = lat.chrom_corrector.aim_dQx if dQx is None else dQx
            dqy = lat.chrom_corrector.aim_dQy if dQy is None else dQy
            lat.correctChrom( dqx, dqy )

    
    def compute_off_momentum_twiss(self, list dp_range, double dp_step=1e-4, bint local_twiss=True):
        """
        compute_off_momentum_twiss(self, list dp_range, double dp_step=1e-4, bint local_twiss=True)
            dp_range: (min dp, max dp)
            dp_step : the step of dp to compute the off-momentum twiss functions at the end of beamline
            local_twiss : whether to calculate the  twiss functions of local elements
        returns:
            dps, twiss_mat
            dps: 1-d array of N dp elements to calculate off-momentum twisses
            twiss_mat: NxNUM_TWS array of the twisses function
        """
            
        cdef:
            CppBeamLine* lat=self.lat
            size_t i
        range_size = len(dp_range)
        def range_error():
            raise ValueError(f"dp_range should be list of [dpmin/p,dpmax/p], but {dp_range} is given!")

        if range_size != 2: range_error()
        if isinstance(dp_range[0],float) and isinstance(dp_range[1],float) and dp_range[0] < dp_range[1]:
            pass
        else:
            range_error()
        if dp_step == 0: raise ValueError(f"dp_step should be non-zero float, but {dp_step} is given!")
        dpmin,dpmax = dp_range
        total_point = math.ceil( (dpmax-dpmin)/dp_step ) + 1
        dp_lst=np.linspace(dpmin,dpmax, total_point, endpoint=True)
        cdef double[:,:] twiss_mat = np.zeros((total_point, TWS_NUM ), dtype=float)
        for i, dp in enumerate(dp_lst):
            lat.compute_off_momentum_twiss(&twiss_mat[i,0], dp, local_twiss)
        return dp_lst, np.asarray(twiss_mat)
    

    def track(self, double[:,::1] beam0, int start_pos=0, int end_pos=-1, int nturn0=1,int nturn1=1):
        """
        track(self, double[:,::1] beam0, int start_pos=0, int end_pos=-1, int nturn0=1,int nturn1=1)
            track(self, double[:,:] beam, start_pos=0, end_pos=-1, nturn0=1,nturn1=1)
            beam: Nx7([X ,PX ,Y ,PY ,Z ,DP ,LOSS) array, loss is zero if particle is not lost
            start_pos: int,the start element to track beams
            end_pos: int, the end element to end the tracking
            nturn0: int, the nturn0-th turn to start tracking
            nturn1: int, the nturn1-th turn to end the tracking
        return 
            N x [X ,PX ,Y ,PY ,Z ,DP ,LOSS ,   LOSSTURN    ,LOSSPOS   ]
        """
        cdef:
            int nbegin, nend, nturn_begin, nturn_end, num_col=beam0.shape[1] , num_particles=beam0.shape[0]
            int length = self.lat.length
            int max_turns= 1000000
            int i
        
        def check_position(position0, arg_name ):
            if position0<0 and position0+length+1 >= 0:
                return position0+length+1
            elif position0>=0 and position0 < length+1:
                return position0
            else:
                raise ValueError(f"{arg_name} should be int between [{-length-1},{length}], but {position0} is given!  ")
        
        nbegin = check_position(start_pos, "start_pos")
        nend = check_position(end_pos, "end_pos")

        if nturn0>nturn1: raise ValueError(f"nturn0 should be no more than nturn1!")

        if nturn0 >0 and nturn0 <max_turns:
            nturn_begin=nturn0
        else:
            raise ValueError(f"nturn0 should be int between [1,{max_turns}], but {nturn0} is given!  ")
        if nturn1 >0 and nturn1 <max_turns:
            nturn_end = nturn1
        else:
            raise ValueError(f"nturn1 should be int between [1,{max_turns}], but {nturn1} is given!  ")
        
        # 构造返回beam
        if num_col == TRK_NUM:
            beam1 = np.array(beam0) 
        elif num_col == 7:
            beam1 = np.zeros((num_particles, TRK_NUM),dtype=float)
            for i in range(num_particles):
                beam1[i,:7]=beam0[i,:7]
        else:
            raise ValueError(f"Shape of beam0 should be (N,7) or (N,{TRK_NUM})!")
        
        cdef double[:,::1] res = beam1
        cdef double *pbeam1 = &res[0,0]

        for i in prange(num_particles,nogil=True):
            self.lat.track(pbeam1+i*TRK_NUM, nbegin, nend,nturn_begin, nturn_end )
        
        return beam1

            

    def display(self,str token,bint detail=False):
        cdef:
            CppBeamLine* lat=NULL
            bytes name
            size_t index
        lat=self.lat
        if token=="KWD":
            self.lat.display(KWD,detail)
        elif token=="TWS":
            self.lat.display(TWS,detail)
        elif token=="DVTs":
            self.lat.display(DVTs,detail)
        elif token=="ID":
            num_id=lat.id_table.num_id
            print(f"{'ID':30}{'total:':10}{num_id:>20}")
            print(60*"=")
            index=0
            for name in lat.id_table.id_table:
                value=lat.id_table.id_dict[name].value
                print(f"{index:<4}:{name.decode('utf8'):25}{value:>30.e6}")
            print(60*"=")
        elif token=="VAR":
            num_indep_vars=lat.vars.num_independent_vars
            num_dep_vars=lat.vars.num_dependent_vars
            print(f"{'VAR':30}{'total:':10}{num_indep_vars+num_dep_vars:<5}{'CoVar:':10}{num_dep_vars:>5}")
            print(60*"=")
            for index in range(num_indep_vars):
                name=lat.vars.ordered_independ_var_names[index]
                lb=(<Var*>lat.vars.independent_vars[name]).lb
                ub=(<Var*>lat.vars.independent_vars[name]).ub
                print(f"{index:<4}:{name.decode('utf8'):25}{lb:15.6}{ub:15.6}")
            for index in range(num_dep_vars):
                name=lat.vars.ordered_depend_var_names[index]
                print(f"{index+num_indep_vars:<4}:{name.decode('utf8'):36}{'CoVar':20}")
            print(60*"=")
        elif token=="CONSTRAINT":
            num_constraint =lat.constraints.num_constraint
            print(f"{'CONSTRAINT':30}{'total:':10}{num_constraint:>20}")
            print(60*"=")
            for index in range(num_constraint):
                name=lat.constraints.ordered_constraint_names[index]
                print(f"{index:<4}:{name.decode('utf8'):56}")
            print(60*"=")
        elif token=="OPTIMIZE":
            num_optima =lat.optima.num_optima
            print(f"{'OPTIMIZE':30}{'total:':10}{num_optima:>20}")
            print(60*"=")
            for index in range(num_optima):
                name=lat.optima.ordered_optima_names[index]
                min_or_max="MINIMIZE" if lat.optima.minormax[index] >0 else "MAXIMIZE"
                print(f"{index:<4}:{name.decode('utf8'):36}{min_or_max:>20}")
            print(60*"=")


    def export(self,str filename, str filetype="atpy"):
        """
        export(self,str filename, str filetype="atpy")
            export the lattice to the form of other simulation programs ( some detail of the transformation may be not correct)
            filename:  names of file to store the output lattice
            filetype: str, one of the ["atpy","elegant", "opa","sad","madx"]
        """
        cdef:
            CppBeamLine* lat=NULL
        lat=self.lat
        if filetype =="atpy":
            with open(filename,"w+") as fn:
                fn.write('lat=r"""\nfrom atpy import*\n')
            self.pyline.export(filename,filetype )
            with open(filename,"a+") as fn:
                fn.write('"""\n')
                fn.write('from atpy import*\n')
                fn.write('translate(lat)\n')
                fn.write('from atpy.tools.translate import*\n')
            self.stat.export(filename, filetype)
            param_name=["betax", "betay", "alphax", "alphay", "etax", "etapx" ]
            tws0_str=", ".join( [ f"'{name}':{lat.line[0].tws(-1,TWS_INDEX[name] ):13.6e}"  for name in  param_name if abs(lat.line[0].tws(-1,TWS_INDEX[name] ) )>1e-10 ] )
            with open(filename,"a+") as fn:
                fn.write(f"\ntws0={{{tws0_str}}}\n")
                fn.write(self.str(filetype) )
        elif filetype =="opa":
            param_name=["betax", "betay", "alphax", "alphay", "etax", "etapx" ]
            tws0_str="; ".join( [ f"{name}={lat.line[0].tws(-1,TWS_INDEX[name] ):13.6e}"  for name in  param_name if abs(lat.line[0].tws(-1,TWS_INDEX[name] ) )>1e-10 ] )
            with open(filename,"w+") as fn:
                fn.write(f"energy={self.stat.stat.energy*1e-9};\r\n")
                fn.write(f"{tws0_str};\r\n")
            self.pyline.export(filename,filetype )
        elif filetype =="sad":
            self.pyline.export(filename,filetype )
            with open(filename,"a+") as fn:
                fn.write(f"CHARGE=1;\n")
                fn.write(f"MASS={lat.globals[MASS0] };\n")
                fn.write(f"MOMENTUM={self.stat.stat.energy};\n")
                fn.write(f"FFS;\nUSE {self.pyline.name};\n{'CELL' if self.stat.stat.period else 'INS'};\nCALC;\n")
        elif filetype =="madx":
            self.pyline.export(filename,filetype )
            pass
        elif filetype =="elegant":
            self.pyline.export(filename,filetype )
            
    

    def str(self,str filetype="atpy"):
        if filetype =="atpy":
            return f'{self.name}={self.kind}("{self.name}",stat,{self.pyline.name},**tws0)\n'





    cdef object _getitem(self,name, int position=0):
        cdef:
            int kind
            CppBeamLine* lat=NULL
            double value
        lat=self.lat
        if name in KWD_INDEX.keys():
            kind = KWD_INDEX[name]
            if cpp_count(lat.line[position].elem.keywords.begin(),lat.line[position].elem.keywords.end(), kind ):
                return lat.line[position].elem.values[kind]
            else:
                return None
        elif name in TWS_INDEX.keys():
            kind =TWS_INDEX[name]
            return lat.line[position ].tws(-1,kind)
        elif name in LOC_INDEX.keys():
            kind = LOC_INDEX[name]
            return lat.line[position ].local[kind]
        elif name in GLB_INDEX.keys():
            kind = GLB_INDEX[name]
            return lat.globals[kind]
    

    cdef void _setitem(self,name, int position, double value):
        cdef:
            int kind
            CppBeamLine* lat=NULL
        lat=self.lat
        if name in KWD_INDEX.keys():
            kind = KWD_INDEX[name]
            if cpp_count(lat.line[position].elem.keywords.begin(),lat.line[position].elem.keywords.end(), kind ):
                lat.line[position].elem.values[kind]=value
            else:
                pass
        elif name in TWS_INDEX.keys():
            kind =TWS_INDEX[name]
            (&lat.line[position ].tws[-1])[kind]=value
        elif name in LOC_INDEX.keys():
            kind = LOC_INDEX[name]
            lat.line[position ].local[kind]=value
        elif name in GLB_INDEX.keys():
            kind = GLB_INDEX[name]
            lat.globals[kind]=value

    def optics(self, start=0, stop=None,step=0.01,key=None,*,bint at_s=False):
        """
        optics(self, start=0, stop=None,step=0.01,key=None,*,bint at_s=False)
        return the twiss function at given steps, usually for ploting
        """
        cdef:
            double tmp_tws[TWS_NUM]
            CppBeamLine* lat=NULL 
            size_t nrow=0, ncol=0
            double ubound=0
            size_t _start_index=0,_stop_index=0
        lat=self.lat
        fill( &(tmp_tws[0]), &tmp_tws[0]+TWS_NUM,0);
        # 根据位置范围or元件范围确定上下边界
        if at_s:
            ubound=lat.line[lat.length ].local[S]
            lbound=0
        else:
            ubound=lat.length
            lbound=-ubound-1
        # 判断输入参数是否超出边界
        if start<lbound or start> ubound :
            raise IndexError("Start point {start} is out of range!")
        if stop:
            if stop<lbound or stop>ubound :raise IndexError("Stop point {stop} is out of range!")
        else:
             _stop=lat.line[lat.length].local[S]   #stop没有输入参数，则无论时at_s 或at_elem，_stop都应是最末端
        if abs(step)<1e-4 :  #step不用太小，
            raise IndexError("Step point {step} is smaller than 1e-4!")
        # 根据at_s or at_elem 确定起始点,根据正负号
        if at_s:
            _start=start
            if stop:_stop=stop
        else:#at_elem
            if not isinstance(start,int):
                ValueError("Start point should be int type when at_s is False, while {start} is input!")
            if stop:
                if isinstance(stop,int):
                    ValueError("Stop point should be int type, while {stop} is input!")
                if stop<0:
                    _stop=lat.line[ubound+stop+1].local[S]
                else:
                    _stop=lat.line[stop].local[S]
            #else stop==None时 _stop 已经在前面赋值
            if start<0:
                _start=lat.line[ubound+start+1].local[S]
            else:
                _start=lat.line[start].local[S]
        _start_index=lat.get_position_at_s(_start)
        _stop_index=lat.get_position_at_s(_stop )
        _step = step
        # 
        if (_start>_stop and _step>0) or (_start<_stop and _step<0 ):
            raise ValueError(f"Start point and stop point should coincident with step, while is start={start}, stop={stop}, step={step}!")
        if _start<_stop:
            lb=_start
            coordinates=[]
            for i in range(_start_index,_stop_index+1 ):
                if i<_stop_index:
                    ub=lat.line[i].local[S]
                else:
                    ub=_stop
                npoints=round( abs((ub-lb)/_step) )+1
                tmp_positions = np.linspace(lb,ub, npoints ).tolist()
                _tmp_start_point = 1 if npoints>1 else 0
                coordinates = coordinates+tmp_positions[_tmp_start_point:]
                lb=ub 
        else: #逆序数据 _start>_stop
            ub=_start
            coordinates=[]
            for i in range(_start_index,_stop_index-1,-1 ):
                if i>_stop_index:
                    lb=lat.line[i-1].local[S]
                else:
                    lb=_stop
                npoints=round( abs((ub-lb)/_step) )+1
                tmp_positions = np.linspace(ub,lb, npoints ).tolist()
                _tmp_start_point = 1 if npoints>1 else 0
                coordinates = coordinates+tmp_positions[_tmp_start_point:]
                ub=lb 
        nrow=len( coordinates )
        ## key关键字参数解析
        arg2=key 
        if arg2 is None:
            param = list(INDEX_TWS.keys() )
        elif isinstance(arg2,str):
            if arg2 == "ALL":
                param = list(INDEX_TWS.keys() )
            elif arg2 in TWS_INDEX.keys():
                param = [ TWS_INDEX[ arg2 ] ]
            else:
                raise TypeError(f"Unknown parameter {arg2} for the 2-nd argument!")
        elif isinstance(arg2,list):
            if arg2==[]: raise ValueError("Empty list is not expected for the 2-nd argument!")
            tws_names=list(TWS_INDEX.keys() )
            unique_names=set(arg2 )
            if len(unique_names)<len(arg2 ): raise ValueError(f"Repeated twiss word in the second argument!")
            if all( [True if name in tws_names else False for name in arg2 ] ):
                param = [TWS_INDEX[name ] for name in arg2 ]
            else:
                raise ValueError(f"Unexpected parameter word in 2-nd argument!")
        else:
            raise TypeError(f"Unknown {arg2} for the key argument!")
        ncol=len(param)
        
        # print("positions: ",coordinates )
        # print("param: ", param )
        ##set return values
        if nrow >1 and ncol>1:
            res = np.zeros((nrow,ncol),dtype="float")
            for i,pos in enumerate(coordinates):
                lat.get_optics_at_s(pos,tmp_tws )
                for j,para in enumerate(param):
                    res[i,j]=tmp_tws[para]
        elif nrow==1 and ncol==1:
            lat.get_optics_at_s(coordinates[0],tmp_tws )
            return tmp_tws[param[0] ]
        elif nrow==1:
            res = np.zeros(ncol)
            lat.get_optics_at_s(coordinates[0],tmp_tws )
            for j,para in enumerate(param ):
                res[j]=tmp_tws[para ]
        else:
            res=np.zeros(nrow)
            for i,pos in enumerate(coordinates):
                lat.get_optics_at_s(pos,tmp_tws )
                res[i]=tmp_tws[param[0] ]
        return res

    def get_DA_area(self):
        """
        get_DA_area
        get DA bounds with the form numpy.array([x_list,y_list])
        """        
        cdef double[:,::1] area_datas=np.zeros((2,self.lat.stat.track_lines) )
        self.lat.get_DA_area( &area_datas[0,0] )
        return np.array(area_datas.base )

    def __getitem__(self,args):
        """
        
        """
        cdef:
            CppBeamLine* lat=self.lat
            Var* var =NULL
            bytes bytes_name
            string string_name

        if isinstance(args,str):
            if args=="VAR":
                #bounds=[values,lb,ub,step]
                bounds=[ [],[],[],[] ]
                for i in range(lat.vars.num_independent_vars):
                    var = <Var*>lat.vars.independent_vars[lat.vars.ordered_independ_var_names[i] ] 
                    bounds[0].append(var.left.data[0] )
                    bounds[1].append(var.lb)
                    bounds[2].append(var.ub)
                    bounds[3].append(var.step)
                return bounds
            elif args=="CONSTRAINT":
                values=[]
                for i in range(lat.constraints.num_constraint):
                    values.append(lat.constraints.values[i].value )
                return values
            elif args=="OPTIMIZE":
                values=[]
                for i in range(lat.optima.num_optima):
                    values.append(lat.optima.values[i].value )
                return values
            elif args=="ID":
                values={}
                for string_name in lat.id_table.id_table:
                    bytes_name = string_name
                    values[bytes_name.decode("utf8") ]=lat.id_table.id_dict[string_name].value
                return values
            elif args=="STAT":
                return lat.stat
            elif args =='KWD':
                return KWD_INDEX.copy()
            elif args =='TWS':
                return TWS_INDEX.copy()
            elif args =='LOC':
                return LOC_INDEX.copy()
            elif args =='GLB':
                return GLB_INDEX.copy()
            elif args in GLB_INDEX.keys():
                return self._getitem(args)
            elif args in self.pyline.elems.keys():
                return self.pyline.elems[args ]
            #elif args in :
            else: # 返回匹配args的元素，返回Element列表
                pattern=re.compile(rf"^{args}$")
                names= [name for name in self.pyline.elems.keys()  if pattern.fullmatch(name) ]
                if names==[]: raise ValueError(rf"Unexpected arg {args}!")
                names.sort(key=lambda name: self.elems_index[name][0] )
                return [self.pyline.elems[name] for name in names ]
        elif isinstance(args,int):
            if args< -self.length or args>self.length-1:
                raise IndexError(f"Index is out of range in the first argument, {args}!")
            positions=[args if args>=0 else self.length+args ]
            return self.pyline.expand[args ]

        elif isinstance(args,list):
            if args==[]: raise ValueError("Empty list is not expected for list type argument!")
            unique_names = set(args)
            if len(unique_names) <len(args) :
                raise ValueError(f"Repeated word in argument!")
            if all([isinstance(argi,int) for argi in unique_names]):
                return [self.pyline.expand[argi] for argi in args ]
            elif all([isinstance(argi,str) for argi in unique_names]):
                glb_names = GLB_INDEX.keys()
                # res=[]
                for name in args:
                    if name not in glb_names: raise ValueError(f"Invalid word {name} in arg list!")
                    # res.append(self._getitem(name))
                return [self._getitem(name) for name in args ]
            else:
                raise ValueError("Input argument should be list of int as index or list of str within global parameters!")
        #切片返回line中元素列表
        elif isinstance(args,slice):
            return self.pyline.expand[args ]
        #返回参数矩阵
        elif isinstance(args,tuple):
            len_args=len(args)
            if len_args==2:
                arg1,arg2=args
                ## 第一个参数解析
                if isinstance(arg1,list):
                    if arg1==[]: raise ValueError("Empty list is not expected for the 1-st argument!")
                    if all( [isinstance(value,int) for value in arg1 ] ):
                        positions=arg1
                        lb=min(positions)
                        ub=max(positions)
                        positions=[value if value>=0 else value+self.length for value in positions ]
                        if lb< -self.length or ub>self.length-1:
                            raise IndexError(f"Index is out of range in the first argument, lb={lb},ub={ub}!")
                    elif all( [isinstance(value,str) for value in arg1 ] ):
                        names = [self.elems_index[name] for name in arg1]
                        positions=sum(names,[])
                        positions.sort()
                    else:
                        raise ValueError(f"The values should either all be int type or all be str type, when list as the first argument!")
                    nrow=len(positions)
                elif isinstance(arg1,int):
                    if arg1< -self.length or arg1>self.length-1:
                        raise IndexError(f"Index is out of range in the first argument, {arg1}!")
                    positions=[arg1 if arg1>=0 else self.length+arg1 ]
                    nrow=1
                elif isinstance(arg1,str):
                    if arg1 in self.elems_index.keys():
                        positions=self.elems_index[arg1 ]
                    else:
                        pattern=re.compile(rf"^{arg1}$")
                        names=[ self.elems_index[name] for name in self.elems_index.keys() if pattern.fullmatch(name) ]
                        positions = sum(names,[]) #展开两层list
                        positions.sort()
                        if positions==[]: raise ValueError(f"Invalid name {arg1}, should be element name or regex represent element names!")
                    nrow=len(positions )
                elif isinstance(arg1,slice):
                    total_index=list(range(self.length ) )
                    positions=total_index[arg1]
                    nrow=len(positions )
                else:
                    raise TypeError("1-st arg should be int,str,slice or list!")
                ## 第二个参数解析
                if isinstance(arg2,str):
                    if arg2 in POSITION_DEPEND_NAMES:
                        names=[arg2 ]
                        ncol =1
                    else:
                        raise ValueError(f"Unknown argument {arg2} as 2-nd argument!")
                elif isinstance(arg2,list):
                    if arg2==[]: raise ValueError("Empty list is not expected for the 2-nd argument!")
                    unique_names = set(arg2)
                    if len(unique_names) <len(arg2) :
                        raise ValueError(f"Repeated word in 2-nd argument!")
                    for name in arg2:
                        if name not in POSITION_DEPEND_NAMES: raise ValueError(f"Invalid word {name} in 2-nd argument!")
                    names=arg2
                    ncol=len(names )
                else:
                    raise TypeError(f"Unexpected arg {arg2}, 2-nd should be str or list!")
                ## 返回参数
                if nrow >1 and ncol>1:
                    res = np.zeros((nrow,ncol),dtype="float")
                    for i,pos in enumerate(positions):
                        for j,name in enumerate(names):
                            res[i,j]=self._getitem(name,pos )
                elif nrow==1 and ncol==1:
                    return self._getitem(names[0],positions[0] )
                elif nrow==1:
                    res = np.zeros(ncol)
                    for j,name in enumerate(names ):
                        res[j]=self._getitem(name,positions[0] )
                else:
                    res=np.zeros(nrow)
                    for i,pos in enumerate(positions):
                        res[i]=self._getitem(names[0], pos )
                return positions,res
            else:
                raise ValueError(f"Expect 2 args!")
            #return args
        else:
            raise ValueError("Invalid args!")

    def __setitem__(self,args,values):
        cdef:
            string string_name
            CppBeamLine* lat=NULL
        lat=self.lat
        if isinstance(args,tuple): 
            len_args=len(args)
            if len_args!=2: raise ValueError(f"Expect 2 args when input is tuple, while {args} given!")
            arg1,arg2=args
            # if isinstance(arg1,list):
            #     raise ValueError(f"Not implement function!")
            #     for py_int in arg1:
            #         if not isinstance(py_int,int) or py_int< -lat.length-1 or py_int>lat.length :
            #             raise ValueError(f"Expect int type for element position in range [{-lat.length-1}, {lat.length}], {py_int} is given!")
                ##### 尚未完成，主要用于设置误差及校正
            
            if isinstance(arg1,list):
                if arg1==[]: raise ValueError("Empty list is not expected for the 1-st argument!")
                if all( [isinstance(value,int) for value in arg1 ] ):
                    positions=arg1
                    lb=min(positions)
                    ub=max(positions)
                    positions=[value if value>=0 else value+self.length for value in positions ]
                    if lb< -self.length or ub>self.length-1:
                        raise IndexError(f"Index is out of range in the first argument, lb={lb},ub={ub}!")
                elif all( [isinstance(value,str) for value in arg1 ] ):
                    names = [self.elems_index[name] for name in arg1]
                    positions=sum(names,[])
                    positions.sort()
                else:
                    raise ValueError(f"The values should either all be int type or all be str type, when list as the first argument!")
                nrow=len(positions)
            elif isinstance(arg1,int):
                if arg1< -self.length or arg1>self.length-1:
                    raise IndexError(f"Index is out of range in the first argument, {arg1}!")
                positions=[arg1 if arg1>=0 else self.length+arg1 ]
                nrow=1
            elif isinstance(arg1,str):
                if arg1 in self.elems_index.keys():
                    positions=self.elems_index[arg1 ]
                else:
                    pattern=re.compile(rf"^{arg1}$")
                    names=[ self.elems_index[name] for name in self.elems_index.keys() if pattern.fullmatch(name) ]
                    positions = sum(names,[]) #展开两层list
                    positions.sort()
                    if positions==[]: raise ValueError(f"Invalid name {arg1}, should be element name or regex represent element names!")
                nrow=len(positions )
            elif isinstance(arg1,slice):
                total_index=list(range(self.length ) )
                positions=total_index[arg1]
                nrow=len(positions )
            else:
                raise TypeError("1-st arg should be int,str,slice or list!")
            ## 第二个参数解析
            if isinstance(arg2,str):
                if arg2 in POSITION_DEPEND_NAMES:
                    names=[arg2 ]
                    ncol =1
                else:
                    raise ValueError(f"Unknown argument {arg2} as 2-nd argument!")
            elif isinstance(arg2,list):
                if arg1==[]: raise ValueError("Empty list is not expected for the 2-nd argument!")
                unique_names = set(arg2)
                if len(unique_names) <len(arg2) :
                    raise ValueError(f"Repeated word in 2-nd argument!")
                for name in arg2:
                    if name not in POSITION_DEPEND_NAMES: raise ValueError(f"Invalid word {name} in 2-nd argument!")
                names=arg2
                ncol=len(names )
            else:
                raise TypeError(f"Unexpected arg {arg2}, 2-nd should be str or list!")
            #解析values
            
            def sizeerror():
                raise ValueError("Size of right values does not cosist with left args!")
            if isinstance(values,float):
                res=np.full((nrow,ncol),values,dtype="float")
            elif isinstance(values,int):
                res=np.full((nrow,ncol),values,dtype="float")
            elif isinstance(values,np.ndarray):
                dims=values.shape
                if values.ndim ==2:
                    if values.shape[0] != nrow: sizeerror()
                elif values.ndim ==1:
                    if (nrow==1 and values.shape[0] != ncol) or (nrow!=1 and values.shape[0] !=nrow ): sizeerror()
                else:
                    sizeerror()
                if values.size !=nrow*ncol: sizeerror()
                res=values.reshape(nrow,ncol)
                # if dims[0] !=nrow or dims[1] !=ncol: raise ValueError("Size of right values does not cosist with left args!")
            elif isinstance(values,list):
                if len(values) !=nrow: sizeerror()
                res=np.array(values,dtype="float").reshape(nrow,ncol)
                if res.size !=nrow*ncol: sizeerror()
            else:
                sizeerror()
            #############################################################
            if nrow >1 and ncol>1:
                for i,pos in enumerate(positions):
                    for j,name in enumerate(names):
                        self._setitem(name,pos, res[i,j] )
            elif nrow==1 and ncol==1:
                self._setitem(names[0],positions[0], values )
            elif nrow==1:
                for j,name in enumerate(names ):
                    self._setitem(name,positions[0], res[0,j] )
            else:
                for i,pos in enumerate(positions):
                    self._setitem(names[0], pos, res[i,0] )

        elif isinstance(args,list):
            for pystr_name in args:
                if not isinstance(pystr_name,str) or pystr_name not in STAT_INDEPEND_FIELDS: 
                    raise ValueError(f"Expect status field names when input is list, while {pystr_name} is given!")
            self.stat[args]=values
            lat.stat=self.stat.stat
        elif isinstance(args,str):
            if args=="VAR":
                if not isinstance(values,list):
                    raise ValueError(f"Expect list, while {values} given!")
                if len(values) != lat.vars.num_independent_vars:
                    raise ValueError(f"Expect {lat.vars.num_independent_vars} {{lb,ub}} list, while {len(values)} given!")
                for index,value in enumerate(values):
                    string_name=lat.vars.ordered_independ_var_names[index]
                    (<Var*?>lat.vars.independent_vars[string_name]).lb=<double?>value[0]
                    (<Var*?>lat.vars.independent_vars[string_name]).ub=<double?>value[1]
            elif args in STAT_INDEPEND_FIELDS:
                self.stat[args]=values
                lat.stat=self.stat.stat
                if self.nkernel>1:
                    for i in range(self.nkernel):
                        self.lat_pool[i].stat = lat.stat
            else:
                raise ValueError(f"Unexpected arg {args}!")
            # return args
        else:
            raise ValueError("Invalid args!")


    def __dealloc__(self):
        cdef int i
        if self.lat is not NULL:
            #print("BeamLine.__dealloc__")
            # self.parser.line=NULL
            # self.cache_parser.line=NULL
            
            # del self.parser.line
            # del self.cache_parser.line
            # print(self.lat.line[1].elem.values[0] )
            # deallocate(self.lat )
            del self.lat
            self.lat=NULL
        if self.lat_pool.size()>0:
            for i in range(self.lat_pool.size() ):
                del self.lat_pool[i]
                self.lat_pool[i]=NULL
            #pass
        


            

    
