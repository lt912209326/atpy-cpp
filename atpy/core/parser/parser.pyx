import re
import warnings

value2enum={"+":ADD, "-": SUB, "*": MUL, "/":DIV, "**":POW, "%":MOD, "//":FLOOR,
            "ABS":ABS, "SQRT":SQRT, "DIM":DIM, "MAX":MAX, "MIN":MIN,  "MAXABS":MAXABS, "MINABS":MINABS, 
            "NUMBER":NUMBER, "TWS":TWS, "KWD":KWD, "LOC":LOC, "GLB":GLB,
            "DELAY":DELAY
             }

cdef:
    tuple chromatic_terms = ("dQx", "dQy", "d2Qx", "d2Qy", "d3Qx", "d3Qy", "detax" , "detapx", "dbetax", "dalphax","dalphay"
    "dbetay", "ddetax", "ddbetax",  "ddbetay","Wx","Wy")
    
    tuple driving_terms = ("dnux_dJx", "dnux_dJy", "dnuy_dJy", 
        "H11110", "H22000", "H00220","H20110", "H31000", "H11200", "H00310", "H20200", "H20020", "H40000", "H00400",
                            "dnux_dJx@0", "dnux_dJy@0", "dnuy_dJy@0",
        "H11110@0", "H22000@0", "H00220@0", "H20110@0", "H31000@0", "H11200@0", 
                            "H00310@0", "H20200@0", "H20020@0", "H40000@0", "H00400@0" )
    tuple monitor_off_momentum_terms = ("lower_Qx", "high_Qy","lower_Qx", "high_Qy")
    tuple off_momentum_sum_terms = ( "sum_sqr_Qx","sum_sqr_Qy" )
    tuple da_tracking_terms = ("DA" ,"DA_SIGMA")
    tuple ma_tracking_terms = ("MA" ,"MA_SIGMA")

cdef class Parser:  # 定义语法分析器的类
    def __cinit__(self, dict eids):
        self.lexer = Lexer()  # 接收词法分析器对象
        self.eids = eids
        self.elems_name = list(eids.keys() )
        self.isdatabase = False
        self.vars_name = []
        self.vars_elem_name = []
        self.id_table=[]
        self.tokens = []
        self.term_names=[]
        self.ordered_var_positions=[]

    
    def __dealloc__(self):
        #print("Parser.__dealloc__")
        self.line=NULL
    
    
    cdef void _set_database(self,CppBeamLine* line)except*:
        self.line = line
        self.isdatabase = True
        

        
    cdef void eat(self,  str kind)except*:
        if self.current_token.kind == kind:
            self.term_names.append(self.current_token.value)
            self.current_token = self.lexer.get_next_token()
        else:
            raise ValueError(f'Parser.eat:Expect {kind }, but got {self.current_token.value} in line {self.current_token.line} :\n "{self.lexer.codes[self.current_token.line-1] }"\n at column {self.current_token.column}!')
        

    cdef int position(self)except*:
        cdef:
            Token token
            position=0,n_th=0
        token=self.current_token
        if token.kind != "ID" or token.value not in self.elems_name:
            raise SyntaxError(f"Unknown element name {token.value} at line {self.current_token.line} column {self.current_token.column}!")
        self.eat("ID")
        n_th=self.get_eid()
        if n_th>len(self.eids[token.value] )-1 or n_th<-len(self.eids[token.value]):
            raise IndexError(f"{n_th} out or range of element {token.value}!")
        position=self.eids[token.value][n_th]
        return position



    cdef tuple parameter(self)except*:
        cdef:
            Token token
            int kind,index
        token=self.current_token
        if self.current_token.kind=="KWD":
            kind=KWD
            index = KWD_INDEX[ token.value ]
            self.eat("KWD")
        elif self.current_token.kind=="TWS":
            kind=TWS
            self.eat("TWS")
            index = TWS_INDEX[ token.value ]
        elif self.current_token.kind=="LOC":
            self.eat("LOC")
            kind=LOC
            index = LOC_INDEX[ token.value ]
        elif self.current_token.kind=="GLB":
            self.eat("GLB")
            kind=GLB
            index = GLB_INDEX[ token.value ]
        else:
            raise SyntaxError(f"Unknown syntax {token.value} at line {self.current_token.line} column {self.current_token.column}!")
        return kind, index



    cdef AST* function(self, str func)except*:
        cdef:
            int start,end, kind, index, n_th,position
            list positions=[]
            AST* node
            Token   token
        node =NULL
        self.eat("(")

        if func in ("MIN", "MAX", "MINABS", "MAXABS"):
            for i in range(2):
                position=self.position()
                positions.append( position )
                self.eat(",")
            kind,index=self.parameter()
            if kind==GLB:
                raise SyntaxError(f"Global property is not allowed in {func} function!")
            node = new RangeFunction(value2enum[func], index,positions[0], positions[1],kind, self.line.line )
        else:
            node = self.expr()
            token = self.current_token
            if token.kind == "," and func =="DIM" :
                self.eat(",")
                node = new BiFunction(node, value2enum[func], self.expr() )
            elif token.kind == ")" and func in ("ABS", "SQRT"):
                node = new MonoFunction(value2enum[func], node)
            else:
                raise ValueError(f'Number args for {func} function  is not correct!')
        self.eat(")")
        return node

    
    
    cdef AST* property(self)except*:
        cdef: 
            AST*    node
            int     position=0, n_th=0 ,kind, index, parameter,data_kind
            bint    isglb=True
            Token   token
        if self.current_token.kind=="ID":
            position=self.position()
            self.eat(".")
            token=self.current_token
            kind,index=self.parameter()
            if kind == KWD:
                if cpp_count(self.line.line[position].elem.keywords.begin(), self.line.line[position].elem.keywords.end(),index)==0:
                    raise KeyError(f"{token.value} is invalid at line {token.line} column {token.column}!")
            if kind == GLB:
                raise SyntaxError(f"Global property  {token.value}  is not allowed with position at line {token.line} column {token.column}!")
        elif self.current_token.kind=="GLB":
            kind=GLB
            index=GLB_INDEX[self.current_token.value ]
            self.eat("GLB")
        else:
            raise SyntaxError(f"Unknown syntax {self.current_token.value}  at line {self.current_token.line} column {self.current_token.column}!")
        node =new Property(kind, position,index, self.line.line, self.line.globals)
        return node



    cdef AST* factor(self)except*:
        cdef: 
            AST* node
            Token token,current_token
            int pos_or_neg=1
            bytes name
            list tmp_idtable=[]
        current_token = self.current_token
        name = current_token.value.encode("utf8")
        cdef string name0 = current_token.value.encode("utf8")
        if current_token.kind == "NUMBER":
            self.eat("NUMBER")
            value = float(current_token.value)
            node = new Number(value)
        elif current_token.kind in ("-","+"):
            self.eat(current_token.kind)
            if current_token.kind =="-":
                pos_or_neg=-1
            if self.current_token.kind=="NUMBER":
                value = pos_or_neg*float(self.current_token.value)
                node = new Number(value)
                self.eat("NUMBER")
            elif self.current_token.kind=="(":
                if current_token.kind =="-":
                    node = new Number(0.0)
                    node = new Node( node, SUB, self.expr() )
                else:
                    node = self.expr()
            else:
                if current_token.kind =="-":
                    node = new Number(0.0)
                    node = new Node( node, SUB, self.factor() )
                else:
                    node = self.factor()

        # property and defined ID
        elif current_token.kind == "ID":
            if current_token.value in self.elems_name:
                node = self.property()
            elif current_token.value in self.id_table or name in self.id_table :
                self.eat("ID")
                # name=current_token.value.encode("utf8")
                node = new Refer(name, <Identity*>self.line.id_table.id_dict[ name ])
            elif cpp_algorithm.find(self.line.id_table.id_table.begin(),self.line.id_table.id_table.end(),name0) != self.line.id_table.id_table.end():
                self.eat("ID")
                name=current_token.value.encode("utf8")
                tmp_idtable = self.line.id_table.id_table
                self.id_table= list( set(self.id_table+tmp_idtable) ) 
                # .append( current_token.value )
                node = new Refer(name, <Identity*>self.line.id_table.id_dict[ name ])
            else:
                raise NameError(f"Undefined ID {current_token.value} at line {current_token.line} column {current_token.column}!")
        elif current_token.kind == "GLB":
            node = self.property()
        # FUNCTION 
        elif current_token.kind =="FUN":
            # (MIN, MAX, DIM, ABS, SQRT):
            self.eat(current_token.kind)
            node = self.function(current_token.value)
        elif current_token.kind in ("(",")"):
            self.eat("(")
            node = self.expr()  # # 创建运算符节点对象
            self.eat(")")
        else:
            raise SyntaxError(f"Unknown syntax {current_token.value} at {current_token.line} column {current_token.column}!")
        # ** // %
        while self.current_token.value in ("**", "//","%"):
            token=self.current_token
            self.eat(token.kind )
            if token.kind == "%":
                node = new Node(node,MOD,self.factor() )
            elif token.kind == "//":
                node = new Node(node,FLOOR,self.factor() )
            elif token.kind == "**":
                node = new Node(node,POW,self.factor() )
        return node  # 返回运算符节点对象
    
    
    cdef AST* term(self)except*:
        cdef: 
            AST* node
            Token token
        node = self.factor()  # 左侧节点对象
        while self.current_token.kind in ("*", "/"):
            token = self.current_token
            self.eat(token.kind)
            node = new Node(node, value2enum[token.kind], self.factor())  # 创建运算符节点对象
        return node  # 返回节点对象


    cdef AST* expr(self)except*:
        cdef: 
            AST* node
            Token token
        node = self.term()  # 左侧节点对象
        while self.current_token.kind in ("+", "-"):
            token = self.current_token
            self.eat(token.kind)
            node = new Node(node, value2enum[token.value], self.term())  #self.term 返回二元运算表达式的节点
        return node  # 返回运算符节点对象


    cdef tuple get_bounds(self)except*:
        cdef:
            AST* node
            double lb=-1e20,ub=1e20,step=1e-10
        if self.current_token.kind=="," and self.lexer.check_next_token()=="LOWER":
            self.eat(",")
            self.eat("LOWER")
            self.eat("ASSIGN")
            node = self.expr()
            lb=node.calc()
            del node
            node=NULL
        if self.current_token.kind=="," and self.lexer.check_next_token()=="UPPER":
            self.eat(",")
            self.eat("UPPER")
            self.eat("ASSIGN")
            node = self.expr()
            ub=node.calc()
            del node
            node=NULL
        if self.current_token.kind=="," and self.lexer.check_next_token()=="STEP":
            self.eat(",")
            self.eat("STEP")
            self.eat("ASSIGN")
            node = self.expr()
            step=node.calc()
            del node
            node=NULL
        return lb,ub,step

    cdef int get_eid(self)except*:
        cdef:
            int eid=0
            int neg_or_pos=1
        self.eat("[")
        if self.current_token.kind=="-":
            neg_or_pos=-1
            self.eat("-")
        if self.current_token.kind=="NUMBER" and "." not in self.current_token.value:
            eid=int(self.current_token.value)
            self.eat("NUMBER")
        else:
            raise SyntaxError(f"Unknown syntax {self.current_token.value}  at line {self.current_token.line} column {self.current_token.column}!")
        
        self.eat("]")
        return neg_or_pos*eid



    cdef void set_Identity(self)except*:
        cdef Token current_token=self.current_token
        cdef bytes name =self.current_token.value.encode("utf8")
        cdef str str_name = self.current_token.value
        if self.current_token.kind =="ID" and (self.lexer.check_next_token()=="DELAY" or self.lexer.check_next_token()=="ASSIGN") :
            if self.current_token.value not in self.id_table:
                self.id_table.append(self.current_token.value )
                self.line.id_table.id_table.push_back(name)
            else:
                del self.line.id_table.id_dict[name]
                warnings.warn(f"Identity {str_name} is redefined at line {current_token.line}, column {current_token.column}, which might cause error dependences!")
            self.eat("ID")
            if self.current_token.kind=="DELAY":
                self.eat("DELAY")
                self.line.id_table.id_dict[name] = new Identity(name,True, self.expr() )
            else:
                self.eat("ASSIGN")
                self.line.id_table.id_dict[name] = new Identity(name,False, self.expr() )
        else:
            raise SyntaxError(f"Unknown syntax {str_name} with {self.lexer.check_next_token()} at {current_token.line}, column {current_token.column}!")


    cdef void set_variable(self)except*:
        cdef:
            Token token, param_token
            AST* node=NULL
            double lb=-1e10,ub=1e10,step=1e-14
            str names, param_name
            string string_vars_name
            list ordered_var_positions=[]
            int eid=0,position=0,index=0,kind=0

        self.eat("VAR")
        self.eat(",")
        if self.current_token.kind == "NAME":
            self.eat("NAME")
            self.eat("ASSIGN" )

            token=self.current_token
            if token.kind=="REGEX":
                pattern=re.compile(rf"^{self.current_token.value[1:] }" )
                self.eat("REGEX")
                eid=self.get_eid()
                self.eat(".")
                param_name=self.current_token.value
                kind,index=self.parameter()
                lb,ub,step =self.get_bounds()
                for name in self.elems_name:
                    if pattern.fullmatch(name):
                        if eid>len(self.eids[name ] )-1 or eid<-len(self.eids[name]):
                            raise IndexError(f"{eid} out or range of element {name} at line {token.line} column {token.column}!")
                        position=self.eids[name][eid]
                        
                        if kind == KWD:
                            if cpp_count(self.line.line[position].elem.keywords.begin(), self.line.line[position].elem.keywords.end(),index)==0:
                                raise KeyError(f"{param_name} is invalid keyword for {name} at line {token.line} column {token.column}!")

                        self.ordered_var_positions.append(position )
                        var_name=f"{name}[{eid}].{param_name}"
                        if var_name in self.vars_name:
                            raise SyntaxError(f"Redefined Var {var_name} at line {token.line} column {token.column}!")
                        string_vars_name=var_name.encode("utf8")
                        self.vars_name.append(var_name) 
                        self.vars_elem_name.append(name) 
                        
                        node = new Property(kind, position,index, self.line.line, self.line.globals)
                        node=new Var(<Property*?>node,lb,ub,step)
                        self.line.vars.append(string_vars_name, node,True)
            elif token.kind=="ID":
                self.term_names=[]
                name = self.current_token.value
                position=self.position()
                self.eat(".")
                param_token=self.current_token
                kind,index=self.parameter()
                
                if kind == KWD:
                    if cpp_count(self.line.line[position].elem.keywords.begin(), self.line.line[position].elem.keywords.end(),index)==0:
                        raise KeyError(f"{param_token.value} is invalid keyword at line {param_token.line} column {param_token.column}!")

                var_name="".join(self.term_names )
                lb,ub,step =self.get_bounds()
                self.ordered_var_positions.append(position )
                if var_name in self.vars_name:
                    raise SyntaxError(f"Redefined Var {var_name} at line {token.line} column {token.column}!")
                self.vars_name.append(var_name)    
                string_vars_name=var_name.encode("utf8")

                node = new Property(kind, position,index, self.line.line, self.line.globals)
                node=new Var(<Property*?>node,lb,ub,step)
                self.line.vars.append(string_vars_name, node,True)
                self.vars_elem_name.append(name) 
            else:
                raise SyntaxError(f"Unknown syntax {token.value}  at line {token.line} column {token.column}!")
        elif self.current_token.kind=="ID":
            token=self.current_token
            name = self.current_token.value
            self.term_names=[]
            position=self.position()
            self.eat(".")   
            param_token=self.current_token
            kind,index=self.parameter()

            if kind == KWD:
                if cpp_count(self.line.line[position].elem.keywords.begin(), self.line.line[position].elem.keywords.end(),index)==0:
                    raise KeyError(f"{param_token.value} is invalid keyword at line {param_token.line} column {param_token.column}!")

            var_name="".join(self.term_names )
            self.ordered_var_positions.append(position )
            if var_name in self.vars_name:
                raise SyntaxError(f"Redefined Var {var_name} at line {token.line} column {token.column}!")
            self.vars_name.append(var_name)   
            string_vars_name=var_name.encode("utf8")
            self.eat("DELAY")

            node= new Property(kind, position,index, self.line.line, self.line.globals)
            node =new Var(<Property*?>node,self.expr())
            self.line.vars.append(string_vars_name,node,False)
            self.vars_elem_name.append(name) 

            # self.line.vars.ordered_depend_var_names.push_back( string_vars_name )
            # self.line.vars.num_dependent_vars+=1
        else:
            raise SyntaxError(f"Unknown syntax {self.current_token.value} at line {self.current_token.line} column {self.current_token.column}!")


    cdef void set_constraint(self)except*:
        
        cdef:
            Token token
            AST* node=NULL
            string name
        self.eat("CONSTRAINT")
        self.eat(",")
        self.eat("EXPR")
        self.eat("DELAY" )
        self.term_names=[]
        node=self.expr()
        name = "".join(self.term_names).encode("utf8")
        self.line.constraints.append(name,node)

        for value in self.term_names:
            if value in chromatic_terms:
                self.line.constraints.time_consuming_terms[CHROMATIC_TERMS]=True
            elif  value in driving_terms:
                self.line.constraints.time_consuming_terms[DRIVING_TERMS]=True
            elif  value in da_tracking_terms:
                self.line.constraints.time_consuming_terms[DA_TRACKING_TERMS]=True
            elif  value in ma_tracking_terms:
                self.line.constraints.time_consuming_terms[MA_TRACKING_TERMS]=True
            elif  value in monitor_off_momentum_terms:
                self.line.constraints.time_consuming_terms[MONITOR_OFF_MOMENTUM_TERMS]=True
            elif  value in off_momentum_sum_terms:
                self.line.constraints.time_consuming_terms[ OFF_MOMENTUM_SUM_TERMS ]=True
        # self.line.constraints.values.push_back(node)
        # self.line.constraints.num_constraint+=1



    cdef void set_optimize(self)except*:
        cdef:
            Token token
            AST* node=NULL
            string name
            bint ismin=True
            
        self.eat("OPTIMIZE")
        self.eat(",")
        self.eat("EXPR")
        self.eat("DELAY" )
        self.term_names=[]
        node=self.expr()
        name = "".join(self.term_names).encode("utf8")
        if self.current_token.kind==",":
            self.eat(",")
            if self.current_token.kind=="MINIMIZE":
                ismin=True
                self.eat("MINIMIZE")
                #self.line.optima.minormax.push_back(1.0)
            elif self.current_token.kind=="MAXIMIZE":
                ismin=False
                self.eat("MAXIMIZE")
                #self.line.optima.minormax.push_back(-1.0)
        else:
            ismin=True
            #self.line.optima.minormax.push_back(1.0)
        self.line.optima.append(name, node,ismin )  
        for value in self.term_names:
            if value in chromatic_terms:
                self.line.optima.time_consuming_terms[CHROMATIC_TERMS]=True
            elif  value in driving_terms:
                self.line.optima.time_consuming_terms[DRIVING_TERMS]=True
            elif  value in da_tracking_terms:
                self.line.optima.time_consuming_terms[DA_TRACKING_TERMS]=True
            elif  value in ma_tracking_terms:
                self.line.optima.time_consuming_terms[MA_TRACKING_TERMS]=True
            elif  value in monitor_off_momentum_terms:
                self.line.optima.time_consuming_terms[MONITOR_OFF_MOMENTUM_TERMS]=True
            elif  value in off_momentum_sum_terms:
                self.line.optima.time_consuming_terms[ OFF_MOMENTUM_SUM_TERMS ]=True
        #self.line.optima.values.push_back(node)
        #self.line.optima.num_optima+=1


    
    cdef void set_chromaticity(self)except*:
        cdef:
            Token token
            AST* node=NULL
            double value=0,coeff=1, pos_or_neg=1
            bytes name
            vector[size_t] position
            vector[double] factor
            str coknob_name
            string corr_name
            int kind=SEXTUPOLE
        self.eat("CHROM")
        self.eat(",")
        if self.current_token.kind=="AIM_DQX":
            self.eat("AIM_DQX")
            self.eat("ASSIGN" )
            node=self.expr()
            value=node.calc()
            del node
            node=NULL
            self.line.chrom_corrector.aim_dQx=value
            self.eat(",")
            self.eat("KNOB")
            self.eat("ASSIGN")
            token=self.current_token
            self.eat("ID")
            name =token.value.encode("utf8")
            if token.value not in self.elems_name or (<CppElement*>self.line.elems.at(name).at(0)).kind!=kind:
                raise TypeError(f"KNOB must be name of Sextupole element at line {token.line} column {token.column}!")
            
            self.line.chrom_corrector.corrector1=name
            if self.line.chrom_corrector.iscorr1 ==True:
                raise SyntaxError(f"Redefined AIM_DQX at line {token.line} column {token.column}!")
            self.line.chrom_corrector.iscorr1=True
            position=self.eids[token.value]
            factor = len(self.eids[token.value])*[1]
            self.line.chrom_corrector.position1.insert(self.line.chrom_corrector.position1.end(),position.begin(),position.end() )
            self.line.chrom_corrector.coeff1.insert(self.line.chrom_corrector.coeff1.end(),factor.begin(),factor.end() )
            self.line.chrom_corrector.unique_position1.push_back(self.eids[token.value][0] )
            self.line.chrom_corrector.unique_coeff1.push_back(1.0)

        elif self.current_token.kind=="AIM_DQY":
            self.eat("AIM_DQY")
            self.eat("ASSIGN" )
            node=self.expr()
            value=node.calc()
            del node
            node=NULL
            self.line.chrom_corrector.aim_dQy=value
            self.eat(",")
            self.eat("KNOB")
            self.eat("ASSIGN")
            token=self.current_token
            self.eat("ID")
            name =token.value.encode("utf8")
            if token.value not in self.elems_name or (<CppElement*>self.line.elems.at(name).at(0)).kind != kind:
                raise TypeError(f"KNOB should be name of Sextupole element at line {token.line} column {token.column}!")
            
            self.line.chrom_corrector.corrector2=name
            if self.line.chrom_corrector.iscorr2 ==True:
                raise SyntaxError(f"Redefined AIM_DQY at line {token.line} column {token.column}!")
            self.line.chrom_corrector.iscorr2=True
            position=self.eids[token.value]
            factor = len(self.eids[token.value])*[1]
            self.line.chrom_corrector.position2.insert(self.line.chrom_corrector.position2.end(),position.begin(),position.end() )
            self.line.chrom_corrector.coeff2.insert(self.line.chrom_corrector.coeff2.end(),factor.begin(),factor.end() )
            self.line.chrom_corrector.unique_position2.push_back(self.eids[token.value][0] )
            self.line.chrom_corrector.unique_coeff2.push_back(1.0)

        elif self.current_token.kind=="ID":
            token=self.current_token
            name=token.value.encode("utf8")
            if token.value not in self.elems_name or (<CppElement*>self.line.elems.at(name).at(0)).kind !=kind:
                raise TypeError(f"KNOB should be name of Sextupole element at line {token.line} column {token.column}!")
            coknob_name=token.value
            self.eat("ID")
            self.eat("DELAY")

            if self.current_token.kind in ("+","-"):
                token=self.current_token
                self.eat(self.current_token.kind)
                pos_or_neg=float(f"{token.value}1")
            if self.current_token.kind=="NUMBER":
                coeff=float(f"{self.current_token.value}")
                self.eat("NUMBER")
                self.eat("*")
            coeff=pos_or_neg*coeff
            position=self.eids[coknob_name]
            factor = [coeff for i in self.eids[coknob_name ]]
            token=self.current_token
            self.eat("ID")
            corr_name=token.value.encode("utf8")
            
            if <string>coknob_name.encode("utf8")== self.line.chrom_corrector.corrector1 or <string>coknob_name.encode("utf8")== self.line.chrom_corrector.corrector2:
                raise SyntaxError(f"Co-KNOB {coknob_name} shouldn't be the same as {token.value} at  line {token.line} column {token.column}!")
            if corr_name== self.line.chrom_corrector.corrector1:
                self.line.chrom_corrector.position1.insert(self.line.chrom_corrector.position1.end(),position.begin(),position.end() )
                self.line.chrom_corrector.coeff1.insert(self.line.chrom_corrector.coeff1.end(),factor.begin(),factor.end() )
                self.line.chrom_corrector.unique_position1.push_back(self.eids[coknob_name][0] )
                self.line.chrom_corrector.unique_coeff1.push_back(coeff )
                
            elif corr_name== self.line.chrom_corrector.corrector2:
                #if coknob_name == token.value:
                    #raise SyntaxError(f"Co-KNOB {coknob_name} shouldn't be the same as {token.value} at  line {token.line} column {token.column}!")
                self.line.chrom_corrector.position2.insert(self.line.chrom_corrector.position2.end(),position.begin(),position.end() )
                self.line.chrom_corrector.coeff2.insert(self.line.chrom_corrector.coeff2.end(),factor.begin(),factor.end() )
                self.line.chrom_corrector.unique_position2.push_back(self.eids[coknob_name][0] )
                self.line.chrom_corrector.unique_coeff2.push_back(coeff )
            else:
                raise SyntaxError(f"Co-KNOB should be bind to KNOB at  line {token.line} column {token.column}!")
                



    cdef double eval(self,str code)except*:
        cdef AST*   node
        cdef double value
        if not self.isdatabase:
            raise RuntimeError('Parser is not linked to database!')
        self.lexer.tokenize(code)
        self.current_token = self.lexer.get_current_token()
        node =self.expr()
        value=node.calc()
        del node
        node=NULL
        return value



    def parse(self):
        cdef AST*   node
        cdef bytes name
        value=None
        values=[]
        # print("Parser.parse:0")
        if not self.isdatabase:
            raise RuntimeError('Parser is not linked to database!')

        # self.lexer.tokenize(code)
        # self.tokens.append(self.lexer.tokens)

        # print("Parser.parse:1")
        self.current_token = self.lexer.get_current_token()
        # node = self.expr()
        # print("Parser.parse:2")
        while self.current_token.kind != "EOF":
            # print(self.current_token.kind,":",self.current_token.value)
            # try:
            self.term_names=[]
            if self.current_token.kind =="ID" and self.lexer.check_next_token() in ("DELAY","ASSIGN"):
                self.set_Identity()
                self.eat("END")
            elif self.current_token.kind =="VAR":
                self.set_variable()
                self.eat("END")
            elif self.current_token.kind =="CONSTRAINT":
                self.set_constraint()
                self.eat("END")
            elif self.current_token.kind =="OPTIMIZE":
                self.set_optimize()
                self.eat("END")
            elif self.current_token.kind =="CHROM":
                self.set_chromaticity()
                self.eat("END")
            else:
                # raise SyntaxError(f"Unknown syntax at line {self.current_token.line} column {self.current_token.column}!")
                node=self.expr()
                value=node.calc()
                values.append(value )
                del node
                if(self.current_token.kind=="END"):
                    self.eat("END")
        if self.ordered_var_positions!=[]:
            self.line.ordered_changed_position=list(set(self.ordered_var_positions))
            #print(self.line.ordered_changed_position)
        if len(values)>1:
            return values
        else:
            return value

        # self.eat("EOF")
    
