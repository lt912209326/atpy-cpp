import re

cdef class Token:
    def __init__(self, str kind0, str value, int line , int column):
        self.kind = kind0 
        #TWS,KWD, LOC, GLB,FUN,ID,ASSIGN,DELAY,
        self.value = value
        self.line = line
        self.column = column


cdef class Lexer:
    def __init__(self):
        self.count = 0
        self.line_num = 0
        # self.elems_index = elems_index
    
    
    cdef Token get_current_token(self):
        return self.tokens[self.count]
    
    
    cdef get_next_token(self):
        self.count += 1
        if self.count <self.token_num:
            return self.tokens[self.count]
        else:
            raise ValueError('No end charactor for expression!')
    
    
    cdef str check_next_token(self):
        if self.count+1<self.token_num:
            return (<Token>self.tokens[self.count+1]).kind
        else:
            raise ValueError('Token error, out of memmory when indexing tokens !')
           
    

    cdef int tokenize(self, str code):
        
        self.tokens=[]
        self.count=0
        reserved_words = ["VAR",  "CONSTRAINT", "OPTIMIZE", "CHROM","NAME","EXPR", "LOWER", "UPPER", "STEP","MINIMIZE", "MAXIMIZE", "AIM_DQX","AIM_DQY", "KNOB"  ]
        keywords = KWD_INDEX.keys()
        twiss = TWS_INDEX.keys()
        local = LOC_INDEX.keys()
        glbs = GLB_INDEX.keys()

        funs = ["ABS", "MIN", "MAX", "SQRT", "DIM", "SIN",    "COS",    "SINH",   "COSH", "EXP"]
        operators= ["+","-","*", "/", "//", "%", "**" ]
        bracket = ["(", ")" ]
        sign = [",","@","."]
    
    
        
        token_specification = [
            ("COMMENT",  r"^[ \t]*#.+"),   # Comment 
            ('DELAY',    r':='),           # Delayed Assignment operator
            ('ASSIGN',   r'='),           # Assignment operator
            ('END',      r';'),            # Statement terminator
            ("REGEX",    r'\$[A-Za-z0-9_.^*\]\[+?|\\}{-]+\$'), 
            ('ID',       r'[A-Za-z][A-Za-z0-9_]*@?[A-Za-z0-9]*'),    # Identifiers
            ('NUMBER',   r'(\d+(\.\d+)?|\.\d+)([eE][-+]?\d+)?'),  # Integer or decimal number
            ('OP',       r'[*/]{2}|[+\-*/%]'),      # Arithmetic operators
            ('NEWLINE',  r'\n'),           # Line endings
            ('SKIP',     r'[ \t]+'),       # Skip over spaces and tabs
            ("PAREN",    r"[]()[]"),
            ("SIGN",     r"[,@.]"),
            ('MISMATCH', r'.'),            # Any other character
        ]
        tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
        line_num = 1
        line_start = 0

        self.codes=code.split("\n")
        
        for mo in re.finditer(tok_regex, code,re.M):
            kind = mo.lastgroup
            value = mo.group()
            column = mo.start() - line_start
            if kind == 'NUMBER':
                value = value
            elif kind == 'ID' and value in reserved_words:
                kind = value
            elif kind == 'ID' and value in keywords:
                kind = "KWD"
            elif kind == 'ID' and value in twiss:
                kind = "TWS"
            elif kind == 'ID' and value in local:
                kind = "LOC"
            elif kind == 'ID' and value in glbs:
                kind = "GLB"
            elif kind == 'ID' and value in funs:
                kind = "FUN"
            elif kind == 'PAREN':
                kind = value
            elif kind == 'OP':
                kind = value
            elif kind == "SIGN":
                kind=value
            elif kind == 'NEWLINE':
                line_start = mo.end()
                line_num += 1
                continue
            elif kind == 'SKIP' or kind == "COMMENT":
                continue
            elif kind == 'MISMATCH':
                raise RuntimeError(f'{value!r} unexpected on line {line_num}')
            # yield Token(kind, value, line_num, column)
            self.tokens.append( Token(kind,  value, line_num, column))
        self.tokens.append( Token("EOF", "EOF", line_num+1, 0))
        self.token_num = len(self.tokens)

statements = '''
VAR , NAME=$ABE.+[12]E$[2].betax;
CONSTRAINT, EXPR=ABS(MK0[0].alhpax);
CHROMX:=SD1[0].betax*SD1[0].etax;

'''
