import matplotlib.pyplot as plt 
import numpy as np

def plot_tune(nux0,nuy0, max_order=5):
    int_nux0,int_nuy0 = int(nux0//1),int(nuy0//1)
    
    x1 = [[int_nux0,int_nuy0 ],
          [int_nux0+1,int_nuy0 ],
          [int_nux0,int_nuy0+1 ],
          [int_nux0+1,int_nuy0+1 ]]
    coeffs = []
    lines=[]
    for m in range(-max_order, max_order+1):
        for n in range(-(max_order-abs(m) ), max_order - abs(m)+1):
            p=[ m*x[0]+n*x[1] for x in x1]
            pmin=min( p )
            pmax=max( p )
            for k in range(pmin,pmax+1 ):
                if n == 0 :
                    if m == 0: continue
                    lines.append( [[k/m, k/m],[ int_nuy0,int_nuy0+1 ],abs(m)+abs(n) ] )
                else:
                    lines.append( [[int_nux0,int_nux0+1 ],[(k-m*int_nux0)/n,(k-m*int_nux0-m)/n ],abs(m)+abs(n) ]  )
    fig,ax=plt.subplots(1)
    colors = [ ["k",1/0.2],
              ["r",1/0.4],["b",1/0.6], ["g",1/0.8],['c',1/1],["y",1/2]]
    lines = sorted(lines, key=lambda x: x[2], reverse=True )
    lns = []
    orders = []
    for line in lines:
        order=line[2]
        ln, = ax.plot(line[0],line[1],"-",color=colors[order-1][0], linewidth=colors[order-1][1],label=f"{order}" )
        if order not in orders:
            orders.append(order)
            lns.append(ln)
        
    ax.set_xlim(int_nux0,int_nux0+1)
    ax.set_ylim(int_nuy0,int_nuy0+1)
    lns.reverse()
    kargs = dict(mode='expand', ncol=min(max_order,6), borderaxespad=1,  loc=3,  bbox_to_anchor=(0, 1,1,0.2), frameon=False)
    ax.legend(handles=lns,  **kargs)
    fig.tight_layout()
    return fig,ax,lines,kargs