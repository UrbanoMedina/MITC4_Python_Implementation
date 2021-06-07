from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def print_Geometry(Element_Connect,Nodes_Coord,ax,c):
        # fig = plt.figure(Fig_Label)
        # #Create a 3d cart. axis
        # ax = fig.add_subplot(111, projection='3d')
    #Plot the lines which compose the element:
        for i in range(Element_Connect.shape[0]): #Number of elements
            El_Con=Element_Connect[i]
            # for j in range(Element_Connect.shape[1]): #Number of nodes per element    
            #     if j==3:
            #         ax.plot([Nodes_Coord[El_Con-1][3][0],Nodes_Coord[El_Con-1][0][0]],
            #                 [Nodes_Coord[El_Con-1][3][1],Nodes_Coord[El_Con-1][0][1]],
            #                 [Nodes_Coord[El_Con-1][3][2],Nodes_Coord[El_Con-1][0][2]],color="k",
            #                 linewidth=0.6)
            #     else:
            #         ax.plot([Nodes_Coord[El_Con-1][j][0],Nodes_Coord[El_Con-1][j+1][0]],
            #                 [Nodes_Coord[El_Con-1][j][1],Nodes_Coord[El_Con-1][j+1][1]],
            #                 [Nodes_Coord[El_Con-1][j][2],Nodes_Coord[El_Con-1][j+1][2]],color="k",
            #                 linewidth=0.6)
            x=Nodes_Coord[Element_Connect[i]-1,0]
            y=Nodes_Coord[Element_Connect[i]-1,1]
            z=Nodes_Coord[Element_Connect[i]-1,2]
            verts = [list(zip(x,y,z))]
            ax.add_collection3d(Poly3DCollection(verts,color=c,alpha=0.2))

def node_labels(ax,Nodes_Coord):
    c1=0
    for i in Nodes_Coord:
        ax.text(i[0],i[1],i[2],'U{}'.format(c1+1))
        c1=c1+1
    
def grafvector (ax,Nodes_Coord,node,component):
    vec_comp=[0,0,0]
    
    if component < 0 :
        m=-1
    else:
        m=1
    
    component = np.abs(component)
    if component == 3 or component == 0:
        if np.max(Nodes_Coord[:,0]) == 0:
            vec_comp=[0.015*m,0,0]
        else:
            vec_comp=[np.max(Nodes_Coord[:,0])*0.15*m,0,0]
    
    elif component == 4 or component == 1:
        if np.max(Nodes_Coord[:,1]) == 0:
            vec_comp=[0,0.015*m,0]
        else:
            vec_comp=[0,np.max(Nodes_Coord[:,1])*0.15*m,0]
            
    elif component == 5 or component == 2:
        if np.max(Nodes_Coord[:,2]) == 0:
            vec_comp=[0,0,0.015*m]
        else:
            vec_comp=[0,0,np.max(Nodes_Coord[:,2])*0.15*m]           
    
    if component == 3 or component == 4 or component == 5:
        co='g'
    else:
        co='tab:red'
    ax.quiver(Nodes_Coord[node-1][0],Nodes_Coord[node-1][1],Nodes_Coord[node-1][2],
              vec_comp[0],vec_comp[1],vec_comp[2],pivot='tip',color=co)
    
#def BC_Node_Plot()