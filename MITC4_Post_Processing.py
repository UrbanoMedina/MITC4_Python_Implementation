import Plotting
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#MITC4 POSTPROCESSING# 

#This function is defined in order to review all the solution variables and
#graph the deformed mesh. 

#Therefore, it receives as an input the output of the processing function. 

    
def Graph(Mag,U,Nodes_Coord,Element_Connect):
    #--------------------------------------------#
    #--SUM THE DISPLACEMENTS TO THE NODAL COORD--#
    #--------------------------------------------#
    
    #Create an array were the displaced nodal coord will be stored:
    nnode= np.shape(Nodes_Coord)[0]
    Nodes_Coord_Disp=np.zeros([nnode,3])
    
    for i in range(len(U)):
        if U[i]<1e-10 and U[i]>-1e-10:
            U[i]=0
        
    
    for i in range(nnode):
        disp_comp=np.zeros(3)
        for j in range(3): #Cycle to obtain for the node the 3 displacement components
            disp_comp[j]=U[i*6+j]*Mag
        Nodes_Coord_Disp[i]=Nodes_Coord[i]+disp_comp        
        
    # #--------------------------------------------#
    # #-----PRINT ORIGINAL VS DISPLACED FIGURE-----#
    # #--------------------------------------------#
    
    #Create the figure to plot the geometry:
    fig = plt.figure(1)
     #Create a 3d cart. axis
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim(np.min(Nodes_Coord_Disp[:,1]),np.max(Nodes_Coord_Disp[:,1]))
    ax.set_xlim(np.min(Nodes_Coord_Disp[:,0]),np.max(Nodes_Coord_Disp[:,0]))
    ax.set_zlim(np.min(Nodes_Coord_Disp[:,2]),np.max(Nodes_Coord_Disp[:,2]))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #Plot each element
    Plotting.print_Geometry(Element_Connect,Nodes_Coord,ax,'b')
    Plotting.print_Geometry(Element_Connect,Nodes_Coord_Disp,ax,'r')
    
    return(fig)