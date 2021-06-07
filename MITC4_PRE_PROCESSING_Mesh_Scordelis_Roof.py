from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math 
from Pre_Processing_Functions import*
import Plotting 
import Mesh_Rect_Cyl


# from dX_and_dU import*
# from Coordinate_Bases import*
# from Matrices_C_B_H import*

############## ------------------------------  #########################
##############  NODAL_POINT/ELEMENT READ IN    #########################
############## ------------------------------  #########################
'''
In this script the following information regarding nodal points and 
elements is defined:
    
    Material properties: Elasticity Modulus - E
                         Poisson Ratio - v
    
    Nodal Coordinates: Matrix with nodal information. Each row correspond to 
                       a node and each column to a x y z component. 
        Example:
            MatrizCoordenadasNodales=np.array([[7,8,3],
                                               [2,9,1],
                                               [1,1,1],
                                               [6,2,2]])
            
    Thickness Direction Vector: Vn 
        Define a matrix in which each row corresponds to a node normal vector and
        each column to the x y z components respectively 
        Example: 
            MatrizVectoresNodales=np.array([[0,0,1], #Node 1 Vn vector components
                                            [0,1,1],
                                            [1,0,1],
                                            [0,-1,1]])  #Node 4 Vn vector components
            
            NOTE: This vectors MUST have a magnitude equal to 1. If not, use the function 
                Normalize_Vector(Matrix_Vectors) in which the input is the matrix of nodal vectors.
    
    Element Connectivity: Each row correspond to an element.
                          Each column corresponds to a node. 
       Example: 
           conec2D =
                      [3,     1,     2,     4; #Element 1 composed of: node 3,1,2,4
                       7,     1,     3,     5;
                       5,     3,     4,     6;
                       6,     4,     2,     8;
                       7,     5,     6,     8] #Element 4
           
    Define which degrees of freedom are restriced and on which nodes.
        Use BC function:
            Input: (node,component,bcdof)
                Node: node in which the BC is defined. 
                Component: define which components are to be restricted.
                    [u,v,w,thetax,thetay,thetaz] are the dof.
                    [0,1,2,3,4,5] each number corresponds to a degree of freedom, so
                                  input the numbers which corresponds to the restricted dofs.
                bcdof: array in which every element corresponds to the position of 
                       restricted dofs in the nodal displacement vector. 
                EXAMPLE:
                     Node 2, components u,W restricted 
                     BC(2,[0,2],bcdofs)
                     Output: bcdofs= [6,8]
                     Means elements in the positions 6 and 8 of the nodal displacement vector
                        are restricted.

    External Forces:
        Define the R vector on the KU=R Matrix Equation
            
        First: create the R vector with all entris equal to zero. Number of elements 
            equal the number of degrees of freedom in the system. 
                Use the R_Vector function. Input: Number of nodes in the system.
            
        Second: input the forces in the R_Vector. Use the Force() function
                Forces(node,component,vect_R,mag)
                    node: node in which the force is applied. 
                    component: [Fx,Fy,Fz,Mx,My,Mz] == [0,1,2,3,4,5]
                    To set a force component then input the number equivalent to the force. 
                    vect_R: R nodal forces vector.
                    mag: magnitude of the force 
        
'''
def Node_Element_Info():
                 
#MATERIAL PROPERTIES 
    E=3.0e6
    v=0
    Material=[E,v]


    #GEOMETRY DEFINITION#
    #In cylindrical_mesh(radius,angle(in degrees),height,
    #number of vertical divisions,number of horizontal divisions)
    Nodes_Coord,Element_Connect,Nodes_Vn,Nodes_Thickness=Mesh_Rect_Cyl.cylindrical_mesh(300,[90,50],300,20,20,3)
    #Normalize the Vn vector
    Nodes_Vn=Vector_Normalize(Nodes_Vn)
    
    #Check the geometry by plotting it. 
    fig = plt.figure(0)
    #Create a 3d cart. axis
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim(np.min(Nodes_Coord[:,1]),np.max(Nodes_Coord[:,1]))
    ax.set_xlim(np.min(Nodes_Coord[:,0]),np.max(Nodes_Coord[:,0]))
    ax.set_zlim(np.min(Nodes_Coord[:,2]),np.max(Nodes_Coord[:,2]))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #Plot each element
    Plotting.print_Geometry(Element_Connect,Nodes_Coord,ax,'b')
    #Add to the graph node labels at nodal points:
    Plotting.node_labels(ax,Nodes_Coord)
    
    
    #ESTABLISH THE BCS ON A LINE#
    #Using the functions BC_Rect_Line(const_coord,const_value,restric_comp,Nodes_Coord,bcdof,nodes):
        #Each input variable is explained in the script Pre_Processing_Functions
        
    #bcdof:array in which the poisition of restrictied dofs in the U system vector is specified. 
    #nodes:array in which the nodes were BCs are applied are specified.
    BC_nodes=[]
    bcdof=[]
    #The inputs of the function: BC_Rect_Line(const_coord,const_value,restric_comp,Nodes_Coord,bcdof,nodes):
    bcdof,BC_nodes=BC_Cylin_Line([0,2],[300,0],[0,2],Nodes_Coord,bcdof,BC_nodes) #Corresponds to DA
    bcdof,BC_nodes=BC_Cylin_Line([0,2],[300,300],[1,3,5],Nodes_Coord,bcdof,BC_nodes) #Corresponds to BC
    bcdof,BC_nodes=BC_Cylin_Line([0,1],[300,90],[0,4,5],Nodes_Coord,bcdof,BC_nodes) #Corresponds to CD


    #Always afert finishing applying BCs:
    bcdof=np.unique(bcdof)
    BC_nodes=np.unique(BC_nodes)

        
    #ESTABLISH THE FORCES
    #Can be defined either by specifying a node or a coordinate:
    
    #First, define the R vector, an array to save the nodes in which the forces are applied,
    #and an array were each force component is stored  
    R=R_Vector(np.shape(Nodes_Coord)[0])
    F_nodes=[]
    F_components=[]
    
    #Defining a force via a node (if using this option, append the node to the 
    #array F_nodes):
    # R=Forces(1,-2,R,613.6) 
    # F_nodes.append(1)
    # F_components.append(-2)

    #Define a force via a coordinate:
    R,F_nodes,F_components=Force_Cylin_Mesh([300,90,300],-2,1000,Nodes_Coord,R,
                                           F_nodes,F_components)
        
    # R=Forces(2,-2,R,1227.2)
    # F_nodes.append(2)
    # F_components.append(-2)
    
    # R=Forces(3,-2,R,1227.2)
    # F_nodes.append(3)
    # F_components.append(-2)
    
    # R=Forces(4,-2,R,1227.2)
    # F_nodes.append(4)
    # F_components.append(-2)
    
    # R=Forces(5,-2,R,613.6)
    # F_nodes.append(5)
    # F_components.append(-2)
    
    # R=Forces(6,-2,R,1227.2)
    # F_nodes.append(6)
    # F_components.append(-2)
    
    # R=Forces(7,-2,R,2454.4)
    # F_nodes.append(7)
    # F_components.append(-2)
    
    # R=Forces(8,-2,R,2454.4)
    # F_nodes.append(8)
    # F_components.append(-2)
    
    # R=Forces(9,-2,R,2454.4)
    # F_nodes.append(9)
    # F_components.append(-2)
    
    # R=Forces(10,-2,R,1227.2)
    # F_nodes.append(10)
    # F_components.append(-2)
    
    # R=Forces(11,-2,R,1227.2)
    # F_nodes.append(11)
    # F_components.append(-2)
    
    # R=Forces(12,-2,R,2454.4)
    # F_nodes.append(12)
    # F_components.append(-2)
    
    # R=Forces(13,-2,R,2454.4)
    # F_nodes.append(13)
    # F_components.append(-2)
    
    # R=Forces(14,-2,R,2454.4)
    # F_nodes.append(14)
    # F_components.append(-2)
    
    # R=Forces(15,-2,R,1227.2)
    # F_nodes.append(15)
    # F_components.append(-2)
    
    # R=Forces(16,-2,R,1227.2)
    # F_nodes.append(16)
    # F_components.append(-2)
    
    # R=Forces(17,-2,R,2454.4)
    # F_nodes.append(17)
    # F_components.append(-2)
    
    # R=Forces(18,-2,R,2454.4)
    # F_nodes.append(18)
    # F_components.append(-2)
    
    # R=Forces(19,-2,R,2454.4)
    # F_nodes.append(19)
    # F_components.append(-2)
    
    # R=Forces(20,-2,R,1227.2)
    # F_nodes.append(20)
    # F_components.append(-2)
    
    # R=Forces(21,-2,R,613.6)
    # F_nodes.append(21)
    # F_components.append(-2)
    
    # R=Forces(22,-2,R,1227.2)
    # F_nodes.append(22)
    # F_components.append(-2)
    
    # R=Forces(23,-2,R,1227.2)
    # F_nodes.append(23)
    # F_components.append(-2)
    
    # R=Forces(24,-2,R,1227.2)
    # F_nodes.append(24)
    # F_components.append(-2)
    
    # R=Forces(25,-2,R,613.6)
    # F_nodes.append(25)
    # F_components.append(-2)
    

    
    
    #Define a force via a coordinate:
    # R,F_nodes,F_components=Force_Cylin_Mesh([300,90,300],-2,0.25,Nodes_Coord,R,
    #                                        F_nodes,F_components)
    
    #Plot the forces vectors
    #Use the grafvector function: inputs are the axis in which the figure is plotted,
    #the nodal coordinates, the node in which the force is applied, and the force component
    #which can be 0,1,2,3,4,5 as stated previously. 
    for i in range(len(F_nodes)):
        Plotting.grafvector(ax,Nodes_Coord,F_nodes[i],F_components[i])

    return(Material,Nodes_Coord,Nodes_Thickness,Nodes_Vn,Element_Connect,bcdof,R)
    
