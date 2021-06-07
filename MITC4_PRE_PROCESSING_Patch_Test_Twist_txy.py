from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math 
from Pre_Processing_Functions import*
import Plotting 



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
    E=2.1e6
    v=0.3
    Material=[E,v]

#Nodes Coordinates
    Nodes_Coord=np.array([[0, 10, 0],
                      [0, 0, 0],
                      [4, 7, 0],
                      [2,2,0],
                      [8,7,0],
                      [8,3,0],
                      [10,10,0],
                      [10,0,0]])
     
#Nodes Thickness
    Nodes_Thickness=np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])

#Nodes Directional Vector:
    Nodes_Vn=np.array([[0,0,1], 
                       [0,0,1],
                       [0,0,1],
                       [0,0,1],
                       [0,0,1],
                       [0,0,1],
                       [0,0,1],
                       [0,0,1]])

#Element Connectivity:
    Element_Connect = np.array([[3,1,2,4],
                                [7,1,3,5],
                                [5,3,4,6],
                                [6,4,2,8],
                                [7,5,6,8]])
    
    
    
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

    
    #Degrees of Freedom: 
    bcdof=[]
    BC_nodes=[]
    
    bcdof=BC(1,[0,1,2],bcdof)
    BC_nodes.append(1)
    
    bcdof=BC(2,[0,1,2],bcdof)
    BC_nodes.append(2)
    
    bcdof=BC(8,[0,1,2],bcdof)
    BC_nodes.append(8)
    
    bcdof=BC(3,[0,1],bcdof)
    BC_nodes.append(3)
    
    bcdof=BC(4,[0,1],bcdof)
    BC_nodes.append(4)
    
    bcdof=BC(5,[0,1],bcdof)
    BC_nodes.append(5)
    
    bcdof=BC(6,[0,1],bcdof)
    BC_nodes.append(6)
    
    bcdof=BC(7,[0,1],bcdof)
    BC_nodes.append(7)
    
    #Set all thetaz degrees of freedom to zero. 
    for i in range(np.shape(Nodes_Coord)[0]):
        bcdof=BC(i+1,[5],bcdof)
        BC_nodes.append(i+1)
        
    #Always afert finishing applying BCs:
    bcdof=np.unique(bcdof)
    BC_nodes=np.unique(BC_nodes)
    
    
    #Forces:
    #Create R Vector:
    # R=R_Vector(np.shape(Nodes_Coord)[0])
    # #Fill the vector with the forces:
    # R=Forces(7,[0],R,[1000])
    # R=Forces(8,[0],R,[1000])
    
    #Forces:
    #First, define the R vector, an array to save the nodes in which the forces are applied,
    #and an array were each force component is stored  
    R=R_Vector(np.shape(Nodes_Coord)[0])
    F_nodes=[]
    F_components=[]
    
    #Defining a force via a node (if using this option, append the node to the 
    #array F_nodes):
    R=Forces(7,2,R,[-1000])
    F_nodes.append(7)
    F_components.append(-2)


    #Plot the forces vectors
    #Use the grafvector function: inputs are the axis in which the figure is plotted,
    #the nodal coordinates, the node in which the force is applied, and the force component
    #which can be 0,1,2,3,4,5 as stated previously. 
    for i in range(len(F_nodes)):
        Plotting.grafvector(ax,Nodes_Coord,F_nodes[i],F_components[i])

    return(Material,Nodes_Coord,Nodes_Thickness,Nodes_Vn,Element_Connect,bcdof,R)
    
    
'''
#Nodes Rotation Vector Vn1 and Vn2:
    #ESTA PARTE PASA AL PROCESAMIENTO
Vec_Vn1_Vn2=np.zeros([np.shape(Nodes_Coord)[0],2,3])



c1=0
for i in Nodes_Vn:
    Vec_Vn1_Vn2[c1,0]=vec_rot(i)[0] #Vn1
    Vec_Vn1_Vn2[c1,1]=vec_rot(i)[1] #Vn2 
    c1=c1+1
H=H_Matrix(0,0,0,Vec_Vn1_Vn2,Nodes_Thickness)   
    



dX_dr,dX_ds,dX_dt=dX(0,0,0,Nodes_Coord,Nodes_Vn,Nodes_Thickness)    
dU_dr,dU_ds,dU_dt=dU(0,0,0,Vec_Vn1_Vn2,Nodes_Thickness)    

gi=gi_covariant(dX_dr,dX_ds,dX_dt)
gicontra=gi_contra(gi)
'''