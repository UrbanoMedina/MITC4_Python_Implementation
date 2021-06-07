import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import Mesh_Rect_Cyl

#Function to set the BC giving a node, the restricted degree of freedom.
#The output will be the bcdof vector which will contain the positions of the 
#restricted degrees of freedom at the Unodal vector for the system.
def BC(nodo,componente,bcdof):
    for i in componente:
        bcdof.append((nodo-1)*6+i)
    return bcdof

#Creates the R vector given the number of points(nodes) in the system.
def R_Vector(num_puntos):
    vec_forces=np.zeros([num_puntos*6,1])
    return vec_forces

#Function to set a force component in the system R vector. Given a node and its,
#force component, the function set its magnitude in the corresponding position of the 
#R vector.

def Forces(node,component,vect_R,mag):
    if component < 0:
        m=-1
    else:
        m=1
    vect_R[(node-1)*6+np.absolute(component)]=mag*m

    return vect_R

#Function to normalize vectors.
#The input is an array in which every row is a vector and columns are the components. 
def Vector_Normalize(Vect_Array):
    VA=np.zeros([Vect_Array.shape[0],Vect_Array.shape[1]])
    for i in range(Vect_Array.shape[0]):
        Norm=np.linalg.norm(Vect_Array[i])
        for j in range(Vect_Array.shape[1]):
            VA[i,j]=Vect_Array[i,j]/Norm
    return(VA)


def BC_Rect_Line(const_coord,const_value,restric_comp,Nodes_Coord,bcdof,nodes):
    #A line in a rectangular mesh is given by a constant x or y value. 
    #To establish the line, it must be stated which component will be left as constant.
    #This will be specified in the const_coord variable in which only a 0 or 1 value must 
    #be entered.
    #0 will represent a constant x coordinate value. 
    #1 will represent a constant y coordinate value. 
    #const_value will be the value which will define the line in question. 
    #restric_comp will be an array with numbers from 0 to 5.
    #The input numbers will represent restricted degrees of freedom as follows:
    #0,1,2,3,4,5=Tx,Ty,Tz,Rx,Ry,Rz
    #So if numbers 1,3 are entered, Ty and Rx will be restricted.
    
    #In this loop all of the nodal coordinates are iterated. The coordinates that meet the
    #criteria defined will be sent to the BC function to apply the restrictions. 
    c=0
    for i in range(Nodes_Coord.shape[0]): 
        if np.abs(Nodes_Coord[:,const_coord][i]-const_value) < 1e-6:
            bcdof=BC(c+1,restric_comp,bcdof)
            nodes.append(c+1)
        c+=1
    return (bcdof,nodes)


def BC_Cylin_Line(const_coord,const_value,restric_comp,Nodes_Coord,bcdof,nodes):
    #For a cylindrical mesh a line will be defined by 2 of the 3 variables
    #(radius, height, angle) left as constant.
    #So const_coord will be an array were the 2 variables components will be specified.
    #Index 0 corresponds to radius
    #Index 1 corresponds to theta
    #Index 2 corresponds to height 
    #So const_coord=[0,2] means that the line will be defined by constant radius and height components.
    #const_value will be an array which will correspond to the constant values of the variables defined in
    #const_coord.
    #restric_comp defined in the BC_Rect_Line function
    
    #To start, Nodes_Coord coordinates will be transformed to cylindrical coordinates.
    Nodes_Coord_Cyl=Mesh_Rect_Cyl.Cart_to_Cyl(Nodes_Coord)
    
    #Then all nodal coordinates are compared with the criteria specified, if it 
    #satisfied, then the node will be sent to the BC function to establish the 
    #restrictions.
    c=0
    for i in range(Nodes_Coord.shape[0]):
        if np.abs(Nodes_Coord_Cyl[:,const_coord[0]][i]-const_value[0]) < 1e-6 \
        and np.abs(Nodes_Coord_Cyl[:,const_coord[1]][i]-const_value[1]) < 1e-6:
            bcdof=BC(c+1,restric_comp,bcdof)
            nodes.append(c+1)
        c+=1
    return (bcdof,nodes)


def Force_Rect_Mesh(F_coord,F_comp,Mag,Nodes_Coord,Rvector,F_nodes,F_components):
    #F_coord: array in which the coordinates of the location in which the force is applied is specified
        #F_coord=[xcoord,ycoord,zcoord]
    #F_comp: array in which the force component is specified. The component is specified 
    #with a number from 0 to 5 which have the following correspondance:
        #[Fx,Fy,Fz,Mx,My,Mz] == [0,1,2,3,4,5]
    #Mag: magnitude of the force vector. 
    #Nodes_Coord: coordinates of the nodes in the mesh.
    #Rvector: vector which corresponds to the system forces vector in the formula:
        #KU=R 
    c=0
    c2=0
    for i in range(Nodes_Coord.shape[0]): 
        if np.abs(Nodes_Coord[i][0] - F_coord[0]) < 1e-6 and \
            np.abs(Nodes_Coord[i][1] - F_coord[1]) < 1e-6 and \
            np.abs(Nodes_Coord[i][2] - F_coord[2]) < 1e-6:
                Rvector = Forces(c+1,F_comp,Rvector,Mag)
                F_nodes.append(c+1)
                F_components.append(F_comp)
                c2+=1
        c+=1
    if c2==0:
        print ('WARNING: SPECIFIED LOCATION DONT CORRESPOND TO ANY NODE')
    return (Rvector,F_nodes,F_components)

def Force_Cylin_Mesh(F_coord,F_comp,Mag,Nodes_Coord,Rvector,F_nodes,F_components):
    #F_coord: array in which the coordinates of the location in which the force is applied is specified.
    #The input of the coordinates will be in cylindrical coordinates:
        #F_coord=[radius,theta,height]
    #F_comp: array in which the force component is specified. The component is specified 
    #with a number from 0 to 5 which have the following correspondance:
        #[Fx,Fy,Fz,Mx,My,Mz] == [0,1,2,3,4,5]
    #Mag: magnitude of the force vector. 
    #Nodes_Coord: coordinates of the nodes in the mesh.
    #Rvector: vector which corresponds to the system forces vector in the formula:
        #KU=R 
      
        
    #So first, Nodes_Coord coordinates will be transformed to cylindrical coordinates.
    Nodes_Coord_Cyl=Mesh_Rect_Cyl.Cart_to_Cyl(Nodes_Coord)
    #Then iterate to find the nodes coordinates which equal F_coord.
    
    c=0
    c2=0
    for i in range(Nodes_Coord.shape[0]): 
        if np.abs(Nodes_Coord_Cyl[i][0] - F_coord[0]) < 1e-6 and \
            np.abs(Nodes_Coord_Cyl[i][1] - F_coord[1]) < 1e-6 and \
            np.abs(Nodes_Coord_Cyl[i][2] - F_coord[2]) < 1e-6:
                Rvector = Forces(c+1,F_comp,Rvector,Mag)
                F_nodes.append(c+1)
                F_components.append(F_comp)
                c2+=1
        c+=1
    if c2==0:
        print ('WARNING: SPECIFIED LOCATION DONT CORRESPOND TO ANY NODE')
    return (Rvector,F_nodes,F_components)