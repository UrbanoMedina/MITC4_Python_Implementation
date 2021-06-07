import numpy as np 
import math
#FUNCTIONS TO CREATE RECTANGULAR AND CYLINDRICAL MESH STRUCTURES.


#Function to create a rectangular 2d mesh. 
    #All z coordinates are set as 0.
    #All Vn vectors are of the form [0,0,1]
def rectangular_mesh(base,height,division1,division2,thickness):
    #division1: number of vertical divisions in the rectangle
    #division2: number of vertical divisions in the rectangle
    dbase= base/division2 #Length of an element base
    dheight= height/division1 #Length of an element height 
    #Array were nodal coordinates and Vn components will be stored.
    Nodes_Coord=np.zeros([(division1+1)*(division2+1),3])
    Nodes_Vn=np.zeros([(division1+1)*(division2+1),3])
    #Array were element connectivity will be stored.
    Element_Connect=np.zeros([division1*division2,4]) #4: number of nodes per element 
    #Array to store each node thickness
    Nodes_Thickness=np.zeros((division1+1)*(division2+1))
    
    #Loop to store nodal coordinates:
    cont=0
    for i in range(division1+1):
        y= height-dheight*i #y coordinates of nodes in row i, were row 0 is the top one.
        for j in range(division2+1):
            Nodes_Coord[cont]=[dbase*j,y,0]
            Nodes_Thickness[cont]=thickness
            cont+=1
            
            
    #Loop to store Vn nodes components:
    cont=0
    for i in range(division1+1):
        for j in range(division2+1):
            Nodes_Vn[cont]=[0,0,1]
            cont+=1
    #Loop to store element connectivity:
    c=0
    for i in range(division1):
        for j in range(division2):
            Element_Connect[c][1]=i*(division2+1)+1+j
            Element_Connect[c][0]=i*(division2+1)+2+j 
            Element_Connect[c][3]=(i+1)*(division2+1)+2+j   
            Element_Connect[c][2]=(i+1)*(division2+1)+1+j
            c+=1
    Element_Connect = Element_Connect.astype('int32')
    return(Nodes_Coord,Element_Connect,Nodes_Vn,Nodes_Thickness)


#Function to create a cylindrical mesh. 
    #R:radius of the cylinder
    #Theta:angle spanned by the cylinder(in degrees)
    #L:Height of the cylinder 
def cylindrical_mesh(R,theta,L,division1,division2,thickness):
    #R:radius of the cylinder 
    #theta:angles spanned by the cylinder. Array form:[theta2,theta1]
        #theta2:maximun angle spanned by the cylinder
        #theta1: starting angle of the cylinder.
        #Example: theta=[90,30] will give a cylinder which spans from 90 to 30 degrees. 
    #division1: number of angle divisions in the cylinder.
    #division2: number of height divisions in the cylinder
    dL= L/division2 #Length of an element´s base
    if theta[0] < theta[1]:
        print ('WARNING: FIRST ANGLE INPUT MUST BE LARGER THAN THE SECOND ONE')
        return ('Please enter again the angles in theta variable')
    
    dtheta= (theta[0]-theta[1])/division1 #Length of an elementa´s height 
    
    #Array were nodal coordinates and Vn components will be stored.
    Nodes_Coord=np.zeros([(division1+1)*(division2+1),3])
    Nodes_Vn=np.zeros([(division1+1)*(division2+1),3])
    
    #Array were element connectivity will be stored.
    Element_Connect=np.zeros([division1*division2,4]) #4: number of nodes per element 
    #Array to store each node thickness
    Nodes_Thickness=np.zeros((division1+1)*(division2+1))
    
    #Loop to store nodal coordinates in cylindrical coordinates
    cont=0
    for i in range(division1+1):
        the= theta[0]-dtheta*i #y coordinates of nodes in row i, were row 0 is the top one.
        for j in range(division2+1):
            Nodes_Coord[cont]=[dL*j,the,R] #The coordinates will be stored as [L,theta,radiud]
            Nodes_Thickness[cont]=thickness
            cont+=1
            

    #Obtain the Vn vector for each nodal coordinate
    Nodes_Vn=Vn_Cylindrical(Nodes_Coord)
    
    #Transform the nodal coordinates in cylindrical coordinates to cartesian coordinates
    #using the Cyl_to_Cart(L,theta,R) function
    Nodes_Coord=Cyl_to_Cart(Nodes_Coord)
    
    #Loop to store element connectivity:
    c=0
    for i in range(division1):
        for j in range(division2):
            Element_Connect[c][3]=i*(division2+1)+1+j   #3
            Element_Connect[c][2]=i*(division2+1)+2+j   #2
            Element_Connect[c][1]=(i+1)*(division2+1)+2+j #1  
            Element_Connect[c][0]=(i+1)*(division2+1)+1+j #0
            c+=1
        
    Element_Connect = Element_Connect.astype('int32')
    return(Nodes_Coord,Element_Connect,Nodes_Vn,Nodes_Thickness)


#Function to transform cylindrical coordinates to cartesian coordinates
def Cyl_to_Cart(Nodes_Coord_Cylindrical): #(L,theta,R)
    Nodes_Coord_Cart=np.zeros([Nodes_Coord_Cylindrical.shape[0],
                               Nodes_Coord_Cylindrical.shape[1]])
    for i in range(Nodes_Coord_Cylindrical.shape[0]):
        
        Nodes_Coord_Cart[i][1]=Nodes_Coord_Cylindrical[i][0]
        Nodes_Coord_Cart[i][0]=math.cos(math.radians(Nodes_Coord_Cylindrical[i][1]))\
                               *Nodes_Coord_Cylindrical[i][2]
        Nodes_Coord_Cart[i][2]=math.sin(math.radians(Nodes_Coord_Cylindrical[i][1]))\
                               *Nodes_Coord_Cylindrical[i][2]
    return Nodes_Coord_Cart


#Function to transform cartesian coordinates to cylindrical coordinates
def Cart_to_Cyl(Nodes_Coord_Cart): #(L,theta,R)
    Nodes_Coord_Cyl=np.zeros([Nodes_Coord_Cart.shape[0],
                               Nodes_Coord_Cart.shape[1]])
    for i in range(Nodes_Coord_Cart.shape[0]):
        
        Nodes_Coord_Cyl[i][0]=math.sqrt(Nodes_Coord_Cart[i][0]**2+Nodes_Coord_Cart[i][2]**2)
        if Nodes_Coord_Cart[i][0] < 1e-6:
            Nodes_Coord_Cyl[i][1]=90
        else:
            Nodes_Coord_Cyl[i][1]=math.degrees(np.arctan(Nodes_Coord_Cart[i][2]/Nodes_Coord_Cart[i][0]))
        
        Nodes_Coord_Cyl[i][2]=Nodes_Coord_Cart[i][1]

    return Nodes_Coord_Cyl

#Functions to obtain the Vn vector at a point given by L,theta,R coordinates. 
#This vector is equal to dP/dR which is the total derivative of a point with respect
#to the R variable. 
def Vn_Cylindrical(Nodes_Coord_Cylindrical): #(L,theta,R)
    Nodes_Vn=np.zeros([Nodes_Coord_Cylindrical.shape[0],
                               Nodes_Coord_Cylindrical.shape[1]])    

    for i in range(Nodes_Coord_Cylindrical.shape[0]):
        Nodes_Vn[i][1]=0
        Nodes_Vn[i][0]=math.cos(math.radians(Nodes_Coord_Cylindrical[i][1]))
        Nodes_Vn[i][2]=math.sin(math.radians(Nodes_Coord_Cylindrical[i][1]))
    return Nodes_Vn




