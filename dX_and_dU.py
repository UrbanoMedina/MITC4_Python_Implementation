import numpy as np
#This function is used to calculate the derivatives dX_dr dX_dr dX_dr and
#dU_dr dU_dr dU_dr
#In which dX_dr=[dx/dr dy/dr dz/dr] and dU_dr=[du/dr du/dr du/dr] 

def dX(r,s,t,Nodal_Coord,Nodes_Vect,Nodes_Thick): 
    Vn1x,Vn1y,Vn1z=Nodes_Vect[0][0], Nodes_Vect[0][1], Nodes_Vect[0][2]
    Vn2x,Vn2y,Vn2z=Nodes_Vect[1][0], Nodes_Vect[1][1], Nodes_Vect[1][2]
    Vn3x,Vn3y,Vn3z=Nodes_Vect[2][0],Nodes_Vect[2][1],Nodes_Vect[2][2]   
    Vn4x,Vn4y,Vn4z=Nodes_Vect[3][0],Nodes_Vect[3][1],Nodes_Vect[3][2]
    
    x1x,x1y,x1z= Nodal_Coord[0][0], Nodal_Coord[0][1], Nodal_Coord[0][2]
    x2x,x2y,x2z= Nodal_Coord[1][0], Nodal_Coord[1][1], Nodal_Coord[1][2] 
    x3x,x3y,x3z= Nodal_Coord[2][0], Nodal_Coord[2][1], Nodal_Coord[2][2]
    x4x,x4y,x4z= Nodal_Coord[3][0], Nodal_Coord[3][1], Nodal_Coord[3][2]
    
    a1,a2,a3,a4=Nodes_Thick[0],Nodes_Thick[1],Nodes_Thick[2],Nodes_Thick[3]
    
    dX_dr=[(x1x*(s + 1))/4 - (x2x*(s + 1))/4 + (x3x*(s - 1))/4 - (x4x*(s - 1))/4 + 
        (Vn1x*a1*t*(s + 1))/8 - (Vn2x*a2*t*(s + 1))/8 + (Vn3x*a3*t*(s - 1))/8 - 
        (Vn4x*a4*t*(s - 1))/8, #dx/dr
        (x1y*(s + 1))/4 - (x2y*(s + 1))/4 + (x3y*(s - 1))/4 - (x4y*(s - 1))/4 + 
        (Vn1y*a1*t*(s + 1))/8 - (Vn2y*a2*t*(s + 1))/8 + (Vn3y*a3*t*(s - 1))/8 - 
        (Vn4y*a4*t*(s - 1))/8, #dy/dr
        (x1z*(s + 1))/4 - (x2z*(s + 1))/4 + (x3z*(s - 1))/4 - (x4z*(s - 1))/4 + 
        (Vn1z*a1*t*(s + 1))/8 - (Vn2z*a2*t*(s + 1))/8 + (Vn3z*a3*t*(s - 1))/8 - 
        (Vn4z*a4*t*(s - 1))/8] #dz/dr
    
    dX_ds=[x1x*(r/4 + 1/4) - x2x*(r/4 - 1/4) + x3x*(r/4 - 1/4) - x4x*(r/4 + 1/4) + 
        (Vn1x*a1*t*(r/4 + 1/4))/2 - (Vn2x*a2*t*(r/4 - 1/4))/2 + (Vn3x*a3*t*(r/4 - 1/4))/2 - 
        (Vn4x*a4*t*(r/4 + 1/4))/2,#dx/ds
        x1y*(r/4 + 1/4) - x2y*(r/4 - 1/4) + x3y*(r/4 - 1/4) - x4y*(r/4 + 1/4) + 
        (Vn1y*a1*t*(r/4 + 1/4))/2 - (Vn2y*a2*t*(r/4 - 1/4))/2 + (Vn3y*a3*t*(r/4 - 1/4))/2 - 
        (Vn4y*a4*t*(r/4 + 1/4))/2,#dy/ds
        x1z*(r/4 + 1/4) - x2z*(r/4 - 1/4) + x3z*(r/4 - 1/4) - x4z*(r/4 + 1/4) + 
        (Vn1z*a1*t*(r/4 + 1/4))/2 - (Vn2z*a2*t*(r/4 - 1/4))/2 + (Vn3z*a3*t*(r/4 - 1/4))/2 - 
        (Vn4z*a4*t*(r/4 + 1/4))/2]#dz/ds
    
    dX_dt=[(Vn1x*a1*(r/4 + 1/4)*(s + 1))/2 - (Vn2x*a2*(r/4 - 1/4)*(s + 1))/2 + 
        (Vn3x*a3*(r/4 - 1/4)*(s - 1))/2 - (Vn4x*a4*(r/4 + 1/4)*(s - 1))/2,#du/dt
        (Vn1y*a1*(r/4 + 1/4)*(s + 1))/2 - (Vn2y*a2*(r/4 - 1/4)*(s + 1))/2 +
        (Vn3y*a3*(r/4 - 1/4)*(s - 1))/2 - (Vn4y*a4*(r/4 + 1/4)*(s - 1))/2,#dy/dt
        (Vn1z*a1*(r/4 + 1/4)*(s + 1))/2 - (Vn2z*a2*(r/4 - 1/4)*(s + 1))/2 + 
        (Vn3z*a3*(r/4 - 1/4)*(s - 1))/2 - (Vn4z*a4*(r/4 + 1/4)*(s - 1))/2] #dz/dt
   
    
    return dX_dr,dX_ds,dX_dt

def dU(r,s,t,Vec_Vn1_Vn2,Nodes_Thick):
#Asigning every component of the Vn1 vector (for eaach node) to a variable: 
    V1_n1x,V1_n1y,V1_n1z= Vec_Vn1_Vn2[0][0][0], Vec_Vn1_Vn2[0][0][1], Vec_Vn1_Vn2[0][0][2]
    V1_n2x,V1_n2y,V1_n2z= Vec_Vn1_Vn2[1][0][0], Vec_Vn1_Vn2[1][0][1], Vec_Vn1_Vn2[1][0][2]
    V1_n3x,V1_n3y,V1_n3z= Vec_Vn1_Vn2[2][0][0], Vec_Vn1_Vn2[2][0][1], Vec_Vn1_Vn2[2][0][2]
    V1_n4x,V1_n4y,V1_n4z= Vec_Vn1_Vn2[3][0][0], Vec_Vn1_Vn2[3][0][1], Vec_Vn1_Vn2[3][0][2]

#Asigning every component of the Vn2 vector (for each node) to a variable:
    V2_n1x,V2_n1y,V2_n1z= Vec_Vn1_Vn2[0][1][0], Vec_Vn1_Vn2[0][1][1], Vec_Vn1_Vn2[0][1][2]
    V2_n2x,V2_n2y,V2_n2z= Vec_Vn1_Vn2[1][1][0], Vec_Vn1_Vn2[1][1][1], Vec_Vn1_Vn2[1][1][2]
    V2_n3x,V2_n3y,V2_n3z= Vec_Vn1_Vn2[2][1][0], Vec_Vn1_Vn2[2][1][1], Vec_Vn1_Vn2[2][1][2]
    V2_n4x,V2_n4y,V2_n4z= Vec_Vn1_Vn2[3][1][0], Vec_Vn1_Vn2[3][1][1], Vec_Vn1_Vn2[3][1][2]
    
    a1,a2,a3,a4=Nodes_Thick[0],Nodes_Thick[1],Nodes_Thick[2],Nodes_Thick[3]
    
    dU_dr=np.zeros([3,20]) #Row corresponds to du_dr,dv_dr,dw_dr
    dU_ds=np.zeros([3,20]) #Row corresponds to du_ds,dv_ds,dw_ds
    dU_dt=np.zeros([3,20]) #Row corresponds to du_dt,dv_dt,dw_dt
    
    dU_dr[0,:]=[ s/4 + 1/4, 0, 0, -(V2_n1x*a1*t*(s + 1))/8, (V1_n1x*a1*t*(s + 1))/8, 
                - s/4 - 1/4, 0, 0, (V2_n2x*a2*t*(s + 1))/8, -(V1_n2x*a2*t*(s + 1))/8, 
                s/4 - 1/4, 0, 0, -(V2_n3x*a3*t*(s - 1))/8, (V1_n3x*a3*t*(s - 1))/8, 1/4 
                - s/4, 0, 0, (V2_n4x*a4*t*(s - 1))/8, -(V1_n4x*a4*t*(s - 1))/8]
    
    dU_dr[1,:]=[ 0, s/4 + 1/4, 0, -(V2_n1y*a1*t*(s + 1))/8, (V1_n1y*a1*t*(s + 1))/8, 
                0, - s/4 - 1/4, 0, (V2_n2y*a2*t*(s + 1))/8, -(V1_n2y*a2*t*(s + 1))/8,
                0, s/4 - 1/4, 0, -(V2_n3y*a3*t*(s - 1))/8, (V1_n3y*a3*t*(s - 1))/8,
                0, 1/4 - s/4, 0, (V2_n4y*a4*t*(s - 1))/8, -(V1_n4y*a4*t*(s - 1))/8]
    
    dU_dr[2,:]=[ 0, 0, s/4 + 1/4, -(V2_n1z*a1*t*(s + 1))/8, (V1_n1z*a1*t*(s + 1))/8,
                0, 0, - s/4 - 1/4, (V2_n2z*a2*t*(s + 1))/8, -(V1_n2z*a2*t*(s + 1))/8, 
                0, 0, s/4 - 1/4, -(V2_n3z*a3*t*(s - 1))/8, (V1_n3z*a3*t*(s - 1))/8, 
                0, 0, 1/4 - s/4, (V2_n4z*a4*t*(s - 1))/8, -(V1_n4z*a4*t*(s - 1))/8]

    
    dU_ds[0,:]=[ r/4 + 1/4, 0, 0, -(V2_n1x*a1*t*(r/4 + 1/4))/2, (V1_n1x*a1*t*(r/4 + 1/4))/2,
                1/4 - r/4, 0, 0, (V2_n2x*a2*t*(r/4 - 1/4))/2, -(V1_n2x*a2*t*(r/4 - 1/4))/2,
                r/4 - 1/4, 0, 0, -(V2_n3x*a3*t*(r/4 - 1/4))/2, (V1_n3x*a3*t*(r/4 - 1/4))/2,
                - r/4 - 1/4, 0, 0, (V2_n4x*a4*t*(r/4 + 1/4))/2, -(V1_n4x*a4*t*(r/4 + 1/4))/2]

    dU_ds[1,:]=[ 0, r/4 + 1/4, 0, -(V2_n1y*a1*t*(r/4 + 1/4))/2, (V1_n1y*a1*t*(r/4 + 1/4))/2, 
                0, 1/4 - r/4, 0, (V2_n2y*a2*t*(r/4 - 1/4))/2, -(V1_n2y*a2*t*(r/4 - 1/4))/2,
                0, r/4 - 1/4, 0, -(V2_n3y*a3*t*(r/4 - 1/4))/2, (V1_n3y*a3*t*(r/4 - 1/4))/2,
                0, - r/4 - 1/4, 0, (V2_n4y*a4*t*(r/4 + 1/4))/2, -(V1_n4y*a4*t*(r/4 + 1/4))/2]

    dU_ds[2,:]=[ 0, 0, r/4 + 1/4, -(V2_n1z*a1*t*(r/4 + 1/4))/2, (V1_n1z*a1*t*(r/4 + 1/4))/2,
                0, 0, 1/4 - r/4, (V2_n2z*a2*t*(r/4 - 1/4))/2, -(V1_n2z*a2*t*(r/4 - 1/4))/2,
                0, 0, r/4 - 1/4, -(V2_n3z*a3*t*(r/4 - 1/4))/2, (V1_n3z*a3*t*(r/4 - 1/4))/2,
                0, 0, - r/4 - 1/4, (V2_n4z*a4*t*(r/4 + 1/4))/2, -(V1_n4z*a4*t*(r/4 + 1/4))/2]

    
    dU_dt[0,:]=[ 0, 0, 0, -(V2_n1x*a1*(r/4 + 1/4)*(s + 1))/2, (V1_n1x*a1*(r/4 + 1/4)*(s + 1))/2,
                0, 0, 0, (V2_n2x*a2*(r/4 - 1/4)*(s + 1))/2, -(V1_n2x*a2*(r/4 - 1/4)*(s + 1))/2,
                0, 0, 0, -(V2_n3x*a3*(r/4 - 1/4)*(s - 1))/2, (V1_n3x*a3*(r/4 - 1/4)*(s - 1))/2, 
                0, 0, 0, (V2_n4x*a4*(r/4 + 1/4)*(s - 1))/2, -(V1_n4x*a4*(r/4 + 1/4)*(s - 1))/2]

    dU_dt[1,:]=[ 0, 0, 0, -(V2_n1y*a1*(r/4 + 1/4)*(s + 1))/2, (V1_n1y*a1*(r/4 + 1/4)*(s + 1))/2,
                0, 0, 0, (V2_n2y*a2*(r/4 - 1/4)*(s + 1))/2, -(V1_n2y*a2*(r/4 - 1/4)*(s + 1))/2,
                0, 0, 0, -(V2_n3y*a3*(r/4 - 1/4)*(s - 1))/2, (V1_n3y*a3*(r/4 - 1/4)*(s - 1))/2,
                0, 0, 0, (V2_n4y*a4*(r/4 + 1/4)*(s - 1))/2, -(V1_n4y*a4*(r/4 + 1/4)*(s - 1))/2]

    dU_dt[2,:]=[ 0, 0, 0, -(V2_n1z*a1*(r/4 + 1/4)*(s + 1))/2, (V1_n1z*a1*(r/4 + 1/4)*(s + 1))/2,
                0, 0, 0, (V2_n2z*a2*(r/4 - 1/4)*(s + 1))/2, -(V1_n2z*a2*(r/4 - 1/4)*(s + 1))/2,
                0, 0, 0, -(V2_n3z*a3*(r/4 - 1/4)*(s - 1))/2, (V1_n3z*a3*(r/4 - 1/4)*(s - 1))/2,
                0, 0, 0, (V2_n4z*a4*(r/4 + 1/4)*(s - 1))/2, -(V1_n4z*a4*(r/4 + 1/4)*(s - 1))/2]
    return dU_dr,dU_ds,dU_dt #Output each one as a matrix. 

