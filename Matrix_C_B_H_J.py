import numpy as np
from dX_and_dU import* 

def Inter_Conv_To_Cartesian(r,s,t,Nodal_Coord,Nodes_Vect,Nodes_Thick):
    h1 = (1/4)*(1+r)*(1+s)
    h2 = (1/4)*(1-r)*(1+s)
    h3 = (1/4)*(1-r)*(1-s)
    h4 = (1/4)*(1+r)*(1-s)
    h_fun=np.array([h1,h2,h3,h4])
    x= np.sum(h_fun*Nodal_Coord[:,0])+(t/2)*(np.sum(Nodes_Thick*Nodes_Vect[:,0]*h_fun))
    y= np.sum(h_fun*Nodal_Coord[:,1])+(t/2)*(np.sum(Nodes_Thick*Nodes_Vect[:,1]*h_fun))
    z= np.sum(h_fun*Nodal_Coord[:,2])+(t/2)*(np.sum(Nodes_Thick*Nodes_Vect[:,2]*h_fun))
    global_coord=np.array([x,y,z])
    return global_coord

def H_Matrix(r,s,t,Vec_Vn1_Vn2,Nodes_Thick):
    h1 = (1/4)*(1+r)*(1+s)
    h2 = (1/4)*(1-r)*(1+s)
    h3 = (1/4)*(1-r)*(1-s)
    h4 = (1/4)*(1+r)*(1-s)
    
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
    
    H_Matrix=np.array([
        [h1, 0, 0, (-t/2)*h1*a1*V2_n1x, (t/2)*h1*a1*V1_n1x, h2, 0, 0, (-t/2)*h2*a2*V2_n2x,
 (t/2)*h2*a2*V1_n2x, h3, 0, 0, (-t/2)*h3*a3*V2_n3x, (t/2)*h3*a3*V1_n3x, h4, 0, 0, 
 (-t/2)*h4*a4*V2_n4x, (t/2)*h4*a4*V1_n4x],
        
        [0, h1, 0, (-t/2)*h1*a1*V2_n1y, (t/2)*h1*a1*V1_n1y, 0, h2, 0, (-t/2)*h2*a2*V2_n2y,
 (t/2)*h2*a2*V1_n2y, 0, h3, 0, (-t/2)*h3*a3*V2_n3y, (t/2)*h3*a3*V1_n3y, 0, h4, 0,
 (-t/2)*h4*a4*V2_n4y, (t/2)*h4*a4*V1_n4y],
    
        [0, 0, h1, (-t/2)*h1*a1*V2_n1z, (t/2)*h1*a1*V1_n1z, 0, 0, h2, (-t/2)*h2*a2*V2_n2z,
 (t/2)*h2*a2*V1_n2z, 0, 0, h3, (-t/2)*h3*a3*V2_n3z, (t/2)*h3*a3*V1_n3z, 0, 0, h4,
 (-t/2)*h4*a4*V2_n4z, (t/2)*h4*a4*V1_n4z]])    
    
    return (H_Matrix)


def B_ID(dX_dr,dX_ds,dX_dt,dU_dr,dU_ds,dU_dt):
    dx_dr,dy_dr,dz_dr=dX_dr[0],dX_dr[1],dX_dr[2]
    dx_ds,dy_ds,dz_ds=dX_ds[0],dX_ds[1],dX_ds[2]
    dx_dt,dy_dt,dz_dt=dX_dt[0],dX_dt[1],dX_dt[2]
    
    du_dr,dv_dr,dw_dr=dU_dr[0],dU_dr[1],dU_dr[2]
    du_ds,dv_ds,dw_ds=dU_ds[0],dU_ds[1],dU_ds[2]
    du_dt,dv_dt,dw_dt=dU_dt[0],dU_dt[1],dU_dt[2]
    
    B_err_id = dx_dr*du_dr + dy_dr*dv_dr + dz_dr*dw_dr;
    B_ess_id = dx_ds*du_ds + dy_ds*dv_ds + dz_ds*dw_ds;
#Deformaciones Cortantes
    B_ers_id = 0.5*(dx_dr*du_ds + dy_dr*dv_ds + dz_dr*dw_ds + du_dr*dx_ds + dv_dr*dy_ds + dw_dr*dz_ds);
    # B_est_id = 0.5*(dx_ds*du_dt + dy_ds*dv_dt + dz_ds*dw_dt + du_ds*dx_dt + dv_ds*dy_dt + dw_ds*dz_dt);
    # B_ert_id = 0.5*(dx_dr*du_dt + dy_dr*dv_dt + dz_dr*dw_dt + du_dr*dx_dt + dv_dr*dy_dt + dw_dr*dz_dt);
    
    #B_ID=np.array([[B_err_id],[B_ess_id],[B_ers_id]])
    return(B_err_id,B_ess_id,B_ers_id)

def B_MX(r,s,Nodal_Coord,Nodes_Vect,Vec_Vn1_Vn2,Nodes_Thick):
    #Evaluation points: 
    A= [0,1,0]
    B= [-1,0,0]
    C= [0,-1,0]
    D= [1,0,0]
    
    #Evaluate for each point the partial derivatives dX and dU using the functions in the 
    #dX_and_dU module 
    #POINT A:
    dX_dr_A,dX_ds_A,dX_dt_A=dX(A[0],A[1],A[2],Nodal_Coord,Nodes_Vect,Nodes_Thick)
    dU_dr_A,dU_ds_A,dU_dt_A=dU(A[0],A[1],A[2],Vec_Vn1_Vn2,Nodes_Thick)
    #POINT B:
    dX_dr_B,dX_ds_B,dX_dt_B=dX(B[0],B[1],B[2],Nodal_Coord,Nodes_Vect,Nodes_Thick)
    dU_dr_B,dU_ds_B,dU_dt_B=dU(B[0],B[1],B[2],Vec_Vn1_Vn2,Nodes_Thick)
    #POINT C:
    dX_dr_C,dX_ds_C,dX_dt_C=dX(C[0],C[1],C[2],Nodal_Coord,Nodes_Vect,Nodes_Thick)
    dU_dr_C,dU_ds_C,dU_dt_C=dU(C[0],C[1],C[2],Vec_Vn1_Vn2,Nodes_Thick)
    #POINT D:
    dX_dr_D,dX_ds_D,dX_dt_D=dX(D[0],D[1],D[2],Nodal_Coord,Nodes_Vect,Nodes_Thick)
    dU_dr_D,dU_ds_D,dU_dt_D=dU(D[0],D[1],D[2],Vec_Vn1_Vn2,Nodes_Thick)
    
    #Components obtained by direct interpolation:
    ert_A=0.5*(dX_dr_A[0]*dU_dt_A[0] + dX_dr_A[1]*dU_dt_A[1] + dX_dr_A[2]*dU_dt_A[2] 
               + dU_dr_A[0]*dX_dt_A[0] + dU_dr_A[1]*dX_dt_A[1] + dU_dr_A[2]*dX_dt_A[2]) 
    
    ert_C=0.5*(dX_dr_C[0]*dU_dt_C[0] + dX_dr_C[1]*dU_dt_C[1] + dX_dr_C[2]*dU_dt_C[2] 
               + dU_dr_C[0]*dX_dt_C[0] + dU_dr_C[1]*dX_dt_C[1] + dU_dr_C[2]*dX_dt_C[2]) 
    
    est_D=0.5*(dX_ds_D[0]*dU_dt_D[0] + dX_ds_D[1]*dU_dt_D[1] + dX_ds_D[2]*dU_dt_D[2] 
               + dU_ds_D[0]*dX_dt_D[0] + dU_ds_D[1]*dX_dt_D[1] + dU_ds_D[2]*dX_dt_D[2])
    
    est_B=0.5*(dX_ds_B[0]*dU_dt_B[0] + dX_ds_B[1]*dU_dt_B[1] + dX_ds_B[2]*dU_dt_B[2] 
               + dU_ds_B[0]*dX_dt_B[0] + dU_ds_B[1]*dX_dt_B[1] + dU_ds_B[2]*dX_dt_B[2])
    
    # B_est_id = 0.5*(dx_ds*du_dt + dy_ds*dv_dt + dz_ds*dw_dt + du_ds*dx_dt + dv_ds*dy_dt + dw_ds*dz_dt);
    # B_ert_id = 0.5*(dx_dr*du_dt + dy_dr*dv_dt + dz_dr*dw_dt + du_dr*dx_dt + dv_dr*dy_dt + dw_dr*dz_dt);
    
    
   #Applying mix interpolation 
    Bert_mx= 0.5*(1+s)*ert_A + 0.5*(1-s)*ert_C
    Best_mx= 0.5*(1+r)*est_D + 0.5*(1-r)*est_B
    
    return (Bert_mx,Best_mx)

def C_local_Matrix(E,v):
    k=5/6
    C_local=(E/(1-v**2))*np.array([[1,v,0,0,0,0],
                [v,1,0,0,0,0],
                [0,0,0,0,0,0],
                [0,0,0,(1-v)/2,0,0],
                [0,0,0,0,k*(1-v)/2,0],
                [0,0,0,0,0,k*(1-v)/2]])
    return C_local

def Jacobian(dX_dr,dX_ds,dX_dt):
    J=np.array([[dX_dr[0],dX_dr[1],dX_dr[2]],
        [dX_ds[0],dX_ds[1],dX_ds[2]],
        [dX_dt[0],dX_dt[1],dX_dt[2]]])
    J_Det=np.linalg.det(J)
    return (J,J_Det)






