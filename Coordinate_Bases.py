import numpy as np
#BASES MODULES#
'''
Functions to calculate the following bases:
    Rotating Vectors Vn1 and Vn2
    gi covariant 
    gi contravariant
    ei local
'''

def gi_covariant(dX_dr,dX_ds,dX_dt):
    g1=[dX_dr[0], #e1
        dX_dr[1], #e2
        dX_dr[2]] #e3
    
    g2=[dX_ds[0], #e1
        dX_ds[1], #e2
        dX_ds[2]] #e3
    
    g3=[dX_dt[0], #e1
        dX_dt[1], #e2
        dX_dt[2]] #e3
    
    return [g1,g2,g3]

def gi_contra(g_cov):
    g_cov=np.array(g_cov)
    g1=g_cov[0]
    g2=g_cov[1]
    g3=g_cov[2]
    gij_cov=np.array([[np.dot(g1,g1),np.dot(g1,g2),np.dot(g1,g3)],
                      [np.dot(g2,g1),np.dot(g2,g2),np.dot(g2,g3)],
                      [np.dot(g3,g1),np.dot(g3,g2),np.dot(g3,g3)]])
    gij_contra=np.linalg.inv(gij_cov)
    g_contra=np.zeros([3,3]) 
    
    for i in range(len(g_cov)): 
        g_i=0
        for j in range(len(g_cov)):
            g_i= g_i + gij_contra[i,j]*g_cov[j]
        g_contra[i]=g_i
    
    return(g_contra)
        
def localei (r,s,t,gi):
    e3=gi[2]/np.linalg.norm(gi[2])
    e1=np.cross(gi[1],e3)/np.linalg.norm(np.cross(gi[1],e3))
    e2=np.cross(e3,e1)
    return [e1,e2,e3]

def vec_rot(Vn):
    e2=[0,1,0]
    if np.linalg.norm(np.cross(e2,Vn)) == 0:
        V1n= [1,0,0]
        V2n=np.cross(Vn,V1n)
    else:
        V1n=np.cross(e2,Vn)/np.linalg.norm(np.cross(e2,Vn))
        V2n=np.cross(Vn,V1n)
    return np.array([V1n,V2n])
    