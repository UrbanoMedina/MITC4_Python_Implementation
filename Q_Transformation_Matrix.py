import numpy as np

def Q_CartA_to_Q_CartB(exA,eyA,ezA,exB,eyB,ezB):
    #Function to transform the B cartesian basis to A cartesian basis. 
        #E_A = Q . E_B
        
    #The B basis will correspond to global cartesian basis.
        #exB=[1,0,0]
        #eyB=[0,1,0]
        #ezB=[0,0,1]
        
    #The A basis corresponds to local cartesian basis.
    
    l1=np.dot(exA,exB)
    l2=np.dot(eyA,exB)
    l3=np.dot(ezA,exB)
    
    m1=np.dot(exA,eyB)
    m2=np.dot(eyA,eyB)
    m3=np.dot(ezA,eyB)
    
    n1=np.dot(exA,ezB)
    n2=np.dot(eyA,ezB)
    n3=np.dot(ezA,ezB)
    
    Q=np.array([[l1**2,m1**2,n1**2,l1*m1,m1*n1,n1*l1],
                [l2**2,m2**2,n2**2,l2*m2,m2*n2,n2*l2],
                [l3**2,m3**2,n3**2,l3*m3,m3*n3,n3*l3],
                [2*l1*l2,2*m1*m2,2*n1*n2,l1*m2+m1*l2,m1*n2+n1*m2,n1*l2+l1*n2],
                [2*l2*l3,2*m2*m3,2*n2*n3,l2*m3+m2*l3,m2*n3+n2*m3,n2*l3+l2*n3],
                [2*l1*l3,2*m1*m3,2*n1*n3,m1*l3+l1*m3,n1*m3+m1*n3,l1*n3+n1*l3]])
    return Q 

def Q_Cart_to_Cconv(ex,ey,ez,g1_contr,g2_contr,g3_contr):
    #Function to transform from convective to global coordinates
        #E_global = Q . E_conv
        
    
    l1=np.dot(g1_contr,ex)
    l2=np.dot(g1_contr,ey)
    l3=np.dot(g1_contr,ez)
    
    m1=np.dot(g2_contr,ex)
    m2=np.dot(g2_contr,ey)
    m3=np.dot(g2_contr,ez)
    
    n1=np.dot(g3_contr,ex)
    n2=np.dot(g3_contr,ey)
    n3=np.dot(g3_contr,ez)
    
    Q=np.array([[l1**2,m1**2,n1**2,l1*m1,m1*n1,n1*l1],
                [l2**2,m2**2,n2**2,l2*m2,m2*n2,n2*l2],
                [l3**2,m3**2,n3**2,l3*m3,m3*n3,n3*l3],
                [2*l1*l2,2*m1*m2,2*n1*n2,l1*m2+m1*l2,m1*n2+n1*m2,n1*l2+l1*n2],
                [2*l2*l3,2*m2*m3,2*n2*n3,l2*m3+m2*l3,m2*n3+n2*m3,n2*l3+l2*n3],
                [2*l1*l3,2*m1*m3,2*n1*n3,m1*l3+l1*m3,n1*m3+m1*n3,l1*n3+n1*l3]])
    return Q 
