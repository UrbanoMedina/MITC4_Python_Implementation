import numpy as np
#Function to obtain from the element node array the corresponding dofs for
#each node in the global displacement matrix 
def Ksys_C_Array(Elem_Nodes,nnel,ndofs):
    Connec_Array=np.zeros(nnel*ndofs)
    k=0
    for i in Elem_Nodes:
        for j in range (ndofs):
            Connec_Array[k]=((i-1)*6+j)
            k=k+1
    return Connec_Array

def Ksys_Assem(Kelem,Connec_Array,sysdof):
    Ksys=np.zeros([sysdof,sysdof])
    for i in range(np.shape(Kelem)[0]):
        for j in range(np.shape(Kelem)[1]):
            Ksys[int(Connec_Array[i]),int(Connec_Array[j])]=Kelem[i,j]
    #Ksys=Ksys+Ksys
    return Ksys

def U_Elem(Element_Connect,U,nel,ndofs,nnel):
    U_Elem=np.zeros([nel,4*ndofs])
    for i in range(nel):
        Conectivity_Array= Ksys_C_Array(Element_Connect[i],nnel,ndofs)
        c=0
        for j in Conectivity_Array:
            U_Elem[i][c]=U[int(j)]
            c=c+1
    return(U_Elem)
    
def Solve_System_K_R(Ksys,R,bcdof):
    #In this function K matrix and R vector will be modified in order to solve
    #the linear system of equations. In Ksys, all the rows corresponding to a 
    #degree of freedom will be set to zero except the element corresponding to 
    #the matrix diagonal which will be set to one.
    #For the R vector, all the elements corresponding to a degree of freedom 
    #will be set to zero. 
    
    #Iterate for each degree of freedom restricted.
    for i, bc in enumerate(bcdof):
        #Iterate for each row in K
        for j in range (Ksys.shape[0]):
            Ksys[bc,j]=0
        Ksys[bc,bc]=1
        R[bc]=0
         
    return(Ksys,R)

            
            
