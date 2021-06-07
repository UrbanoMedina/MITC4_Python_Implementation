import numpy as np
import math 
import dX_and_dU 
import Coordinate_Bases 
import Matrix_C_B_H_J
import MITC4_PRE_PROCESSING_Patch_Test_Bending_txx
import Coordinate_Bases 
import Q_Transformation_Matrix
import K_Matrix_Functions
import Plotting
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
############## ------------------------------------- ########################
############## MITC4 PROCESSING: DISP./STRAIN/STRESS #########################
############## ------------------------------------- #########################

'''
This script will have as an input the nodal point and element information from 
the MITC4_PRE_PROCESSING.py script. The information is extracted with the following
line of code:
    
    Material,Nodes_Coord,Nodes_Thickness
    ,Nodes_Vn,Element_Connect,bcdof,R=MITC4_Pre_Processing()
    
    Each variable is described in the MITC4_PRE_PROCESSING.py script.
    
After reading the variables, the processing will start.
    The steps are the following:
        1. Obtain each element stiffness matrix by numerical integration. 
        2. Assemble the system stiffness matrix. 
        3. Remove from the stiffness matrix the rows and columns corresponding 
        to the known displacements (given by the boundary conditions).
        4. Modify the force vector to only leave the BCs forces and zero forces. 
            NOTE: Forces corresponding to zero displacements are non-zero. 
        5. Solve the system for displacements: U=K/R
        6. Obtain strains: e=BU
        7.Obtain stress: tau=Ce
        8. Output displacement/stress/strain to the POST_PROCESSING script. 

Numerical Integration Steps: 
    1. Initiate a loop to integrate each element and obtain for each one the 
    stiffness matrix. 
    2. After initiating the loop, obtain for each element the following information:
        -From the Element_Connect array obtain the nodes which compose the element. 
            Note: the order of the nodes will be the same as in the array. 
        -Given the 4 nodes:
            -Create a Nodal coordinates array: each row corresponding to a node
            and each column to the x y z components.
            -Obtain the Vn vector for each node. Each row corresponding to a node and 
            columns to the components. 
            -Obtain the V1 and V2 vectors for each node. This information will be stored
            in a multidimensional array  of shape [4,2,3], interpreted as:
            4 matrices / Each Matrix will have 2 rows / Each row will have 3 columns. 
            Each matrix will correspond to a node.
            Row 1 will correspond to the V1k vector.
            Row 2 will correspond to the V2k vector.
            Each column will correspond to the x y z components of the vectors. 
    3. Obtain the Rotational_Matrix for each node, which is used to express the vectors 
    Vn, V1 and V2  in cartesian coordinates. 
        -First, obtain the rotational matrix for each node. This rotational matrix 
        will be 6x6. This will be done in a loop.
        -Second, after each iteration, assemble the rotational matrix in a
        multiidimensional array of shap[4,6,6]. There will be the 4 matrices and 
        each matrix will correspond to a node. 
        Finally, assemble the whole element rotation matrix, which will be 24x24
    4. Start the numerical integration. 
        - Start 3 loops, for each r s t variable, evaluating Gauss Points. 
        - 8 sets of r s t points will be evaluated. 
        - For each set of points evaluate:
            -Coordinate basis: gi_covariant,gi_contravariant,eilocal
            -Bconv(6x20), then transform it to Bglobal(6x20)
            -Then transform Bglobal(6x20) to Bglobal(6x24). When doing this transformation,
            the theta3 degree of freedom is included. 
            -Transform Clocal to Cglobal.
            -Obtain dvol=determinan(J)
            -Obtain the element stifness matrix by:
                Kelem = Kelem + B_GLOBAL6x24' * C_GLOBAL * B_GLOBAL6x24 * dVOL
            -Once the 8 point cycle finishes, the Kelem will be full evaluated. 
            -Assemble Kelem in the Ksyst.
    
'''
def MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
Nodes_Vn,Element_Connect,bcdof,R): 

    # Material,Nodes_Coord,Nodes_Thickness,\
    # Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_Bending_txx.Node_Element_Info()
    
    #General Information of the given geometry: 
    nel= np.shape(Element_Connect)[0] #Number of elements
    nnel=4                            #Nodes por element.
    ndofs=6                           #DOFs per node.
    nnode= np.shape(Nodes_Coord)[0]   #Number of nodes in the system
    sysdof= nnode*ndofs               #Number of degrees of freedom in the system
    elemdof=nnel*ndofs                #DOFS per element
    num_gp=8                          #To numerically integrate, 8 sets of GP will be evaluated
    gauss_points=[-1/math.sqrt(3),1/math.sqrt(3)]
    U_unknown=np.zeros(sysdof)     #All nodal displacements components will be stored in this array. 
    B_Global_GP=np.zeros([nel,num_gp,ndofs,nnel*ndofs]) #For each element, 8 B global matrices will be stored
                                                        #corresponding to each set of gauss points. 
    C_Global_GP=np.zeros([nel,num_gp,ndofs,ndofs]) #For each element, 8 B global matrices will be stored
                                                        #corresponding to each set of gauss points. 
                                                        #This n*8 matrices will be used to calculate stress at GP. 
                                                        #This n*8 matrices will be used to calculate strains at GP.
    GP_Global_Coord=np.zeros([nel,num_gp,3]) #This array will be used to store for each element the gauss points´
                                             #global coordinates, which will be needed to graph. 
    Ksys_Final=np.zeros([sysdof,sysdof]) #System Stifness Matrix 
    for q in range(nel):
        #FOR THE ELEMENT i:
        Elem_Nodes=Element_Connect[q] #These nodes compose the element. 
        Elem_Nodes_Coord=np.zeros([nnel,3]) #Coordinates of the elemental nodes.
        Elem_Nodes_Vn=np.zeros([nnel,3])    #Components of each Vn node vector. 
        Elem_Nodes_V1_V2=np.zeros([nnel,2,3]) #Components of each V1 and V2 node vector. 
        Elem_Nodes_Thickness=np.zeros([nnel])
        Kelem=np.zeros([elemdof,elemdof])
        
        c1=0
        #Fill Nodes and Vn coordinates/components 
        for j in Elem_Nodes:
            Elem_Nodes_Coord[c1]=Nodes_Coord[j-1] #Reason [j-1]: nodes are numbered starting from 1 
            Elem_Nodes_Vn[c1]=Nodes_Vn[j-1]
            Elem_Nodes_Thickness[c1]=Nodes_Thickness[j-1]
            c1=c1+1
    
        #Fill V1 and V2 components 
        c2=0
        for i in Elem_Nodes_Vn:
            Elem_Nodes_V1_V2[c2,0]=Coordinate_Bases.vec_rot(i)[0] #Vn1
            Elem_Nodes_V1_V2[c2,1]=Coordinate_Bases.vec_rot(i)[1] #Vn2 
            c2=c2+1    
    
        #----------------------------------------------#
        #-OBTAIN THE ROTATIONAL MATRIX FOR THE ELEMENT-#
        #----------------------------------------------#    
        Trot_k=np.zeros([4,6,6]) #Matrix in which each node rotational matrix will be stored.
        T_k=np.zeros([4,6,6]) #Matrix to store Q for each node 
        Trot_Elem=np.zeros([24,24])#Matrix in which the element rotational matrix will be stored 
        for i in range(nnel):
            T_k[i,:,:]=Q_Transformation_Matrix.Q_CartA_to_Q_CartB(
                [Elem_Nodes_V1_V2[i][0][0],Elem_Nodes_V1_V2[i][0][1],Elem_Nodes_V1_V2[i][0][2]],
                [Elem_Nodes_V1_V2[i][1][0],Elem_Nodes_V1_V2[i][1][1],Elem_Nodes_V1_V2[i][1][2]],
                [Elem_Nodes_Vn[i][0],Elem_Nodes_Vn[i][1],Elem_Nodes_Vn[i][2]],
                [1,0,0],[0,1,0],[0,0,1])
            Trot_k[i,0,0]=1
            Trot_k[i,1,1]=1
            Trot_k[i,2,2]=1
            Trot_k[i,3,3]=math.sqrt(T_k[i,0,0])
            Trot_k[i,3,4]=math.sqrt(T_k[i,1,0])
            Trot_k[i,3,5]=math.sqrt(T_k[i,2,0])
            Trot_k[i,4,3]=math.sqrt(T_k[i,0,1])
            Trot_k[i,4,4]=math.sqrt(T_k[i,1,1])
            Trot_k[i,4,5]=math.sqrt(T_k[i,2,1])
            Trot_k[i,5,3]=math.sqrt(T_k[i,0,2])
            Trot_k[i,5,4]=math.sqrt(T_k[i,1,2])
            Trot_k[i,5,5]=math.sqrt(T_k[i,2,2])
            
        Trot_Elem[0:6,0:6]=Trot_k[0,:,:]
        Trot_Elem[6:12,6:12]=Trot_k[1,:,:]
        Trot_Elem[12:18,12:18]=Trot_k[2,:,:]
        Trot_Elem[18:24,18:24]=Trot_k[3,:,:]
    
                #--------------------------------------------#
                #-#INITIATE CYCLES OF NUMERICAL INTEGRATION--#
                #--------------------------------------------#
        Ksys_El=np.zeros([sysdof,sysdof]) #Used to express the element K in th system K. 
        gp_counter=0
        for t in gauss_points:
            for s in gauss_points:
                for r in gauss_points:
                #------------------------------------------------------------#    
                #-EVALUATE COORDINATE BASIS AT THE (r,s,t) POINT IN QUESTION-#
                #------------------------------------------------------------# 
                    dX_dr,dX_ds,dX_dt=dX_and_dU.dX(r,s,t,Elem_Nodes_Coord,
                                                    Elem_Nodes_Vn,Elem_Nodes_Thickness) 
                    gi_cov=Coordinate_Bases.gi_covariant(dX_dr,dX_ds,dX_dt)
                    gi_contr=Coordinate_Bases.gi_contra(gi_cov)
                    ei_local=Coordinate_Bases.localei(r,s,t,gi_cov)
                    #print(Matrix_C_B_H_J.Inter_Conv_To_Cartesian(r,s,t,Elem_Nodes_Coord,
                                                                 #Elem_Nodes_Vn,Elem_Nodes_Thickness))
    
                #------------------------------------------------------------#    
                #----STORE THE GLOBAL CART. COORDINATES OF THE R,S,T POINT---#
                #------------------------------------------------------------# 
                    GP_Global_Coord[q,gp_counter]=Matrix_C_B_H_J.Inter_Conv_To_Cartesian(r,s,t,Elem_Nodes_Coord,Elem_Nodes_Vn,
                                                           Elem_Nodes_Thickness)
                                        
                #-----------------------#    
                #-ASSEMBLE THE B MATRIX-#
                #-----------------------#
                    B_conv=np.zeros([6,20])
                #Obtain the dU_dr,dU_ds,dU_dt derivatives:
                    dU_dr,dU_ds,dU_dt=dX_and_dU.dU(r,s,t,Elem_Nodes_V1_V2,
                                                   Elem_Nodes_Thickness)
                #Obtain the Direct Interpolation components:
                    B_err_id,B_ess_id,B_ers_id=Matrix_C_B_H_J.B_ID(dX_dr,
                                                                 dX_ds,dX_dt,dU_dr,
                                                                 dU_ds,dU_dt)
                #Obtain the Mix Interpolation Components:
                    B_ert_mx,B_est_mx=Matrix_C_B_H_J.B_MX(r,s,Elem_Nodes_Coord,Elem_Nodes_Vn,
                                                      Elem_Nodes_V1_V2,Elem_Nodes_Thickness)
                #Fill the B_conv matrix:
                    B_conv[0,:]=B_err_id
                    B_conv[1,:]=B_ess_id
                    B_conv[3,:]=2*B_ers_id
                    B_conv[4,:]=2*B_est_mx
                    B_conv[5,:]=2*B_ert_mx
                #Transform Bconv to Bglobal 
                    Q_conv_to_global=Q_Transformation_Matrix.Q_Cart_to_Cconv([1,0,0],
                                                                             [0,1,0],[0,0,1],
                                                                             gi_contr[0],
                                                                             gi_contr[1],
                                                                             gi_contr[2])
                    B_Global_6x20=np.matmul(Q_conv_to_global,B_conv)
                 #Transform Bglobal6x20 to Bglobal6x24
                    B_Global_6x24=np.zeros([6,24]) #Here the third degree of freedom is being added
                    #The columns corresponding to theta3 degree of freedom are left 0
                    B_Global_6x24[:,0:5]=B_Global_6x20[:,0:5]
                    B_Global_6x24[:,6:11]=B_Global_6x20[:,5:10]
                    B_Global_6x24[:,12:17]=B_Global_6x20[:,10:15]
                    B_Global_6x24[:,18:23]=B_Global_6x20[:,15:20]
                    
                    #Store the B matrix for the (r,s,t) point evaluated:
                    B_Global_GP[q,gp_counter]=B_Global_6x24
                #-----------------------#    
                #-ASSEMBLE THE C MATRIX-#
                #-----------------------#
                
                    #Obtain the local C matrix:
                    C_Matrix=Matrix_C_B_H_J.C_local_Matrix(Material[0],Material[1])
                    #Obtain the transformational matrix
                    Q_C_Matrix=Q_Transformation_Matrix.Q_CartA_to_Q_CartB(ei_local[0],ei_local[1],
                                                                          ei_local[2],[1,0,0],[0,1,0],
                                                                          [0,0,1])
                    #Obtain the C_Global 
                    C_Global=np.matmul(np.matmul(np.transpose(Q_C_Matrix),C_Matrix),Q_C_Matrix)
                    #Store the C matrix for the (r,s,t) point evaluated:
                    C_Global_GP[q,gp_counter]=C_Global
                #-----------------------------#    
                #-OBTAIN JACOBIAN DETERMINANT-#
                #-----------------------------#
                    J,J_Det=Matrix_C_B_H_J.Jacobian(dX_dr,dX_ds,dX_dt)
                
                #------------------------------#    
                #-SUM FOR THE K ELEMENT MATRIX-#
                #------------------------------#
                    Kelem=Kelem+np.matmul(np.matmul(np.transpose(B_Global_6x24),C_Global),B_Global_6x24)*J_Det
                    gp_counter=gp_counter+1
        #---------------------------------------------------#    
        #-ASSEMBLE Kelem in Ksys AFTER NUM. IN. IS FINISHED-#
        #---------------------------------------------------#
        #First rotate Kelem 
        #print (np.diagonal(Kelem).min())
        Kelem=np.matmul(np.matmul(np.transpose(Trot_Elem),Kelem),Trot_Elem)
        #Obtain the connectivity array for the element.
        Connec_Array=K_Matrix_Functions.Ksys_C_Array(Elem_Nodes,nnel,ndofs)
        #Assemble the local Kelem [24,24] in the global Ksys[48,48]
        Ksys_El=K_Matrix_Functions.Ksys_Assem(Kelem,Connec_Array,sysdof)
        #The system stiffness matrix will be the sum of Ksys_El
        Ksys_Final=Ksys_Final+Ksys_El
    
    #---------------------------------------------------#    
    #--------APPLICATION OF BOUNDARY CONDITIONS---------#
    #---------------------------------------------------#
    #The goal is to solve the linear system KU=R for U. This way U will be equal to
    #U=inv(K)*R. But, K and R must be modified to only account for the unkown displacments, 
    #meaning that the columns and rows from Ksys which correspond to BCs must be removed. The 
    #same with the force vector but only the rows clearly.
    #Example: If a 24x24 matrix is considered and the elements of the nodal displacement array
    #[6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 5, 11, 17, 23] are equal to zero, then in Ksys all of this
    #columns and rows must be removed, leaving Ksys with a shape of [10x10].
    #For the force vector, this rows must be removed. 
    
    #Eliminate for the in the force vectors elements corresponding to the dofs. 
    #R_red=np.zeros([sysdof-len(bcdof)])
    # Rred=np.delete(R,bcdof,axis=0) 
    
    # #Eliminate the rows and columns for the system stiffness matrix 
    # Kred=np.delete(Ksys_Final,bcdof,axis=0) 
    # Kred=np.delete(Kred,bcdof,axis=1)
    #---------------------------------------------------#    
    #------------SOLVER THE LINEAR SYSTEM---------------#
    #---------------------------------------------------#
    
    #Obtain the modified stiffness matrix and load vector:
    Kmod,Rmod=K_Matrix_Functions.Solve_System_K_R(Ksys_Final,R,bcdof)
    
    U_unknown=np.matmul(np.linalg.inv(Kmod),Rmod)
    # EVALUATE UNKNOWN DISPLACEMENTS:
    # U_unknown=np.matmul(np.linalg.inv(Kred),Rred)
    
    #ASSEMBLE SYSTEM NODAL DISPLACEMENT VECTOR in U_unknown:
    #Find the indices for which U_unkown account to
    # u_indices=np.linspace(0,sysdof-1,sysdof) #All of the dofs indices 
    # u_indices=np.delete(u_indices,bcdof,axis=0) 
    
    # c4=0 
    # for i in u_indices:
    #     U_unknown[int(i)]=U_unknown[c4] #Store all the nodal displacement components
    #     c4=c4+1
    
    #Store each node displacement component
    U_Nodal=np.zeros([nnode,ndofs])
    # U_Nodal2=np.zeros([nnode,ndofs])
    for i in range(nnode):
        for j in range(ndofs):
            U_Nodal[i][j]= U_unknown[(i)*6+j]
    
        
    #FIND THE REACTIONS OF THE SYSTEM 
    System_Nodal_Reactions=np.matmul(Ksys_Final,U_unknown)
    
    #--------------------------------------------#
    #--SUM THE DISPLACEMENTS TO THE NODAL COORD--#
    #--------------------------------------------#
    
    #Create an array were the displaced nodal coord will be stored:
    Nodes_Coord_Disp=np.zeros([nnode,3])
    Magni=1e6 #Displacement magnification of the displacements to be able to
    #visualize them better.
    
    for i in range(nnode):
        disp_comp=np.zeros(3)
        for j in range(3): #Cycle to obtain for the node the 3 displacement components
            disp_comp[j]=U_unknown[i*6+j]*Magni
        Nodes_Coord_Disp[i]=Nodes_Coord[i]+disp_comp        
        
    # #--------------------------------------------#
    # #-----PRINT ORIGINAL VS DISPLACED FIGURE-----#
    # #--------------------------------------------#
    
    # #Create the figure to plot the geometry:
    # fig = plt.figure(1)
    #  #Create a 3d cart. axis
    # ax = fig.add_subplot(111, projection='3d')
    # ax.set_ylim(np.min(Nodes_Coord_Disp[:,1]),np.max(Nodes_Coord_Disp[:,1]))
    # ax.set_xlim(np.min(Nodes_Coord_Disp[:,0]),np.max(Nodes_Coord_Disp[:,0]))
    # ax.set_zlim(np.min(Nodes_Coord_Disp[:,2]),np.max(Nodes_Coord_Disp[:,2]))
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # #Plot each element
    # Plotting.print_Geometry(Element_Connect,Nodes_Coord,ax,'b')
    # Plotting.print_Geometry(Element_Connect,Nodes_Coord_Disp,ax,'r')
    
    
    
    
    #---------------------------------------------------#    
    #-----GAUSS POINTS STRESS/STRAIN CALCULATION--------#
    #---------------------------------------------------#
    
    #First organize the nodal displacements per element:
    U_Elem=K_Matrix_Functions.U_Elem(Element_Connect,U_unknown, nel, ndofs,nnel)
    
    
    #---------------------------------------------------#    
    #----------OBTAIN STRESS/STRAIN FOR NODES-----------#
    #---------------------------------------------------#
    
    Nodes_Natural_Coordinates=[[1,1,0],
                               [-1,1,0],
                               [-1,-1,0],
                               [1,-1,0]]
    
    #Evaluate B and C matrix for each node.
        #Create zero arrays to store them
    B_Nodal_Global=np.zeros([nel,nnel,ndofs,nnel*ndofs])
    C_Nodal_Global=np.zeros([nel,nnel,ndofs,ndofs])
    
    #Iterate each element for each of the 4 nodal natural coordinates.
    for q in range (nel):
        #Obtain the element information
        Elem_Nodes=Element_Connect[q] #These nodes compose the element. 
        Elem_Nodes_Coord=np.zeros([nnel,3]) #Coordinates of the elemental nodes.
        Elem_Nodes_Vn=np.zeros([nnel,3])    #Components of each Vn node vector. 
        Elem_Nodes_V1_V2=np.zeros([nnel,2,3]) #Components of each V1 and V2 node vector. 
        Elem_Nodes_Thickness=np.zeros([nnel])
        
        c1=0
        #Fill Nodes and Vn coordinates/components 
        for j in Elem_Nodes:
            Elem_Nodes_Coord[c1]=Nodes_Coord[j-1] #Reason [j-1]: nodes are numbered starting from 1 
            Elem_Nodes_Vn[c1]=Nodes_Vn[j-1]
            Elem_Nodes_Thickness[c1]=Nodes_Thickness[j-1]
            c1=c1+1
    
        #Fill V1 and V2 components 
        c2=0
        for i in Elem_Nodes_Vn:
            Elem_Nodes_V1_V2[c2,0]=Coordinate_Bases.vec_rot(i)[0] #Vn1
            Elem_Nodes_V1_V2[c2,1]=Coordinate_Bases.vec_rot(i)[1] #Vn2 
            c2=c2+1    
        gp_counter=0
        for j in Nodes_Natural_Coordinates:
            r,s,t=j[0],j[1],j[2]
            #Obtain the coordinate basis information for the node. 
            dX_dr,dX_ds,dX_dt=dX_and_dU.dX(r,s,t,Elem_Nodes_Coord,
                                                    Elem_Nodes_Vn,Elem_Nodes_Thickness) 
            gi_cov=Coordinate_Bases.gi_covariant(dX_dr,dX_ds,dX_dt)
            gi_contr=Coordinate_Bases.gi_contra(gi_cov)
            ei_local=Coordinate_Bases.localei(r,s,t,gi_cov)
            #-----------------------#    
            #-ASSEMBLE THE B MATRIX-#
            #-----------------------#
            B_conv=np.zeros([6,20])
            #Obtain the dU_dr,dU_ds,dU_dt derivatives:
            dU_dr,dU_ds,dU_dt=dX_and_dU.dU(r,s,t,Elem_Nodes_V1_V2,Elem_Nodes_Thickness)
            #Obtain the Direct Interpolation components:
            B_err_id,B_ess_id,B_ers_id=Matrix_C_B_H_J.B_ID(dX_dr,dX_ds,dX_dt,dU_dr,
                                                           dU_ds,dU_dt)
            #Obtain the Mix Interpolation Components:
            B_ert_mx,B_est_mx=Matrix_C_B_H_J.B_MX(r,s,Elem_Nodes_Coord,Elem_Nodes_Vn,
                                                      Elem_Nodes_V1_V2,Elem_Nodes_Thickness)
            #Fill the B_conv matrix:
            B_conv[0,:]=B_err_id
            B_conv[1,:]=B_ess_id
            B_conv[3,:]=2*B_ers_id
            B_conv[4,:]=2*B_est_mx
            B_conv[5,:]=2*B_ert_mx
            #Transform Bconv to Bglobal 
            Q_conv_to_global=Q_Transformation_Matrix.Q_Cart_to_Cconv([1,0,0],
                                                                     [0,1,0],[0,0,1],
                                                                     gi_contr[0],
                                                                     gi_contr[1],
                                                                     gi_contr[2])
            B_Global_6x20=np.matmul(Q_conv_to_global,B_conv)
            #Transform Bglobal6x20 to Bglobal6x24
            B_Global_6x24=np.zeros([6,24]) #Here the third degree of freedom is being added
            #The columns corresponding to theta3 degree of freedom are left 0
            B_Global_6x24[:,0:5]=B_Global_6x20[:,0:5]
            B_Global_6x24[:,6:11]=B_Global_6x20[:,5:10]
            B_Global_6x24[:,12:17]=B_Global_6x20[:,10:15]
            B_Global_6x24[:,18:23]=B_Global_6x20[:,15:20]
                    
            #Store the B matrix for the (r,s,t) point evaluated:
            B_Nodal_Global[q,gp_counter]=B_Global_6x24
            #-----------------------#    
            #-ASSEMBLE THE C MATRIX-#
            #-----------------------#
                
            #Obtain the local C matrix:
            C_Matrix=Matrix_C_B_H_J.C_local_Matrix(Material[0],Material[1])
            #Obtain the transformational matrix
            Q_C_Matrix=Q_Transformation_Matrix.Q_CartA_to_Q_CartB(ei_local[0],ei_local[1],
                                                                          ei_local[2],[1,0,0],[0,1,0],
                                                                          [0,0,1])
            #Obtain the C_Global 
            C_Global=np.matmul(np.matmul(np.transpose(Q_C_Matrix),C_Matrix),Q_C_Matrix)
            #Store the C matrix for the (r,s,t) point evaluated:
            C_Nodal_Global[q,gp_counter]=C_Global        
            gp_counter=gp_counter+1
    
    #Generate the arrays to store the strain/stress for each element,for each node. 
    Nodes_Elemental_Strains=np.zeros([nel,nnel,6])
    Nodes_Elemental_Stress=np.zeros([nel,nnel,6])    
    Nodes_Elemental_VM=np.zeros([nel,nnel]) #Store the Von Misses stress
    
    # #Iterate to obtain Nodal stress/strain for each point:
    for i in range(nel):
        for j in range(nnel):
            Nodes_Elemental_Strains[i,j]=np.matmul(B_Nodal_Global[i][j],U_Elem[i])
            Nodes_Elemental_Stress[i,j]=np.matmul(C_Nodal_Global[i][j],Nodes_Elemental_Strains[i][j])
            Nodes_Elemental_VM[i,j]=(math.sqrt(0.5*((Nodes_Elemental_Stress[i][j][0]-
                                                      Nodes_Elemental_Stress[i][j][1])**2
                                                      +(Nodes_Elemental_Stress[i][j][1]-
                                                      Nodes_Elemental_Stress[i][j][2])**2
                                                      +(Nodes_Elemental_Stress[i][j][2]-
                                                      Nodes_Elemental_Stress[i][j][0])**2
                                                      +6*(Nodes_Elemental_Stress[i][j][4]**2
                                                          +Nodes_Elemental_Stress[i][j][5]**2
                                                          +Nodes_Elemental_Stress[i][j][3]**2))))
    
    #        
    #Create an array in which each row corresponds a node and 4 columns. 
    #First 3 columns will correspond to the displaced nodal coordinates
    #The last column will correspond to the VM stress at that point
    Nodal_Disp_VM=np.zeros([Nodes_Coord_Disp.shape[0],4])
    #Assign the displaced nodal coordinats to the first 3 columns
    Nodal_Disp_VM[:,0:3]=Nodes_Coord_Disp
    #Assign to the last column VM stress
    for i in range(nel):
        Conec=Element_Connect[i] #Save the nodes which correspond the element i
        for j in range (nnel):
            Nodal_Disp_VM[Conec[j]-1][3]=Nodes_Elemental_VM[i][j]
    
    
    
    
    #
    #--------------------------------------------#
    #----OBTAIN STRESS/STRAIN FOR GAUSS POINTS----
    #--------------------------------------------#
    
    #Generate the arrays to store the strain/stress for each element,for each node. 
    GP_Elemental_Strains=np.zeros([nel,num_gp,6])
    GP_Elemental_Stress=np.zeros([nel,num_gp,6])    
    GP_Elemental_VM=np.zeros([nel,num_gp]) #Store the Von Misses stress
    
    # #Iterate to obtain GP stress/strain for each point:
    for i in range(nel):
        for j in range(num_gp):
            GP_Elemental_Strains[i,j]=np.matmul(B_Global_GP[i][j],U_Elem[i])
            GP_Elemental_Stress[i,j]=np.matmul(C_Global_GP[i][j],GP_Elemental_Strains[i][j])
            GP_Elemental_VM[i,j]=(math.sqrt(0.5*((GP_Elemental_Stress[i][j][0]-
                                                      GP_Elemental_Stress[i][j][1])**2
                                                      +(GP_Elemental_Stress[i][j][1]-
                                                      GP_Elemental_Stress[i][j][2])**2
                                                      +(GP_Elemental_Stress[i][j][2]-
                                                      GP_Elemental_Stress[i][j][0])**2
                                                      +6*(GP_Elemental_Stress[i][j][4]**2
                                                          +GP_Elemental_Stress[i][j][5]**2
                                                          +GP_Elemental_Stress[i][j][3]**2))))

    return(Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,
           GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,Nodes_Elemental_Strains,
           Nodes_Elemental_Stress,Nodes_Elemental_VM)
    
