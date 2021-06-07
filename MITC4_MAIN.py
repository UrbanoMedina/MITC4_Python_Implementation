
import MITC4_PRE_PROCESSING_Patch_Test_txx
import Processing
import MITC4_Post_Processing
import numpy as np
import MITC4_PRE_PROCESSING_Patch_Test_txy
import MITC4_PRE_PROCESSING_Patch_Test_tyx
import MITC4_PRE_PROCESSING_Patch_Test_Bending_txx
import MITC4_PRE_PROCESSING_Patch_Test_Shear_txz
import MITC4_PRE_PROCESSING_Patch_Test_Twist_txy
import MITC4_PRE_PROCESSING_Twist_16x16_Rect_Mesh

#MAIN MITC4#


# #Patch Test membrane txx#

# #FIRST: PRE-PROCESSING
# #Obtain from the case associated pre-processing file all the
# # information requiered for the analysis. Note: This function has to be imported.
# #NOTE: this file must be imported

# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_txx.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)
 

# xcoord=np.zeros([4,Element_Connect.shape[0]])
# ycoord=np.zeros([4,Element_Connect.shape[0]])
# zcoord=np.zeros([4,Element_Connect.shape[0]])

# for count,i in enumerate(Element_Connect):
#     for j in range(4):
#         xcoord[j,count]=Nodes_Coord[i[j]-1][0]
#         ycoord[j,count]=Nodes_Coord[i[j]-1][1]
#         zcoord[j,count]=Nodes_Coord[i[j]-1][2]
        
# #Patch Test shear txy#

# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_txy.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)
 

# #Patch Test shear txx bending#
    
# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_Bending_txx.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)
 
# #Patch Test shear txz bending#

# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_Shear_txz.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)
        

# #Patch Test twist txz bending#
    
# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_Shear_txz.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)


# #Patch Test twist txy bending#
    
# Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Patch_Test_Twist_txy.Node_Element_Info()
# #SECOND: PROCESSING 
# #Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
# Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
#            GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
#            Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
#            Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
#     Nodes_Vn,Element_Connect,bcdof,R)
# #THIRD: POST-PROCESSING 
# #Graph the deformed surface
# Mag=100 #Note: displacement magnification in order to visualize the deformed surface
# figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)
    
#Patch Test twist txy bending refined mesh#
 

Material,Nodes_Coord,Nodes_Thickness,\
     Nodes_Vn,Element_Connect,bcdof,R=MITC4_PRE_PROCESSING_Twist_16x16_Rect_Mesh.Node_Element_Info()
#SECOND: PROCESSING 
#Given the PRE-PROCESSING INFORMATION, SOLVE THE CASE.
Nodes_Coord,Element_Connect,nnode,U_unknown,U_Nodal,System_Nodal_Reactions,\
           GP_Elemental_Strains,GP_Elemental_Stress,GP_Elemental_VM,\
           Nodes_Elemental_Strains,Nodes_Elemental_Stress,\
           Nodes_Elemental_VM= Processing.MITC4_Processing(Material,Nodes_Coord,Nodes_Thickness,\
    Nodes_Vn,Element_Connect,bcdof,R)
#THIRD: POST-PROCESSING 
#Graph the deformed surface
Mag=100 #Note: displacement magnification in order to visualize the deformed surface
figure=MITC4_Post_Processing.Graph(Mag,U_unknown,Nodes_Coord,Element_Connect)






        