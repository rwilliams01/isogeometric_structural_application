##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your actual configuration                #####
##################################################################
import sys
import os
#kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
#kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
#kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
#sys.path.append(kratos_libs_path)
#sys.path.append(kratos_applications_path)

##################################################################
##################################################################

sys.path.append('./square_KL_shell')
import square_KL_shell_include

from square_KL_shell_include import *

# Initialize model
model = square_KL_shell_include.Model('square_KL_shell', os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF infinite_plate.gid |
# vvvvvvvvvvvvvvvvvvvvv

time = 0.0
model.Solve(time, 0, 0, 0, 0)
model.WriteOutput(time)

#boundary condition
R0 = 100.0
tol = 1.0e-6
for node in model.model_part.Nodes:
    if abs(node.X0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )

    if abs(node.X0-100) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )

    if abs(node.Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )

    if abs(node.Y0-100) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )




for cond in model.model_part.Conditions:
    cond.SetValue(PRESSURE, -1)

##solve the model
time = 1.0
model.Solve(time, 0, 0, 0, 0 )
model.WriteOutput(time )

