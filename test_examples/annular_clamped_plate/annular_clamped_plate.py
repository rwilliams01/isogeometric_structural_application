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
import math
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

sys.path.append('./annular_clamped_plate')
import annular_clamped_plate_include

from annular_clamped_plate_include import *

# Initialize model
model = annular_clamped_plate_include.Model('annular_clamped_plate', os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# user-defined script is used (will be appended automatically)
# =====================
# | USER SCRIPT FOR CALCULATION OF infinite_plate.gid |
# vvvvvvvvvvvvvvvvvvvvv

def get_displacement(p0, x, y, R0, R, B, nu):
    r = (x**2 + y**2)**0.5
    c0 = (p0 * R0**2 * (-6 * R**4 - 5 * R**2 * R0**2 + R0**4 - 2 * R**4 * nu + 3 * R**2 * R0**2 * nu - R0**4 * nu + 4 * R**2 * (R0**2 * (1 + nu) + R**2 * (3 + nu)) * math.log(R0) - 16 * R**4 * (1 + nu) * math.log(R0)**2 + 8 * R**4 * (1 + nu) * math.log(R) * (-1 + 2 * math.log(R0)))) / (64 * B * (-R0**2 * (-1 + nu) + R**2 * (1 + nu)))
    c1 = -p0 * R**2 * R0**2 * (R**2 + R0**2 - R**2 * nu + R0**2 * nu + 4 * R**2 * (1 + nu) * math.log(R) - 4 * R**2 * (1 + nu) * math.log(R0)) / (16 * B * (-R0**2 * (-1 + nu) + R**2 * (1 + nu)))
    c2 = (p0 * (3 * R**4 + 2 * R**2 * R0**2 - R0**4 + R**4 * nu - 2 * R**2 * R0**2 * nu + R0**4 * nu + 4 * R**4 * (1 + nu) * math.log(R) - 4 * R**2 * R0**2 * (-1 + nu) * math.log(R0))) / (32 * B * (-R0**2 * (-1 + nu) + R**2 * (1 + nu)))
    c3 = -p0 * R**2 / (8 * B)
    w = c0 + c1*math.log(r) + c2 * r**2 + c3 * r**2 * math.log(r) + p0 * r**4 / (64*B)
    return w

time = 0.0
model.Solve(time, 0, 0, 0, 0)
#model.WriteOutput(time)

#boundary condition
R0 = 100.0
tol = 1.0e-5
for node in model.model_part.Nodes:
    if abs(node.X0-100) < tol and abs(node.Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-108.333333) < tol and abs(node.Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-100) < tol and abs(node.Y0-9.174737) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-108.333333) < tol and abs(node.Y0-9.939298) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-96.560154) < tol and abs(node.Y0-27.763231) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-104.606833) < tol and abs(node.Y0-30.076834) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-89.247117) < tol and abs(node.Y0-46.235586) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-96.684377 ) < tol and abs(node.Y0-50.088551) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-78.033009) < tol and abs(node.Y0-63.388348) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-84.535759) < tol and abs(node.Y0-68.67071) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-63.388348) < tol and abs(node.Y0-78.033009) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-68.67071) < tol and abs(node.Y0-84.535759) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-46.235586) < tol and abs(node.Y0-89.247117) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-50.088551) < tol and abs(node.Y0-96.684377) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-27.763231) < tol and abs(node.Y0-96.560154) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-30.076834) < tol and abs(node.Y0-104.606833) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-9.174737) < tol and abs(node.Y0-100) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-9.939298) < tol and abs(node.Y0-108.333333) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-0.0) < tol and abs(node.Y0-100) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0-0.0) < tol and abs(node.Y0-108.333333) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )



load = -1.0
for cond in model.model_part.Conditions:
    cond.SetValue(PRESSURE, load)

##solve the model
time = 1.0
model.Solve(time, 0, 0, 0, 0 )
#model.WriteOutput(time )

E = model.model_part.Properties[1].GetValue(YOUNG_MODULUS)
nu = model.model_part.Properties[1].GetValue(POISSON_RATIO)
h = model.model_part.Properties[1].GetValue(THICKNESS)
B = E*h**3 / (12.0*(1.0-nu**2))
R0 = 100.0
R = 200.0

tol = 1.0e-6
for node in model.model_part.Nodes:
    if abs(node.X0) < tol and abs(node.Y0-R) < tol:
        center_node = node
w0 = center_node.GetSolutionStepValue(DISPLACEMENT_Z)
print("displacement at edge: " + str(w0))
w0_ana = get_displacement(load, center_node.X0, center_node.Y0, R0, R, B, nu)
print("analytical displacement at edge: " + str(w0_ana))
print("relative error: " + str(abs(w0-w0_ana)/abs(w0_ana)))


#for node in model.model_part.Nodes:
 #    print(node.GetSolutionStepValue(DISPLACEMENT_Z))
