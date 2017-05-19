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
    if abs(node.X0- model.model_part.Nodes[1].X0) < tol and abs(node.Y0- model.model_part.Nodes[1].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[2].X0) < tol and abs(node.Y0- model.model_part.Nodes[2].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[15].X0) < tol and abs(node.Y0- model.model_part.Nodes[15].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[16].X0) < tol and abs(node.Y0- model.model_part.Nodes[16].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[29].X0) < tol and abs(node.Y0- model.model_part.Nodes[29].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[30].X0) < tol and abs(node.Y0- model.model_part.Nodes[30].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[43].X0) < tol and abs(node.Y0- model.model_part.Nodes[43].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[44].X0) < tol and abs(node.Y0- model.model_part.Nodes[44].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[57].X0) < tol and abs(node.Y0- model.model_part.Nodes[57].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[58].X0) < tol and abs(node.Y0- model.model_part.Nodes[58].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[71].X0) < tol and abs(node.Y0- model.model_part.Nodes[71].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[72].X0) < tol and abs(node.Y0- model.model_part.Nodes[72].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[85].X0) < tol and abs(node.Y0- model.model_part.Nodes[85].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[86].X0) < tol and abs(node.Y0- model.model_part.Nodes[86].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[99].X0) < tol and abs(node.Y0- model.model_part.Nodes[99].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[100].X0) < tol and abs(node.Y0- model.model_part.Nodes[100].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[113].X0) < tol and abs(node.Y0- model.model_part.Nodes[113].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[114].X0) < tol and abs(node.Y0- model.model_part.Nodes[114].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[127].X0) < tol and abs(node.Y0- model.model_part.Nodes[127].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[128].X0) < tol and abs(node.Y0- model.model_part.Nodes[128].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[141].X0) < tol and abs(node.Y0- model.model_part.Nodes[141].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[142].X0) < tol and abs(node.Y0- model.model_part.Nodes[142].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[155].X0) < tol and abs(node.Y0- model.model_part.Nodes[155].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[156].X0) < tol and abs(node.Y0- model.model_part.Nodes[156].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[169].X0) < tol and abs(node.Y0- model.model_part.Nodes[169].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[170].X0) < tol and abs(node.Y0- model.model_part.Nodes[170].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[183].X0) < tol and abs(node.Y0- model.model_part.Nodes[183].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[184].X0) < tol and abs(node.Y0- model.model_part.Nodes[184].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[197].X0) < tol and abs(node.Y0- model.model_part.Nodes[197].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[198].X0) < tol and abs(node.Y0- model.model_part.Nodes[198].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[211].X0) < tol and abs(node.Y0- model.model_part.Nodes[211].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[212].X0) < tol and abs(node.Y0- model.model_part.Nodes[212].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[225].X0) < tol and abs(node.Y0- model.model_part.Nodes[225].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[226].X0) < tol and abs(node.Y0- model.model_part.Nodes[226].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[239].X0) < tol and abs(node.Y0- model.model_part.Nodes[239].Y0) < tol:
        node.Fix( DISPLACEMENT_X )
        node.Fix( DISPLACEMENT_Y )
        node.Fix( DISPLACEMENT_Z )
    if abs(node.X0- model.model_part.Nodes[240].X0) < tol and abs(node.Y0- model.model_part.Nodes[240].Y0) < tol:
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
