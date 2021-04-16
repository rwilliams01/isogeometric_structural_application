##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.BRepApplication import *
#from KratosMultiphysics.MultigridSolversApplication import *
#from KratosMultiphysics.FiniteCellApplication import *
kernel = Kernel()   #defining kernel

import model_iga_include
from model_iga_include import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_export2 = MultiNURBSPatchGLVisExporter()

import geometry_factory
'''
sys.path.append("/home/hbui/Benchmark_kratos/structural_application/std_problems/internal_pressurized_cylinder")
import analytical_solution
'''
E = 2.1e11
nu = 0.3
P = 100.0
r1 = 100.0
r2 = 200.0
'''
ana_sol = analytical_solution.Solution(r1, r2, P, 0.0, E, nu)
'''
def CreateMultiPatch():
    ## create arc 1
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 90.0)
    arc1 = arc1_ptr.GetReference()
    arc1.Id = 1

    ## create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r2, 0.0, 90.0)
    arc2 = arc2_ptr.GetReference()
    arc2.Id = 2

    ## create ring patch by connect the two arcs
    ring_patch_ptr = bsplines_patch_util.CreateLoftPatch(arc2, arc1)
    ring_patch = ring_patch_ptr.GetReference()
    ring_patch.Id = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(ring_patch_ptr)

    #### elevate the degree
    multipatch_refine_util.DegreeElevate(mpatch[1], [1, 1])
    # multipatch_refine_util.DegreeElevate(mpatch[1], [2, 3])

    return mpatch

def Refine(mpatch, ins_knots):
    print("###############REFINEMENT###############")
    multipatch_refine_util.InsertKnots(mpatch[1], [ins_knots, ins_knots])

    return mpatch

def CreateModelPart(mpatch):
    mpatch_util = MultiPatchUtility()
    element_name = "KinematicLinearBezier2D"
    load_condition_name = "LinePressureBezier2D"

    mpatch_mp = MultiPatchModelPart2D(mpatch)

    mpatch_mp.BeginModelPart()
    model_part = mpatch_mp.GetModelPart()
    import structural_solver_advanced
    structural_solver_advanced.AddVariables( model_part )
    model_part.AddNodalSolutionStepVariable(THREED_STRESSES)

    mpatch_mp.CreateNodes()

    #problem data
    body_force = ZeroVector(2)
    gravity = Vector(2)
    gravity[0] = 0.0
    gravity[1] = 0.0
    prop = model_part.Properties[1]
    prop.SetValue(GRAVITY, gravity )
    prop.SetValue(BODY_FORCE, body_force )
    prop.SetValue(NUM_IGA_INTEGRATION_METHOD, 2)
    prop.SetValue(INTEGRATION_ORDER, 2)
    prop.SetValue(DENSITY,            0 )
    prop.SetValue(YOUNG_MODULUS, E )
    prop.SetValue(POISSON_RATIO, nu )
    prop.SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
    prop.SetValue(THICKNESS, 1)

    patch_ids = [1]
    for sid in patch_ids:
#        print("sid", sid)
        patch_ptr = mpatch[sid]
        #        print(patch_ptr)
        patch = patch_ptr.GetReference()
        #        print(patch)

        ## add volume elements
        last_elem_id = mpatch_util.GetLastElementId(model_part)
        elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)
    
        ## add loading conditions
        last_cond_id = mpatch_util.GetLastConditionId(model_part)
        load_conds = mpatch_mp.AddConditions(patch, BoundarySide2D.V1, load_condition_name, last_cond_id+1, prop)
        for cond in load_conds:
            cond.SetValue(PRESSURE, -P)
    print(elems)
    mpatch_mp.EndModelPart()
    print(mpatch_mp)

    # fix displacement on the bottom
    tol = 1.0e-6
    for node in model_part.Nodes:
        if abs(node.Y0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_Y)
            node.SetSolutionStepValue(DISPLACEMENT_Y, 0.0)


    # fix displacement on the left
    for node in model_part.Nodes:
        if abs(node.X0 - 0.0) < tol:
            node.Fix(DISPLACEMENT_X)
            node.SetSolutionStepValue(DISPLACEMENT_X, 0.0)

    return mpatch_mp

def CreateModel():
    mpatch = CreateMultiPatch()

    ins_knots = []
    nsampling = 1
    for i in range(1, nsampling):
        ins_knots.append(float(i)/nsampling)
    mpatch = Refine(mpatch, ins_knots)

    mpatch.Enumerate()
    print(mpatch)
    mpatch_export.Export(mpatch, "internal_pressurized_cylinder.m")
    mpatch_export2.Export(mpatch, "internal_pressurized_cylinder.mesh")

    mpatch_mp = CreateModelPart(mpatch)

    #############ANALYSIS MODEL#######################################
    model_part = mpatch_mp.GetModelPart()
    params = model_iga_include.StaticParameters()
    params["builder_and_solver_type"] = "residual-based block"
    model = model_iga_include.Model('internal_pressurized_cylinder', os.getcwd()+"/", model_part, params)
    model.InitializeModel()

    return [mpatch_mp, model]

def main():
    #############ANALYSIS MODEL#######################################
    [mpatch_mp, model] = CreateModel()
    mpatch = mpatch_mp.GetMultiPatch()

    time = 0.0
    model.Solve(time, 0, 0, 0, 0)

    transfer_util = BezierPostUtility()
    transfer_util.TransferVariablesToNodes(THREED_STRESSES, model.model_part, SuperLUSolver())

    ######Synchronize back the results to multipatch
    mpatch_mp.SynchronizeBackward(DISPLACEMENT)
    mpatch_mp.SynchronizeBackward(THREED_STRESSES)
    ##################################################################
    params_post = {}
    params_post['name'] = "internal_pressurized_cylinder"
    params_post['division mode'] = "uniform"
    params_post['uniform division number'] = 20
#    params_post['division mode'] = "non-uniform"
#    params_post['division number u'] = 10
#    params_post['division number v'] = 10
#    params_post['division number w'] = 1
    params_post['variables list'] = [DISPLACEMENT, THREED_STRESSES]
    params_post['backend'] = ["GiD", "Glvis"]
    dim = 2
    model_iga_include.PostMultiPatch(mpatch, dim, time, params_post)

#    ## post processing
#    r1 = 1.0
#    b1 = 4.0
#    disp_grid_function_1 = mpatch[1].GetReference().GridFunction(DISPLACEMENT)
#    top_disp = disp_grid_function_1.GetValue([r1, 0.0, 0.0])
#    print("displacement at (1.0, 0.0, 0.0): " + str(top_disp[2]))
#    disp_grid_function_2 = mpatch[2].GetReference().GridFunction(DISPLACEMENT)
#    left_side_disp = disp_grid_function_2.GetValue([0.0, r1, 0.0])
#    print("displacement at (0.0, 1.0, 0.0): " + str(left_side_disp[2]))

if __name__ == "__main__":
    main()


