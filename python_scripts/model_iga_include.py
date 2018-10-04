##################################################################
################### model_iga_include.py   #######################
##################################################################
##### supplementary module for isogeometric analysis         #####
#####                                            with KRATOS #####
##### copyright Hoang-Giang BUI                              #####
#####              Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
import math
#kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################

#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MKLSolversApplication import *
## REMARK: external applications must be imported separately in the other code
#from KratosMultiphysics.MortarApplication import *
#from KratosMultiphysics.IsogeometricMortarApplication import *
kernel = Kernel()   #defining kernel

##################################################################
### Give the basic parameters for static analysis. This is used for initialization of Model
def StaticParameters():
    analysis_parameters = {}
    # content of analysis_parameters:
    # perform_contact_analysis_flag
    # penalty value for normal contact
    # maximum number of uzawa iterations
    # friction coefficient
    # penalty value for frictional contact
    # contact_double_check_flag
    # contact_ramp_penalties_flag
    # maximum penalty value for normal contact
    # ramp criterion for normal contact
    # ramp factor for normal contact
    # maximum penalty value for frictional contact
    # ramp criterion for frictional contact
    # ramp factor for frictional contact
    # analysis type: static (0), quasi-static (1) or dynamic (2)
    perform_contact_analysis_flag = False
    penalty = 0.0
    maxuzawa = 0.0
    friction = 0.0
    frictionpenalty = 0.0
    contact_double_check_flag = False
    contact_ramp_penalties_flag = False
    maxpenalty = 0.0
    rampcriterion = 0.0
    rampfactor = 0.0
    fricmaxpenalty = 0.0
    fricrampcriterion = 0.0
    fricrampfactor = 0.0
    analysis_parameters['perform_contact_analysis_flag'] = perform_contact_analysis_flag
    analysis_parameters['penalty'] = penalty
    analysis_parameters['maxuzawa'] = maxuzawa
    analysis_parameters['friction'] = friction
    analysis_parameters['frictionpenalty'] = frictionpenalty
    analysis_parameters['contact_double_check_flag'] = contact_double_check_flag
    analysis_parameters['contact_ramp_penalties_flag'] = contact_ramp_penalties_flag
    analysis_parameters['maxpenalty'] = maxpenalty
    analysis_parameters['rampcriterion'] = rampcriterion
    analysis_parameters['rampfactor'] = rampfactor
    analysis_parameters['fricmaxpenalty'] = fricmaxpenalty
    analysis_parameters['fricrampcriterion'] = fricrampcriterion
    analysis_parameters['fricrampfactor'] = fricrampfactor
    analysis_parameters['print_sparsity_info_flag'] = False
    analysis_parameters['analysis_type'] = 0
    analysis_parameters['dissipation_radius'] = 0.1
    analysis_parameters['decouple_build_and_solve'] = False
    analysis_parameters['builder_and_solver_type'] = "residual-based elimination deactivation"
    analysis_parameters['solving_scheme'] = 'monolithic'
    analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
    analysis_parameters['abs_tol'] = 1.0e-10
    analysis_parameters['rel_tol'] = 1.0e-13
    analysis_parameters['echo_level'] = 2
    analysis_parameters['max_iter'] = 10
    return analysis_parameters

### Give the basic parameters for static analysis. This is used for initialization of Model
def QuasiStaticParameters():
    analysis_parameters = StaticParameters()
    analysis_parameters['analysis_type'] = 1
    analysis_parameters['dissipation_radius'] = 0.1
    return analysis_parameters

##################################################################
### Interface class with KRATOS to perform analysis
class Model:
    def __init__( self, problem_name, path, model_part, analysis_parameters ):
        #setting the domain size for the problem to be solved
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = model_part
        self.model_part.Name = problem_name
        self.model_part.SetBufferSize(2)
        self.path = path
        self.problem_name = problem_name
        self.domain_size = 2
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
        number_of_time_steps = 1
        self.analysis_parameters = analysis_parameters

        self.abs_tol = self.analysis_parameters['abs_tol']
        self.rel_tol = self.analysis_parameters['rel_tol']

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, self.abs_tol, self.rel_tol )
        ##################################################################
        ## POST_PROCESSING DEFINITIONS ###################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostAscii
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.isogeometric_classtical_post_utility = BezierClassicalPostUtility(self.model_part)
        self.generate_post_model_part = False

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        structural_solver_advanced.AddDofs( self.model_part )

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = SuperLUSolver() # MKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(self.analysis_parameters['echo_level'])
        (self.solver.solver).max_iter = self.analysis_parameters['max_iter'] #control the maximum iterations of Newton Raphson loop
        (self.solver.solver).CalculateReactionsFlag = False
        (self.solver.solver).MoveMeshFlag = False

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )

    def WriteOutput( self, time ):
        if self.generate_post_model_part == False:
            print ('Before GenerateModelPart')
            self.isogeometric_post_utility.GenerateModelPart(self.model_part_post, PostElementType.Quadrilateral)
            self.isogeometric_post_utility.GenerateModelPart2(self.model_part_post)
            self.generate_post_model_part = True
            print ('Generate PostModelPart completed')
        if self.generate_post_model_part == True:
            self.isogeometric_post_utility.TransferNodalResults(DISPLACEMENT, self.model_part_post)
            self.isogeometric_post_utility.TransferNodalResults(REACTION, self.model_part_post)
#            self.isogeometric_post_utility.TransferIntegrationPointResults(STRESSES, self.model_part_post, self.solver_post)
            print ('Synchronize PostModelPart completed')
        self.gid_io.InitializeMesh( time )
        post_mesh = self.model_part_post.GetMesh()
        print(post_mesh)
        self.gid_io.WriteNodeMesh( post_mesh )
        self.gid_io.WriteMesh( post_mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, post_mesh )
        print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part_post.Nodes, time, 0)
#        self.gid_io.WriteNodalResults(STRESSES, self.model_part_post.Nodes, time, 0)
        self.gid_io.FinalizeResults()

    def InitializeModel( self ):
        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        ## ELEMENTS on layers ############################################

        ## NODES on layers ###############################################

        ## CONTACT MASTER NODES ##########################################

        ## CONTACT SLAVE NODES ###########################################

        ## INNER BOUNDARY NODES ##########################################

        ##################################################################
        #print ("layer sets stored" )
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################

        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        self.model_part.Check(self.model_part.ProcessInfo)
        print ("activation utility initialized")
        ##################################################################
        print ("model successfully initialized")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()

##################################################################
### Create FEM mesh for post
def CreatePostMesh(mpatch, dim, params):
    #################DEFAULT SETTINGS#################################
    if 'name' not in params:
        params['name'] = "iga model"
    if 'base element name' not in params:
        params['base element name'] = "KinematicLinear"
    if 'last node id' not in params:
        params['last node id'] = 1
    if 'last element id' not in params:
        params['last element id'] = 1
    if 'last condition id' not in params:
        params['last condition id'] = 1
    if 'division mode' not in params:
        params['division mode'] = "uniform"
    if 'uniform division number' not in params:
        params['uniform division number'] = 10
    # default variables list
    if 'variables list' not in params:
        params['variables list'] = []
        params['variables list'].append(DISPLACEMENT)

    #################POST PROCESSING##################################
    if dim == 2:
        fem_mesh = NonConformingMultipatchLagrangeMesh2D(mpatch)
    elif dim == 3:
        fem_mesh = NonConformingMultipatchLagrangeMesh3D(mpatch)
    fem_mesh.SetBaseElementName(params['base element name'])
    fem_mesh.SetLastNodeId(params['last node id'])
    fem_mesh.SetLastElemId(params['last element id'])
    fem_mesh.SetLastPropId(params['last condition id'])
    if params['division mode'] == "uniform":
        fem_mesh.SetUniformDivision(params['uniform division number'])
    elif params['division mode'] == "non-uniform":
        for patch_ptr in mpatch.Patches():
            patch = patch_ptr.GetReference()
            fem_mesh.SetDivision(patch.Id, 0, params['division number u'])
            fem_mesh.SetDivision(patch.Id, 1, params['division number v'])
            fem_mesh.SetDivision(patch.Id, 2, params['division number w'])

    return fem_mesh

### Write a post mesh to GiD
def WriteGiD(post_model_part, time, params):
    #######WRITE TO GID
    write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
    write_elements = WriteConditionsFlag.WriteConditions
    #write_elements = WriteConditionsFlag.WriteElementsOnly
    post_mode = GiDPostMode.GiD_PostBinary
    multi_file_flag = MultiFileFlag.MultipleFiles
    gid_io = StructuralGidIO(params['name'], post_mode, multi_file_flag, write_deformed_flag, write_elements)
    gid_io.InitializeMesh( time )
    post_mesh = post_model_part.GetMesh()
    gid_io.WriteMesh( post_mesh )
    print("mesh written...")
    gid_io.FinalizeMesh()
    gid_io.InitializeResults( time, post_mesh )
    print("write nodal results")
    for var in params['variables list']:
        gid_io.WriteNodalResults(var, post_model_part.Nodes, time, 0)
    gid_io.FinalizeResults()

### Create a post_model_part out from mpatch
def CreatePostModelPart(mpatch, dim, params):
    fem_mesh = CreatePostMesh(mpatch, dim, params)

    post_model_part = ModelPart("iga-fem mesh " + params['name'])
    for var in params['variables list']:
        post_model_part.AddNodalSolutionStepVariable(var)
    fem_mesh.WriteModelPart(post_model_part)

    return post_model_part

### Post-process a multipatch
def PostMultiPatch(mpatch, dim, time, params):
    post_model_part = CreatePostModelPart(mpatch, dim, params)
    print(post_model_part)
    WriteGiD(post_model_part, time, params)

##################################################################
### Compute the strain energy of the model_part
def ComputeStrainEnergy(model_part):
    senergy = 0.0
    for element in model_part.Elements:
        if element.GetValue(IS_INACTIVE) == False:
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
            SE = element.GetValuesOnIntegrationPoints(STRAIN_ENERGY, model_part.ProcessInfo)
            for i in range(0, len(W)):
                senergy = senergy + SE[i][0] * W[i][0] * J0[i][0]
    return senergy

### Compute the absolute L2-error on element
def ComputeL2errorOnElement(element, analytical_solution, process_info):
    error = 0.0
    u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, process_info)
    J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
    Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
    W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
    for i in range(0, len(u)):
        ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], Q[i][2])
        error = error + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
    return error

### Compute the relative L2-error
def ComputeL2error(model_part, analytical_solution):
    nom = 0.0
    denom = 0.0
    for element in model_part.Elements:
        if element.GetValue(IS_INACTIVE) == False:
            u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for i in range(0, len(u)):
                ana_u = analytical_solution.get_displacement(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2)) * W[i][0] * J0[i][0]
    error = math.sqrt(nom / denom)
    print("Global displacement (L2) error:", error)
    return error

### Compute the absolute H1-error on element
def ComputeH1error(element, analytical_solution):
    error = 0.0
    o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model_part.ProcessInfo)
    J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
    Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
    W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
    for i in range(0, len(o)):
        ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], Q[i][2])
        error = error + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
    return error

### Compute the relative H1-error
def ComputeH1error(model_part, analytical_solution):
    nom = 0.0
    denom = 0.0
    for element in model_part.Elements:
        if element.GetValue(IS_INACTIVE) == False:
            o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model_part.ProcessInfo)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for i in range(0, len(o)):
                ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], Q[i][2])
                nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
    error = math.sqrt(nom / denom)
    print("Global stress (H1) error:", math.sqrt(nom / denom))
    return error

### Compute the absolute H1-error on element
def ComputeH1errorOnElement(element, analytical_solution, process_info):
    error = 0.0
    o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, process_info)
    J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
    Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
    W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
    for i in range(0, len(o)):
        ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], Q[i][2])
        error = error + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
    return error

### Compute energy norm error
### REF: https://www.sharcnet.ca/Software/Ansys/17.0/en-us/help/ans_vm/Hlp_V_CH2_5.html
def ComputeEnergyError(model_part, analytical_solution):
    senergy = 0.0
    ana_senergy = 0.0
    error = 0.0
    for element in model_part.Elements:
        if element.GetValue(IS_INACTIVE) == False:
            o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, model_part.ProcessInfo)
            e = element.GetValuesOnIntegrationPoints(THREED_STRAIN, model_part.ProcessInfo)
            J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
            Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
            W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
            for i in range(0, len(o)):
                senergy = senergy + 0.5 * (e[i][0]*o[i][0] + e[i][1]*o[i][1] + e[i][2]*o[i][2] + e[i][3]*o[i][3] + e[i][4]*o[i][4] + e[i][5]*o[i][5]) * W[i][0] * J0[i][0]
                ana_o = analytical_solution.get_stress_3d(Q[i][0], Q[i][1], Q[i][2])
                ana_e = analytical_solution.get_strain_3d(Q[i][0], Q[i][1], Q[i][2])
                ana_senergy = ana_senergy + 0.5 * (ana_e[0]*ana_o[0] + ana_e[1]*ana_o[1] + ana_e[2]*ana_o[2] + ana_e[3]*ana_o[3] + ana_e[4]*ana_o[4] + ana_e[5]*ana_o[5]) * W[i][0] * J0[i][0]
                error = error + 0.5 * ((e[i][0] - ana_e[0])*(o[i][0] - ana_o[0]) + (e[i][1] - ana_e[1])*(o[i][1] - ana_o[1]) + (e[i][2] - ana_e[2])*(o[i][2] - ana_o[2]) + (e[i][3] - ana_e[3])*(o[i][3] - ana_o[3]) + (e[i][4] - ana_e[4])*(o[i][4] - ana_o[4]) + (e[i][5] - ana_e[5])*(o[i][5] - ana_o[5])) * W[i][0] * J0[i][0]
    error = math.sqrt(error/ana_senergy)
    print("Global energy norm error:", error)
    return [error, senergy]


#############################################################


