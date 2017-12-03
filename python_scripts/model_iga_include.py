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
#kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
kernel = Kernel()   #defining kernel

##################################################################
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
    analysis_parameters['solving_scheme'] = 'monolithic'
    analysis_parameters['stop_Newton_Raphson_if_not_converge'] = True
    analysis_parameters['abs_tol'] = 1.0e-10
    analysis_parameters['rel_tol'] = 1.0e-13
    analysis_parameters['echo_level'] = 2
    analysis_parameters['max_iter'] = 10
    return analysis_parameters

##################################################################
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

        abs_tol = self.analysis_parameters['abs_tol']
        rel_tol = self.analysis_parameters['rel_tol']

        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        self.solver.CalculateReactionFlag = False
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
        plinear_solver = MKLPardisoSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(self.analysis_parameters['echo_level'])
        (self.solver.solver).max_iter = self.analysis_parameters['max_iter'] #control the maximum iterations of Newton Raphson loop

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
        print ("model successfully initialized")
    
    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        if( mpi.rank == 0 ):
            self.mergefile.close()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        #self.deac.Reactivate( self.model_part, from_reac, to_reac )
        #self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
        
##################################################################
