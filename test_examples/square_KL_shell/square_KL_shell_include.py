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
#kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path ):
        #setting the domain size for the problem to be solved
        self.domain_size = 3
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("isogeometric_simulation")
        self.model_part_post = ModelPart("isogeometric_mesh")
        self.path = path
        self.problem_name = problem_name
        ##################################################################
        ## DEFINE SOLVER #################################################
        ##################################################################
        # reading simulation parameters
        number_of_time_steps = 1
        self.analysis_parameters = []
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
        self.analysis_parameters.append(perform_contact_analysis_flag)
        self.analysis_parameters.append(penalty)
        self.analysis_parameters.append(maxuzawa)
        self.analysis_parameters.append(friction)
        self.analysis_parameters.append(frictionpenalty)
        self.analysis_parameters.append(contact_double_check_flag)
        self.analysis_parameters.append(contact_ramp_penalties_flag)
        self.analysis_parameters.append(maxpenalty)
        self.analysis_parameters.append(rampcriterion)
        self.analysis_parameters.append(rampfactor)
        self.analysis_parameters.append(fricmaxpenalty)
        self.analysis_parameters.append(fricrampcriterion)
        self.analysis_parameters.append(fricrampfactor)
        #PrintSparsityInfoFlag
        self.analysis_parameters.append(False)
        self.analysis_parameters.append(0)
        
        abs_tol =     1.0e-10
        rel_tol =   1e-6
        
        ## generating solver
        import structural_solver_advanced
        self.solver = structural_solver_advanced.SolverAdvanced( self.model_part, self.domain_size, number_of_time_steps, self.analysis_parameters, abs_tol, rel_tol )
        structural_solver_advanced.AddVariables( self.model_part )
        
        self.model_part.AddNodalSolutionStepVariable(STRESSES)
        self.model_part_post.AddNodalSolutionStepVariable(DISPLACEMENT)
        self.model_part_post.AddNodalSolutionStepVariable(REACTION)
        self.model_part_post.AddNodalSolutionStepVariable(STRESSES)
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostAscii
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = BezierModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.isogeometric_post_utility = IsogeometricClassicalPostUtility(self.model_part)
        self.generate_post_model_part = False
        self.meshWritten = True
        self.solver.CalculateReactionFlag = True
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        print ("input data read OK")
        #print "+++++++++++++++++++++++++++++++++++++++"
        #for node in self.model_part.Nodes:
            #print node
        #print "+++++++++++++++++++++++++++++++++++++++"
        
        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################        
        structural_solver_advanced.AddDofs( self.model_part )
        structural_solver_advanced.AddDofs( self.model_part_post )
        #ekate_solver_parallel.AddDofs( self.model_part )

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        #defining linear solver
        plinear_solver = SkylineLUFactorizationSolver()
        self.solver.structure_linear_solver = plinear_solver
        self.solver.Initialize()
        (self.solver.solver).SetEchoLevel(2);
        
        #defined linear solver for post-processing
#        self.solver_post = SkylineLUFactorizationSolver()
        self.solver_post = SkylineLUFactorizationSolver()

        ##################################################################
        ## INITIALISE RESTART UTILITY ####################################
        ##################################################################
        #restart_utility= RestartUtility( self.problem_name )
        
    #def SetUpActivationLevels( self, model_part, activation_list, cond_activation_list ):
     #   for element in self.model_part.Elements:
      #      element.SetValue(ACTIVATION_LEVEL, activation_list[element.Id])
       # for condition in self.model_part.Conditions:
        #    if( not (condition.GetValue(IS_TYING_MASTER) or condition.GetValue(IS_CONTACT_MASTER) ) ):
         #       condition.SetValue(ACTIVATION_LEVEL, activation_list[cond_activation_list[condition.Id-1]])

    #def write_restart_file( self, time ):
     #   print("------------> restart file written for time step: "+str(time))
      #  self.restart_utility.ChangeFileName(problem_name+str(time))
       # self.restart_utility.StoreNodalVariables(model_part)
        #self.restart_utility.StoreInSituStress(model_part)
        #self.restart_utility.StoreConstitutiveLawVariables(model_part)

#    def restart_time_step( self, time, Dt ):
#        print("############ time step solution has to be restarted ############")
#        time = time-Dt
#        model_part.CloneTimeStep(time)
#        for step in range(1,11):
#            time = time+ Dt/10.0
#            model_part.CloneTimeStep(time)
#            #####################################################################################################
#            model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
#            model_part.ProcessInfo.SetValue( FIRST_TIME_STEP, False )
#            #####################################################################################################
#            solver.Solve()
#            print("~~~~~~~~~~~~~~ RESTARTED STEP ( DT= "+str(Dt/10.0)+" / Step= "+str(step)+" ) ~~~~~~~~~~~~~~")
#        print("############ restart finished ############")
#
#    def write_to_file( self, time ):
#        for i in range(0, len(self.layer_nodes_sets['top'])):
#            settlements.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['top'][i]].GetSolutionStepValue(DISPLACEMENT_Z))+"\n")
#    for i in range(0, len(layer_nodes_sets['side'])):
#        pressure_air.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(AIR_PRESSURE))+"\n")
#        pressure_water.write(str(time)+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetZ())+"/"+str(model_part.Nodes[layer_nodes_sets['side'][i]].GetSolutionStepValue(WATER_PRESSURE))+"\n")
#
            
    def WriteOutput( self, time ):
        if self.generate_post_model_part == False:
            print ('Before GenerateModelPart')
#            self.isogeometric_post_utility.GenerateModelPart(self.model_part_post, PostElementType.Quadrilateral)
            self.isogeometric_post_utility.GenerateModelPart2(self.model_part_post)
            self.generate_post_model_part = True
            print ('Generate PostModelPart completed')
        if self.generate_post_model_part == True:
            self.isogeometric_post_utility.TransferNodalResults(DISPLACEMENT, self.model_part_post)
            self.isogeometric_post_utility.TransferNodalResults(REACTION, self.model_part_post)
            self.isogeometric_post_utility.TransferIntegrationPointResults(STRESSES, self.model_part_post, self.solver_post)
            print ('Synchronize PostModelPart completed')
        self.gid_io.InitializeMesh( time )
        post_mesh = self.model_part_post.GetMesh()
        print(post_mesh)
        #self.gid_io.WriteNodeMesh( post_mesh )
        self.gid_io.WriteMesh( post_mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, post_mesh )
        print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part_post.Nodes, time, 0)
        self.gid_io.WriteNodalResults(STRESSES, self.model_part_post.Nodes, time, 0)
        self.gid_io.FinalizeResults()
                
    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
        
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStress() )
        self.model_part.Properties[1].SetValue(THICKNESS, 1.0 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,1.0000e+7)
        self.model_part.Properties[1].SetValue(POISSON_RATIO,0.0 )
        self.model_part.Properties[1].SetValue(PRESSURE,-1)

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
        
        #self.deac = DeactivationUtility()
        
        #self.SetUpActivationLevels( self.model_part, self.activation_flags, self.cond_activation_flags )
        #self.deac.Initialize( self.model_part )
        #self.model_part.Check(self.model_part.ProcessInfo)
        #sys.exit()
        #print ("activation utility initialized")
        #print ("model successfully initialized")
    
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
