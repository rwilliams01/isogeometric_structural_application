#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.IsogeometricApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import sys

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(REFERENCE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ELASTIC_BEDDING_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_LAGRANGE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD)
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS)
    model_part.AddNodalSolutionStepVariable(PRESTRESS)
    model_part.AddNodalSolutionStepVariable(STRESSES)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXCESS_PORE_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    #auxiliary variables misused for mesh rezoning ;-)
    model_part.AddNodalSolutionStepVariable(IS_VISITED)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_MULTIPLIER)
    model_part.AddNodalSolutionStepVariable(GAP)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(LAGRANGE_WATER_PRESSURE)
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES)
    model_part.AddNodalSolutionStepVariable(MOMENTUM)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    print "variables for the dynamic structural solution added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(WATER_PRESSURE)
        node.AddDof(AIR_PRESSURE)
        node.AddDof(ROTATION_X)
        node.AddDof(ROTATION_Y)
        node.AddDof(ROTATION_Z)
    print "dofs for the dynamic structural solution added correctly"

#######################################################################
class SolverAdvanced():
    def __init__( self, model_part, domain_size, analysis_parameters, abs_tol, rel_tol ):
        self.model_part = model_part
        self.MaxNewtonRapshonIterations = 100
        self.analysis_parameters = self.CheckAndConvertParameters(analysis_parameters)
        self.echo_level = 0
        self.dissipation_radius = self.analysis_parameters['dissipation_radius']
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        #definition of the solvers
        self.structure_linear_solver = SuperLUSolver()
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.000001,1e-9)
        self.CalculateReactionFlag = False
        #######################################################################

    def CheckAndConvertParameters(self, analysis_parameters):
        if( type( analysis_parameters ) == dict ):
            return analysis_parameters
        elif( type( analysis_parameters ) == list ):
            new_analysis_parameters = {}
            new_analysis_parameters['perform_contact_analysis_flag'] = analysis_parameters[0]
            new_analysis_parameters['penalty'] = analysis_parameters[1]
            new_analysis_parameters['maxuzawa'] = analysis_parameters[2]
            new_analysis_parameters['friction'] = analysis_parameters[3]
            new_analysis_parameters['frictionpenalty'] = analysis_parameters[4]
            new_analysis_parameters['contact_double_check_flag'] = analysis_parameters[5]
            new_analysis_parameters['contact_ramp_penalties_flag'] = analysis_parameters[6]
            new_analysis_parameters['maxpenalty'] = analysis_parameters[7]
            new_analysis_parameters['rampcriterion'] = analysis_parameters[8]
            new_analysis_parameters['rampfactor'] = analysis_parameters[9]
            new_analysis_parameters['fricmaxpenalty'] = analysis_parameters[10]
            new_analysis_parameters['fricrampcriterion'] = analysis_parameters[11]
            new_analysis_parameters['fricrampfactor'] = analysis_parameters[12]
            new_analysis_parameters['print_sparsity_info_flag'] = analysis_parameters[13]
            new_analysis_parameters['analysis_type'] = analysis_parameters[14]
            if(len(analysis_parameters) > 15):
                new_analysis_parameters['dissipation_radius'] = analysis_parameters[15]
            else:
                if new_analysis_parameters['analysis_type'] == 2:
                    new_analysis_parameters['dissipation_radius'] = 0.1
                else:
                    new_analysis_parameters['dissipation_radius'] = 1.0
            new_analysis_parameters['decouple_build_and_solve'] = False
            return new_analysis_parameters
        else:
            print 'unsupported type of analysis parameters'
            sys.exit(0)

        #######################################################################

    def Initialize(self):
        #definition of time integration scheme
        if( self.analysis_parameters['analysis_type'] == 0 ):
            print("using static scheme")
            self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
            #self.time_scheme = ParallelResidualBasedIncrementalUpdateStaticScheme()
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 1 ):
            print("using newmark quasi-static scheme")
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, True )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0):
                self.time_scheme = ResidualBasedNewmarkScheme(self.dissipation_radius)
            else:
                self.time_scheme = ResidualBasedNewmarkScheme() #pure Newmarkscheme
            self.MoveMeshFlag = True
        elif( self.analysis_parameters['analysis_type'] == 2 ):
            print("using newmark dynamic scheme")
            self.model_part.ProcessInfo.SetValue( QUASI_STATIC_ANALYSIS, False )
            if(self.dissipation_radius >= 0.0 and self.dissipation_radius <= 1.0):
                self.time_scheme = ResidualBasedNewmarkScheme(self.dissipation_radius)
            else:
                self.time_scheme = ResidualBasedNewmarkScheme() #pure Newmarkscheme
            #self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.dissipation_radius)
            #self.time_scheme.Check(self.model_part)
        else:
            print("analysis type is not defined! Define in analysis_parameters['analysis_type']:")
            print("   'using static scheme': static analysis")
            print("   'using newmark quasi-static scheme': quasi-static analysis")
            print("   'using newmark dynamic scheme': dynamic analysis")
            sys.exit(0)
        #definition of the convergence criteria
        self.conv_criteria = MultiPhaseFlowCriteria(self.toll,self.absolute_tol)
        if(self.analysis_parameters['decouple_build_and_solve'] == False):
            builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
        else:
            builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(LinearSolver())

        #creating the solution strategy
        self.ReformDofSetAtEachStep = True
        #KLUDGE: this has to be True!
        self.MoveMeshFlag = True
        self.space_utils = UblasSparseSpace()
        #self.space_utils = ParallelUblasSparseSpace()
        #importing strategy
        #import ekate_strategy
        self.solver = UzawaSolvingStrategy( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )

    #######################################################################
    def Solve(self):
        (self.solver).Solve()

    #######################################################################
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    #######################################################################
    def SolveLagrange(self):
        (self.solver).SolveLagrange()

    #######################################################################
    def SolveLagrangeLocal(self):
        (self.solver).SolveLagrangeLocal()

class UzawaSolvingStrategy:
    #######################################################################
    def __init__( self, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Parameters, space_utils, builder_and_solver ):
        #save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.Parameters = Parameters
        self.PerformContactAnalysis = self.Parameters['perform_contact_analysis_flag']
        self.PrintSparsity = self.Parameters['print_sparsity_info_flag']
        self.space_utils = space_utils
        #contact utility
        self.cu = IsogeometricContactUtility( 3 )
        #default values for some variables
        self.max_iter = 30
        self.echo_level = 1
        self.builder_and_solver = builder_and_solver

        #local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()
        ##local matrices and vectors
        #self.A = CompressedMatrix()
        #self.Dx = Vector()
        #self.b = Vector()

        #initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False
        #provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag);
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep);

        self.solveCounter = 0; #hbui added this variable


    #######################################################################
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)
        if (self.scheme.ElementsAreInitialized() == False):
            self.scheme.InitializeElements(self.model_part)

    #######################################################################
    def SolveLagrange( self ):
        #print self.model_part
        ## - storing original condition size before adding virtual conditions.
        ## - performing contact search
        ## - creating virtual link conditions for the assembling
        if( self.PerformContactAnalysis == False ):
            self.PerformNewtonRaphsonIteration()
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory
            if(self.ReformDofSetAtEachStep == True):
                self.Clear();
            return
        print "setting up contact conditions"
        last_real_node = len(self.model_part.Nodes)
        originalPosition = self.cu.SetUpContactConditionsLagrangeTying(self.model_part )
        self.PerformNewtonRaphsonIteration()
        self.cu.CleanLagrangeTying( self.model_part, originalPosition, last_real_node )
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();


    #######################################################################
    def Solve(self):
        #print self.model_part
        ## - storing original condition size before adding virtual conditions.
        ## - performing contact search
        ## - creating virtual link conditions for the assembling
        self.solveCounter = self.solveCounter + 1
        if( self.PerformContactAnalysis == False ):
            self.PerformNewtonRaphsonIteration()
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory
            if(self.ReformDofSetAtEachStep == True):
                self.Clear();
            return
        print "setting up contact conditions"
        originalPosition =  self.cu.SetUpContactConditions(self.model_part, self.Parameters['penalty'], self.Parameters['frictionpenalty'], self.Parameters['contact_double_check_flag'] )
        uzawaConverged = False
        ##  First step: reform DOF set and check if uzawa iteration is necessary
        self.PerformNewtonRaphsonIteration()
        self.cu.Update( self.model_part, originalPosition, self.Parameters['friction'], self.Parameters['contact_ramp_penalties_flag'], self.Parameters['rampcriterion'], self.Parameters['fricrampcriterion'], self.Parameters['rampfactor'], self.Parameters['fricrampfactor'], self.Parameters['maxpenalty'], self.Parameters['fricmaxpenalty']  )
        if( self.cu.IsConverged( self.model_part, 0,  originalPosition, self.Parameters['friction'] ) == True ):
            uzawaConverged = True
            (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
            self.cu.Clean( self.model_part, originalPosition );
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory
            if(self.ReformDofSetAtEachStep == True):
                self.Clear()
            return
        ## beginning of UZAWA loop
        (self.builder_and_solver).SetReshapeMatrixFlag(False)
        for uzawaStep in range(1, self.Parameters['maxuzawa'] ):
            print "I am inside the uzawa loop, iteration no. " + str(uzawaStep)
            ## solving the standard newton-raphson iteration
            self.PerformNewtonRaphsonIteration()
            ## updating the lagrange multipliers
            self.cu.Update( self.model_part, originalPosition, self.Parameters['friction'], self.Parameters['contact_ramp_penalties_flag'], self.Parameters['rampcriterion'], self.Parameters['fricrampcriterion'], self.Parameters['rampfactor'], self.Parameters['fricrampfactor'], self.Parameters['maxpenalty'], self.Parameters['fricmaxpenalty']  )
            ## checking convergence
            if( self.cu.IsConverged( self.model_part, uzawaStep, originalPosition, self.Parameters['friction'] ) == True ):
                uzawaConverged = True
                break
        if(self.Parameters['maxuzawa'] == 1):
            print "Congratulations. Newton-Raphson loop has converged."
        else:
            if( uzawaConverged == False ):
                print "That's bad. Uzawa algorithm failed to converge within maximum number of iterations."
            else:
                print 'Congratulations. Uzawa loop has converged.'
        ### end of UZAWA loop
        ### cleaning up the conditions
        self.cu.Clean( self.model_part, originalPosition )
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();

    def PerformNewtonRaphsonIteration( self ):
        print("time = " + str(self.model_part.ProcessInfo[TIME]))
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True
        #perform initializations for the current step
        #this operation implies:
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            self.InitializeSolutionStep()
            self.SolutionStepIsInitialized = True
        #perform prediction
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        self.iterationCounter = 0 #hbui added this variable
        self.iterationCounter = self.iterationCounter + 1
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

        original_penalty = 0.0
        if( self.PerformContactAnalysis == True ):
            for cond in self.model_part.Conditions:
                if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                    original_penalty = cond.GetValue( PENALTY )[0]
                    break

        #non linear loop
        converged = False
        it = 0
        while(it < self.max_iter and converged == False):
            #increase penalty...
            if( self.PerformContactAnalysis == True ):
                for cond in self.model_part.Conditions:
                    if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                        penalty = cond.GetValue( PENALTY )
                        for i in range(0,len(penalty)):
                            penalty[i] = 1.0*penalty[i]
                        cond.SetValue( PENALTY, penalty )
            #end of increase penalty
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            self.iterationCounter = self.iterationCounter + 1
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

            #verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #update iteration count
            it = it + 1

        if( it == self.max_iter and converged == False):
            print("Iteration did not converge")
            sys.exit("Stop, the time step did not converge at time step " + str(self.model_part.ProcessInfo[TIME]))

        if( self.PerformContactAnalysis == True ):
            for cond in self.model_part.Conditions:
                if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                    penalty = cond.GetValue( PENALTY )
                    for i in range(0,len(penalty)):
                        penalty[i] = original_penalty
                    cond.SetValue( PENALTY, penalty )

    #######################################################################
    #######################################################################

    #######################################################################
    def Predict(self):
        self.scheme.Predict(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

    #######################################################################
    def InitializeSolutionStep(self):
        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            #initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.scheme,self.model_part);
            #reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)
            #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA,self.pDx,self.pb,self.model_part.Elements,self.model_part.Conditions,self.model_part.ProcessInfo)
            #updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()
        if(self.SolutionStepIsInitialized == False):
            self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
            self.scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)

    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag):
        #reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)
        self.space_utils.SetToZeroVector(self.b)

        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        #provide data for the preconditioner and linear solver
        self.linear_solver.ProvideAdditionalData(self.A,self.Dx,self.b,self.builder_and_solver.GetDofSet(),self.model_part)

        #build and solve the problem
        if(self.Parameters['decouple_build_and_solve'] == False):
            self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
        else:
            self.builder_and_solver.Build(self.scheme,self.model_part,self.A,self.b)
            self.linear_solver.Solve(self.A,self.Dx,self.b)

        #full output if needed
        if( self.PrintSparsity ):
            #hbui edited
            #self.PlotSparsityScheme( self.A )
#            wr = UblasMatrixIO()
#            wr.WriteHB(self.A, self.b, "matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".hb.dat")
#            self.space_utils.WriteMatrixMarketMatrix("matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".mm",self.A,False)
            petsc_utils.DumpUblasCompressedMatrixVector("tempAb", self.A, self.b, False)

        if(echo_level >= 3):
            print "SystemMatrix = ", self.A
        #printA = []
        #printdx = []
        #printb = []
        #for i in range(0,len(self.Dx)):
        #    if( abs(self.Dx[i]) < 1.0e-10 ):
         #       printdx.append(0.0)
         #   else:
         #       printdx.append(self.Dx[i])
         #   if( abs(self.b[i]) < 1.0e-6 ):
         #       printb.append(0.0)
         #   else:
         #       printb.append(self.b[i])
         #   row = []
         #   for j in range(0,len(self.Dx)):
         #       if( abs(self.A[(i,j)]) < 1.0 ):
         #           row.append( 0.0 )
         #       else:
         #           row.append(self.A[(i,j)])
         #   printA.append(row)
            print "solution obtained = ", self.Dx
            #formatted_printdx = [ '%.6f' % elem for elem in printdx ]
            #print formatted_printdx
            #formatted_printb = [ '%.4f' % elem for elem in printb ]
            print "RHS = ", self.b
        #print formatted_printb
        #print "Matrix: "
        #for i in range(0,len(self.Dx)):
        #    formatted_printA = [ '%.1f' % elem for elem in printA[i] ]
        #    print(formatted_printA)
        self.AnalyseSystemMatrix(self.A)

        #perform update
        self.scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

        #move the mesh as needed
        if(MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);

        #to account for prescribed displacement, the displacement at prescribed nodes need to be updated
        for node in self.model_part.Nodes:
            if node.IsFixed(DISPLACEMENT_X):
                curr_disp = node.GetSolutionStepValue(DISPLACEMENT_X)
                delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, curr_disp + delta_disp)
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, 0.0) # set the prescribed displacement to zero to avoid update in the second step
            if node.IsFixed(DISPLACEMENT_Y):
                curr_disp = node.GetSolutionStepValue(DISPLACEMENT_Y)
                delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, curr_disp + delta_disp)
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, 0.0) # set the prescribed displacement to zero to avoid update in the second step
            if node.IsFixed(DISPLACEMENT_Z):
                curr_disp = node.GetSolutionStepValue(DISPLACEMENT_Z)
                delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, curr_disp + delta_disp)
                node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z, 0.0) # set the prescribed displacement to zero to avoid update in the second step

        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        #calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag == True):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0

        return normDx

    #######################################################################
    def FinalizeSolutionStep(self,CalculateReactionsFlag):
        if(CalculateReactionsFlag == True):
            self.builder_and_solver.CalculateReactions(self.scheme,self.model_part,self.A,self.Dx,self.b)

        #Finalisation of the solution step,
        self.scheme.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.scheme.Clean()
        #reset flags for the next step
        self.SolutionStepIsInitialized = False

    #######################################################################
    def Clear(self):
        self.space_utils.ClearMatrix(self.pA)
        self.space_utils.ResizeMatrix(self.A,0,0)

        self.space_utils.ClearVector(self.pDx)
        self.space_utils.ResizeVector(self.Dx,0)

        self.space_utils.ClearVector(self.pb)
        self.space_utils.ResizeVector(self.b,0)

        #updating references
        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)
        self.builder_and_solver.Clear()

    #######################################################################
    def SetEchoLevel(self,level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

#######################################################################
    def AnalyseSystemMatrix(self,  A):
        max = 0.0
        for i in range(0,  A.Size1()):
           if( abs(A[(i, i)]) > max ):
               max = A[(i, i)]

#        nonzeros = 0
#        for i in range(0,  A.Size1()):
#            for j in range(0,  A.Size2()):
#                if( abs(A[(i, j)]) > 1e-16 ):
#                    nonzeros = nonzeros + 1

        print("#############################")
        print("Number of rows: " +str(A.Size1()) )
        print("Number of columns: " +str(A.Size2()) )
#        print("Number of entries: " +str(nonzeros) )
        print("Max in Diagonal: " +str(max) )
        print("#############################")
    #######################################################################
    def PlotSparsityScheme(self, A):
        try:
            import Gnuplot, Gnuplot.PlotItems, Gnuplot.funcutils
        except ImportError:
            # kludge in case Gnuplot hasn't been installed as a module yet:
            import __init__
            Gnuplot = __init__
            import PlotItems
            Gnuplot.PlotItems = PlotItems
            import funcutils
            Gnuplot.funcutils = funcutils
        print("gnuplot-python imported")
        g = Gnuplot.Gnuplot(debug=1)
        g.clear()
        #g.plot(Gnuplot.Func('sin(x)'))
        #self.wait('hit Return to continue')
        file = open("matrix.dat",'w')
        for i in range(0, A.Size1()):
            for j in range(0, A.Size2()):
                tmp = A[(i,j)]
                if( (tmp > 1.0e-9) or (tmp < -1.0e-9) ):
                   #file.write( str(tmp) +"\t" )
                   file.write( "1.0 " )
                else:
                   file.write("0.0 ")
            file.write("\n")
        file.close()
        g("set term postscript")
        g("set size square")
        g("set output 'matrix.ps'")
        g("set zrange [0.5:1.5]")
        g("set pm3d map")
        g("splot 'matrix.dat' matrix with dots")


    def wait(self,str=None, prompt='Press return to show results...\n'):
        if str is not None:
            print str
        raw_input(prompt)

