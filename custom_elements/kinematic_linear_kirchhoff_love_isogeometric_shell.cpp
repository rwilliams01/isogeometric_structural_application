//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Apr 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell.h"
#include "isogeometric_structural_application/isogeometric_structural_application.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"
#include "structural_application/structural_application.h"

namespace Kratos
{

//************************************************************************************
//***** Constructor and Destructor ***************************************************
//************************************************************************************
KinematicLinearKirchhoffLoveIsogeometricShell::KinematicLinearKirchhoffLoveIsogeometricShell(IndexType NewId,
        GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
{
    mIsInitialized = false;
//    mpIsogeometricGeometry =
//        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
    /*
    Important remarks:
        + GetGeometry() and (*mpIsogeometricGeometry) refer to the same instance of IsogeometricGeometryType in the memory. However, GetGeometry() only provides access to the functions wrapped by Geometry interface, whereby (*mpIsogeometricGeometry) provides access to functions exclusive to IsogeometricGeometryType. It is ok to replace every instances of GetGeometry() by (*mpIsogeometricGeometry) but to keep the code looks compatible (especiall for comparison with old code), GetGeometry() can still be kept, but take note to the wrapped functions.
     */
}

KinematicLinearKirchhoffLoveIsogeometricShell::KinematicLinearKirchhoffLoveIsogeometricShell(IndexType NewId,
        GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
{
    mIsInitialized = false;
//    mpIsogeometricGeometry =
//        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

/**
 * Destructor. Never to be called manually
 */
KinematicLinearKirchhoffLoveIsogeometricShell::~KinematicLinearKirchhoffLoveIsogeometricShell()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer KinematicLinearKirchhoffLoveIsogeometricShell::Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
            new KinematicLinearKirchhoffLoveIsogeometricShell(NewId,
                    GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer KinematicLinearKirchhoffLoveIsogeometricShell::Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
            new KinematicLinearKirchhoffLoveIsogeometricShell(NewId,
                    pGeom, pProperties));
}

void KinematicLinearKirchhoffLoveIsogeometricShell::Initialize()
{
    KRATOS_TRY        //EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

        if (mIsInitialized)
        {
            //dimension of the problem
            unsigned int dim = 3;

            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize(mpIsogeometricGeometry->size(), dim, false);

            for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
                for (unsigned int i = 0; i < dim; ++i) // hbui edited
                    mInitialDisp(node, i) =
                        (*mpIsogeometricGeometry)[node].GetSolutionStepValue(
                            DISPLACEMENT)[i];

            return;
        }

        ///////////////////////////////////////////////////////////////
        // One time initialisation
        ///////////////////////////////////////////////////////////////

        ////////////////////Initialize geometry_data/////////////////////////////
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_1))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_2))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_3))
//        KRATOS_WATCH(GetValue(NURBS_WEIGHT))
//        KRATOS_WATCH(GetValue(EXTRACTION_OPERATOR))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_1))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_2))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_3))

        // try to read the extraction operator from the elemental data
        Matrix ExtractionOperator;
        bool manual_initilization = false;
        if( this->Has( EXTRACTION_OPERATOR ) )
        {
            ExtractionOperator = this->GetValue( EXTRACTION_OPERATOR );
            manual_initilization = true;
        }
        else if( this->Has( EXTRACTION_OPERATOR_MCSR ) )
        {
            Matrix Temp = this->GetValue( EXTRACTION_OPERATOR_MCSR );

            // make a simple check
            if(Temp.size1() != 2)
                KRATOS_THROW_ERROR(std::logic_error, "Invalid MCSR matrix for extraction operator found at element", this->Id())

            // choose the best storage scheme based ratio between number of nonzeros and the full size of the matrix
            unsigned int size_ex_n = (unsigned int)(Temp(0, 0) - 1);
            unsigned int size_ex_nz = Temp.size2() - 1;
            if( ( (double)(size_ex_nz) ) / (size_ex_n * size_ex_n) < 0.2 )
                ExtractionOperator = IsogeometricMathUtils::MCSR2CSR(Temp);
            else
                ExtractionOperator = IsogeometricMathUtils::MCSR2MAT(Temp);

            manual_initilization = true;
        }
        else if( this->Has( EXTRACTION_OPERATOR_CSR_ROWPTR )
             && this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
             && this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
        {
            Vector rowPtr = this->GetValue( EXTRACTION_OPERATOR_CSR_ROWPTR ); // must be 0-base
            Vector colInd = this->GetValue( EXTRACTION_OPERATOR_CSR_COLIND ); // must be 0-base
            Vector values = this->GetValue( EXTRACTION_OPERATOR_CSR_VALUES );
//            int m = mpIsogeometricGeometry->size();
//            int n;
//            unsigned int dim = mpIsogeometricGeometry->NURBS_WorkingSpaceDimension();
//            if(dim == 2)
//                n = (1 + this->GetValue(NURBS_DEGREE_1)) * (1 + this->GetValue(NURBS_DEGREE_2));
//            else if(dim == 3)
//                n = (1 + this->GetValue(NURBS_DEGREE_1)) * (1 + this->GetValue(NURBS_DEGREE_2)) * (1 + this->GetValue(NURBS_DEGREE_3));
//            KRATOS_WATCH(m)
//            KRATOS_WATCH(n)
//            KRATOS_WATCH(rowPtr)
//            KRATOS_WATCH(colInd)
//            KRATOS_WATCH(values)
            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
//            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(m, n, rowPtr, colInd, values);

            manual_initilization = true;
        }
//        else
//            KRATOS_THROW_ERROR(std::logic_error, "The extraction operator was not given for element", Id())

        // initialize the geometry
        if(manual_initilization)
            mpIsogeometricGeometry->AssignGeometryData(
                this->GetValue(NURBS_KNOTS_1),
                this->GetValue(NURBS_KNOTS_2),
                this->GetValue(NURBS_KNOTS_3),
                this->GetValue(NURBS_WEIGHT),
                ExtractionOperator,
                this->GetValue(NURBS_DEGREE_1),
                this->GetValue(NURBS_DEGREE_2),
                this->GetValue(NURBS_DEGREE_3),
                2 // only need to compute 2 integration rules
            );

        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

        InitializeJacobian();
        ////////////////////End Initialize geometry_data/////////////////////////////

        //Set Up Initial displacement for StressFreeActivation of Elements
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        mInitialDisp.resize(mpIsogeometricGeometry->size(), dim, false);

        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
            for (unsigned int i = 0; i < dim; ++i)
                mInitialDisp(node, i) =
                    (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT)[i];

        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // law was at each integration point
        const GeometryType::IntegrationPointsArrayType& integration_points =
                mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);
        if (mConstitutiveLawVector.size() != integration_points.size())
        {
            mConstitutiveLawVector.resize(integration_points.size());
        }

        InitializeMaterial();

        mIsInitialized = true;

    KRATOS_CATCH( "" )
}


void KinematicLinearKirchhoffLoveIsogeometricShell::InitializeJacobian()
{
    unsigned int dim = 3;

    //number of integration points used, mThisIntegrationMethod refers to the
    //integration method defined in the constructor
    const GeometryType::IntegrationPointsArrayType& integration_points =
            mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    //initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference
    // configuration
    GeometryType::JacobiansType J0(integration_points.size());

    mTotalDomainInitialSize = 0.00;

    //calculating the Jacobian
    J0 = mpIsogeometricGeometry->Jacobian(J0, mThisIntegrationMethod);

    //calculating the inverse Jacobian
    double DetJ0;
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();
        //calculating and storing inverse of the jacobian and the parameters needed
        Matrix JtJ = prod(trans(J0[PointNumber]), J0[PointNumber]);
        DetJ0 = sqrt(MathUtils<double>::Det(JtJ));
        //calculating the total area
        mTotalDomainInitialSize += DetJ0 * IntegrationWeight;
    }
}



/**
 * Initialization of the Material law at each integration point
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::InitializeMaterial()
{
    KRATOS_TRY

        //calculating shape functions values
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
//                mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
//                mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
//            std::cout << "consitutive law vector " << i << " received clone" << std::endl;
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(),
                (*mpIsogeometricGeometry),
                row(Ncontainer, i)
            );
//            std::cout << "consitutive law vector " << i << " is initialized" << std::endl;
        }

    KRATOS_CATCH( "" )


}



void KinematicLinearKirchhoffLoveIsogeometricShell::ResetConstitutiveLaw()
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    KRATOS_CATCH( "" )
}


void KinematicLinearKirchhoffLoveIsogeometricShell::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    // reset all resistant forces at node
    for ( unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i )
    {
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue( REACTION_X ) = 0.0;
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue( REACTION_Z ) = 0.0;

    }
}
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp = Matrix();

    bool need_calculate_stiffness = false;
    for ( unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node )
    {
        if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_X)
        || (*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Y)
        || (*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Z))
        {
            need_calculate_stiffness = true;
            break;
        }
    }

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, need_calculate_stiffness, true, false);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    bool MaterialUpdateFlag = false;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, MaterialUpdateFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * This function calculates all system contributions
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag,
                                  bool MaterialUpdateFlag)
{
    KRATOS_TRY

    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = 3;

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    mpIsogeometricGeometry->Initialize(mThisIntegrationMethod);
    #endif

    // get integration points
    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

     unsigned int strain_size = 3;
     unsigned int mat_size = number_of_ctrl_points*dim;

     if(CalculateStiffnessMatrixFlag == true)
     {
         if(rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
             rLeftHandSideMatrix.resize(mat_size, mat_size, false);

         noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
     }

     if(CalculateResidualVectorFlag == true)
     {
         if(rRightHandSideVector.size() != mat_size)
             rRightHandSideVector.resize(mat_size, false);

         noalias(rRightHandSideVector) = ZeroVector(mat_size);
     }

     // Current displacements
     Matrix CurrentDisp(number_of_ctrl_points, dim);
     for(unsigned int node =0; node < mpIsogeometricGeometry->size() ;++node)
         noalias(row(CurrentDisp, node)) = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

    // calculate the Jacobian
    GeometryType::JacobiansType J0(integration_points.size());
    J0 = mpIsogeometricGeometry->Jacobian0(J0, mThisIntegrationMethod);

    // material parameters
    const double& E = GetProperties()[YOUNG_MODULUS];
    const double& NU = GetProperties()[POISSON_RATIO];
    const double& t = GetProperties()[THICKNESS];

    // loop over integration points
    for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        // Matrix JtJ = prod(trans(J0[PointNumber]), J0[PointNumber]);
        // DetJ0 = sqrt(MathUtils<double>::Det(JtJ));

        // declare variables
        Vector MembraneStressVector(strain_size);
        Vector BendingStressVector(strain_size);

        Matrix DN_De;
        Vector Ncontainer;

        // local gradient values at an integration point
        DN_De = mpIsogeometricGeometry->ShapeFunctionsLocalGradients( DN_De, integration_points[PointNumber]);

        // shape function values at an integration point
        Ncontainer = mpIsogeometricGeometry->ShapeFunctionsValues( Ncontainer , integration_points[PointNumber]);

        // calculate covariant base vector
        array_1d<double,3> G1;
        array_1d<double,3> G2;
        array_1d<double,3> G3;

        for(unsigned int i = 0; i < dim; ++i)
        {
            G1(i) = J0[PointNumber](i, 0);
        }

        for(unsigned int i = 0; i < dim; ++i)
        {
            G2(i) = J0[PointNumber](i, 1);
        }

        noalias(G3) = MathUtils<double>::CrossProduct(G1, G2);
        double dA = norm_2(G3); // note that this is equivalent to sqrt(det(JtJ)) above

        /////////Bm
        Matrix Bm(strain_size, mat_size);
        CalculateMembraneBOperator(Bm, G1, G2, DN_De);

        //////// get local second local gradients
        /* D2N_De2[i](D2N_Dxi2, D2N_Dxieta; D2N_Detaxi, D2N_Deta2]*/
        ShapeFunctionsSecondDerivativesType D2N_De2;
        D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);

        //std::cout << "at point " << PointNumber << ":" << std::endl;
        //KRATOS_WATCH(integration_points[PointNumber])
        //KRATOS_WATCH(Ncontainer)
        //KRATOS_WATCH(DN_De)
        //KRATOS_WATCH(D2N_De2)
        //std::cout << "----------------------------" << std::endl;

        //////// Bb
        Matrix Bb(strain_size, mat_size);
        CalculateBendingBOperator(Bb, G1, G2, DN_De, D2N_De2);

        //////// Hookean Matrix D
        Matrix D(strain_size, strain_size);
        CalculateHookeanMatrix(D, G1, G2, E, NU);

        /////// membrane strains and curvature changes
        Vector MembraneStrainVector(strain_size);
        Vector CurvatureChangeVector(strain_size);

        CalculateStrain(MembraneStrainVector, Bm, CurrentDisp);
        CalculateStrain(CurvatureChangeVector, Bb, CurrentDisp);

        noalias(MembraneStressVector) = t * prod(D, MembraneStrainVector);
        noalias(BendingStressVector) = (pow(t,3)/12) * prod(D, CurvatureChangeVector);

        //calculating weights for integration on the reference configuration
        double Weight = integration_points[PointNumber].Weight();

        // if (dim == 2) Weight *= t;

        if (CalculateStiffnessMatrixFlag == true) //compute the contribution to LHS
        {
            noalias(rLeftHandSideMatrix) += t*Weight * dA * prod(trans(Bm), Matrix(prod(D, Bm)));
            noalias(rLeftHandSideMatrix) += pow(t,3)/12*Weight * dA * prod(trans(Bb), Matrix(prod(D, Bb)));
        } // compute the contribution to LHS

        // compute the contribution to RHS
        if(CalculateResidualVectorFlag == true)
        {
            //contribution to external forces
            if(GetProperties().Has(BODY_FORCE))
                CalculateAndAddExtForceContribution(rRightHandSideVector, Ncontainer, Weight, dA );

            if(GetProperties().Has(DENSITY) && GetProperties().Has(GRAVITY))
                AddBodyForcesToRHS(rRightHandSideVector, Ncontainer, Weight, dA);

            // contribution of internal forces
            AddInternalForcesToRHS( rRightHandSideVector, Bm, MembraneStressVector, Weight, dA );
            AddInternalForcesToRHS( rRightHandSideVector, Bb, BendingStressVector, Weight, dA );
        } // compute the contribution to RHS
    }//loop over integration points

    if(CalculateResidualVectorFlag == true)
    {
        // modify the right hand side to account for prescribed displacement
        // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
        // // However, I have to temporarily disable it to keep the consistency.
        for ( unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node )
        {
            if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_X))
            {
                double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                for( unsigned int i = 0; i < mat_size; ++i )
                    rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
            }
            if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Y))
            {
                double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                for( unsigned int i = 0; i < mat_size; ++i )
                    rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
            }
            if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Z))
            {
                double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                for( unsigned int i = 0; i < mat_size; ++i )
                    rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
            }
        }
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    //clean the internal data of the geometry
    mpIsogeometricGeometry->Clean();
    #endif

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
 * This function calculates the damping matrix
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateDampingMatrix( MatrixType& rDampMatrix,
                                    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
//TODO
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
 * This function calculates the mass matrix
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateMassMatrix( MatrixType& rMassMatrix,
                                    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
//TODO
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
 * returns the used integration method
 */
KinematicLinearKirchhoffLoveIsogeometricShell::IntegrationMethod KinematicLinearKirchhoffLoveIsogeometricShell::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector
*/
void KinematicLinearKirchhoffLoveIsogeometricShell::EquationIdVector( EquationIdVectorType& rResult,
                                     ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = number_of_nodes * 3;

    if ( rResult.size() != dim )
        rResult.resize( dim );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list
 */
void KinematicLinearKirchhoffLoveIsogeometricShell::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//////////////////////////////////////////////////////////////////////
////////////////// private subroutines ///////////////////////////////
//////////////////////////////////////////////////////////////////////
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateStrain(Vector& StrainVector, const Matrix& B, const Matrix& Displacements)
{
     unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
     unsigned int strain_size = 3;

     noalias(StrainVector) = ZeroVector(strain_size);
     for(unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
         for(unsigned int j = 0; j < strain_size; ++j)
             for(unsigned int k = 0; k < dim; ++k)
                 StrainVector[j] += B(j, i*dim + k) * Displacements(i, k);
}

void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateAndAddExtForceContribution(
    VectorType& rRightHandSideVector, const Vector& N, const double& Weight, const double& DetJ
)
{
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

    const Vector& BodyForce = GetProperties()[BODY_FORCE];

    for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            rRightHandSideVector(i*dim + j) += Weight * DetJ * N(i) * BodyForce(j);
}

void KinematicLinearKirchhoffLoveIsogeometricShell::AddBodyForcesToRHS(
    VectorType& rRightHandSideVector, const Vector& N, const double& Weight, const double& DetJ
)
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

    const Vector& gravity = GetProperties()[GRAVITY];
    const double& density = GetProperties()[DENSITY];

    for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
        {
            rRightHandSideVector(i * dim + j) += N(i) * density * gravity(j) * DetJ * Weight;
        }
    }
}

void KinematicLinearKirchhoffLoveIsogeometricShell::AddInternalForcesToRHS(
    Vector& rRightHandSideVector,
    const Matrix& B,
    const Vector& StressVector,
    const double& Weight,
    const double& DetJ
)
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int strain_size = 3;
    Vector InternalForces(dim);

    for(unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
    {
       for( unsigned int j = 0; j < dim; ++j)
       {
              InternalForces(j) = 0.0;
              for(unsigned int k = 0; k < strain_size; ++k)
                InternalForces(j) += B(k, dim * i + j) * StressVector(k) * Weight * DetJ;

              rRightHandSideVector(dim * i + j) -= InternalForces(j);
       }
     }
}

void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateBendingBOperator(
    Matrix& Bb,
    const array_1d<double,3>& G1,
    const array_1d<double,3>& G2,
    const Matrix& DN_De,
    const ShapeFunctionsSecondDerivativesType& D2N_De2 )
{
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = 3;
    unsigned int mat_size = number_of_ctrl_points*dim;
    unsigned int strain_size = 3;

    if( Bb.size1() != strain_size || Bb.size2() != mat_size)
        Bb.resize(strain_size, mat_size, false);

    noalias(Bb) = ZeroMatrix(strain_size, mat_size);

     // bending Bb operator
     //////// covariant of unit normal base vector
     array_1d<double,3> G3;
     array_1d<double,3> G1xG2;
     MathUtils<double>::CrossProduct(G1xG2,G1,G2);
     G3 = G1xG2/norm_2(G1xG2);
     //////// compute G1,1 ; G1,2 ; G2,2
     /*array_1d<double,3> G11 ;
     array_1d<double,3> G12 ;
     array_1d<double,3> G22 ;
     G11[2] = G12[2] = G22[2] = 0;
     */
     Vector G11 = ZeroVector(3);
     Vector G12 = ZeroVector(3);
     Vector G22 = ZeroVector(3);

     for(unsigned int j=0; j< number_of_ctrl_points; j++)
     {
         G11(0) += D2N_De2[j](0,0)*GetGeometry()[j].X0();
         G11(1) += D2N_De2[j](0,0)*GetGeometry()[j].Y0();
         G11(2) += D2N_De2[j](0,0)*GetGeometry()[j].Z0();

         G22(0) += D2N_De2[j](1,1)*GetGeometry()[j].X0();
         G22(1) += D2N_De2[j](1,1)*GetGeometry()[j].Y0();
         G22(2) += D2N_De2[j](1,1)*GetGeometry()[j].Z0();

         G12(0) += D2N_De2[j](0,1)*GetGeometry()[j].X0();
         G12(1) += D2N_De2[j](0,1)*GetGeometry()[j].Y0();
         G12(2) += D2N_De2[j](0,1)*GetGeometry()[j].Z0();
     }

     //////// cross product of two vectors G1,1 x G2 ; G1 x G1,1; ...
     array_1d<double,3> G11xG2, G22xG2, G12xG2;
     array_1d<double,3> G1xG11, G1xG22, G1xG12 ;
     array_1d<double,3> G2xG3, G3xG1 ;

     MathUtils<double>::CrossProduct(G11xG2 ,G11 ,G2);
     MathUtils<double>::CrossProduct(G22xG2 ,G22 ,G2);
     MathUtils<double>::CrossProduct(G12xG2 ,G12 ,G2);
     MathUtils<double>::CrossProduct(G1xG11 ,G1 ,G11);
     MathUtils<double>::CrossProduct(G1xG22 ,G1 ,G22);
     MathUtils<double>::CrossProduct(G1xG12 ,G1 ,G12);
     MathUtils<double>::CrossProduct(G2xG3 ,G2 ,G3);
     MathUtils<double>::CrossProduct(G3xG1 ,G3 ,G1);

     /////// compute Bb
     Bb(strain_size, mat_size);
     for(unsigned int i=0; i< number_of_ctrl_points; i++)
     {
         for(unsigned int j=0; j< dim; j++)
         {
             Bb(0, i*dim+j) = -D2N_De2[i](0,0)*G3(j) + ( 1/norm_2(G1xG2) )*( DN_De(i,0)*G11xG2(j) + DN_De(i,1)*G1xG11(j)
                                + inner_prod(G3, G11)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) ) );
             Bb(1, i*dim+j) = -D2N_De2[i](1,1)*G3(j) + ( 1/norm_2(G1xG2) )*( DN_De(i,0)*G22xG2(j) + DN_De(i,1)*G1xG22(j)
                                + inner_prod(G3, G22)*(  DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) ) );
             Bb(2, i*dim+j) = -2*D2N_De2[i](0,1)*G3(j) + 2*( 1/norm_2(G1xG2) )*( DN_De(i,0)*G12xG2(j) + DN_De(i,1)*G1xG12(j)
                                + 2*inner_prod(G3, G12)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) ) );
         }
     }
}// Bending B operator

void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateMembraneBOperator(
    Matrix& Bm,
    const array_1d<double,3>& G1,
    const array_1d<double,3>& G2,
    const Matrix& DN_De)
{
    unsigned int dim = 3;
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int strain_size = 3;
    unsigned int mat_size = dim*number_of_ctrl_points;


    if( Bm.size1() != strain_size || Bm.size2() != mat_size)
        Bm.resize(strain_size, mat_size, false);

    for(unsigned int i=0 ; i< number_of_ctrl_points; i++)
    {
        for(unsigned int j=0; j< dim; j++)
        {
            Bm(0,i*dim + j) = DN_De(i,0)*G1(j);
            Bm(1,i*dim + j) = DN_De(i,1)*G2(j);
            Bm(2,i*dim + j) = DN_De(i,1)*G1(j) + DN_De(i,0)*G2(j);
        }
    }
}

void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateHookeanMatrix(Matrix& D,
    const Vector& G1, const Vector& G2, const double& E, const double& NU)
{
    /////////////////////////////////////////////////////////////
    ////////////// elastic matrix D///////////////////////////////
    double G_11 = inner_prod(G1, G1); double G_12 = inner_prod(G1, G2);
    double G_21 = inner_prod(G2, G1); double G_22 = inner_prod(G2, G2);
    Matrix G_ab(2,2);
    G_ab(0, 0) = G_11; G_ab(0, 1) = G_12;
    G_ab(1, 0) = G_21; G_ab(1, 1) = G_22;

    double temp_ab;
    Matrix InvG_ab(2,2); // note that this is the contravariant matrix
    MathUtils<double>::InvertMatrix(G_ab, InvG_ab, temp_ab);

    // material parameters
    double lambda = E*NU /(1.0 +NU)/(1.0 - 2.0*NU);
    double mu = 0.5*E/(1.0 + NU);

    double aux = E/(1-pow(NU,2) );

    if (D.size1() != 3 || D.size2() != 3)
        D.resize(3, 3, false);

    D(0,0) = aux*pow(InvG_ab(0,0), 2);
    D(0,1) = aux*( NU*InvG_ab(0,0)*InvG_ab(1,1) + (1-NU)*pow(InvG_ab(0,1), 2) );
    D(0,2) = aux*InvG_ab(0,0)*InvG_ab(0,1);
    D(1,0) = D(0,1);
    D(1,1) = aux*pow(InvG_ab(1,1), 2);
    D(1,2) = aux*InvG_ab(1,1)*InvG_ab(0,1);
    D(2,0) = D(0,2);
    D(2,1) = D(1,2);
    D(2,2) = aux*0.5*( (1-NU)*InvG_ab(0,0)*InvG_ab(1,1) + (1+NU)*pow(InvG_ab(0,1), 2) ) ;
    ////////////// elastic matrix D///////////////////////////////
    //////////////////////////////////////////////////////////////
}






//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
    {
        rValues[Point] = mConstitutiveLawVector[Point]->GetValue(rVariable, rValues[Point]);
    }
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if( rVariable == DISPLACEMENT )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = ZeroVector(3);
            for(std::size_t i = 0; i < GetGeometry().size(); ++i)
            {
                const array_1d<double, 3>& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                noalias(rValues[point]) += Ncontainer(point, i) * displacement;
            }
        }

        return;
    }

    /*if( rVariable == INTEGRATION_POINT_GLOBAL )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
        }

        return;
    }*/

    /*if( rVariable == INTEGRATION_POINT_LOCAL )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = integration_points[point];
        }

        return;
    }*/
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if( rVariable == K0 )
    {
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = GetValue( K0 );
        }
        return;
    }
    else if( rVariable == STRAIN_ENERGY )
    {
        std::vector<Vector> StrainList(rValues.size());
        GetValueOnIntegrationPoints(STRAIN, StrainList, rCurrentProcessInfo);

        std::vector<Vector> StressList(rValues.size());
        GetValueOnIntegrationPoints(STRESSES, StressList, rCurrentProcessInfo);

        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            // calculate strain energy as C = 0.5 * <epsilon, sigma>
            rValues[i] = 0.5 * inner_prod(StrainList[i], StressList[i]);
        }
    }
    /*else if( rVariable == JACOBIAN_0 )
    {
        // initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        // compute the Jacobian
        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            Matrix JtJ = prod(trans(J0[i]), J0[i]); // TODO check this
            rValues[i] = sqrt(MathUtils<double>::Det(JtJ));
        }
    }*/
    /*else if( rVariable == INTEGRATION_WEIGHT )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            rValues[i] = integration_points[i].Weight();
        }
    }
    else if( rVariable == MATERIAL_DENSITY )
    {
        for( unsigned int i = 0; i < rValues.size(); ++i )
        {
            rValues[i] = this->GetValue(MATERIAL_DENSITY);
        }
    }*/
    else
    {
        //reading integration points and local gradients
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
        }
    }
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
        std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchhoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/*int KinematicLinearKirchhoffLoveIsogeometricShell::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        KRATOS_WATCH("At Check")

        unsigned int dimension = this->mpIsogeometricGeometry->WorkingSpaceDimension();

        //verify valid id
        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Invalid element ->", this->Id());
        }

        //verify valid domain size
        #ifndef IGNORE_NEGATIVE_JACOBIAN
        if (mTotalDomainInitialSize < 0)
        {
            std::stringstream ss;
            ss << "error on element -> " << this->Id() << ": ";
            ss << "Domain size can not be less than 0. Please check Jacobian. mTotalDomainInitialSize = " << mTotalDomainInitialSize;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        #endif

        //verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "constitutive law not provided for property ",
                    this->GetProperties().Id());
        }

        //verify that the constitutive law has the correct dimension
        if (dimension == 2)
        {
            if (this->GetProperties().Has(THICKNESS) == false)
                KRATOS_THROW_ERROR(std::logic_error,
                        "THICKNESS not provided for element", this->Id());
        }

        //check constitutive law
        int ok = 0;
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            ok = mConstitutiveLawVector[i]->Check(GetProperties(),
                    (*mpIsogeometricGeometry), rCurrentProcessInfo);
            if (ok != 0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the consitutive law", i)
            }

//            if( mConstitutiveLawVector[i]->IsIncremental() )
//                KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//            if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//                KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//            if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//                KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
        }

        //check Jacobian (just for debugging)
        //check Jacobian  should be detected by Area() ot Volume()
//        #ifdef CHECK_JACOBIAN
//        GeometryType::CoordinatesArrayType P;
//
//        P[0] = 0.0;
//        P[1] = 0.0;
//        P[2] = 0.0;
//
//        double J0 = mpIsogeometricGeometry->DeterminantOfJacobian( P );
//
//        if(J0 < 0.0)
//        {
//            KRATOS_THROW_ERROR(std::logic_error, "Negative Jacobian is detected", __FUNCTION__)
//        }
//        #endif

//        KRATOS_WATCH("Check completed")

        return ok;

    KRATOS_CATCH( "" );

}*/


} // Namespace Kratos
