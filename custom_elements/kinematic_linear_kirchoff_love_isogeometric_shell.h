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
#include "includes/c2c_variables.h"
#include "custom_elements/kinematic_linear_kirchoff_love_isogeometric_shell.h"
#include "structural_application/structural_application.h"
#include "isogeometric_structural_application/isogeometric_structural_application.h"

//#define DEBUG_DKGQ

namespace Kratos
{

//************************************************************************************
//***** Constructor and Destructor ***************************************************
//************************************************************************************
KinematicLinearKirchoffLoveIsogeometricShell::KinematicLinearKirchoffLoveIsogeometricShell()
{
}

KinematicLinearKirchoffLoveIsogeometricShell::KinematicLinearKirchoffLoveIsogeometricShell( IndexType NewId, 
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

KinematicLinearKirchoffLoveIsogeometricShell::KinematicLinearKirchoffLoveIsogeometricShell( IndexType NewId, 
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
KinematicLinearKirchoffLoveIsogeometricShell::~KinematicLinearKirchoffLoveIsogeometricShell()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer KinematicLinearKirchoffLoveIsogeometricShell::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new KinematicLinearKirchoffLoveIsogeometricShell(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer KinematicLinearKirchoffLoveIsogeometricShell::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new KinematicLinearKirchoffLoveIsogeometricShell(NewId, pGeom, pProperties));
}

void KinematicLinearKirchoffLoveIsogeometricShell::Initialize()
{
    KRATOS_TRY

    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    const GeometryType::IntegrationPointsArrayType& integration_points
            = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
    }

    InitializeMaterial();

    KRATOS_CATCH("")
}

void KinematicLinearKirchoffLoveIsogeometricShell::InitializeMaterial()
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
        mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
        mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

         //verify that the constitutive law has the correct dimension
//                if ( dimension == 2 )
//                {
//                    if ( this->GetProperties().Has( THICKNESS ) == false )
//                        KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

//                    if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 3 )
//                        KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
//                }
//                else if(dimension == 3)
//                {
//                    if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 6 )
//                        KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
//                }

        //check constitutive law
        mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), *(ProcessInfo*)0 );
//                if( mConstitutiveLawVector[i]->IsIncremental() )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//                if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//                if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//                    KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
    }

    KRATOS_CATCH( "" )
}

void KinematicLinearKirchoffLoveIsogeometricShell::ResetConstitutiveLaw()
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************ 
//************************************************************************************
/**
 * calculates only the RHS vector
 */
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateRightHandSide( VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    bool MaterialUpdateFlag = false;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector, 
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, 
                  CalculateResidualVectorFlag,
                  MaterialUpdateFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
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
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                  VectorType& rRightHandSideVector, 
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag,
                                  bool MaterialUpdateFlag)
{
    KRATOS_TRY

    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    mpIsogeometricGeometry->Initialize(mThisIntegrationMethod);
    #endif

    // get integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = 
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

     unsigned int strain_size = dim*(dim+1)/2;
     unsigned int mat_size = number_of_ctrl_points*dim;

     if(CalculateStiffnessMatrixFlag == true)
     {
         if(rLeftHandSideMatrix.size1() != mat_size)
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

    // loop over integration points
    for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        // declare variables
        Vector MembraneStressVector = ZeroVector(strain_size);
        Vector BendingStressVector = ZeroVector(strain_size);

        Matrix DN_De;
        Vector Ncontainer;

        // material parameters
        double E = GetProperties()[YOUNG_MODULUS];
        double NU = GetProperties()[POISSON_RATIO];
        double lambda = E*NU /(1.0 +NU)/(1.0 - 2.0*NU);
        double mu = 0.5*E/(1.0 + NU);
        double kappa = GetProperties()[KAPPA];
	double t = GetProperties()[THICKNESS];

        // local gradient values at an integration point
        DN_De = mpIsogeometricGeometry->ShapeFunctionsLocalGradients( DN_De, integration_points[PointNumber]);

        // shape function values at an integration point
        Ncontainer = mpIsogeometricGeometry->ShapeFunctionsValues( Ncontainer , integration_points[PointNumber]);

        // calculate covariant base vector
        array_1d<double,3> G1; G1[2] = 0.0;
        array_1d<double,3> G2; G2[2] = 0.0;

        for(unsigned int i = 0; i < dim; ++i)
        {
            G1(i) = J0[PointNumber](i, 0);
        }

        for(unsigned int i = 0; i < dim; ++i)
        {
            G2(i) = J0[PointNumber](i, 1);
        }

        /////////Bm
        Matrix Bm(strain_size, mat_size);
        CalculateMembraneBOperator(Bm, G1, G2, DN_De);

        //////// get local second local gradients
        /* D2N_De2[i](D2N_Dxi2, D2N_Dxieta; D2N_Detaxi, D2N_Deta2]*/
        ShapeFunctionsSecondDerivativesType D2N_De2;
        D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);

        //////// Bb
        Matrix Bb(strain_size, mat_size);
        CalculateBendingBOperator(Bb, G1, G2, DN_De, D2N_De2);

        //////// Hookean Matrix D
        Matrix D;
        CalculateHookeanMatrix(D, G1, G2);

        /////// membrane strains and curvature changes
        Vector MembraneStrainVector(strain_size);
        Vector CurvatureChangeVector(strain_size);

        CalculateStrain(MembraneStrainVector, Bm, CurrentDisp);
        CalculateStrain(CurvatureChangeVector, Bb, CurrentDisp);

	noalias(MembrainStressVector) = prod(D, t*MembraneStrainVector);
	noalias(BendingStressVector) = prod(D, pow(t,3)/12*CurvatureChangeVector);



        //calculating weights for integration on the reference configuration
        double Weight = integration_points[PointNumber].Weight();

        if (dim == 2)
            Weight *= t;

        if (CalculateStiffnessMatrixFlag == true) //compute the contribution to LHS
        {
            noalias(rLeftHandSideMatrix) += prod(trans(Bm), (t*Weight * mDetJ0[PointNumber]) * Matrix(prod(D, Bm)));
            noalias(rLeftHandSideMatrix) += prod(trans(Bb), (pow(t,3)/12*Weight * mDetJ0[PointNumber]) * Matrix(prod(D, Bb)));
        } // compute the contribution to LHS

        // compute the contribution to RHS 
        if(CalculateResidualVectorFlag == true)
        {
            //contribution to external forces
            if(GetProperties().Has(BODY_FORCE))
                CalculateAndAddExtForceContribution(rRightHandSideVector, Ncontainer, Weight, mDetJ0[PointNumber] );  


            if(GetProperties().Has(DENSITY) && GetProperties().Has(GRAVITY))
                AddBodyForcesToRHS(rRightHandSideVector, Ncontainer, Weight, mDetJ0[PointNumber]);

            // contribution of internal forces 
            AddInternalForcesToRHS( rRightHandSideVector, Bm, MembraneStressVector, Weight, mDetJ0[PointNumber] );
            AddInternalForcesToRHS( rRightHandSideVector, Bb, BendingStressVector, Weight, mDetJ0[PointNumber] );
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
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateDampingMatrix( MatrixType& rDampMatrix,
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
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateMassMatrix( MatrixType& rMassMatrix,
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
KinematicLinearKirchoffLoveIsogeometricShell::IntegrationMethod KinematicLinearKirchoffLoveIsogeometricShell::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector
*/
void KinematicLinearKirchoffLoveIsogeometricShell::EquationIdVector( EquationIdVectorType& rResult, 
                                      ProcessInfo& CurrentProcessInfo)
{
    unsigned int dofs_per_node = 6;
    unsigned int mat_size = GetGeometry().size() * dofs_per_node;

    if ( rResult.size() != mat_size )
        rResult.resize( mat_size, false );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
    {
        int index = i * dofs_per_node;
        rResult[index    ] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list
 */
void KinematicLinearKirchoffLoveIsogeometricShell::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateStrain(Vector& StrainVector, const Matrix& B, const Matrix& Displacements)
{
     unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
     unsigned int strain_size = dim * (dim + 1) / 2;

     noalias(StrainVector) = ZeroVector(strain_size);
     for(unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
         for(unsigned int j = 0; j < strain_size; ++j)
             for(unsigned int k = 0; k < dim; ++k)
                 StrainVector[j] += B(j, i*dim + k) * Displacements(i, k);
}

void KinematicLinearKirchoffLoveIsogeometricShell::CalculateAndAddExtForceContribution(
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

void KinematicLinearKirchoffLoveIsogeometricShell::AddBodyForcesToRHS(
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

void KinematicLinearKirchoffLoveIsogeometricShell::AddInternalForcesToRHS(
    Vector& rRightHandSideVector,
    const Matrix& B,
    const Vector& StressVector,
    const double& Weight,
    const double& DetJ
)
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim + 1) / 2;
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

void KinematicLinearKirchoffLoveIsogeometricShell::CalculateBendingBOperator(
    Matrix& Bb,
    const array_1d<double,3>& G1,
    const array_1d<double,3>& G2,
    const Matrix& DN_De,
    const ShapeFunctionsSecondDerivativesType& D2N_De2 )
{
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int mat_size = number_of_ctrl_points*dim;
    unsigned int strain_size = dim * (dim + 1) / 2;

    if( Bb.size1() != 3|| Bb.size2() != mat_size)
        Bb.resize(3, mat_size);

    noalias(Bb) = ZeroMatrix(strain_size, mat_size);
    
     // bending Bb operator
     //////// covariant of unit normal base vector
     array_1d<double,3> G3; 
     MathUtils<double>::CrossProduct(G3,G1,G2);
     G3 /= norm_2(G3);

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

         G22(0) += D2N_De2[j](1,1)*GetGeometry()[j].X0();
         G22(1) += D2N_De2[j](1,1)*GetGeometry()[j].Y0();

         G12(0) += D2N_De2[j](0,1)*GetGeometry()[j].X0();
         G12(1) += D2N_De2[j](0,1)*GetGeometry()[j].Y0();
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
             Bb(0, i*dim+j) = -D2N_De2[i](0,0)*G3(j) + ( 1/norm_2(G3) )*( DN_De(i,0)*G11xG2(j) + DN_De(i,1)*G1xG11(j) ) 
                                + inner_prod(G3, G11)/norm_2(G3)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) );
             Bb(1, i*dim+j) = -D2N_De2[i](1,1)*G3(j) + ( 1/norm_2(G3) )*( DN_De(i,0)*G22xG2(j) + DN_De(i,1)*G1xG22(j) ) 
                                + inner_prod(G3, G22)/norm_2(G3)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) );
             Bb(2, i*dim+j) = -2*D2N_De2[i](0,1)*G3(j) + 2*( 1/norm_2(G3) )*( DN_De(i,0)*G12xG2(j) + DN_De(i,1)*G1xG12(j) ) 
                                + 2*inner_prod(G3, G12)/norm_2(G3)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG1(j) );
         }
     }
}// Bending B operator

void KinematicLinearKirchoffLoveIsogeometricShell::CalculateMembraneBOperator(
    Matrix& Bm,
    const array_1d<double,3>& G1,
    const array_1d<double,3>& G2,
    const Matrix& DN_De)
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int strain_size = dim*(dim+1)/2;
    unsigned int mat_size = dim*number_of_ctrl_points;

    if( Bm.size1() != 3 || Bm.size2() != mat_size)
        Bm.resize(3, mat_size);

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

void KinematicLinearKirchoffLoveIsogeometricShell::CalculateHookeanMatrix(Matrix& D, const Vector& G1, const Vector& G2)
{
    /////////////////////////////////////////////////////////////
    ////////////// elastic matrix D///////////////////////////////
    // compute curvature coefficient
    double G_11 = inner_prod(G1, G1); double G_12 = inner_prod(G1, G2);
    double G_21 = inner_prod(G2, G1); double G_22 = inner_prod(G2, G2);
    Matrix G_ab(2,2);
    G_ab(0, 0) = G_11; G_ab(0, 1) = G_12;
    G_ab(1, 0) = G_21; G_ab(1, 1) = G_22;

    double temp_ab;
    Matrix InvG_ab(2,2);
    MathUtils<double>::InvertMatrix(G_ab, InvG_ab, temp_ab);

    // material parameters
    double E = GetProperties()[YOUNG_MODULUS];
    double NU = GetProperties()[POISSON_RATIO];
    double lambda = E*NU /(1.0 +NU)/(1.0 - 2.0*NU);
    double mu = 0.5*E/(1.0 + NU);

    double aux = E/(1-pow(NU,2) );

    if (D.size1() != 3||D.size2() != 3)
        D.resize(3,3);

    
    D(0,0) = aux*pow(InvG_ab(0,0), 2);
    D(0,1) = aux*( NU*InvG_ab(0,0)*InvG_ab(1,1) + (1-NU)*pow(InvG_ab(0,1), 2) ); 
    D(0,2) = aux*InvG_ab(0,0)*InvG_ab(0,1);
    D(1,0) = D(0,1);
    D(1,1) = aux*pow(InvG_ab(1,1), 2); 
    D(1,2) = aux*InvG_ab(1,1)*InvG_ab(0,1);
    D(2,0) = D(0,2);
    D(2,1) = D(1,2);
    D(2,2) = aux*0.5*( (1-NU)*InvG_ab(0,0)*InvG_ab(1,1) + (1+NU)*pow(InvG_ab(0,1), 2) ) ;
    

    /* D(0, 0) = aux*pow(G_11, 2);
    D(0, 1) = aux*( NU*G_11*G_22 + (1-NU)*pow(G_12, 2) ); 
    D(0, 2) = aux*G_11*G_12;
    D(1, 1) = aux*pow(G_22, 2); 
    D(1, 2) = aux*G_22*G_12;
    D(2, 2) = aux*0.5*( (1-NU)*G_11*G_22 + (1+NU)*pow(G_12, 2) ) ;

    D(1, 0) = D(0, 1);
    D(2, 0) = D(0, 2);
    D(2, 1) = D(1, 2);
    */

    //KRATOS_WATCH(aux)
    //KRATOS_WATCH(D)
    //char key; cin >> key;         
    ////////////// elastic matrix D///////////////////////////////
    //////////////////////////////////////////////////////////////
}






//************************************************************************************
//************************************************************************************
void KinematicLinearKirchoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
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
void KinematicLinearKirchoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
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

    if( rVariable == INTEGRATION_POINT_GLOBAL )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            rValues[point] = GetGeometry().GlobalCoordinates(rValues[point], integration_points[point]);
        }

        return;
    }

    if( rVariable == INTEGRATION_POINT_LOCAL )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            noalias(rValues[point]) = integration_points[point];
        }

        return;
    }
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchoffLoveIsogeometricShell::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
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
    else if( rVariable == JACOBIAN_0 )
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
    }
    else if( rVariable == INTEGRATION_WEIGHT )
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
    }
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
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
        std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void KinematicLinearKirchoffLoveIsogeometricShell::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    this->GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

} // Namespace Kratos
