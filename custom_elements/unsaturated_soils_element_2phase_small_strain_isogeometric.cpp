/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 23 Oct 2014 $
 *   Revision:            $Revision: 1.1 $
 *
 * *********************************************************/

// System includes
#include <ctime>

// External includes
#include "boost/timer.hpp"

// Project includes
#include "includes/define.h"
#include "includes/legacy_structural_app_vars.h"
#include "custom_elements/unsaturated_soils_element_2phase_small_strain_isogeometric.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

extern Variable<double> STABILISATION_FACTOR;
extern Variable<double> SHEAR_MODULUS;

UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric( IndexType NewId,
        GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    mIsInitialized = false;
    mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

//************************************************************************************
//************************************************************************************
UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric( IndexType NewId,
        GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mNodesPressMin = 1;
    mNodesPressMax = GetGeometry().size();
    mNodesDispMin = 1;
    mNodesDispMax = GetGeometry().size();

    mIsInitialized = false;
    mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

Element::Pointer UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric( NewId, GetGeometry().Create( ThisNodes ),
                             pProperties ) );
}

UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::~UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric()
{
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        if (mIsInitialized)
        {
            //dimension of the problem
            unsigned int dim = GetGeometry().WorkingSpaceDimension();

            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize(GetGeometry().size(), dim, false);

            for (unsigned int node = 0; node < GetGeometry().size(); ++node)
                for (unsigned int i = 0; i < dim; ++i) // hbui edited
                    mInitialDisp(node, i) =
                        GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT)[i];

//            for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; ++i )
//            {
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL ) =
//                    ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS ) =
//                    ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL_DT ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS_DT ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_NULL ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_EINS ) = ZeroVector( 3 );
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_OLD ) = ZeroVector( 3 );
//
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL ) =
//                    ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS ) =
//                    ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_DT ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_ACCELERATION ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) = 0;
//                ( GetGeometry()[i] ).GetSolutionStepValue( REACTION_WATER_PRESSURE ) = 0;
//            }

            return;
        }

        ///////////////////////////////////////////////////////////////
        // One time initialisation
        ///////////////////////////////////////////////////////////////

        ////////////////////Initialize geometry_data/////////////////////////////
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_1))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_2))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_3))
//        KRATOS_WATCH(GetValue(NURBS_WEIGHTS))
//        KRATOS_WATCH(GetValue(EXTRACTION_OPERATOR))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_1))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_2))
//        KRATOS_WATCH(GetValue(NURBS_DEGREE_3))

        // try to read the extraction operator from the elemental data
        Matrix ExtractionOperator;
        if( this->Has( EXTRACTION_OPERATOR ) )
        {
            ExtractionOperator = this->GetValue( EXTRACTION_OPERATOR );
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
        }
        else if( this->Has( EXTRACTION_OPERATOR_CSR_ROWPTR )
             && this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
             && this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
        {
            Vector rowPtr = this->GetValue( EXTRACTION_OPERATOR_CSR_ROWPTR ); // must be 0-base
            Vector colInd = this->GetValue( EXTRACTION_OPERATOR_CSR_COLIND ); // must be 0-base
            Vector values = this->GetValue( EXTRACTION_OPERATOR_CSR_VALUES );
            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
        }

//        KRATOS_WATCH(ExtractionOperator)

        // initialize the geometry
        mpIsogeometricGeometry->AssignGeometryData(
            this->GetValue(NURBS_KNOTS_1),
            this->GetValue(NURBS_KNOTS_2),
            this->GetValue(NURBS_KNOTS_3),
            this->GetValue(NURBS_WEIGHTS),
            ExtractionOperator,
            this->GetValue(NURBS_DEGREE_1),
            this->GetValue(NURBS_DEGREE_2),
            this->GetValue(NURBS_DEGREE_3),
            2 // only need to compute 2 integration rules
        );

        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

        InitializeJacobian();
        ////////////////////End Initialize geometry_data/////////////////////////////

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize( GetGeometry().size(), dim, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            for ( unsigned int i = 0; i < 3; ++i )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

        //Set up initial value for time integration scheme
//        for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; ++i )
//        {
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL ) =
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS ) =
//                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL_DT ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS_DT ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_NULL ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_EINS ) = ZeroVector( 3 );
//            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_OLD ) = ZeroVector( 3 );
//
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL ) =
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS ) =
//                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_DT ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_ACCELERATION ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) = 0;
//            ( GetGeometry()[i] ).GetSolutionStepValue( REACTION_WATER_PRESSURE ) = 0;
//        }

        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
            mReferencePressures.resize( integration_points.size() );
            InitializeMaterial();
        }

        mIsInitialized = true;

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::InitializeJacobian()
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

    //number of integration points used, mThisIntegrationMethod refers to the
    //integration method defined in the constructor
    const GeometryType::IntegrationPointsArrayType& integration_points =
            mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    //initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference
    // configuration
    GeometryType::JacobiansType J0(integration_points.size());

    mInvJ0.resize(integration_points.size());

    mTotalDomainInitialSize = 0.00;

    for (unsigned int i = 0; i < integration_points.size(); ++i)
    {
        mInvJ0[i].resize(dim, dim, false);
        noalias(mInvJ0[i]) = ZeroMatrix(dim, dim);
    }

    mDetJ0.resize(integration_points.size(), false);

    noalias(mDetJ0) = ZeroVector(integration_points.size());

    //calculating the Jacobian
    J0 = mpIsogeometricGeometry->Jacobian(J0, mThisIntegrationMethod);

    //calculating the inverse Jacobian
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix(J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);
        //calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
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
            mConstitutiveLawVector[i]->ResetMaterial(
                GetProperties(),
                (*mpIsogeometricGeometry),
                row(Ncontainer, i)
            );
        }
    }

    KRATOS_CATCH( "" )
}


/**
* THIS method is called from the scheme at the start of each solution step
* @param rCurrentProcessInfo
*/
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    //calculating shape functions values
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        mConstitutiveLawVector[PointNumber]->InitializeSolutionStep(
            GetProperties(),
            (*mpIsogeometricGeometry),
            row(Ncontainer, PointNumber),
            CurrentProcessInfo
        );
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
    unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize1 = ( number_of_nodes_disp * dim + number_of_nodes_press );
    unsigned int MatSizeU = number_of_nodes_disp * dim;
    unsigned int MatSizeP = number_of_nodes_press;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize1 )
            rLeftHandSideMatrix.resize( MatSize1, MatSize1, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize1, MatSize1 ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize1 )
            rRightHandSideVector.resize( MatSize1, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize1 ); //resetting RHS
    }

    //calculating shape functions values and local gradients
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    double Weight;

    double capillaryPressure;

    double capillaryPressure_Dt;

    double waterPressure;

    double porosity;

    double density;

    double DetJ = 0.0;

    Matrix Help_K_UU( MatSizeU, MatSizeU );

    Matrix Help_K_UW( MatSizeU, MatSizeP );

    Matrix Help_K_WU( MatSizeP, MatSizeU );

    Matrix Help_K_WW( MatSizeP, MatSizeP );

    Vector Help_R_U( MatSizeU );

    Vector Help_R_W( MatSizeP );

    Matrix DN_DX_PRESS( number_of_nodes_press, 3 );

    Vector N_DISP( number_of_nodes_disp );

    Vector N_PRESS( number_of_nodes_press );

    Matrix tanC_W( 3, 3 );

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_K_UU ) = ZeroMatrix( MatSizeU, MatSizeU );
        noalias( Help_K_UW ) = ZeroMatrix( MatSizeU, MatSizeP );

        noalias( Help_K_WU ) = ZeroMatrix( MatSizeP, MatSizeU );
        noalias( Help_K_WW ) = ZeroMatrix( MatSizeP, MatSizeP );
    }

    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_R_U ) = ZeroVector( MatSizeU );
        noalias( Help_R_W ) = ZeroVector( MatSizeP );
    }

    //Initialize local variables
    Matrix B( 6, MatSizeU );

    Matrix TanC_U( 6, 6 );

    Vector StrainVector( 6 );

    Vector StressVector( 6 );

    Matrix DN_DX_DISP( number_of_nodes_disp, dim );

    Matrix CurrentDisp( number_of_nodes_disp, dim );

    //Current displacements
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        CurrentDisp( node, 0 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisp( node, 1 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisp( node, 2 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Z );
    }

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space to compute the average of pressure shape function
    /////////////////////////////////////////////////////////////////////////
    Vector N_PRESS_averaged( number_of_nodes_press );
    N_PRESS_averaged = ZeroVector( number_of_nodes_press );
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        Weight = integration_points[PointNumber].Weight();
        DetJ   = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer, PointNumber );

        noalias( N_PRESS_averaged ) += N_PRESS * Weight * DetJ;
    }
    N_PRESS_averaged /= mTotalDomainInitialSize;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space to compute the average of current capillaryPressure
    /////////////////////////////////////////////////////////////////////////
    double averageCapillaryPressure_Dt = 0.0;
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        Weight = integration_points[PointNumber].Weight();
        DetJ   = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer, PointNumber );
        capillaryPressure_Dt = GetDerivativeDCapillaryPressureDt(N_PRESS);

        averageCapillaryPressure_Dt += capillaryPressure_Dt * Weight * DetJ;
     }
     averageCapillaryPressure_Dt /= mTotalDomainInitialSize;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space to compute lhs and rhs contribution
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        noalias( DN_DX_PRESS ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        noalias( DN_DX_DISP ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        Weight = integration_points[PointNumber].Weight();

        DetJ = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer, PointNumber );
        noalias( N_DISP ) = row( Ncontainer, PointNumber );

        //Initializing B_Operator at the current integration point
        CalculateBoperator( B, DN_DX_DISP );

        //Calculate the current strain vector using the B-Operator
        CalculateStrain( B, CurrentDisp, StrainVector );

        GetPressures( N_PRESS, capillaryPressure, waterPressure );
        // Remarks: comment to maintain the consistent linearisation
//        if( waterPressure < 0.0 )
//        {
//            waterPressure = 0.0;
//        }

        porosity = GetPorosity( DN_DX_DISP );

        if ( porosity > 1.0 || porosity < 0.0 )
        {
            KRATOS_THROW_ERROR(std::logic_error, "Error porosity = ", porosity)
        }

        density = GetAveragedDensity( capillaryPressure, porosity );

        CalculateStressAndTangentialStiffnessUnsaturatedSoils( StressVector, TanC_U, tanC_W, StrainVector, waterPressure,  PointNumber, rCurrentProcessInfo );
        if ( CalculateStiffnessMatrixFlag == true )
        {
            //Calculation of spatial Stiffnes and Mass Matrix
            CalculateStiffnesMatrixUU( Help_K_UU, TanC_U, B, DN_DX_DISP, N_DISP, density,  capillaryPressure, Weight, DetJ );
            CalculateStiffnesMatrixUW( Help_K_UW, tanC_W, DN_DX_DISP, N_DISP, N_PRESS, capillaryPressure, Weight, DetJ );
            CalculateStiffnesMatrixWU( Help_K_WU, DN_DX_DISP, DN_DX_PRESS, N_PRESS, capillaryPressure,  Weight, DetJ );
            CalculateStiffnesMatrixWW( Help_K_WW, DN_DX_DISP, DN_DX_PRESS, N_PRESS, capillaryPressure, Weight, DetJ );
            //Calculate stabilisation term
            CalculateStiffnesMatrixWWs( Help_K_WW, N_PRESS, N_PRESS_averaged, Weight, DetJ );
        }

        if ( CalculateResidualVectorFlag == true )
        {
            //Calculation of spatial Loadvector
            AddBodyForcesToRHSVectorU( Help_R_U, N_DISP, density, Weight, DetJ );
            AddInternalForcesToRHSU( Help_R_U, B, StressVector, Weight, DetJ );
            AddInternalForcesToRHSW( Help_R_W, DN_DX_DISP, DN_DX_PRESS, N_PRESS, capillaryPressure, Weight, DetJ );
            //Calculate stabilisation term
            capillaryPressure_Dt = GetDerivativeDCapillaryPressureDt(N_PRESS);
            AddInternalForcesToRHSWs( Help_R_W, N_PRESS, N_PRESS_averaged, capillaryPressure_Dt, averageCapillaryPressure_Dt, Weight, DetJ );
        }
    }
    ///////////////////////////////////////////////////////////////////////
    // END Integration in space sum_(beta=0)^(number of quadrature points)
    ///////////////////////////////////////////////////////////////////////

//    KRATOS_WATCH(Help_K_UU)
//    KRATOS_WATCH(Help_K_UW)
//    KRATOS_WATCH(Help_K_WU)
//    KRATOS_WATCH(Help_K_WW)

    if ( CalculateStiffnessMatrixFlag == true )
    {
        AssembleTimeSpaceStiffnessFromStiffSubMatrices( rLeftHandSideMatrix, Help_K_UU, Help_K_UW, Help_K_WU, Help_K_WW );
    }

    if ( CalculateResidualVectorFlag == true )
    {
        AssembleTimeSpaceRHSFromSubVectors( rRightHandSideVector, Help_R_U, Help_R_W );
    }

//    KRATOS_WATCH(rLeftHandSideMatrix)
//    KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    DampMatrix(rMassMatrix, rCurrentProcessInfo);
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
    unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
//    unsigned int number_of_nodes = number_of_nodes_disp+number_of_nodes_press;
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize1 = ( number_of_nodes_disp * dim + number_of_nodes_press );
    unsigned int MatSizeU = number_of_nodes_disp * dim;
    unsigned int MatSizeP = number_of_nodes_press;

    if ( rDampMatrix.size1() != MatSize1 )
        rDampMatrix.resize( MatSize1, MatSize1 );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize1, MatSize1 ); //resetting LHS

    //calculating shape functions values and local gradients
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    //auxiliary terms
//    InitializeGalerkinScheme(rCurrentProcessInfo);

    double Weight;

    double capillaryPressure;

    double waterPressure;

    Matrix Help_D_UU( MatSizeU, MatSizeU );

    Matrix Help_D_UW( MatSizeU, MatSizeP );

    Matrix Help_D_WU( MatSizeP, MatSizeU );

    Matrix Help_D_WW( MatSizeP, MatSizeP );

    Matrix DN_DX_DISP( number_of_nodes_disp, 3 );

    Vector N_PRESS( number_of_nodes_press );

    noalias( Help_D_UU ) = ZeroMatrix( MatSizeU, MatSizeU );

    noalias( Help_D_UW ) = ZeroMatrix( MatSizeU, MatSizeP );

    noalias( Help_D_WU ) = ZeroMatrix( MatSizeP, MatSizeU );

    noalias( Help_D_WW ) = ZeroMatrix( MatSizeP, MatSizeP );

    double DetJ = 0.0;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space to compute the average of pressure shape function
    /////////////////////////////////////////////////////////////////////////
    Vector N_PRESS_averaged( number_of_nodes_press );
    N_PRESS_averaged = ZeroVector( number_of_nodes_press );
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        Weight = integration_points[PointNumber].Weight();
        DetJ   = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer, PointNumber );

        noalias( N_PRESS_averaged ) += N_PRESS * Weight * DetJ;
    }
    N_PRESS_averaged /= mTotalDomainInitialSize;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        noalias( DN_DX_DISP ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        //double DetDef= Determinant_DeformationTensor(DN_DX_DISP);

        Weight = integration_points[PointNumber].Weight();
        // Jacobian on current quadrature point

        DetJ = mDetJ0[PointNumber];
        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer, PointNumber );

        GetPressures( N_PRESS, capillaryPressure, waterPressure );

        // Calculation of spatial Stiffnes and Mass Matrix
        CalculateDampingMatrixWU( Help_D_WU, DN_DX_DISP, N_PRESS,
                                  capillaryPressure, Weight, DetJ );
        CalculateDampingMatrixWW( Help_D_WW, DN_DX_DISP, N_PRESS,
                                  capillaryPressure, Weight, DetJ );
        CalculateDampingMatrixWWs( Help_D_WW, DN_DX_DISP, N_PRESS,
                                  N_PRESS_averaged, Weight, DetJ );
//        CalculateMassMatrix(HelpMassMatrix, N, Weight,DetJ,density);
    }
    ///////////////////////////////////////////////////////////////////////
    // END Integration in space sum_(beta=0)^(number of quadrature points)
    ///////////////////////////////////////////////////////////////////////

    AssembleTimeSpaceStiffnessFromDampSubMatrices( rDampMatrix, Help_D_UU, Help_D_UW, Help_D_WU, Help_D_WW );

//    KRATOS_WATCH(Help_D_UU)
//    KRATOS_WATCH(Help_D_UW)
//    KRATOS_WATCH(Help_D_WU)
//    KRATOS_WATCH(Help_D_WW)

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
}

////************************************************************************************
////************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    //calculating shape functions values
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    unsigned int NumberOfIntegrationPoints =
        mpIsogeometricGeometry->IntegrationPointsNumber(mThisIntegrationMethod);

    for (unsigned int Point = 0; Point < NumberOfIntegrationPoints; ++Point)
    {
        mConstitutiveLawVector[Point]->FinalizeSolutionStep(
            GetProperties(),
            (*mpIsogeometricGeometry),
            row(Ncontainer, Point),
            CurrentProcessInfo
        );
    }

    Vector Dummy_Vector( 9 );

    noalias( Dummy_Vector ) = mConstitutiveLawVector[0]->GetValue( INTERNAL_VARIABLES, Dummy_Vector );

    for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_X ) = Dummy_Vector( 0 );
        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Y ) = Dummy_Vector( 1 );
        GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Z ) = Dummy_Vector( 2 );
        GetGeometry()[i].GetSolutionStepValue( PRESSURE ) = Dummy_Vector( 3 );
        GetGeometry()[i].GetSolutionStepValue( ERROR_RATIO ) = Dummy_Vector( 4 );
        GetGeometry()[i].GetSolutionStepValue( EXCESS_PORE_WATER_PRESSURE ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) - 9.81 * 1000.0 * ( 20.0 - GetGeometry()[i].Z() );
    }
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateOnIntegrationPoints( const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

    //calculating shape functions values and local gradients
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size() );

    Vector N_PRESS( number_of_nodes_press );

    double capillaryPressure;

    double waterPressure;

    double saturation;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points) ////
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        // Shape Functions on current spatial quadrature point
        if ( N_PRESS.size() != number_of_nodes_press )
            N_PRESS.resize( number_of_nodes_press );

        noalias( N_PRESS ) = row( Ncontainer, PointNumber );
        GeometryType::CoordinatesArrayType gp_position;
        gp_position = GetGeometry().GlobalCoordinates( gp_position, GetGeometry().IntegrationPoints()[PointNumber] );

        GetPressures( N_PRESS, capillaryPressure, waterPressure );

        saturation = GetSaturation( capillaryPressure );

        if ( rVariable == SATURATION )
        {
            Output[PointNumber] = saturation;
        }

        if ( rVariable == WATER_PRESSURE )
        {
            Output[PointNumber] = waterPressure;
        }

        if ( rVariable == EXCESS_PORE_WATER_PRESSURE )
        {
            Output[PointNumber] = waterPressure - mReferencePressures[PointNumber];
        }

        if ( rVariable == AIR_PRESSURE )
        {
            Output[PointNumber] = 0.0;
        }
    }

    KRATOS_CATCH( "" )
}

/**
* Calculate Vector Variables at each integration point, used for postprocessing etc.
* @param rVariable Global name of the variable to be calculated
* @param output Vector to store the values on the qudrature points, output of the method
* @param rCurrentProcessInfo
*/
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
{
    GetValueOnIntegrationPoints( rVariable, Output, rCurrentProcessInfo );
}


//************************************************************************************
//************************************************************************************

inline void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateAndAddExtForceContribution( const Vector& N, const ProcessInfo& CurrentProcessInfo, Vector& BodyForce, VectorType& rRightHandSideVector, double weight )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::EquationIdVector( EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo ) const
{
    unsigned int dim_press = 1; //one pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize =
        ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press +
        ( mNodesDispMax - mNodesDispMin + 1 ) * dim_disp;

    if ( rResult.size() != MatSize )
        rResult.resize( MatSize );

    unsigned int adddisp = ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press;

    for ( unsigned int i = ( mNodesPressMin - 1 ); i < mNodesPressMax; ++i )
    {
        int index = ( i - mNodesPressMin + 1 ) * dim_press;
        rResult[index] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
    }

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; ++i )
    {
        int index = adddisp + ( i - mNodesDispMin + 1 ) * dim_disp;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo&
        CurrentProcessInfo ) const
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = ( mNodesPressMin - 1 ); i < mNodesPressMax; ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
    }

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetValuesVector( Vector& values, int Step )
{
    unsigned int dim_press = 1;//two pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize =
        ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press + ( mNodesDispMax -
                mNodesDispMin + 1 ) * dim_disp;

    if ( values.size() != MatSize )
        values.resize( MatSize );

    unsigned int adddisp = ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press;

    for ( unsigned int i = ( mNodesPressMin - 1 ); i < mNodesPressMax; ++i )
    {
        int index = ( i - mNodesPressMin + 1 ) * dim_press;
        values( index ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
    }

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; ++i )
    {
        int index = adddisp + ( i - mNodesDispMin + 1 ) * dim_disp;

        values( index )     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values( index + 1 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        values( index + 2 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }

}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AssembleTimeSpaceStiffnessFromDampSubMatrices
( MatrixType& rLeftHandSideMatrix, const Matrix& D_UU, const Matrix& D_UW,
  const Matrix&  D_WU, const Matrix&  D_WW )
{

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;
    unsigned int number_of_nodes_press = mNodesPressMax - mNodesPressMin + 1;
    unsigned int dim_disp = dimension;
    unsigned int dim_press = 1;
    unsigned int index_time_prim;
    unsigned int index_time_sec;
    unsigned int index_space_prim;
    unsigned int index_space_sec;
    unsigned int addIndex_disp = number_of_nodes_press * dim_press;

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim_disp; ++i )
        {
            index_space_prim = prim * dim_disp + i;

            index_time_prim = addIndex_disp + prim * dim_disp + i;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; ++sec )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    += ( -1 ) * D_UU( index_space_prim, index_space_sec );
                }
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim_disp; ++i )
        {
            index_space_prim = prim * dim_disp + i;

            index_time_prim = addIndex_disp + prim * dim_disp + i;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; ++sec )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                += ( -1 ) * D_UW( index_space_prim, index_space_sec );
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_press; ++prim )
    {

        index_space_prim = prim;

        index_time_prim = prim * dim_press;

        for ( unsigned int sec = 0; sec < number_of_nodes_disp; ++sec )
        {
            for ( unsigned int j = 0; j < dim_disp; j++ )
            {
                index_space_sec = sec * dim_disp + j;

                index_time_sec = addIndex_disp + sec * dim_disp + j;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                += ( -1 ) * D_WU( index_space_prim, index_space_sec );
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_press; ++prim )
    {
        index_space_prim = prim;

        index_time_prim = prim * dim_press;

        for ( unsigned int sec = 0; sec < number_of_nodes_press; ++sec )
        {
            index_space_sec = sec;

            index_time_sec = sec * dim_press;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += ( -1 ) * D_WW( index_space_prim, index_space_sec );
        }
    }
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AssembleTimeSpaceStiffnessFromStiffSubMatrices
( MatrixType& rLeftHandSideMatrix, const Matrix& K_UU, const Matrix& K_UW,
  const Matrix& K_WU, const Matrix& K_WW )
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;
    unsigned int number_of_nodes_press = mNodesPressMax - mNodesPressMin + 1;
    unsigned int dim_disp = dimension;
    unsigned int dim_press = 1;
    unsigned int index_time_prim;
    unsigned int index_time_sec;
    unsigned int index_space_prim;
    unsigned int index_space_sec;
    unsigned int addIndex_disp = number_of_nodes_press * dim_press;

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim_disp; ++i )
        {
            index_space_prim = prim * dim_disp + i;

            index_time_prim = addIndex_disp + prim * dim_disp + i;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; ++sec )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    += ( -1 ) * K_UU( index_space_prim, index_space_sec );
                }
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim_disp; ++i )
        {
            index_space_prim = prim * dim_disp + i;

            index_time_prim = addIndex_disp + prim * dim_disp + i;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; ++sec )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                += ( -1 ) * K_UW( index_space_prim, index_space_sec );
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_press; ++prim )
    {

        index_space_prim = prim;

        index_time_prim = prim * dim_press;

        for ( unsigned int sec = 0; sec < number_of_nodes_disp; ++sec )
        {
            for ( unsigned int j = 0; j < dim_disp; j++ )
            {
                index_space_sec = sec * dim_disp + j;

                index_time_sec = addIndex_disp + sec * dim_disp + j;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                += ( -1 ) * K_WU( index_space_prim, index_space_sec );
            }
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_press; ++prim )
    {
        index_space_prim = prim;

        index_time_prim = prim * dim_press;

        for ( unsigned int sec = 0; sec < number_of_nodes_press; ++sec )
        {
            index_space_sec = sec;

            index_time_sec = sec * dim_press;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += ( -1 ) * K_WW( index_space_prim, index_space_sec );
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AssembleTimeSpaceRHSFromSubVectors( VectorType&
        rRightHandSideVector, const Vector& R_U, const Vector& R_W )
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
    unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
    unsigned int dim_disp = dimension;
    unsigned int dim_press = 1;
    unsigned int index_time;
    unsigned int addIndex_disp = number_of_nodes_press * dim_press;

//             noalias(rRightHandSideVector)= ZeroVector(number_of_nodes_press*dim_1+(number_of_nodes_disp-number_of_nodes_press)*dim_2);

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim_disp; ++i )
        {
            index_time = addIndex_disp + prim * dim_disp + i;

            rRightHandSideVector( index_time ) +=
                R_U( prim * dim_disp + i );
        }
    }

    for ( unsigned int prim = 0; prim < number_of_nodes_press; ++prim )
    {
        index_time = prim * dim_press;

        rRightHandSideVector( index_time ) +=
            R_W( prim );
    }

    KRATOS_CATCH( "" )

}

//************************************************************************************
//CALCULATE EXTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AddBodyForcesToRHSVectorU( Vector& R, Vector& N_DISP, double density, double Weight, double detJ )
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector gravity( dim );

    double density = 0.0;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        noalias( gravity ) = GetValue( GRAVITY );
        density = GetValue( DENSITY );
    }
    else
    {
        noalias( gravity ) = GetProperties()[GRAVITY];
        density = GetProperties()[DENSITY];
    }

    for ( unsigned int prim = 0; prim < ( mNodesDispMax - mNodesDispMin + 1 ); ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            R( prim*dim + i ) +=
                N_DISP( prim ) * density * gravity( i ) * Weight * detJ;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE INTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AddInternalForcesToRHSU( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
{
    KRATOS_TRY

    noalias( R ) -= detJ * Weight * prod( trans( B_Operator ), StressVector );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE STIFFNESS MATRICES DISPLACEMENT*******************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStiffnesMatrixUU( Matrix& K, const
        Matrix& C, const Matrix& B_Operator, const Matrix& DN_DX_DISP, Vector& N_DISP,
        double density, double capillaryPressure, double Weight, double detJ )
{
    KRATOS_TRY
// boost::timer timer3;
// clock_t t1;
// clock_t t2;
//
// t1 = clock();    // Referenzwert zum Programmstart speichern
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;


    Vector gravity( dim );

    double density_soil = 0.0;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        noalias( gravity ) = GetValue( GRAVITY );
        density_soil = GetValue( DENSITY );
    }
    else
    {
        noalias( gravity ) = GetProperties()[GRAVITY];
        density_soil = GetProperties()[DENSITY];
    }
    double density_water = GetValue( DENSITY_WATER );
    double density_air = GetValue( DENSITY_AIR );

    double saturation = GetSaturation( capillaryPressure );

    double porosity_divu = 0.0;
//                 double porosity_divu= GetDerivativeDPorosityDDivU(DN_DX_DISP);


    double DrhoDdivU =  porosity_divu * ( -density_soil + ( 1 - saturation ) * density_air + saturation * density_water );
// std::cout << "time in 1071 --> " << timer3.elapsed() << std::endl;
// t2 = clock();
// std::cout <<"line 1077   "<< double(t2-t1)<< std::endl;

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            for ( unsigned int sec = 0; sec < number_of_nodes_disp; ++sec )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {

                    K( prim*dim + i, sec*dim + j ) += N_DISP( prim ) * DrhoDdivU * gravity( i ) * DN_DX_DISP( sec, j ) * Weight * detJ;

//                                                 for(unsigned int alpha=0; alpha<6; alpha++)
//                                                         for(unsigned int beta=0; beta<6; beta++)
//                                                                 K(prim*dim+i,sec*dim+j) +=(-1)*(B_Operator(alpha,prim*3+i)*C(alpha,beta)*B_Operator(beta,sec*3+j))*detJ*Weight;
                }
            }
        }
    }

// std::cout << "time in 1090 --> " << timer3.elapsed() << std::endl;
// t2 = clock();
// std::cout <<"line 1098   "<< double(t2-t1)<< std::endl;
//                 noalias(K)-= prod(trans(B_Operator),(Weight*detJ)*Matrix(prod(C,B_Operator)) );
    noalias( K ) -=
        prod( trans( B_Operator ), ( Weight * detJ ) * Matrix( prod( C, B_Operator ) ) );

// t2 = clock();
// std::cout <<"line 1101   "<< double(t2-t1)<< std::endl;
// std::cout << "time in 1092 --> " << timer3.elapsed() << std::endl;
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStiffnesMatrixUW( Matrix&
        Help_K_UW, Matrix& tanC_W, const Matrix& DN_DX_DISP, Vector& N_DISP, Vector& N_PRESS
        , double capillaryPressure,
        double Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

    Vector gravity( dim );

    // double density_soil = 0.0;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        noalias( gravity ) = GetValue( GRAVITY );
        //density_soil = GetValue( DENSITY );
    }
    else
    {
        noalias( gravity ) = GetProperties()[GRAVITY];
        //density_soil = GetProperties()[DENSITY];
    }
    double density_water = GetValue( DENSITY_WATER );
    double density_air = GetValue( DENSITY_AIR );


    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double porosity = GetPorosity( DN_DX_DISP );

    double DrhoDp_w = porosity * ( density_air - density_water ) * DSDpc;


    for ( unsigned int prim = 0; prim < displacement_size; ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            for ( unsigned int sec = 0; sec < pressure_size; ++sec )
            {
                Help_K_UW( prim*dim + i, sec ) +=
                    N_DISP( prim ) * DrhoDp_w * gravity( i ) * N_PRESS( sec )
                    * Weight * DetJ;

                for ( unsigned int gamma = 0; gamma < 3; gamma++ )
                {
                    Help_K_UW( prim*dim + i, sec ) +=
                        ( -1 ) * ( DN_DX_DISP( prim, gamma ) * tanC_W( i, gamma ) * N_PRESS( sec ) ) * Weight * DetJ;
                }
            }
        }
    }
}

//************************************************************************************
//CALCULATE FORCEVECTORS WATER********************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AddInternalForcesToRHSW( Vector& Help_R_W, const
        Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double  DetJ )
{

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double porosity = GetPorosity( DN_DX_DISP );

    double DS_Dpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double saturation = GetSaturation( capillaryPressure );

    double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

    double div_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

    Vector flow_water( dim );
    noalias( flow_water ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        Help_R_W( prim ) +=
            N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt * Weight * DetJ/* * GetProperties()[SCALE]*/ ;
        Help_R_W( prim ) +=
            N_PRESS( prim ) * saturation * div_Dt * Weight * DetJ/* * GetProperties()[SCALE]*/ ;

        for ( unsigned int gamma = 0; gamma < dim; gamma++ )
        {
            Help_R_W( prim ) +=
                ( -1 ) * ( DN_DX_PRESS( prim, gamma ) * flow_water( gamma ) )
                * Weight * DetJ/* * GetProperties()[SCALE]*/ ;
        }
    }
}

//************************************************************************************
//CALCULATE STABILISATION TERM FOR FORCEVECTORS WATER*********************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::AddInternalForcesToRHSWs(
    Vector& Help_R_W,
    Vector& N_PRESS,
    Vector& N_PRESS_averaged,
    double capillaryPressure_Dt,
    double averageCapillaryPressure_Dt,
    double Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    double t, G;
    if(GetProperties().Has(STABILISATION_FACTOR) && GetProperties().Has(SHEAR_MODULUS))
    {
        t = GetProperties()[STABILISATION_FACTOR];
        G = GetProperties()[SHEAR_MODULUS];
    }
    else
    {
        t = 0.0;
        G = 1.0;
    }

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        Help_R_W( prim ) += (-0.5 * t / G) *
                        ( N_PRESS( prim ) - N_PRESS_averaged( prim ) ) *
                        ( capillaryPressure_Dt - averageCapillaryPressure_Dt ) *
                        Weight * DetJ;
    }
}

//************************************************************************************
//CALCULATE STIFFNESS MATRICES WATER**************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStiffnesMatrixWU
( Matrix& Help_K_WU, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
{

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

    double DnDdivU = GetDerivativeDPorosityDDivU( DN_DX_DISP );

    Vector flow_water( dim );
    noalias( flow_water ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < displacement_size; ++sec )
        {
            for ( unsigned int j = 0; j < dim; ++j )
            {
                Help_K_WU( prim, sec*dim + j ) +=
                    N_PRESS( prim ) * DnDdivU * DSDpc * Dpc_Dt * DN_DX_DISP( sec, j )
                    * Weight * DetJ/* * GetProperties()[SCALE]*/ ;
            }
        }
    }

}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStiffnesMatrixWW
( Matrix& Help_K_WW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS,
  Vector& N_PRESS, double capillaryPressure, double
  Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

    double porosity = GetPorosity( DN_DX_DISP );

    double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

    double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

    Vector Dflow_waterDpw( dim );

    noalias( Dflow_waterDpw ) = GetDerivativeDWaterFlowDpw( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

    double Dflow_waterDgradpw = GetDerivativeDWaterFlowDGradpw( DN_DX_DISP, capillaryPressure );

    double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < pressure_size; ++sec )
        {
            Help_K_WW( prim, sec ) +=
                ( -1 ) * N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                * Weight * DetJ/* * GetProperties()[SCALE]*/;

            Help_K_WW( prim, sec ) +=
                ( -1 ) * N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                * Weight * DetJ/* * GetProperties()[SCALE]*/;

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_K_WW( prim, sec ) +=
                    ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_waterDpw( gamma )
                    * N_PRESS( sec )
                    * Weight * DetJ/* * GetProperties()[SCALE]*/;
                Help_K_WW( prim, sec ) +=
                    ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_waterDgradpw
                    * DN_DX_PRESS( sec, gamma ) * Weight * DetJ/* * GetProperties()[SCALE]*/;
            }
        }
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStiffnesMatrixWWs
( Matrix& Help_K_WW, Vector& N_PRESS, Vector& N_PRESS_averaged, double Weight, double DetJ )
{
}

//************************************************************************************
//CALCULATE DAMPING MATRICES WATER****************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateDampingMatrixWU
( Matrix& Help_D_WU, const Matrix&
  DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

    double saturation = GetSaturation( capillaryPressure );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < displacement_size; ++sec )
        {
            for ( unsigned int j = 0; j < dim; j++ )
            {
                Help_D_WU( prim, sec*dim + j ) +=
                    N_PRESS( prim ) * saturation * DN_DX_DISP( sec, j )
                    * Weight * DetJ/* * GetProperties()[SCALE]*/;
            }
        }
    }
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateDampingMatrixWW( Matrix& Help_D_WW, const Matrix&
        DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
{
    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;
    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
    double porosity = GetPorosity( DN_DX_DISP );

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < pressure_size; ++sec )
        {
            Help_D_WW( prim, sec ) +=
                ( -1 ) * N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                * Weight * DetJ /* * GetProperties()[SCALE]*/;
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateDampingMatrixWWs( Matrix& Help_D_WW, const Matrix&
        DN_DX_DISP, Vector& N_PRESS, Vector& N_PRESS_averaged, double Weight, double DetJ )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

    double t, G;
    if(GetProperties().Has(STABILISATION_FACTOR) && GetProperties().Has(SHEAR_MODULUS))
    {
        t = GetProperties()[STABILISATION_FACTOR];
        G = GetProperties()[SHEAR_MODULUS];
    }
    else
    {
        t = 0.0;
        G = 1.0;
    }

    for ( unsigned int prim = 0; prim < pressure_size; ++prim )
    {
        for ( unsigned int sec = 0; sec < pressure_size; ++sec )
        {
            Help_D_WW( prim, sec ) += (0.5 * t / G) *
                ( N_PRESS( prim ) - N_PRESS_averaged( prim ) ) *
                ( N_PRESS( sec ) - N_PRESS_averaged( sec ) ) *
                Weight * DetJ;
        }
    }
}

//************************************************************************************
//PRIMARY VARIABLES AND THEIR DERIVATIVES
//************************************************************************************
Matrix UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateDisplacementGradient( const Matrix& DN_DX_DISP )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Matrix result( dim, dim );
    noalias( result ) = ZeroMatrix( dim, dim );

    Vector disp_alpha( dim );

    for ( unsigned int point = ( mNodesDispMin - 1 ); point < mNodesDispMax; point++ )
    {
        noalias( disp_alpha ) = GetGeometry()[point].GetSolutionStepValue( DISPLACEMENT );

        for ( unsigned int j = 0; j < 3; j++ )
        {
            for ( unsigned int k = 0; k < dim; k++ )
            {
                result( j, k ) += disp_alpha( j ) * DN_DX_DISP(( point - mNodesDispMin + 1 ), k );
            }
        }
    }

    return result;
}

//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetGradientWater( const Matrix& DN_DX_PRESS )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector result( dim );
    noalias( result ) = ZeroVector( dim );

    double presW_alpha;

    for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ; ++i )
    {
        presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        for ( unsigned int k = 0; k < 3; k++ )
        {
            result( k ) +=
                presW_alpha
                * DN_DX_PRESS(( i - mNodesPressMin + 1 ), k );
        }
    }

    return result;
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetPressures( const Vector& N_PRESS,
        double& capillaryPressure, double& waterPressure )
{
    capillaryPressure = 0.0;

    waterPressure = 0.0;

    double presW_alpha;

    for ( unsigned int i = ( mNodesPressMin - 1 ) ; i < mNodesPressMax ; ++i )
    {
        presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        capillaryPressure +=
            ( -presW_alpha )
            * N_PRESS( i - mNodesPressMin + 1 );

        waterPressure +=
            presW_alpha
            * N_PRESS( i - mNodesPressMin + 1 );
    }
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDPressuresDt( const Vector& N_PRESS, double& capillaryPressure_Dt, double& waterPressure_Dt )
{
    capillaryPressure_Dt = 0.0;

    waterPressure_Dt = 0.0;

    double presW_alpha_Dt;

    for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ; ++i )
    {
        presW_alpha_Dt =
            GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

        capillaryPressure_Dt +=
            ( -presW_alpha_Dt )
            * N_PRESS( i - mNodesPressMin + 1 );

        waterPressure_Dt +=
            presW_alpha_Dt
            * N_PRESS( i - mNodesPressMin + 1 );
    }
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDCapillaryPressureDt( const Vector& N_PRESS )
{
    double capillaryPressure_Dt = 0.0;

    double presW_alpha_Dt;

    for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ; ++i )
    {
        presW_alpha_Dt = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

        capillaryPressure_Dt +=
            ( -presW_alpha_Dt )
            * N_PRESS( i - mNodesPressMin + 1 );
    }

    return capillaryPressure_Dt;
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

//POROSITY AND ITS DERIVATIVES

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetPorosity( const Matrix& DN_DX_DISP )
{
    double porosity = GetValue(POROSITY);
    return porosity;
}

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP )
{
//              double initialPorosity= GetProperties()[POROSITY];
//              double div= GetDivU(DN_DX_DISP);

    double porosity_divu = 0.0;
//              double porosity_divu= (1-initialPorosity)*exp(-div);


    return porosity_divu;
}



//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetGravity()
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    Vector gravity( dim );
    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];
    return gravity;
}

//************************************************************************************
//************************************************************************************
double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDivU( const Matrix& DN_DX_DISP )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha( dim );

    for ( unsigned int i = mNodesDispMin - 1 ; i < mNodesDispMax ; ++i )
    {
//                     //nodal Displacements
        noalias( u_alpha ) =
            GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT );
        //thus strain in vector format

        for ( unsigned int k = 0; k < dim; k++ )
        {
            div += ( u_alpha( k ) )
                   * DN_DX_DISP( i - mNodesDispMin + 1, k );
        }
    }

    return div;
}

//************************************************************************************
//************************************************************************************
double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDDivUDt( const Matrix& DN_DX_DISP )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha_Dt( dim );

    for ( unsigned int i = mNodesDispMin - 1 ; i < mNodesDispMax ; ++i )
    {
//                                //nodal Displacements
        noalias( u_alpha_Dt ) =
            GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT );

        //thus strain in vector format

        for ( unsigned int k = 0; k < dim; k++ )
        {
            div += u_alpha_Dt( k ) * DN_DX_DISP( i - mNodesDispMin + 1, k );
        }
    }

    return div;
}

//AVERAGED DENSITY
//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetAveragedDensity( double capillaryPressure,
        double porosity )
{
    double result = 0.0;

    double density_soil = 0.0;
    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
    {
        density_soil = GetValue(DENSITY);
    }
    else
    {
        density_soil = GetProperties()[DENSITY];
    }
    double density_air = GetValue(DENSITY_AIR);
    double density_water = GetValue(DENSITY_WATER);
    double saturation = GetSaturation( capillaryPressure );

    result = ( 1 - porosity ) * density_soil +
             porosity * ( saturation * density_water + ( 1 - saturation ) * density_air );

    return result;
}


//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

//SATURATION AND ITS DERIVATIVES

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetSaturation( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0.0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double saturation = 0.0;

    if ( capillaryPressure < 0.0 )
        capillaryPressure = 0.0;

    saturation = pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c ) );

// For Liakopolous Benchmark
// saturation =  1.0-1.9722*1e-11*pow(capillaryPressure,2.4279);


    return saturation;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDSaturationDpc( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double result = 0.0;

    if ( capillaryPressure < 0.0 )
    {
        capillaryPressure = 0.0;
    }

    result = ( -c ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * b *
             pow(( capillaryPressure / airEntryPressure ), ( b - 1 ) ) * 1.0 / airEntryPressure;

// For Liakopolous Benchmark
// result =  -1.9722*2.4279*1e-11*pow(capillaryPressure,1.4279);

    return result;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetSecondDerivativeD2SaturationDpc2( double capillaryPressure )
{
    double airEntryPressure = GetValue(AIR_ENTRY_VALUE);

    if ( airEntryPressure <= 0 )
        airEntryPressure = 1.0;

    double b = GetValue(FIRST_SATURATION_PARAM);

    double c = GetValue(SECOND_SATURATION_PARAM);

    double result = 0.0;

    if ( capillaryPressure < 0.0 )
        capillaryPressure = 0.0;

    result = ( -c ) * b / airEntryPressure * (
                 ( -c - 1.0 ) * pow(( 1 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 2.0 ) )
                 * b / airEntryPressure * pow(( capillaryPressure / airEntryPressure ), ( 2.0 * ( b - 1.0 ) ) )
                 + pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * ( b - 1.0 )
                 * pow(( capillaryPressure / airEntryPressure ), ( b - 2.0 ) ) * 1.0 / airEntryPressure );

// For Liakopolous Benschmark
// result =  -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure,0.4279);

    return result;
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetFlowWater( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //Calculation of relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );

    if ( relPerm <= 0.01 )
        relPerm = 0.01;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    Vector gravity( dim );

    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];

    Vector result( dim );

    noalias( result ) = ZeroVector( dim );

    Vector grad_water( dim );

    noalias( grad_water ) = GetGradientWater( DN_DX_PRESS );

    for ( unsigned int i = 0; i < dim; ++i )
    {
        result( i ) = -relPerm * GetValue(PERMEABILITY_WATER) /
                      ( GetValue(DENSITY_WATER) * 9.81 )
                      * ( grad_water( i ) - GetValue(DENSITY_WATER)
                          * gravity( i ) );
    }

    return result;
}

//************************************************************************************
//************************************************************************************

Vector UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDWaterFlowDpw( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
    //Calculation of Derivative relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );
    double relPerm_pw = ( -1 ) * DSDpc;

    if ( relPerm <= 0.0001 )
    {
        relPerm = 0.0001;
        relPerm_pw = 0.0;
    }

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm_pw= -2.207*1.0121*pow((1-saturation),0.0121)*(-DSDpc)*(-1);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    Vector result( dim );

    noalias( result ) = ZeroVector( dim );

    Vector gravity( dim );

    if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        noalias( gravity ) = GetValue(GRAVITY);
    else
        noalias( gravity ) = GetProperties()[GRAVITY];

    Vector grad_water( dim );

    noalias( grad_water ) = GetGradientWater( DN_DX_PRESS );

    for ( unsigned int i = 0; i < dim; ++i )
    {
        result( i ) = -relPerm_pw * GetValue(PERMEABILITY_WATER) /
                      ( GetValue(DENSITY_WATER) * 9.81 )
                      * ( grad_water( i ) - GetValue(DENSITY_WATER)
                          * gravity( i ) );
    }

    return result;
}

//************************************************************************************
//************************************************************************************

double UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetDerivativeDWaterFlowDGradpw( const Matrix& DN_DX_DISP, double capillaryPressure )
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //Calculation of Derivative relative Permeability after Mualem
    double relPerm = GetSaturation( capillaryPressure );

    if ( relPerm <= 0.0001 )
        relPerm = 0.0001;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

    double result( dim );

    result =
        -relPerm * GetValue(PERMEABILITY_WATER) / ( GetValue(DENSITY_WATER) * 9.81 );

    return result;
}

//************************************************************************************
//************************************************************************************
//STRESSES, STRAINS AND CONSTITUTIVE MODELL (UNSATURATED CASE)
//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateEffectiveStress( Vector&
        StressVector, Matrix& tanC_W, const
        double waterPressure )
{
    double capillaryPressure = -waterPressure;

    double saturation = GetSaturation( capillaryPressure );

    double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

    noalias( tanC_W ) = ZeroMatrix( 3, 3 );

    for ( unsigned int i = 0; i < 3; ++i )
    {
        StressVector( i ) =
            StressVector( i ) - ( saturation * waterPressure );
    }

    for ( unsigned int i = 0; i < 3; ++i )
    {
        tanC_W( i, i ) = ( DSDpc * waterPressure - saturation );
    }
}

//************************************************************************************
//************************************************************************************

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStressAndTangentialStiffnessUnsaturatedSoils
( Vector& StressVector, Matrix& tanC_U, Matrix& tanC_W,
  Vector& StrainVector, double waterPressure, int PointNumber, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if ( tanC_W.size1() != dim || tanC_W.size2() != dim )
        tanC_W.resize( dim, dim );

    noalias( tanC_W ) = ZeroMatrix( dim, dim );

    if ( tanC_U.size1() != 6 || tanC_U.size2() != 6 )
        tanC_U.resize( 6, 6 );

    noalias( tanC_U ) = ZeroMatrix( 6, 6 );

    if ( StressVector.size() != 6 )
        StressVector.resize( 6 );

    //Set suction in const. law
    mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, -waterPressure,  rCurrentProcessInfo );

    //retrieve the material response
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
        StrainVector,
        ZeroMatrix( 1 ),
        StressVector,
        tanC_U,
        rCurrentProcessInfo,
        GetProperties(),
        GetGeometry(),
        row( GetGeometry().ShapeFunctionsValues(), PointNumber ),
        true,
        1,
        true );

    CalculateEffectiveStress( StressVector, tanC_W, waterPressure );

    KRATOS_CATCH( "" )
}


/**
 * Computes the strain vector
 */
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
{
    KRATOS_TRY
    noalias( StrainVector ) = ZeroVector( 6 );

    for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
    {
        for ( unsigned int item = 0; item < 6; item++ )
            for ( unsigned int dim = 0; dim < 3; dim++ )
                StrainVector[item] += B( item, 3 * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
    }

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::CalculateBoperator
( Matrix& B_Operator, const Matrix& DN_DX )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    noalias( B_Operator ) = ZeroMatrix( 6, number_of_nodes * 3 );

    for ( unsigned int i = 0; i < number_of_nodes; ++i )
    {
        B_Operator( 0, i*3 ) = DN_DX( i, 0 );
        B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
        B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
        B_Operator( 3, i*3 ) = DN_DX( i, 1 );
        B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
        B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
        B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
        B_Operator( 5, i*3 ) = DN_DX( i, 2 );
        B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
    }

    return;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::InitializeMaterial()
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial(
            GetProperties(), GetGeometry(), ZeroVector( 0 ) );
    }

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            if ( rValues[i].size1() != 3 || rValues[i].size2() != 3 )
                rValues[i].resize( 3, 3 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i] );
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    if ( rVariable == MATERIAL_PARAMETERS )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
            rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
    }

    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PRESTRESS, rValues[i] );
        }
    }


    if ( rVariable == PLASTIC_STRAIN_VECTOR )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
        }
    }


    //To Plot Internal variables
    if ( rVariable == INTERNAL_VARIABLES )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            if ( rValues[i].size() != 9 )
                rValues[i].resize( 9 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );

        }
    }

    //To Plot Stresses
    if ( rVariable == STRESSES )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );

            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );

        }
    }

    //To Plot Fluid Flows
    if ( rVariable == FLUID_FLOWS )
    {
        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

        Vector N_PRESS( number_of_nodes_press );

        Matrix DN_DX_PRESS( number_of_nodes_press, 3 );
        Matrix DN_DX_DISP( 0, 0 );

        double capillaryPressure;

        double waterPressure;

        double saturation;

        Vector waterFlow( 3 );

        //calculating shape functions values and local gradients
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        const GeometryType::IntegrationPointsArrayType& integration_points =
            mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            if ( rValues[PointNumber].size() != 9 )
                rValues[PointNumber].resize( 9 );

            // Shape Functions on current spatial quadrature point
            if ( N_PRESS.size() != number_of_nodes_press )
                N_PRESS.resize( number_of_nodes_press );

            noalias( N_PRESS ) = row( Ncontainer, PointNumber );

            GetPressures( N_PRESS, capillaryPressure, waterPressure );

            saturation = GetSaturation( capillaryPressure );

            rValues[PointNumber][0] = waterPressure;

            rValues[PointNumber][1] = 0.0;

            rValues[PointNumber][2] = saturation;

            noalias( DN_DX_PRESS ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            noalias( DN_DX_DISP ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            noalias( waterFlow ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

            rValues[PointNumber][3] = waterFlow( 0 );

            rValues[PointNumber][4] = waterFlow( 1 );

            rValues[PointNumber][5] = waterFlow( 2 );

            rValues[PointNumber][6] = 0.0;

            rValues[PointNumber][7] = 0.0;

            rValues[PointNumber][8] = 0.0;
        }
    }

    //To Plot Coordinates of Integration Points
    if ( rVariable == COORDINATES )
    {
        for ( unsigned int i = 0; i < integration_points.size(); ++i )
        {
            if ( rValues[i].size() != 3 )
                rValues[i].resize( 3 );

            Geometry<Node<3> >::CoordinatesArrayType dummy;

            GetGeometry().GlobalCoordinates( dummy, integration_points[i] );

            noalias( rValues[i] ) = dummy;
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if( rVariable == PLASTICITY_INDICATOR )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );

        //reading integration points and local gradients
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
        }
        return;
    }

    if( rVariable == K0 )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = GetValue( K0 );
        }

        return;
    }

    unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( rValues.size() != integration_points.size() )
        rValues.resize( integration_points.size() );

    //calculating shape functions values and local gradients
    GeometryType::ShapeFunctionsGradientsType DN_De;
    Matrix Ncontainer;

    mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Ncontainer,
        DN_De,
        mThisIntegrationMethod
    );

    Vector N_PRESS( number_of_nodes_press );

    double capillaryPressure;

    double waterPressure;

    double saturation;


    /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        GeometryType::CoordinatesArrayType gp_position;
        gp_position = GetGeometry().GlobalCoordinates( gp_position, GetGeometry().IntegrationPoints()[PointNumber] );

        // Shape Functions on current spatial quadrature point
        if ( N_PRESS.size() != number_of_nodes_press )
            N_PRESS.resize( number_of_nodes_press );

        noalias( N_PRESS ) = row( Ncontainer, PointNumber );

        GetPressures( N_PRESS, capillaryPressure, waterPressure );

        saturation = GetSaturation( capillaryPressure );

        if ( rVariable == SATURATION )
        {
            rValues[PointNumber] = saturation;
        }

        if ( rVariable == WATER_PRESSURE )
        {
            rValues[PointNumber] = waterPressure;
        }

        if ( rVariable == EXCESS_PORE_WATER_PRESSURE )
        {
            rValues[PointNumber] = waterPressure - mReferencePressures[PointNumber];
        }
    }

    KRATOS_CATCH( "" )
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetValueOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::SetValueOnIntegrationPoints( const Kratos::Variable< ConstitutiveLaw::Pointer >& rVariable, std::vector< ConstitutiveLaw::Pointer >& rValues, const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CONSTITUTIVE_LAW )
    {
        for ( unsigned int i = 0; i < rValues.size(); ++i )
        {
            mConstitutiveLawVector[i] = rValues[i];
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::SetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::SetValueOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::IntegrationMethod UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

int UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Element found with Id 0 or negative", __FUNCTION__);
        }

        if (mTotalDomainInitialSize < 0)
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Domain size can not be less than 0. Please check Jacobian.", __FUNCTION__);
        }

        //verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "constitutive law not provided for property ",
                    this->GetProperties().Id());
        }

        return 1;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos


