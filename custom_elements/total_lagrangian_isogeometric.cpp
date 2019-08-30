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

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 20 Jan 2015 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/total_lagrangian_isogeometric.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

    TotalLagrangianIsogeometric::TotalLagrangianIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
        /*
        Important remarks:
            + GetGeometry() and (*mpIsogeometricGeometry) refer to the same instance of IsogeometricGeometryType in the memory. However, (*mpIsogeometricGeometry) only provides access to the functions wrapped by Geometry interface, whereby (*mpIsogeometricGeometry) provides access to functions exclusive to IsogeometricGeometryType. It is ok to replace every instances of (*mpIsogeometricGeometry) by (*mpIsogeometricGeometry) but to keep the code looks compatible (especiall for comparison with old code), (*mpIsogeometricGeometry) can still be kept, but take note to the wrapped functions.
         */
    }

//************************************************************************************
//************************************************************************************

    TotalLagrangianIsogeometric::TotalLagrangianIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
        /*
        Important remarks:
            + GetGeometry() and (*mpIsogeometricGeometry) refer to the same instance of IsogeometricGeometryType in the memory. However, (*mpIsogeometricGeometry) only provides access to the functions wrapped by Geometry interface, whereby (*mpIsogeometricGeometry) provides access to functions exclusive to IsogeometricGeometryType. It is ok to replace every instances of (*mpIsogeometricGeometry) by (*mpIsogeometricGeometry) but to keep the code looks compatible (especiall for comparison with old code), (*mpIsogeometricGeometry) can still be kept, but take note to the wrapped functions.
         */
    }

    Element::Pointer TotalLagrangianIsogeometric::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangianIsogeometric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    TotalLagrangianIsogeometric::~TotalLagrangianIsogeometric()
    {
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

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

//        mThisIntegrationMethod = mpIsogeometricGeometry->GetDefaultIntegrationMethod();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

        const GeometryType::IntegrationPointsArrayType& integration_points =
                mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        //resizing jacobian inverses containers
        mInvJ0.resize( integration_points.size() );
        mDetJ0.resize( integration_points.size(), false );

        GeometryType::JacobiansType J0;
        J0 = mpIsogeometricGeometry->Jacobian( J0, mThisIntegrationMethod );
        mTotalDomainInitialSize = 0.00;

        //Constitutive Law initialisation

        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
            //InitializeMaterial();
        }

        InitializeMaterial();

        //calculating the inverse J0

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();

            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );

            //calculating the total area
            mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo,
                                        bool CalculateStiffnessMatrixFlag,
                                        bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY
        const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int StrainSize = dim * (dim + 1) / 2;

        Matrix B( StrainSize, number_of_nodes * dim );

        Matrix F( dim, dim );

        Matrix D( StrainSize, StrainSize );

        Matrix C( dim, dim );

        Vector StrainVector( StrainSize );

        Vector StressVector( StrainSize );

        Matrix DN_DX( number_of_nodes, dim );

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize )
                rLeftHandSideMatrix.resize( MatSize, MatSize, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size() != MatSize )
                rRightHandSideVector.resize( MatSize, false );

            rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
        }

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points
            = mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod );

        //calculating shape function values and local gradients
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        //calculating actual jacobian
        GeometryType::JacobiansType J;

        mpIsogeometricGeometry->Jacobian( J );

        //auxiliary terms
        Vector BodyForce(dim);
        noalias(BodyForce) = ZeroVector(dim);

        if(GetProperties().Has(BODY_FORCE))
        {
            noalias(BodyForce) += GetProperties()[BODY_FORCE];
        }

        if(GetProperties().Has(DENSITY) && GetProperties().Has(GRAVITY))
        {
            Vector gravity = GetProperties()[GRAVITY];
            double density = GetProperties()[DENSITY];
            noalias(BodyForce) += density * gravity;
        }

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

            //deformation gradient
            noalias( F ) = prod( J[PointNumber], mInvJ0[PointNumber] );

            //strain calculation
            noalias( C ) = prod( trans( F ), F );

            CalculateStrain( C, StrainVector );

            Comprobate_State_Vector( StrainVector );
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector,
                F,
                StressVector,
                D,
                rCurrentProcessInfo,
                GetProperties(),
                (*mpIsogeometricGeometry),
                row( Ncontainer, PointNumber ),
                true,
                CalculateStiffnessMatrixFlag,
                true );

            //calculating operator B
            CalculateB( B, F, DN_DX, StrainVector.size() );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];

            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //contributions to stiffness matrix calculated on the reference config
                noalias( rLeftHandSideMatrix ) += prod( trans( B ), ( IntToReferenceWeight ) * Matrix( prod( D, B ) ) ); //to be optimized to remove the temporary
                CalculateAndAddKg( rLeftHandSideMatrix, DN_DX, StressVector, IntToReferenceWeight );
            }

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                //contribution of external forces
                // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

                // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), StressVector );
            }
        }


        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );


    }

//************************************************************************************
//************************************************************************************

    double TotalLagrangianIsogeometric::CalculateIntegrationWeight( double GaussPointWeight, double DetJ0 )
    {
        //to permorm the integration over the reference domain we need to include
        // the thickness in 2D
        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();
        double weight = GaussPointWeight;

        weight *= DetJ0;

        if ( dimension == 2 ) weight *= GetProperties()[THICKNESS];

        return weight;
    }

////************************************************************************************
////************************************************************************************

    void TotalLagrangianIsogeometric::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                    (*mpIsogeometricGeometry), row( mpIsogeometricGeometry->ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                    CurrentProcessInfo );
    }

////************************************************************************************
////************************************************************************************

    void TotalLagrangianIsogeometric::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        //         std::cout << "in TL: calling FinalizeSolutionStep" << std::endl;
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                    (*mpIsogeometricGeometry),
                    row( mpIsogeometricGeometry->ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                    CurrentProcessInfo );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::InitializeMaterial()
    {
        KRATOS_TRY

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), (*mpIsogeometricGeometry),
                        row( mpIsogeometricGeometry->ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }
        }
        else
            KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
            KRATOS_CATCH( "" )
        }

    void TotalLagrangianIsogeometric::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
                mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), (*mpIsogeometricGeometry), row( mpIsogeometricGeometry->ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    inline void TotalLagrangianIsogeometric::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& rRightHandSideVector,
        double weight
    )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = mpIsogeometricGeometry->PointsNumber();
        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; ++j ) rRightHandSideVector[index + j] += weight * N[i] * BodyForce[j];
        }

        KRATOS_CATCH( "" )
    }



//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateAndAddKg(
        MatrixType& K,
        Matrix& DN_DX,
        Vector& StressVector,
        double weight )
    {
        KRATOS_TRY
        // unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();
        //Matrix<double> StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
        //Matrix<double> ReducedKg(DN_Dx.RowsNumber(),DN_Dx.RowsNumber());
        //Matrix<double>::MatMulAndAdd_B_D_Btrans(ReducedKg,weight,DN_Dx,StressTensor);
        //MathUtils<double>::ExpandAndAddReducedMatrix(K,ReducedKg,dimension);

        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();
        Matrix StressTensor = MathUtils<double>::StressVectorToTensor( StressVector );
        Matrix ReducedKg = prod( DN_DX, weight * Matrix( prod( StressTensor, trans( DN_DX ) ) ) ); //to be optimized
        MathUtils<double>::ExpandAndAddReducedMatrix( K, ReducedKg, dimension );

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateStrain(
        const Matrix& C,
        Vector& StrainVector )
    {
        KRATOS_TRY
        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();

        if ( dimension == 2 )
        {
            if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

            StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            StrainVector[2] = C( 0, 1 );
        }

        if ( dimension == 3 )
        {
            if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

            StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

            StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

            StrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

            StrainVector[3] = C( 0, 1 ); // xy

            StrainVector[4] = C( 1, 2 ); // yz

            StrainVector[5] = C( 0, 2 ); // xz
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateB(
        Matrix& B,
        Matrix& F,
        Matrix& DN_DX,
        unsigned int StrainSize )
    {
        KRATOS_TRY
        const unsigned int number_of_nodes = mpIsogeometricGeometry->PointsNumber();
        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();
        //
        //unsigned int dim2 = number_of_nodes*dimension;
        //if(B.size1() != StrainSize || B.size2()!=dim2)
        // B.resize(StrainSize,dim2);
        //Matrix Bi;

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = dimension * i;

            if ( dimension == 2 )
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
            }
            else
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 0, index + 2 ) = F( 2, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 1, index + 2 ) = F( 2, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 2 ) * DN_DX( i, 2 );
                B( 2, index + 1 ) = F( 1, 2 ) * DN_DX( i, 2 );
                B( 2, index + 2 ) = F( 2, 2 ) * DN_DX( i, 2 );
                B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
                B( 3, index + 2 ) = F( 2, 0 ) * DN_DX( i, 1 ) + F( 2, 1 ) * DN_DX( i, 0 );
                B( 4, index + 0 ) = F( 0, 1 ) * DN_DX( i, 2 ) + F( 0, 2 ) * DN_DX( i, 1 );
                B( 4, index + 1 ) = F( 1, 1 ) * DN_DX( i, 2 ) + F( 1, 2 ) * DN_DX( i, 1 );
                B( 4, index + 2 ) = F( 2, 1 ) * DN_DX( i, 2 ) + F( 2, 2 ) * DN_DX( i, 1 );
                B( 5, index + 0 ) = F( 0, 2 ) * DN_DX( i, 0 ) + F( 0, 0 ) * DN_DX( i, 2 );
                B( 5, index + 1 ) = F( 1, 2 ) * DN_DX( i, 0 ) + F( 1, 0 ) * DN_DX( i, 2 );
                B( 5, index + 2 ) = F( 2, 2 ) * DN_DX( i, 0 ) + F( 2, 0 ) * DN_DX( i, 2 );
            }

            //CalculateBi(Bi,F,DN_DX,i);
            //MathUtils<double>::WriteMatrix(B,Bi,0,index);
        }

        KRATOS_CATCH( "" )
    }



//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
    {
        int number_of_nodes = mpIsogeometricGeometry->size();
        int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int dim2 = number_of_nodes * dim;

        if ( rResult.size() != dim2 )
            rResult.resize( dim2, false );

        for ( int i = 0; i < number_of_nodes; ++i )
        {
            int index = i * dim;
            rResult[index] = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i )
        {
            ElementalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_Y ) );

            if ( mpIsogeometricGeometry->WorkingSpaceDimension() == 3 )
            {
                ElementalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_Z ) );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
//        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();
//        unsigned int NumberOfNodes = mpIsogeometricGeometry->size();
//        unsigned int MatSize = dimension * NumberOfNodes;
//
//        if ( rMassMatrix.size1() != MatSize )
//            rMassMatrix.resize( MatSize, MatSize, false );
//
//        rMassMatrix = ZeroMatrix( MatSize, MatSize );
//
//        double TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];
//
//        if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];
//
//        Vector LumpFact;
//
//        LumpFact = mpIsogeometricGeometry->LumpingFactors( LumpFact );
//
//        for ( unsigned int i = 0; i < NumberOfNodes; ++i )
//        {
//            double temp = LumpFact[i] * TotalMass;
//
//            for ( unsigned int j = 0; j < dimension; ++j )
//            {
//                unsigned int index = i * dimension + j;
//                rMassMatrix( index, index ) = temp;
//            }
//        }

        //consistent mass matrix
        unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        unsigned int mat_size = number_of_nodes * dim;

        if (rMassMatrix.size1() != mat_size)
        {
            rMassMatrix.resize(mat_size, mat_size, false);
            noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);
        }

        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        const GeometryType::IntegrationPointsArrayType& integration_points =
            mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        for(unsigned int Point = 0; Point < integration_points.size(); ++Point)
        {
            double IntToReferenceWeight = integration_points[Point].Weight();

            if (dim == 2)
                IntToReferenceWeight *= GetProperties()[THICKNESS];

            noalias(rMassMatrix) += outer_prod(row(Ncontainer, Point), row(Ncontainer, Point))
                                    * IntToReferenceWeight * mDetJ0[Point];
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

        if ( rDampMatrix.size1() != MatSize )
            rDampMatrix.resize( MatSize, MatSize, false );

        noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

    TotalLagrangianIsogeometric::IntegrationMethod TotalLagrangianIsogeometric::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( Output.size() != mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size() )
            Output.resize( mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size(), false );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        unsigned int StrainSize;

        if ( mpIsogeometricGeometry->WorkingSpaceDimension() == 2 ) StrainSize = 3;
        else StrainSize = 6;

        Vector StrainVector( StrainSize );

        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            {
                if ( Output[ii].size() != StrainVector.size() )
                    Output[ii].resize( StrainVector.size(), false );

                Output[ii] = mConstitutiveLawVector[ii]->GetValue( INSITU_STRESS, Output[ii] );
            }
        }
        else
        {
            if ( Output.size() != mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size() )
                Output.resize( mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size() );

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
                Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int StrainSize;

        if ( dim == 2 )
            StrainSize = 3;
        else
            StrainSize = 6;

        Matrix F( dim, dim );

        Matrix D( StrainSize, StrainSize );

        Matrix C( dim, dim );

        Vector StrainVector( StrainSize );

        Vector StressVector( StrainSize );

        Matrix DN_DX( number_of_nodes, dim );

        Matrix PlasticStrainVector( mpIsogeometricGeometry->size(), mpIsogeometricGeometry->WorkingSpaceDimension() );

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod );

        const Matrix& Ncontainer = mpIsogeometricGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        //calculating actual jacobian
        GeometryType::JacobiansType J;

        J = mpIsogeometricGeometry->Jacobian( J );

        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {

            //deformation gradient
            noalias( F ) = prod( J[PointNumber], mInvJ0[PointNumber] );

            //strain calculation
            noalias( C ) = prod( trans( F ), F );

            CalculateStrain( C, StrainVector );
            Comprobate_State_Vector( StrainVector );

            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                if ( Output[PointNumber].size2() != StrainVector.size() )
                    Output[PointNumber].resize( 1, StrainVector.size(), false );

                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                    Output[PointNumber]( 0, ii ) = StrainVector[ii];
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                if ( Output[PointNumber].size2() != StrainVector.size() )
                    Output[PointNumber].resize( 1, StrainVector.size(), false );

                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                    StrainVector,
                    F,
                    StressVector,
                    D,
                    rCurrentProcessInfo,
                    GetProperties(),
                    (*mpIsogeometricGeometry),
                    row( Ncontainer, PointNumber ),
                    true,
                    false,
                    false );

                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                {
                    Output[PointNumber]( 0, ii ) = StressVector[ii];
                }
            }
            else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
            {
                double size = StrainVector.size();
                PlasticStrainVector.resize( 1, size );

                if ( Output[PointNumber].size2() != StrainVector.size() )
                    Output[PointNumber].resize( 1, size, false );

                mConstitutiveLawVector[PointNumber]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, PlasticStrainVector );

                Output[PointNumber] = PlasticStrainVector;
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        //   std::cout << mConstitutiveLawVector[0] << std::endl;
        //  if (rVariable==INSITU_STRESS)
        //  {
        //                    for( unsigned int PointNumber = 0; PointNumber<mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size(); ++PointNumber )
        //   {
        //    mConstitutiveLawVector[PointNumber]->SetValue(INSITU_STRESS, rValues[PointNumber],
        //      rCurrentProcessInfo );
        //   }
        //  }
        //                if (rVariable==MATERIAL_PARAMETERS)
        //                {
        //                    for( unsigned int PointNumber = 0; PointNumber<mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size(); ++PointNumber )
        //                    {
        //                        mConstitutiveLawVector[PointNumber]->SetValue( MATERIAL_PARAMETERS,
        //                                rValues[PointNumber], rCurrentProcessInfo );
        //                    }
        //                }

        for ( unsigned int PointNumber = 0; PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }


//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        for ( unsigned int PointNumber = 0; PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size() )
            rValues.resize( mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size(), false );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
    }


//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int& size = mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size();

        if ( rValues.size() != size )
            rValues.resize( size );

        if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size();
                    ++PointNumber )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( PRESTRESS, rValues[PointNumber] );
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
        if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }
        if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( MATERIAL_PARAMETERS, rValues[PointNumber] );
            }
        }

        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size();
                    ++PointNumber )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );

            }
        }


    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        if ( rVariable == PK2_STRESS_TENSOR )
        {

            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }
    }


//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( VELOCITY_X, Step );
            values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dim == 3 )
                values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( VELOCITY_Z, Step );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
        const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            unsigned int index = i * dim;
            values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dim == 3 )
                values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
    {

        double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
        double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
        double c2 = 0.00; // norma de la velocidad actual dentro del elemento
        double c = 0.00;
        double wmax = 0.00;
        Vector Values( mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size() );
        Vector Velocities;

        GetFirstDerivativesVector( Velocities, 0 );

        if ( rVariable == DELTA_TIME )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < mpIsogeometricGeometry->IntegrationPoints( mThisIntegrationMethod ).size();
                    ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
                Values[PointNumber] = c1;
            }
        }

        c1 = ( *std::max_element( Values.begin(), Values.end() ) );

        c2 = norm_2( Velocities );

        c = ( c1 > c2 ) ? c1 : c2;


        double le = mpIsogeometricGeometry->Length();
        //KRATOS_WATCH(le)

        /// maxima frecuencia de un elemento
        wmax = ( lamda * c ) / le;
        Output = 2.0 / wmax;
        //KRATOS_WATCH(Output)

    }

//************************************************************************************
//************************************************************************************

    void TotalLagrangianIsogeometric::Comprobate_State_Vector( Vector& Result )
    {
        for ( unsigned int i = 0.00; i < Result.size(); ++i )
        {
            if ( fabs( Result( i ) ) < 1E-9 )
            {
                Result( i ) = 0.00;
            }
        }
    }


//************************************************************************************
//************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int  TotalLagrangianIsogeometric::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int dimension = this->mpIsogeometricGeometry->WorkingSpaceDimension();

        //verify that the variables are correctly initialized

        if ( VELOCITY.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

        if ( DISPLACEMENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

        if ( ACCELERATION.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

        if ( DENSITY.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );

        if ( BODY_FORCE.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

        if ( THICKNESS.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );

        //verify that the dofs exist
        for ( unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i )
        {
            if ( (*mpIsogeometricGeometry)[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", (*mpIsogeometricGeometry)[i].Id() );

            if ( (*mpIsogeometricGeometry)[i].HasDofFor( DISPLACEMENT_X ) == false || (*mpIsogeometricGeometry)[i].HasDofFor( DISPLACEMENT_Y ) == false || (*mpIsogeometricGeometry)[i].HasDofFor( DISPLACEMENT_Z ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", (*mpIsogeometricGeometry)[i].Id() );
        }

        //verify that the constitutive law exists
        if ( GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", GetProperties().Id() );
        }

        //Verify that the body force is defined
        if ( GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", GetProperties().Id() )
        }

        //verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
            if ( GetProperties().Has( THICKNESS ) == false )
                KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", Id() );

            if ( GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
                KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", Id() );
        }
        else
        {
            if ( GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
                KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", Id() );
        }

        //check constitutive law
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            return mConstitutiveLawVector[i]->Check( GetProperties(), (*mpIsogeometricGeometry), rCurrentProcessInfo );
        }

        //check if it is in the XY plane for 2D case


        return 0;

        KRATOS_CATCH( "" );
    }


    void TotalLagrangianIsogeometric::save( Serializer& rSerializer ) const
    {
//  std::cout << "Saving the TotalLagrangianIsogeometric #" << Id() << std::endl;
        rSerializer.save( "Name", "TotalLagrangianIsogeometric" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void TotalLagrangianIsogeometric::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
//  std::cout << "Loading the TotalLagrangianIsogeometric #" << Id() << std::endl;
        mThisIntegrationMethod = mpIsogeometricGeometry->GetDefaultIntegrationMethod();
    }





} // Namespace Kratos


