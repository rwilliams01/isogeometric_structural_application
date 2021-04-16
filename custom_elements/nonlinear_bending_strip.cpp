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
/* **************************************************************************************
*
*   Last Modified by:    $Author: G.D.Huynh $
*   Date:                $Date: 25.042018 $
*   Revision:            $Revision: 1.0 $
*
* ***************************************************************************************/


// System includes
// External includes
#include <boost/timer.hpp>
// Project includes
#include "includes/define.h"
#include "custom_elements/nonlinear_bending_strip.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

#define ENABLE_PROFILING

#include <iostream>
using namespace std;

namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
NonLinearBendingStrip::NonLinearBendingStrip()
{
}

// Constructor
NonLinearBendingStrip::NonLinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

// Constructor
NonLinearBendingStrip::NonLinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

//***********************************************************************************
//***********************************************************************************
Element::Pointer NonLinearBendingStrip::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new NonLinearBendingStrip( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Element::Pointer NonLinearBendingStrip::Create( IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties ) const
{
return Element::Pointer( new NonLinearBendingStrip( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
NonLinearBendingStrip::~NonLinearBendingStrip()
{
}

//***********************************************************************************
//***********************************************************************************
void NonLinearBendingStrip::EquationIdVector( EquationIdVectorType& rResult,
                                    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
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

    //KRATOS_WATCH(rResult)

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void NonLinearBendingStrip::GetDofList( DofsVectorType& ElementalDofList,
                              const ProcessInfo& rCurrentProcessInfo ) const
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************
void NonLinearBendingStrip::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

        mDim = 3;
        mNumberOfNodes = mpIsogeometricGeometry->size();
        mStrainSize = 3;
        mNumberOfDof = mNumberOfNodes*mDim;
        mThickness = GetProperties()[THICKNESS];

        // material parameters
        mE = GetProperties()[YOUNG_MODULUS];
        mNU = GetProperties()[POISSON_RATIO];
        mLambda = mE*mNU /(1.0 +mNU)/(1.0 - 2.0*mNU);
        mMu = 0.5*mE/(1.0 + mNU);

        ////////////////////////////////////////////////////////////////
        // get nodal coordinates vector R in undeformed configuration //
        if(mNodalCoordinates.size() != mNumberOfNodes)
            mNodalCoordinates.resize(mNumberOfNodes);

        for(unsigned int I=0; I< mNumberOfNodes; ++I)
            mNodalCoordinates[I].resize(mDim);

        for (unsigned int I=0; I< mNumberOfNodes; ++I)
        {
            mNodalCoordinates[I](0) = (*mpIsogeometricGeometry)[I].X0();
            mNodalCoordinates[I](1) = (*mpIsogeometricGeometry)[I].Y0();
            mNodalCoordinates[I](2) = (*mpIsogeometricGeometry)[I].Z0();
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        mpIsogeometricGeometry->Initialize(mThisIntegrationMethod);
        #endif

        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // law was at each integration point
        const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        mNumberOfIntegrationPoint = integration_points.size();

        //////////////////////////////////////////////////////////////
        ///// compute  Jacobian ,Inverse Jacobian, Det Jacobian /////
        mInvJ0.resize(integration_points.size());
        mN.resize(integration_points.size());
        mDN_De.resize(integration_points.size());

        for (unsigned int i = 0; i < integration_points.size(); ++i)
        {
            mInvJ0[i].resize(mDim, mDim, false);
            noalias(mInvJ0[i]) = ZeroMatrix(mDim, mDim);

            mN[i].resize(mNumberOfNodes);
            noalias(mN[i]) = ZeroVector(mNumberOfNodes);

            mDN_De[i].resize(mNumberOfNodes, 2);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes,2);

        }


        mDetJ0.resize(integration_points.size(), false);
        // TODO remove the storage for Jacobian to save memory
        noalias(mDetJ0) = ZeroVector(integration_points.size());

        mIntegrationWeight.resize(integration_points.size());



        // calculate the Jacobian
        mJ0.resize(integration_points.size());
        mJ0 = mpIsogeometricGeometry->Jacobian0(mJ0, mThisIntegrationMethod);
        double DetJ_temp;
        mTotalDomainInitialSize = 0.0;

        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++ PointNumber)
        {
            mN[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsValues( mN[PointNumber] , integration_points[PointNumber]);
            mDN_De[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsLocalGradients( mDN_De[PointNumber], integration_points[PointNumber]);

            MathUtils<double>::InvertMatrix(mJ0[PointNumber], mInvJ0[PointNumber],  DetJ_temp);

            Matrix JtJ = prod(trans(mJ0[PointNumber]), mJ0[PointNumber]);
            mDetJ0[PointNumber] = sqrt(MathUtils<double>::Det(JtJ));

            //getting informations for integration
            mIntegrationWeight[PointNumber] = integration_points[PointNumber].Weight();

            mTotalDomainInitialSize += mDetJ0[PointNumber]* mIntegrationWeight[PointNumber];


        }


        mIsInitialized = true;


        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void NonLinearBendingStrip::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
//***********************************************************************************
void NonLinearBendingStrip::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag,
                                bool CalculateResidualVectorFlag )
{
    KRATOS_TRY


    if (CalculateStiffnessMatrixFlag==true)
    {
        if (rLeftHandSideMatrix.size1() != mNumberOfDof)
            rLeftHandSideMatrix.resize(mNumberOfDof, mNumberOfDof);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mNumberOfDof, mNumberOfDof);
    }

    if (CalculateResidualVectorFlag==true)
    {
        if (rRightHandSideVector.size() != mNumberOfDof)
            rRightHandSideVector.resize(mNumberOfDof);
        noalias(rRightHandSideVector) = ZeroVector(mNumberOfDof);
    }

    ////////////////////////////////////////////////////////////////////////////////
    ////// compute residual vector and stiffness matrix over integration points/////

    #ifdef ENABLE_BEZIER_GEOMETRY
    //initialize the geometry
    mpIsogeometricGeometry->Initialize(mThisIntegrationMethod);
    #endif

    // get integration points
    const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

    // Current displacements
    Matrix CurrentDisplacement(mNumberOfNodes, mDim);
    for(unsigned int node =0; node < mpIsogeometricGeometry->size() ;++node)
        noalias(row(CurrentDisplacement, node)) = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

    ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// loop over integration points
    for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
    {

            //////////// get shape function values and their derivatives
            ShapeFunctionsSecondDerivativesType D2N_De2;
            D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);

            /////// i. covariant base vectors
            std::vector<Vector> u_a;
            std::vector<Vector> a;
            std::vector<Vector> A;  // covarianti base vector of undeformed configuration
            ReferenceCovariantBaseVector(A, mDN_De[PointNumber], mNodalCoordinates);
            FirstDerivativeDisplacement_a( u_a, mDN_De[PointNumber] , CurrentDisplacement);
            DeformedCovariantBaseVector( a, A, u_a );

            // normal directors
            Vector a3Vector, A3Vector, aa3Vector, AA3Vector; double a3, A3;
            NormalDirector(a3Vector,  aa3Vector, a3, a);
            NormalDirector(A3Vector, AA3Vector, A3, A);

            ////ii. derivative of covariant base vectors
            std::vector<std::vector<Vector> > u_ab;
            std::vector<std::vector<Vector> > a_ab;
            std::vector<std::vector<Vector> > A_ab;
            DerivativeReferenceCovariantBaseVector(A_ab, D2N_De2, mNodalCoordinates);
            SecondDerivativeDisplacement_ab( u_ab,  D2N_De2 , CurrentDisplacement);
            DerivativeDeformedCovariantBaseVector( a_ab , A_ab,  u_ab);

            std::vector<Vector> UnitBasisVector;
            this->UnitBaseVectors(UnitBasisVector);

            // a_ar
            std::vector< std::vector<std::vector<Vector>> > a_ar;
            DerivativeCovariantBaseVector_r( a_ar, mDN_De[PointNumber], UnitBasisVector);

            // aa3_rVector
            std::vector<std::vector<Vector>> aa3_rVector;
            DerivativeNonNormalizedDirector_r( aa3_rVector , a_ar, a);

            // a3_r
            std::vector<std::vector<double>> a3_r;
            DerivativeDirectorNorm_r(a3_r , aa3_rVector, a3Vector);

            // a3_rVector
            std::vector<std::vector<Vector>> a3_rVector;
            DerivativeDirector_r(a3_rVector ,  aa3_rVector, a3_r, a3Vector, a3);

            // a_abr
            std::vector< std::vector<std::vector<std::vector<Vector>>> > a_abr;
            DerivativeCovariantBaseVector_abr( a_abr, D2N_De2, UnitBasisVector);

            /*// aa3_rsVector
            std::vector<std::vector<std::vector<std::vector<Vector>>> > aa3_rsVector;
            SecondDerivativeNonNormalizedDirector_rs( aa3_rsVector , a_ar);

            // a3_rs
            std::vector<std::vector<std::vector<std::vector<double>>> > a3_rs;
            SecondDerivativeDirectorNorm_rs(a3_rs , aa3_rsVector,  aa3Vector, a3Vector,  a3, aa3_rVector);

            // a3_rsVector
            std::vector<std::vector<std::vector<std::vector<Vector>>> > a3_rsVector;
            SecondDerivativeDirector_rs( a3_rsVector,  aa3_rsVector, a3_rs, a3Vector, a3,  aa3_rVector, a3_r);
            */

            // metric and curvature coefficents
            Matrix Aab(2,2);
            Matrix aab(2,2);
            Matrix Bab(2,2);
            Matrix bab(2,2);
            CovariantMetricCoefficient(Aab, A);
            CovariantMetricCoefficient(aab, a);
            CovariantCurvatureCoefficient(Bab, A_ab, A3Vector);
            CovariantCurvatureCoefficient(bab, a_ab, a3Vector);

            double detJA = sqrt(MathUtils<double>::Det(Aab));
            // contravariant base vectors
            std::vector<Vector> AA;
            ContravariantBaseVector( AA, A, Aab);

            // local Cartesian basis
            std::vector<Vector> EE;
            LocalCartesianBasisVector(EE,  A, A3Vector);


            // transformation coeff
            std::vector< std::vector<std::vector<std::vector<double>>> > TransformationCoeff;
            LocalTransformationCoefficient(TransformationCoeff, EE, AA);

            // membrane and bending strains
            Matrix kTensor;
            computeCurvatureChange(kTensor,  Bab, bab);

            // transform strains to local form
            Matrix kkTensor = ZeroMatrix(2,2);
            Vector kkVector(3);
            LocalTransformationOfTensor(kkTensor , kTensor, TransformationCoeff);
            kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);

            // material matrix
            Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
            CalculateConstitutiveMatrix(D0);

            Vector moVector= ZeroVector(mStrainSize);
            noalias(moVector) = pow(mThickness,3.0)/12.0*prod(D0,kkVector) ;

            // first derivative of bending strains w.r.t displacement
            std::vector<std::vector<Vector> > kLocalVector_r;
            FirstDerivativeLocalCurvatureChange_r( kLocalVector_r,  a_abr, a_ab, a3Vector , a3_rVector, TransformationCoeff);

            // bending bending B matrix
            Matrix BBb(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBb, kLocalVector_r);

            ///////////////////////////////////////////////
            ///// nonlinear part of stiffness matrix //////
            ///////////////////////////////////////////////

            //std::vector<std::vector<Matrix>> kLocalVector_rs;
            //SecondDerivativeLocalCurvatureChange_rs(kLocalVector_rs , a_abr, a3_rVector,  a_ab, a3_rsVector, TransformationCoeff);


            if(CalculateStiffnessMatrixFlag == true)
            {
                AddLinearStiffnessMatrix(rLeftHandSideMatrix, D0, BBb, BBb, pow(mThickness,3.0)/12.0*detJA , mIntegrationWeight[PointNumber]);
                //AddNonlinearStiffnessMatrix(rLeftHandSideMatrix,moVector, kLocalVector_rs,  detJA , mIntegrationWeight[PointNumber]);
            }

            if(CalculateResidualVectorFlag == true)
            {
                AddInternalForces(rRightHandSideVector , moVector,kLocalVector_r, detJA , mIntegrationWeight[PointNumber]);

            }


    }// loop over integration points


    //KRATOS_WATCH(rRightHandSideVector)


    #ifdef ENABLE_BEZIER_GEOMETRY
    // clean the geometry internal data
    mpIsogeometricGeometry->Clean();
    #endif


    KRATOS_CATCH( "" )
}

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void NonLinearBendingStrip::AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
        std::vector<std::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight)
    {
       noalias(RightHandSideVector) -= StressResultants(0)*StrainVector_r[0][0]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(1)*StrainVector_r[1][1]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(2)*StrainVector_r[0][1]*DetJ*Weight;
    }

    void NonLinearBendingStrip::AddLinearStiffnessMatrix(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += prod( trans(BlhsMatrix), Matrix(prod(Di, BrhsMatrix)) )*DetJ*Weight;
    }

    void NonLinearBendingStrip::AddNonlinearStiffnessMatrix( MatrixType& LeftHandSideMatrix,const Vector& StressResultants
        , const std::vector<std::vector<Matrix>>& StrainVector_rs, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += StressResultants(0)*StrainVector_rs[0][0]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += StressResultants(1)*StrainVector_rs[1][1]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += StressResultants(2)*StrainVector_rs[0][1]*DetJ*Weight;
    }
    ////////////////////////////////////////////////////////
    ///////////////////// strain tensors ///////////////////
    ////////////////////////////////////////////////////////


    void NonLinearBendingStrip::computeCurvatureChange(Matrix& kTensor, Matrix& Bab, Matrix& bab)
    {
        kTensor.resize(2,2);
        kTensor = ZeroMatrix(2,2);

        noalias(kTensor) = Bab - bab;

    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////
    void NonLinearBendingStrip::DeformedCovariantBaseVector(std::vector<Vector>& a, std::vector<Vector>& A
                    , std::vector<Vector>& u_a )
    {

        a.resize(2);
        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            a[alpha].resize(mDim);
            a[alpha] = ZeroVector(mDim);
        }

        for(unsigned int alpha=0; alpha < 2; alpha++)
        {
            noalias(a[alpha]) = A[alpha] + u_a[alpha];
        }
    }


    void NonLinearBendingStrip::ReferenceCovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
                                                                            , const std::vector<Vector>& X)
    {

        // resize A
        A.resize(2);
        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            A[alpha].resize(mDim);
            A[alpha] = ZeroVector(mDim);
        }



        for(unsigned alpha=0; alpha < 2; alpha++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                for(unsigned int I=0; I<mNumberOfNodes; I++)
                {
                    // compute a1, a2
                    A[alpha](i) += DN_De(I,alpha)*X[I](i) ;
                }
            }
        }
    }


    void NonLinearBendingStrip::FirstDerivativeDisplacement_a(std::vector<Vector>& u_a, const Matrix& DN_De , const Matrix& u)
    {
        u_a.resize(2);

        for(unsigned int i=0; i< 2; i++)
        {
            u_a[i].resize(mDim);
            u_a[i] = ZeroVector(mDim);
        }

        for(unsigned alpha=0; alpha < 2; alpha++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                for(unsigned int I=0; I<mNumberOfNodes; I++)
                {
                    // compute a1, a2
                    u_a[alpha](i) += DN_De(I,alpha)*u(I,i) ;
                }
            }
        }
    }

    void NonLinearBendingStrip::SecondDerivativeDisplacement_ab(std::vector<std::vector<Vector> >& u_ab
        , const ShapeFunctionsSecondDerivativesType& D2N_De2 , const Matrix& u)
    {
        u_ab.resize(2);

        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            u_ab[alpha].resize(2);

            for(unsigned int beta=0; beta<2; beta++)
            {

                    u_ab[alpha][beta].resize(mDim);
                    u_ab[alpha][beta] = ZeroVector(mDim);

            }
        }

        for(unsigned alpha=0; alpha<2; alpha++)
        {
            for(unsigned beta=0; beta<2; beta++)
            {
                for(unsigned int i=0; i< mDim; ++i)
                {
                    for(unsigned int I=0; I< mNumberOfNodes; ++I)
                    {
                        // compute a_11
                        u_ab[alpha][beta](i) += D2N_De2[I](alpha,beta)*u(I,i);
                    }
                }
            }
        }
    }

    void NonLinearBendingStrip::NormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3,  std::vector<Vector>& a)
    {
        aa3Vector = MathUtils<double>::CrossProduct(a[0], a[1]);
        a3 = MathUtils<double>::Norm3(aa3Vector);
        a3Vector = aa3Vector/a3;
    }



    void NonLinearBendingStrip::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
    {
        AA.resize(2);
        for(unsigned alpha=0; alpha<2;alpha++)
        {
            AA[alpha].resize(3);
            AA[alpha] = ZeroVector(3);
        }


        double temp_ab;
        Matrix AAab(2,2);
        MathUtils<double>::InvertMatrix(Aab, AAab, temp_ab);

        noalias(AA[0]) = AAab(0,0)*A[0] +  AAab(0,1)*A[1] ;

        noalias(AA[1]) = AAab(1,0)*A[0] + AAab(1,1)*A[1] ;
    }

    void NonLinearBendingStrip::DerivativeReferenceCovariantBaseVector(std::vector<std::vector<Vector> >& A_ab,
        const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector>& X)
    {
        A_ab.resize(2);

        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            A_ab[alpha].resize(2);

            for(unsigned int beta=0; beta<2; beta++)
            {

                    A_ab[alpha][beta].resize(mDim);
                    A_ab[alpha][beta] = ZeroVector(mDim);

            }
        }

        for(unsigned alpha=0; alpha<2; alpha++)
        {
            for(unsigned beta=0; beta<2; beta++)
            {
                for(unsigned int i=0; i< mDim; ++i)
                {
                    for(unsigned int I=0; I< mNumberOfNodes; ++I)
                    {
                        // compute a_11
                        A_ab[alpha][beta](i) += D2N_De2[I](alpha,beta)*X[I](i) ;
                    }
                }
            }
        }

    }

    void NonLinearBendingStrip::DerivativeDeformedCovariantBaseVector(std::vector<std::vector<Vector> >& a_ab
        , std::vector<std::vector<Vector> >& A_ab
        , std::vector<std::vector<Vector> >& u_ab)
    {
        a_ab.resize(2);

        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            a_ab[alpha].resize(2);

            for(unsigned int beta=0; beta<2; beta++)
            {
                a_ab[alpha][beta].resize(mDim);
                a_ab[alpha][beta] = ZeroVector(mDim);
            }
        }

        for(unsigned alpha=0; alpha<2; alpha++)
        {
            for(unsigned beta=0; beta<2; beta++)
            {
                noalias(a_ab[alpha][beta]) = A_ab[alpha][beta] + u_ab[alpha][beta];
            }
        }

    }

    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////

    void NonLinearBendingStrip::CovariantCurvatureCoefficient(Matrix& Bab
        , std::vector<std::vector<Vector> >& A_ab, Vector& A3Vector)
    {
        Bab.resize(2,2);
        Bab = ZeroMatrix(2,2);

        for(unsigned int alpha=0; alpha<2; alpha++)
        {
            for(unsigned int beta=0; beta<2; beta++)
            {
                Bab(alpha,beta)= MathUtils<double>::Dot3(A_ab[alpha][beta], A3Vector);
            }
        }
    }

    void NonLinearBendingStrip::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
    {
        Aab.resize(2,2);
        Aab = ZeroMatrix(2,2);

        for(unsigned int alpha=0; alpha<2; alpha++)
        {
            for(unsigned int beta=0; beta<2; beta++)
            {
                Aab(alpha,beta)= MathUtils<double>::Dot3(A[alpha], A[beta]);
            }
        }
    }

    // second derivative of base vector w.r.t u_r
    void NonLinearBendingStrip::DerivativeCovariantBaseVector_abr( std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
        ,const ShapeFunctionsSecondDerivativesType& D2N_De2, std::vector<Vector>& UnitBasisVector)
    {
        a_abr.resize(2);
        for(unsigned int alpha = 0; alpha < 2; alpha++)
        {
            a_abr[alpha].resize(2);
            for(unsigned int beta = 0; beta < 2; beta++)
            {
                a_abr[alpha][beta].resize(mNumberOfNodes);
                for(unsigned int I = 0; I< mNumberOfNodes; I++)
                {
                    a_abr[alpha][beta][I].resize(mDim);
                    for(unsigned int i=0; i< mDim; i++)
                    {
                        a_abr[alpha][beta][I][i].resize(mDim);
                        a_abr[alpha][beta][I][i] = ZeroVector(mDim);
                    }
                }

            }
        }


        for(unsigned int alpha = 0; alpha < 2; alpha++)
        {
            for(unsigned int beta = 0; beta < 2; beta++)
            {
                for(unsigned int I=0; I < mNumberOfNodes; I++)
                {
                    for(unsigned int i=0; i<mDim; i++)
                    {
                        a_abr[alpha][beta][I][i] = D2N_De2[I](alpha,beta)*UnitBasisVector[i];
                    }
                }
            }
        }


    }

    void NonLinearBendingStrip::DerivativeCovariantBaseVector_r(std::vector< std::vector<std::vector<Vector>> >& a_ar, const Matrix& DN_De, std::vector<Vector>& UnitBasisVector)
    {
        a_ar.resize(2);
        for(unsigned int alpha = 0; alpha<2; alpha++)
        {
            a_ar[alpha].resize(mNumberOfNodes);
            for(unsigned int I=0; I< mNumberOfNodes; I++)
            {
                a_ar[alpha][I].resize(mDim);

                for(unsigned int i=0; i<mDim; i++)
                {
                    a_ar[alpha][I][i].resize(mDim);
                    a_ar[alpha][I][i] = ZeroVector(mDim);
                }
            }

        }


        for(unsigned int alpha = 0; alpha<2; alpha++)
        {
            for(unsigned int I=0; I < mNumberOfNodes; I++)
            {
                for(unsigned int i=0; i<mDim; i++)
                {
                    a_ar[alpha][I][i] = DN_De(I,alpha)*UnitBasisVector[i];
                }
            }
        }

    }
    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void NonLinearBendingStrip::DerivativeNonNormalizedDirector_r( std::vector<std::vector<Vector>>& aa3_rVector
        , std::vector< std::vector<std::vector<Vector>> >& a_ar, std::vector<Vector>& a)
    {
        aa3_rVector.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            aa3_rVector[I].resize(mDim);
            for(unsigned int i=0; i< mDim; i++)
            {
                aa3_rVector[I][i].resize(mDim);
                aa3_rVector[I][i] = ZeroVector(mDim);
            }
        }


        for(unsigned int I=0; I < mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                aa3_rVector[I][i] = MathUtils<double>::CrossProduct(a_ar[0][I][i] ,a[1]) + MathUtils<double>::CrossProduct(a[0], a_ar[1][I][i] );
            }
        }


    }

    void NonLinearBendingStrip::SecondDerivativeNonNormalizedDirector_rs(
          std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
        , std::vector< std::vector<std::vector<Vector>> >& a_ar)
    {
        aa3_rsVector.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            aa3_rsVector[I].resize(mNumberOfNodes);
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                aa3_rsVector[I][J].resize(mDim);
                for(unsigned int i=0; i< mDim; i++)
                {
                    aa3_rsVector[I][J][i].resize(mDim);
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        aa3_rsVector[I][J][i][j].resize(mDim);
                        aa3_rsVector[I][J][i][j] = ZeroVector(mDim);
                    }
                }
            }
        }


        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                for(unsigned int i=0; i< mDim; i++)
                {
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        aa3_rsVector[I][J][i][j] = MathUtils<double>::CrossProduct(a_ar[0][I][i], a_ar[1][J][j]) + MathUtils<double>::CrossProduct(a_ar[0][J][j], a_ar[1][I][i]);
                    }
                }
            }
        }


    }

    void NonLinearBendingStrip::DerivativeDirectorNorm_r(std::vector<std::vector<double>>& a3_r
        , std::vector<std::vector<Vector>>& aa3_rVector, Vector& a3Vector)
    {
        a3_r.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            a3_r[I].resize(mDim);
            for(unsigned int i=0; i< mDim; i++)
            {
                a3_r[I][i] = 0.0;
            }
        }


        for(unsigned int I=0; I < mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                a3_r[I][i] = MathUtils<double>::Dot3(a3Vector ,aa3_rVector[I][i]) ;
            }
        }


    }

    void NonLinearBendingStrip::SecondDerivativeDirectorNorm_rs(
          std::vector<std::vector<std::vector<std::vector<double>>> >& a3_rs
        , std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
        , Vector& aa3Vector, Vector& a3Vector, double& a3,  std::vector<std::vector<Vector>>& aa3_rVector)
    {
        a3_rs.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            a3_rs[I].resize(mNumberOfNodes);
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                a3_rs[I][J].resize(mDim);
                for(unsigned int i=0; i< mDim; i++)
                {
                    a3_rs[I][J][i].resize(mDim);
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        a3_rs[I][J][i][j]= 0.0;
                    }
                }
            }
        }


        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                for(unsigned int i=0; i< mDim; i++)
                {
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        a3_rs[I][J][i][j] = ( MathUtils<double>::Dot3(aa3_rsVector[I][J][i][j], aa3Vector)
                         + MathUtils<double>::Dot3(aa3_rVector[I][i], aa3_rVector[J][j])
                         - MathUtils<double>::Dot3(aa3_rVector[I][i], a3Vector)*MathUtils<double>::Dot3(aa3_rVector[J][j], a3Vector) )/a3;
                    }
                }
            }
        }

    }

    void NonLinearBendingStrip::DerivativeDirector_r( std::vector<std::vector<Vector>>& a3_rVector
        ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r, Vector& a3Vector, double& a3)
    {
        a3_rVector.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            a3_rVector[I].resize(mDim);
            for(unsigned int i=0; i< mDim; i++)
            {
                a3_rVector[I][i].resize(mDim);
                a3_rVector[I][i] = ZeroVector(mDim);
            }
        }


        for(unsigned int I=0; I < mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                a3_rVector[I][i] = (aa3_rVector[I][i] - a3_r[I][i]*a3Vector)/a3;
            }
        }


    }

    void NonLinearBendingStrip::SecondDerivativeDirector_rs(
          std::vector<std::vector<std::vector<std::vector<Vector>>> >& a3_rsVector
        , std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
        , std::vector<std::vector<std::vector<std::vector<double>>> >& a3_rs, Vector& a3Vector, double& a3
        ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r)
    {
        a3_rsVector.resize(mNumberOfNodes);
        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            a3_rsVector[I].resize(mNumberOfNodes);
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                a3_rsVector[I][J].resize(mDim);
                for(unsigned int i=0; i< mDim; i++)
                {
                    a3_rsVector[I][J][i].resize(mDim);
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        a3_rsVector[I][J][i][j].resize(mDim);
                        a3_rsVector[I][J][i][j]= ZeroVector(mDim);;
                    }
                }
            }
        }


        for(unsigned int I = 0; I< mNumberOfNodes; I++)
        {
            for(unsigned int J = 0; J< mNumberOfNodes; J++)
            {
                for(unsigned int i=0; i< mDim; i++)
                {
                    for(unsigned int j=0; j< mDim; j++)
                    {
                        a3_rsVector[I][J][i][j] = (aa3_rsVector[I][J][i][j] - a3_rs[I][J][i][j]*a3Vector)/a3
                            + (2.0*a3_r[I][i]*a3_r[J][j]*a3Vector - a3_r[I][i]*aa3_rVector[J][j] - a3_r[J][j]*aa3_rVector[I][i])*pow(a3,-2.0);
                    }
                }
            }
        }

    }

    void NonLinearBendingStrip::FirstDerivativeLocalCurvatureChange_r(std::vector<std::vector<Vector> >& kLocalVector_r
        , std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
        , std::vector<std::vector<Vector> >& a_ab
        , Vector& a3Vector , std::vector<std::vector<Vector>>& a3_rVector
        , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {
        kLocalVector_r.resize(2);
        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            kLocalVector_r[alpha].resize(2);

            for(unsigned int beta=0; beta<2; beta++)
            {
                kLocalVector_r[alpha][beta].resize(mNumberOfDof);
                kLocalVector_r[alpha][beta] = ZeroVector(mNumberOfDof);
            }
        }

        Vector temp(mNumberOfDof);
        for(unsigned int gamma=0; gamma<2;gamma++)
        {
            for(unsigned int delta=0; delta<2; delta++)
            {
                temp = ZeroVector(mNumberOfDof);

                for(unsigned int alpha=0; alpha<2;alpha++)
                {
                    for(unsigned int beta=0; beta<2; beta++)
                    {
                        for(unsigned int I=0; I< mNumberOfNodes; I++)
                        {
                            for(unsigned int i=0; i< mDim; i++)
                            {
                                temp(I*mDim + i) += ( MathUtils<double>::Dot3(a_abr[alpha][beta][I][i] ,a3Vector) + MathUtils<double>::Dot3(a_ab[alpha][beta], a3_rVector[I][i]) )
                                *TransformationCoeff[gamma][delta][alpha][beta] ;
                            }
                        }
                    }
                }

                noalias(kLocalVector_r[gamma][delta])= temp*(-1.0);
            }
        }

        kLocalVector_r[0][1] *= 2.0;
        kLocalVector_r[1][0] *= 2.0;
    }

    void NonLinearBendingStrip::SecondDerivativeLocalCurvatureChange_rs(std::vector<std::vector<Matrix> >& kLocalVector_rs
        ,std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr, std::vector<std::vector<Vector>>& a3_rVector
        , std::vector<std::vector<Vector> >& a_ab, std::vector<std::vector<std::vector<std::vector<Vector>>> >& a3_rsVector
        , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {
        if (kLocalVector_rs.size() != 2)
            kLocalVector_rs.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(kLocalVector_rs[i].size() !=2)
                kLocalVector_rs[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (kLocalVector_rs[i][j].size1() != mNumberOfDof)
                {
                    kLocalVector_rs[i][j].resize(mNumberOfDof, mNumberOfDof);
                    noalias(kLocalVector_rs[i][j]) = ZeroMatrix(mNumberOfDof, mNumberOfDof);
                }
            }
        }


        Matrix temp(mNumberOfDof, mNumberOfDof);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                temp = ZeroMatrix(mNumberOfDof, mNumberOfDof);

                for(unsigned int alpha=0; alpha<2;++alpha)
                {
                    for(unsigned int beta=0; beta<2; ++beta)
                    {
                        for(unsigned int I=0; I< mNumberOfNodes; ++I)
                        {
                            for(unsigned int J=0; J < mNumberOfNodes; ++J)
                            {
                                for(unsigned int i=0; i< mDim; ++i)
                                {
                                    for(unsigned int j=0; j< mDim; ++j)
                                    {
                                        temp(I*mDim+i, J*mDim+j) += ( MathUtils<double>::Dot3(a_abr[alpha][beta][I][i], a3_rVector[J][j])
                                           + MathUtils<double>::Dot3(a_abr[alpha][beta][J][j], a3_rVector[I][i])
                                           + MathUtils<double>::Dot3(a_ab[alpha][beta], a3_rsVector[I][J][i][j]) )
                                           *TransformationCoeff[gamma][delta][alpha][beta] ;

                                    }
                                }
                            }
                        }
                    }
                }

                noalias(kLocalVector_rs[gamma][delta]) = temp*(-1.0);
            }
        }

        kLocalVector_rs[0][1] *= 2.0;
        kLocalVector_rs[1][0] *= 2.0;
    }





    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////
    void NonLinearBendingStrip::CreatingBmatrix(Matrix& BMatrix, const std::vector<std::vector<Vector> >& LocalStrainVector_r)
    {
        if(BMatrix.size1() != mStrainSize)
            BMatrix.resize(mStrainSize, mNumberOfDof);
        noalias(BMatrix) = ZeroMatrix(mStrainSize, mNumberOfDof);

        for (int i=0; i< mNumberOfDof; i++)
        {
            BMatrix(0,i) = LocalStrainVector_r[0][0](i);
            BMatrix(1,i) = LocalStrainVector_r[1][1](i);
            BMatrix(2,i) = LocalStrainVector_r[0][1](i);
        }

    }


    void NonLinearBendingStrip::LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {
        if(T.size1()!=2)
            T.resize(2, 2);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                double temp=0.0;
                for(unsigned int alpha=0; alpha<2;++alpha)
                {
                    for(unsigned int beta=0; beta<2; ++beta)
                    {
                        temp += M(alpha,beta)*TransformationCoeff[gamma][delta][alpha][beta] ;
                    }
                }
                T(gamma,delta)= temp;
            }
        }
    }

    void NonLinearBendingStrip::UnitBaseVectors(std::vector<Vector>& e)
    {
            if (e.size() != 3)
                e.resize(3);

            e[0]=e[1]=e[2]=ZeroVector(mDim);
            e[0](0)=e[1](1)=e[2](2)= 1.0;
    }



    void NonLinearBendingStrip::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A, Vector& A3Vector)
    {
        EE.resize(3);
        for(unsigned int i=0; i<3; ++i)
        {
            EE[i].resize(3);
        }

        EE[0] = A[0]/MathUtils<double>::Norm3(A[0]);

        Vector EE2_temp = ZeroVector(3);
        noalias(EE2_temp) +=  A[1];
        noalias(EE2_temp) -=  MathUtils<double>::Dot3(A[1], EE[0])*EE[0];
        noalias(EE[1]) = ( EE2_temp )/MathUtils<double>::Norm3(EE2_temp);

        EE[2] = A3Vector;


    }



    void NonLinearBendingStrip::LocalTransformationCoefficient(std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff, std::vector<Vector>& EE, std::vector<Vector>& AA)
    {
        TransformationCoeff.resize(2);
        for(unsigned int gamma=0; gamma <2; gamma++)
        {
            TransformationCoeff[gamma].resize(2);
            for(unsigned int delta = 0; delta<2; delta++)
            {
                TransformationCoeff[gamma][delta].resize(2);
                for(unsigned int alpha=0; alpha< 2; alpha++)
                {
                    TransformationCoeff[gamma][delta][alpha].resize(2);

                    for(unsigned int beta=0; beta<2; beta++)
                    {
                        TransformationCoeff[gamma][delta][alpha][beta] = 0.0;
                    }
                }

            }
        }


        for(unsigned int gamma=0; gamma <2; gamma++)
        {
            for(unsigned int delta = 0; delta<2; delta++)
            {
                for(unsigned int alpha=0; alpha< 2; alpha++)
                {
                    for(unsigned int beta=0; beta<2; beta++)
                    {
                        TransformationCoeff[gamma][delta][alpha][beta] = MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                    }
                }

            }
        }

    }

    ///////////// material matrix
    void NonLinearBendingStrip::CalculateConstitutiveMatrix(Matrix& D)
    {
        if (D.size1() != 3)
            D.resize(3,3);
        D = ZeroMatrix(3,3);

        D(1,1) = pow(10.0 ,4.0)*mE;
    }



} // Namespace Kratos.
