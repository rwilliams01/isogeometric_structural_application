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
#include "custom_elements/linear_bending_strip.h"
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
LinearBendingStrip::LinearBendingStrip()
{
}

// Constructor
LinearBendingStrip::LinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

// Constructor
LinearBendingStrip::LinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

//***********************************************************************************
//***********************************************************************************
Element::Pointer LinearBendingStrip::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new LinearBendingStrip( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Element::Pointer LinearBendingStrip::Create( IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties ) const
{
return Element::Pointer( new LinearBendingStrip( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
LinearBendingStrip::~LinearBendingStrip()
{
}

//***********************************************************************************
//***********************************************************************************
void LinearBendingStrip::EquationIdVector( EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo )
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
void LinearBendingStrip::GetDofList( DofsVectorType& ElementalDofList,
                              ProcessInfo& rCurrentProcessInfo )
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
void LinearBendingStrip::Initialize()
{
    KRATOS_TRY
    

        ////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////
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
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

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
void LinearBendingStrip::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
//***********************************************************************************
void LinearBendingStrip::CalculateAll( MatrixType& rLeftHandSideMatrix,
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

    Vector NodalDisplacement = ZeroVector(mNumberOfDof);
    for(unsigned int I=0; I<mNumberOfNodes; I++)
    {
        for(unsigned int i=0; i<mDim; i++)
        {
            NodalDisplacement(I*mDim +i) = CurrentDisplacement(I,i);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////// loop over integration points
    for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
    {

            //////////// get shape function values and their derivatives
            ShapeFunctionsSecondDerivativesType D2N_De2;
            D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    

            /////// i. covariant base vectors
            std::vector<Vector> A;  // covarianti base vector of undeformed configuration
            CovariantBaseVector( A, mDN_De[PointNumber], mNodalCoordinates);

            ////i. normal vector
            Vector A3Vector, AA3Vector; double A3;
            ReferenceNormalDirector( A3Vector,  AA3Vector, A3, A);

            ////ii. derivative of covariant base vectors
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > A_ab;
            DerivativeCovariantBaseVector( A_ab, D2N_De2, mNodalCoordinates);

            // metric and curvature coefficents
            Matrix Aab(2,2);
            CovariantMetricCoefficient( Aab, A);

            double detJA = sqrt(MathUtils<double>::Det(Aab));

            // contravariant base vectors
            std::vector<Vector> AA;
            ContravariantBaseVector( AA, A, Aab);

            /////////////////////////////////////////////
            ///// local Cartesian basis /////////////////
            std::vector<Vector> EE;
            LocalCartesianBasisVector(EE,  A);

            // first derivative of curvature change
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > kLocalVector_r;
            FirstDerivativeLocalCurvatureChange_r(kLocalVector_r, A, A3, A3Vector, mDN_De[PointNumber], D2N_De2,A_ab, EE, AA);


            // bending B matrix
            Matrix BBb(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBb, kLocalVector_r);


            // material matrix
            Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
            CalculateConstitutiveMatrix(D0);

            //////////////////////////////////
            // compute the contribution to LHS
            Matrix GaussStiffnessMatrix = ZeroMatrix(mNumberOfDof,mNumberOfDof);
            AddStiffnessMatrixComponents(GaussStiffnessMatrix, D0, BBb, BBb,  pow(mThickness,3.0)/12.0*detJA , mIntegrationWeight[PointNumber]);
            noalias(rLeftHandSideMatrix) += GaussStiffnessMatrix;

            //////////////////////////////////
            // compute the contribution to RHS
            noalias(rRightHandSideVector) -= prod(GaussStiffnessMatrix, NodalDisplacement);

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
    void LinearBendingStrip::AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight)
    {
       noalias(RightHandSideVector) -= StressResultants(0)*StrainVector_r[0][0]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(1)*StrainVector_r[1][1]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(2)*StrainVector_r[0][1]*DetJ*Weight;
    }

    void LinearBendingStrip::AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight )
    {
        const Vector BodyForce = GetProperties()[BODY_FORCE];

        for(unsigned int I=0; I < mNumberOfNodes; ++I)
            for(unsigned int i=0; i<mDim; ++i)
                RightHandSideVector(I*mDim + i) += N(I)*BodyForce(i)*DetJ*Weight;
    }

    void LinearBendingStrip::AddStiffnessMatrixComponents(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += prod( trans(BlhsMatrix), Matrix(prod(Di, BrhsMatrix)) )*DetJ*Weight;
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////
    void LinearBendingStrip::CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
                                                                            , const std::vector<Vector>& X)
    {
        A.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            A[i].resize(mDim);
            A[i] = ZeroVector(mDim);
        }

        for(unsigned alpha=0; alpha < 2; ++alpha)
        {
            for(unsigned int i=0; i<mDim; ++i)
            {
                for(unsigned int I=0; I<mNumberOfNodes; ++I)
                {
                    // compute a1, a2
                    A[alpha](i) += DN_De(I,alpha)*X[I](i) ;
                }
            }
        }
        
    }



    void LinearBendingStrip::ReferenceNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A)
    {
        AA3Vector = MathUtils<double>::CrossProduct(A[0], A[1]);
        A3 = MathUtils<double>::Norm3(AA3Vector);
        A3Vector = AA3Vector/A3;
    }

    void LinearBendingStrip::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
    {
        AA.resize(2);
        for(unsigned i=0; i<2;++i)
        {
            AA[i].resize(3);
            AA[i] = ZeroVector(3);
        }


        double temp_ab;
        Matrix AAab(2,2);
        MathUtils<double>::InvertMatrix(Aab, AAab, temp_ab);

        noalias(AA[0]) += AAab(0,0)*A[0] ;
        noalias(AA[0]) += AAab(0,1)*A[1] ;

        noalias(AA[1]) += AAab(1,0)*A[0] ;
        noalias(AA[1]) += AAab(1,1)*A[1] ;
    }


    void LinearBendingStrip::DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab, 
        const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector>& X)
    {
        A_ab.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            A_ab[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                    A_ab[i][j].resize(mDim);
                    A_ab[i][j] = ZeroVector(mDim);
            }
        }

        for(unsigned alpha=0; alpha<2; ++alpha)
        {
            for(unsigned beta=0; beta<2; ++beta)
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

    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void LinearBendingStrip::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
    {
        if(Aab.size1()!=2)
            Aab.resize(2,2);

        for(unsigned int alpha=0; alpha<2; ++alpha)
        {
            for(unsigned int beta=0; beta<2; ++beta)
            {
                Aab(alpha,beta)= MathUtils<double>::Dot3(A[alpha], A[beta]);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void LinearBendingStrip::FirstDerivativeLocalCurvatureChange_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& kLocalVector_r
        , std::vector<Vector>& A, double& A3, Vector& A3Vector, const Matrix& DN_De
        , const ShapeFunctionsSecondDerivativesType& D2N_De2, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
        ,std::vector<Vector>& EE, std::vector<Vector>& AA)
    {
        if (kLocalVector_r.size() != 2)
            kLocalVector_r.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(kLocalVector_r[i].size() !=2)
                kLocalVector_r[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (kLocalVector_r[i][j].size() != mNumberOfDof)
                {
                    kLocalVector_r[i][j].resize(mNumberOfDof);
                    kLocalVector_r[i][j] = ZeroVector(mNumberOfDof);
                }
            }
        }


        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                Vector temp = ZeroVector(mNumberOfDof);

                for(unsigned int alpha=0; alpha<2;++alpha)
                {
                    for(unsigned int beta=0; beta<2; ++beta)
                    {
                        for(unsigned int I=0; I< mNumberOfNodes; ++I)
                        {
                            for(unsigned int i=0; i< mDim; ++i)
                            {
                                Vector A_abxA2(mDim);
                                A_abxA2 = MathUtils<double>::CrossProduct(A_ab[alpha][beta], A[1]);
                
                                Vector A1xA_ab(mDim);
                                A1xA_ab = MathUtils<double>::CrossProduct(A[0], A_ab[alpha][beta]);
                
                                Vector A2xA3(mDim);
                                A2xA3 = MathUtils<double>::CrossProduct(A[1], A3Vector);
                
                                Vector A3xA1(mDim);
                                A3xA1 = MathUtils<double>::CrossProduct(A3Vector, A[0]);
                
                                temp(I*mDim + i) +=  (- A3Vector(i)*D2N_De2[I](alpha,beta) + (1/A3)*( A_abxA2(i)*DN_De(I,0) 
                                           + A1xA_ab(i)*DN_De(I,1) + MathUtils<double>::Dot3(A_ab[alpha][beta], A3Vector)*
                                             ( A2xA3(i)*DN_De(I,0) + A3xA1(i)*DN_De(I,1) ) ) )*
                                             MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]); 
                            }
                        }
                    }
                }

                noalias(kLocalVector_r[gamma][delta])= temp;
            }
        }

        kLocalVector_r[0][1] *= 2.0;
        kLocalVector_r[1][0] *= 2.0;
    }


 
    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////
    void LinearBendingStrip::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A)
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

        EE[2] = MathUtils<double>::CrossProduct(EE[0], EE[1]); 

        
    }


    void LinearBendingStrip::UnitBaseVectors(std::vector<Vector>& e)
    {
        if (e.size() != 3)
            e.resize(3);

        e[0]=e[1]=e[2]=ZeroVector(mDim);
        e[0](0)=e[1](1)=e[2](2)= 1.0;
    }

    void LinearBendingStrip::CreatingBmatrix(Matrix& BMatrix, const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& LocalStrainVector_r)
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

    ///////////// material matrix
    void LinearBendingStrip::CalculateConstitutiveMatrix(Matrix& D)
    {
        if (D.size1() != 3)
            D.resize(3,3);
        D = ZeroMatrix(3,3);

        D(1,1) = pow(10.0 ,4.0)*mE;
    }



} // Namespace Kratos.
