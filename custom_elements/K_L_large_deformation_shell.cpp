//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: DongGiang $
//   Date:                $Date: 25 August 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes 


// External includes 

// Project includes 
#include "custom_elements/K_L_large_deformation_shell.h"
#include "isogeometric_structural_application/isogeometric_structural_application.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "utilities/math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"
#include "phase_field_application/phase_field_application.h"
#include "phase_field_application/custom_utilities/eig/eig3.h"
#include <cmath>
#include <iostream>
using namespace std;



#define ENABLE_BEZIER_GEOMETRY

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_LEVEL3
//#define DEBUG_LEVEL4
//#define DEBUG_LEVEL5
//#define DEBUG_LEVEL6



//#define DEBUG_DKGQ

namespace Kratos
{
    //extern Variable<Vector> STRESSES;
    extern Variable<array_1d<double, 3> > PRESCRIBED_DELTA_DISPLACEMENT;
    //************************************************************************************
    //***** Constructor and Destructor ***************************************************
    //************************************************************************************
    KirchhoffLoveLargeDeformationShell::KirchhoffLoveLargeDeformationShell
            (IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
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

    KirchhoffLoveLargeDeformationShell::KirchhoffLoveLargeDeformationShell(IndexType NewId,
                    GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
                                                    Element(NewId, pGeometry, pProperties)
    {
        mIsInitialized = false;
        //    mpIsogeometricGeometry = 
        //        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
        mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
    }

    /**
    * Destructor. Never to be called manually
    */
    KirchhoffLoveLargeDeformationShell::~KirchhoffLoveLargeDeformationShell()
    {
    }


    //********************************************************
    //**** Operations ****************************************
    //********************************************************

    Element::Pointer KirchhoffLoveLargeDeformationShell::Create(IndexType NewId,
        NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KirchhoffLoveLargeDeformationShell(NewId,
                                GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer KirchhoffLoveLargeDeformationShell::Create(IndexType NewId,
        GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KirchhoffLoveLargeDeformationShell(NewId,
                                                            pGeom, pProperties));
    }


        /////////////////////////////////////////////////////////////////////////////
    //************************************************************************************
    /**
    * returns the used integration method
    */
    KirchhoffLoveLargeDeformationShell::IntegrationMethod KirchhoffLoveLargeDeformationShell::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

        //************************************************************************************
    //************************************************************************************
    /**
    * Setting up the EquationIdVector
    */
    void KirchhoffLoveLargeDeformationShell::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        DofsVectorType ElementalDofList;
        GetDofList(ElementalDofList, rCurrentProcessInfo);

        if (rResult.size() != ElementalDofList.size())
            rResult.resize(ElementalDofList.size(), false);

        for(unsigned int i = 0; i < ElementalDofList.size(); ++i)
            rResult[i] = ElementalDofList[i]->EquationId();
    }

    void KirchhoffLoveLargeDeformationShell::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    


    void KirchhoffLoveLargeDeformationShell::Initialize()
    {
        KRATOS_TRY //EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

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
             and this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
             and this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
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
        mKappa = GetProperties()[KAPPA];
                                
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

            mDN_De[i].resize(mNumberOfNodes, mDim);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes, mDim);

        }
        

        mDetJ0.resize(integration_points.size(), false);        
        // TODO remove the storage for Jacobian to save memory
        noalias(mDetJ0) = ZeroVector(integration_points.size());

        mIntegrationWeight.resize(integration_points.size());


        // calculate the Jacobian
        mJ0.resize(integration_points.size());
        mJ0 = mpIsogeometricGeometry->Jacobian0(mJ0, mThisIntegrationMethod);
        double DetJ_temp;
     
        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++ PointNumber)
        {       
            mN[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsValues( mN[PointNumber] , integration_points[PointNumber]);
            mDN_De[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsLocalGradients( mDN_De[PointNumber], integration_points[PointNumber]);

            MathUtils<double>::InvertMatrix(mJ0[PointNumber], mInvJ0[PointNumber],  DetJ_temp);

            Matrix JtJ = prod(trans(mJ0[PointNumber]), mJ0[PointNumber]);
            mDetJ0[PointNumber] = sqrt(MathUtils<double>::Det(JtJ));

            //getting informations for integration
            mIntegrationWeight[PointNumber] = integration_points[PointNumber].Weight();
        }

        // integration points and weights used for thickness calculation
        mIntegrationPoint1D.resize(4);
        mIntegrationPoint1D(0) = -0.861136311594053 ;
        mIntegrationPoint1D(1) = -0.339981043584856;
        mIntegrationPoint1D(2) = 0.339981043584856;
        mIntegrationPoint1D(3) = 0.861136311594053;
               
        mWeight1D.resize(4);
        mWeight1D(0) = mWeight1D(3) = 0.347854845137454; 
        mWeight1D(1) = mWeight1D(2) = 0.652145154862546; 
        mDetJ1D = 0.5*mThickness;
   


        mIsInitialized = true;

        mIsIsotropicMaterial = true;
   


        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    void KirchhoffLoveLargeDeformationShell::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        // reset all resistant forces at node
        for ( unsigned int i = 0; i < mNumberOfNodes; ++i )
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
    void KirchhoffLoveLargeDeformationShell::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType temp = Matrix();

        bool need_calculate_stiffness = false;
        for ( unsigned int node = 0; node < mNumberOfNodes; ++node )
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
    void KirchhoffLoveLargeDeformationShell::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        bool MaterialUpdateFlag = false;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, MaterialUpdateFlag);

    }


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////// private subroutines /////////////////////////////////

    void KirchhoffLoveLargeDeformationShell::CalculateAll( MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag,
        bool MaterialUpdateFlag)
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

        Vector RightHandSideVector_trial = ZeroVector(mNumberOfDof);

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

        if(mIsIsotropicMaterial == true)
        {
            ///////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////// loop over integration points
            for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
            {
                //////////// get shape function values and their derivatives
                ShapeFunctionsSecondDerivativesType D2N_De2;
                D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    

                /////// i. covariant base vectors
                std::vector<Vector> a;
                std::vector<Vector> A;  // covarianti base vector of undeformed configuration
                CovariantBaseVector( A, mDN_De[PointNumber], mNodalCoordinates);
                CovariantBaseVector(a,mDN_De[PointNumber], mNodalCoordinates, CurrentDisplacement);


                // normal directors
                Vector a3Vector, A3Vector, aa3Vector, AA3Vector; double a3, A3;
                ReferencedNormalDirector( A3Vector,  AA3Vector, A3, A);
                DeformedNormalDirector(a3Vector, aa3Vector, a3, a);

                ////ii. derivative of covariant base vectors
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > A_ab;
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > a_ab;
                DerivativeCovariantBaseVector( A_ab, D2N_De2, mNodalCoordinates);
                DerivativeCovariantBaseVector( a_ab, D2N_De2, mNodalCoordinates, CurrentDisplacement);     

        
                // metric coefficents
                Matrix Aab(2,2);
                CovariantMetricCoefficient( Aab, A);
        
                // contravariant base vectors
                std::vector<Vector> AA;
                ContravariantBaseVector( AA, A, Aab);

                // local Cartesian basis
                std::vector<Vector> EE;
                LocalCartesianBasisVector(EE,  A);

                // local membrane strain
                Matrix eeTensor;
                computeMembraneStrain(eeTensor, a, A,  EE, AA);

                Vector eeVector(3);
                IsotropicTensorUtility<2>::StressTensorToVector(eeTensor, eeVector);
            
                // local curvature changes
                Matrix kkTensor;
                computeCurvatureChange(kkTensor,  A_ab, a_ab,a3Vector,A3Vector,  EE, AA );

                Vector kkVector(3);
                IsotropicTensorUtility<2>::StressTensorToVector(kkTensor, kkVector);


                // first derivative of membrane strains w.r.t displacement
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > eLocalVector_r;
                FirstDerivativeLocalMembraneStrain_r( eLocalVector_r, mDN_De[PointNumber],  a, EE,  AA);

                // membrane B matrix
                Matrix BBm(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBm, eLocalVector_r);
    
                // first derivative of bending strains w.r.t displacement
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > kLocalVector_r;
                FirstDerivativeLocalCurvatureChange_r(kLocalVector_r, a, a3Vector, aa3Vector,  a3 , a_ab, mDN_De[PointNumber],  D2N_De2,  EE,  AA);

                // bending bending B matrix
                Matrix BBb(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBb, kLocalVector_r);

                // material matrix
                Matrix TanC(3,3);
                CalculateElasticMatrix(TanC, mE,  mNU);

                // stress resultant
                Vector nVector(mStrainSize);
                computeNormalForces(nVector, TanC,  eeVector);
                    
                Vector moVector(mStrainSize);
                computeBendingMoments(moVector, TanC,  kkVector); 


                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>> eLocalVector_rs;
                SecondDerivativeLocalMembraneStrain_rs(eLocalVector_rs, mDN_De[PointNumber], EE, AA);

                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>> kLocalVector_rs;
                SecondDerivativeLocalCurvatureChange_rs(kLocalVector_rs, mDN_De[PointNumber], D2N_De2, a_ab, a, a3Vector, aa3Vector,  a3, EE,  AA);

                //KRATOS_WATCH(eLocalVector_rs)
                //KRATOS_WATCH(kLocalVector_rs)


                if(CalculateStiffnessMatrixFlag == true)
                {
                    AddLinearMembraneStiffness(rLeftHandSideMatrix, TanC,  BBm, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddLinearBendingStiffness(rLeftHandSideMatrix, TanC, BBb, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);

                    AddNonlinearMembraneStiffness(rLeftHandSideMatrix, nVector, eLocalVector_rs,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddNonlinearBendingStiffness(rLeftHandSideMatrix, moVector, kLocalVector_rs,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                }

                if(CalculateResidualVectorFlag == true)
                {
                
                    AddInternalForces(rRightHandSideVector , nVector,eLocalVector_r, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddInternalForces(rRightHandSideVector , moVector,kLocalVector_r, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);

                    ////////// add external forces to RHS
                    AddExternalForces(rRightHandSideVector, mN[PointNumber], mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                }

                      

            }// loop over integration points

            //KRATOS_WATCH(rRightHandSideVector)

            //KRATOS_WATCH(rLeftHandSideMatrix)

        }

        KRATOS_CATCH("")
    }

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution

    void KirchhoffLoveLargeDeformationShell::AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
         boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight)
    {
        noalias(RightHandSideVector) -= StressResultants(0)*StrainVector_r[0][0]*DetJ*Weight;
        noalias(RightHandSideVector) -= StressResultants(1)*StrainVector_r[1][1]*DetJ*Weight;
        noalias(RightHandSideVector) -= StressResultants(2)*StrainVector_r[0][1]*DetJ*Weight;
    }

    void KirchhoffLoveLargeDeformationShell::AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight )
    {
        const Vector BodyForce = GetProperties()[BODY_FORCE];

        for(unsigned int I=0; I < mNumberOfNodes; ++I)
            for(unsigned int i=0; i<mDim; ++i)
                RightHandSideVector(I*mDim + i) += N(I)*BodyForce(i)*DetJ*Weight;
    }

    void KirchhoffLoveLargeDeformationShell::AddLinearMembraneStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC
                                                , const Matrix& Bm, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) +=  mThickness*DetJ*Weight*prod( trans(Bm), Matrix(prod(TanC, Bm)) );
    }

    void KirchhoffLoveLargeDeformationShell::AddLinearBendingStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC
                                                , const Matrix& Bb,  const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += pow(mThickness,3)/12*DetJ*Weight*prod( trans(Bb), Matrix(prod(TanC, Bb)) );
    }

    void KirchhoffLoveLargeDeformationShell::AddNonlinearMembraneStiffness( MatrixType& LeftHandSideMatrix,const Vector& nVector
        , const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eVector_rs, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += nVector(0)*eVector_rs[0][0]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += nVector(1)*eVector_rs[1][1]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += nVector(2)*eVector_rs[0][1]*DetJ*Weight;
    }


    void KirchhoffLoveLargeDeformationShell::AddNonlinearBendingStiffness( MatrixType& LeftHandSideMatrix,const Vector& moVector
        , const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kVector_rs, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += moVector(0)*kVector_rs[0][0]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += moVector(1)*kVector_rs[1][1]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += moVector(2)*kVector_rs[0][1]*DetJ*Weight;
    }


    ////////////////////////////////////////////////////////////////////////////////////////// end 
    

    ///////////////////// normal forces
    void KirchhoffLoveLargeDeformationShell::computeNormalForces(Vector& nVector, const Matrix& C, const Vector& eVector)
    {
        nVector = ZeroVector(mStrainSize);
        noalias(nVector) = prod(C, mThickness*eVector);
    }

    void KirchhoffLoveLargeDeformationShell::computeBendingMoments(Vector& moVector, const Matrix& C, const Vector& kVector)
    {
        moVector = ZeroVector(mStrainSize);
        noalias(moVector) = prod(C, pow(mThickness,3)/12*kVector);
    } 


    ///////////////////// strain vectors
    
    void KirchhoffLoveLargeDeformationShell::computeMembraneStrain(Matrix& eTensor,  std::vector<Vector>& a,  std::vector<Vector>& A, std::vector<Vector>& EE, std::vector<Vector>& AA )
    {
        eTensor.resize(2,2);
        eTensor = ZeroMatrix(2,2);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                double temp=0;
                for(int alpha =0 ; alpha<2; alpha++)
                {
                    for(int beta=0; beta<2; beta++)
                    {
                        temp += ( 0.5*( MathUtils<double>::Dot3(a[alpha],a[beta]) - MathUtils<double>::Dot3(A[alpha], A[beta]) ) )
                                *MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                    }
                }

                eTensor(gamma,delta) = temp;
            }
        }

        eTensor(0,1) *=2;
        eTensor(1,0) *=2;

    }


    void KirchhoffLoveLargeDeformationShell::computeCurvatureChange(Matrix& kTensor,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
        ,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& a_ab,  Vector& a3Vector,  Vector& A3Vector, std::vector<Vector>& EE, std::vector<Vector>& AA)
    {
        kTensor.resize(2,2);
        kTensor = ZeroMatrix(2,2);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                double temp=0;
                for(int alpha =0 ; alpha<2; alpha++)
                {
                    for(int beta=0; beta<2; beta++)
                    {
                        temp += ( MathUtils<double>::Dot3(A_ab[alpha][beta], A3Vector) - MathUtils<double>::Dot3(a_ab[alpha][beta], a3Vector) ) 
                                *MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                    }
                }

                kTensor(gamma,delta) = temp;
            }
        }

        kTensor(0,1) *=2;
        kTensor(1,0) *=2;

    }



    /////////////////////////////////////////////////////////////////////////////////////
    double KirchhoffLoveLargeDeformationShell::KroneckerDelta(int i, int j)
    {
        return (i=j) ? 1.0 : 0.0 ;
    }

    void KirchhoffLoveLargeDeformationShell::UnitBaseVectors(std::vector<Vector>& e)
    {
        if (e.size() != 3)
            e.resize(3);

        e[0]=e[1]=e[2]=ZeroVector(mDim);
        e[0](0)=e[1](1)=e[2](2)= 1.0;
    }


    void KirchhoffLoveLargeDeformationShell::CalculateElasticMatrix(Matrix& C, const double& E, const double& NU)
    {
        unsigned int strain_size = 3;
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size);

            double c1 = E * ( 1.00 - pow(NU,2) );
            double c2 = E * ( 1.00 - pow(NU,2) )*NU;
            double c3 = 0.5 * E / ( 1 + NU );

            C( 0, 0 ) = c1;
            C( 0, 1 ) = c2;
            C( 0, 2 ) = 0.0;
            C( 1, 0 ) = c2;
            C( 1, 1 ) = c1;
            C( 1, 2 ) = 0.0;
            C( 2, 0 ) = 0.0;
            C( 2, 1 ) = 0.0;
            C( 2, 2 ) = c3;
     
    }


    ////////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLargeDeformationShell::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A)
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


    void KirchhoffLoveLargeDeformationShell::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
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

    void KirchhoffLoveLargeDeformationShell::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
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

    ////////////////////////////////////////////////////////////////////
    //////////////////////////// membrane part /////////////////////////

    //////////////// first derivatives
    // covariant base vectors
    void KirchhoffLoveLargeDeformationShell::CovariantBaseVector(std::vector<Vector>& a, const Matrix& DN_De
                    , const std::vector<Vector>& X, const Matrix& u )
    {
        if (a.size() != 2)
            a.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if (a[i].size() != mDim)
            {
                a[i].resize(mDim);
                a[i] = ZeroVector(mDim);
            }
        }

        for(unsigned int alpha=0; alpha < 2; ++alpha)
        {
           for(unsigned int i=0; i< mDim; ++i)
            {
                for(unsigned int I=0; I <mNumberOfNodes; ++I)
                {
                    // compute a1, a2
                    a[alpha](i) += DN_De(I,alpha)*( X[I](i) + u(I,i) );
                }
            }
        }
    }
    

    void KirchhoffLoveLargeDeformationShell::CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
                                                                            , const std::vector<Vector>& X)
    {
        if (A.size() != 2)
           A.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if (A[i].size() != mDim)
            {
                A[i].resize(mDim);
                A[i] = ZeroVector(mDim);
            }
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



    void KirchhoffLoveLargeDeformationShell::FirstDerivativeLocalMembraneStrain_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& eLocalVector_r
        , const Matrix& DN_De,  std::vector<Vector>& a, std::vector<Vector>& EE, std::vector<Vector>& AA)
    {
        if (eLocalVector_r.size() != 2)
            eLocalVector_r.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(eLocalVector_r[i].size() !=2)
                eLocalVector_r[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (eLocalVector_r[i][j].size() != mNumberOfDof)
                {
                    eLocalVector_r[i][j].resize(mNumberOfDof);
                    eLocalVector_r[i][j] = ZeroVector(mNumberOfDof);
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
                                // compute a_11
                                temp(I*mDim + i) += 0.5*( DN_De(I,alpha)*a[beta](i) + DN_De(I,beta)*a[alpha](i) )*
                                                MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                            }
                        }
                    }
                }

                noalias(eLocalVector_r[gamma][delta])= temp;
            }
        }

        eLocalVector_r[0][1] *= 2;
        eLocalVector_r[1][0] *= 2;
    }
        
    ////////////////// second derivatives
    void KirchhoffLoveLargeDeformationShell::DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab, 
                                                 const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector>& X, const Matrix& u)
    {
        if (a_ab.size() != 2)
            a_ab.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(a_ab[i].size() !=2)
                a_ab[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (a_ab[i][j].size() != mDim)
                {
                    a_ab[i][j].resize(mDim);
                    a_ab[i][j] = ZeroVector(mDim);
                }
            }
        }

        for(unsigned alpha=0; alpha<2; ++alpha)
        {
            for(unsigned beta=0; beta<2; ++beta)
            {
                for(unsigned int I=0; I< mNumberOfNodes; ++I)
                {
                    for(unsigned int i=0; i< mDim; ++i)
                    {
                        // compute a_11
                        a_ab[alpha][beta](i) += D2N_De2[I](alpha,beta)*( X[I](i) + u(I,i) );
                    }
                }
            }
        }

    }

    void KirchhoffLoveLargeDeformationShell::DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab, 
        const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector>& X)
    {
        if (A_ab.size() != 2)
            A_ab.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(A_ab[i].size() !=2)
                A_ab[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (A_ab[i][j].size() != mDim)
                {
                    A_ab[i][j].resize(mDim);
                    A_ab[i][j] = ZeroVector(mDim);
                }
            }
        }

        for(unsigned alpha=0; alpha<2; ++alpha)
        {
            for(unsigned beta=0; beta<2; ++beta)
            {
                for(unsigned int I=0; I< mNumberOfNodes; ++I)
                {
                    for(unsigned int i=0; i< mDim; ++i)
                    {
                        // compute a_11
                        A_ab[alpha][beta](i) += D2N_De2[I](alpha,beta)*X[I](i) ;
                    }
                }
            }
        }

    }

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeLocalMembraneStrain_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eLocalVector_rs,
        const Matrix& DN_De, std::vector<Vector>& EE, std::vector<Vector>& AA)
    {
        if (eLocalVector_rs.size() != 2)
            eLocalVector_rs.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(eLocalVector_rs[i].size() !=2)
                eLocalVector_rs[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (eLocalVector_rs[i][j].size1() != mNumberOfDof)
                {
                    eLocalVector_rs[i][j].resize(mNumberOfDof, mNumberOfDof);
                    eLocalVector_rs[i][j] = ZeroMatrix(mNumberOfDof, mNumberOfDof);
                }
            }
        }

        std::vector<Vector> e;
        this->UnitBaseVectors(e);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                Matrix temp = ZeroMatrix(mNumberOfDof, mNumberOfDof);

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
                                        // compute a_11
                                        temp(I*mDim+i, J*mDim+j) +=  ( 0.5*( DN_De(I,alpha)*DN_De(J,beta)
                                                                    + DN_De(I,beta)*DN_De(J,alpha)  ) *this->KroneckerDelta(i,j) )
                                                            *MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                                    }
                                }
                            }
                        }
                    }
                }

                noalias(eLocalVector_rs[gamma][delta])= temp;
            }
        }

        eLocalVector_rs[0][1] *= 2;
        eLocalVector_rs[1][0] *= 2;
    }

    ////////////////////////////////////////////////////////////////////
    //////////////////////////// bending part /////////////////////////
    void KirchhoffLoveLargeDeformationShell::DeformedNormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3,  std::vector<Vector>& a)
    {
        aa3Vector = MathUtils<double>::CrossProduct(a[0], a[1]);
        a3 = MathUtils<double>::Norm3(aa3Vector);
        a3Vector = aa3Vector/a3;
    }

    void KirchhoffLoveLargeDeformationShell::ReferencedNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A)
    {
        AA3Vector = MathUtils<double>::CrossProduct(A[0], A[1]);
        A3 = MathUtils<double>::Norm3(AA3Vector);
        A3Vector = AA3Vector/A3;
    }

    void KirchhoffLoveLargeDeformationShell::FirstDerivativeLocalCurvatureChange_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& kLocalVector_r
        , std::vector<Vector>& a
        , Vector& a3Vector, Vector& aa3Vector, double& a3 ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
        , const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2, std::vector<Vector>& EE, std::vector<Vector>& AA)
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

        std::vector<Vector> e;
        this->UnitBaseVectors(e);

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
                                Vector a3_rVector(mDim), aa3_rVector(mDim) ; // need to declare the size of vector
                                double a3_r;
                
                                noalias(aa3_rVector) =  DN_De(I,0)*MathUtils<double>::CrossProduct(e[i],a[1]) + DN_De(I,1)*MathUtils<double>::CrossProduct(a[0],e[i]);
                                a3_r = MathUtils<double>::Dot3(aa3_rVector, a3Vector);   // not use noalias here
                                a3_rVector  = (aa3_rVector*a3 - aa3Vector*a3_r)*pow(a3,-2);
           
                                temp(I*mDim + i) += (- D2N_De2[I](alpha,beta)*a3Vector(i) - MathUtils<double>::Dot3(a_ab[alpha][beta], a3_rVector) )*
                                                            MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]) ;

                            }
                        }
                    }
                }

                noalias(kLocalVector_r[gamma][delta])= temp;
            }
        }

        kLocalVector_r[0][1] *= 2;
        kLocalVector_r[1][0] *= 2;
    }

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeLocalCurvatureChange_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kLocalVector_rs
        ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2
        ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
        ,std::vector<Vector>& a,Vector& a3Vector, Vector& aa3Vector, const double& a3
        , std::vector<Vector>& EE, std::vector<Vector>& AA)
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
                    kLocalVector_rs[i][j] = ZeroMatrix(mNumberOfDof, mNumberOfDof);
                }
            }
        }

        std::vector<Vector> e;
        this->UnitBaseVectors(e);

        for(unsigned int gamma=0; gamma<2;++gamma)
        {
            for(unsigned int delta=0; delta<2; ++delta)
            {
                Matrix temp = ZeroMatrix(mNumberOfDof, mNumberOfDof);

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
                                        /////////////////////// a_abrVector dot a3_sVector
                                        //////////// a_abrVector
                                        Vector a_abrVector= ZeroVector(mDim);
                                        a_abrVector = D2N_De2[I](alpha,beta)*e[i];

                                        //////////// a3_sVector
                                        Vector a3_sVector =ZeroVector(mDim);

                                        // aa3_sVector
                                        Vector aa3_sVector;
                                        aa3_sVector = DN_De(J,0)*MathUtils<double>::CrossProduct(e[j],a[1]) + MathUtils<double>::CrossProduct(a[0],e[j])*DN_De(J,1);

                                        // a3_s
                                        double a3_s = MathUtils<double>::Dot3(aa3_sVector, a3Vector);

                                        // a3_sVector
                                        noalias(a3_sVector) = pow(a3,-2)*(aa3_sVector*a3 - aa3Vector*a3_s);


                                        /////////////////////// a_absVector dot a3_rVector
                                        ////////// a_absVector
                                        Vector a_absVector= ZeroVector(mDim);
                                        a_absVector = D2N_De2[J](alpha,beta)*e[j];

                                        //////////// a3_rVector
                                        Vector a3_rVector =ZeroVector(mDim);

                                        // aa3_rVector
                                        Vector aa3_rVector = ZeroVector(mDim);
                                        noalias(aa3_rVector) = DN_De(I,0)*MathUtils<double>::CrossProduct(e[i],a[1]) + MathUtils<double>::CrossProduct(a[0],e[i])*DN_De(I,1);

                                        // a3_r
                                        double a3_r = MathUtils<double>::Dot3(aa3_rVector, a3Vector);

                                        // a3_rVector
                                        noalias(a3_rVector) = pow(a3,-2)*(aa3_rVector*a3 - aa3Vector*a3_r);

                                        /////////////////////////  a_abVector dot a3_rsVector
                                        ///////// a3_rsVector
                                        // aa3_rsVector
                                        Vector aa3_rsVector = ZeroVector(mDim);
                                        noalias(aa3_rsVector) = DN_De(I,0)*MathUtils<double>::CrossProduct(e[i],e[j])*DN_De(J,1) + 
                                                             DN_De(J,0)*MathUtils<double>::CrossProduct(e[j],e[i])*DN_De(I,1);
                                        // a3_rs
                                        double a3_rs = pow(a3,-1)*( MathUtils<double>::Dot3(aa3_rsVector,aa3Vector) + MathUtils<double>::Dot3(aa3_rVector, aa3_sVector)
                                                        - MathUtils<double>::Dot3(aa3_rVector,a3Vector)*MathUtils<double>::Dot3(aa3_sVector,a3Vector) );

                                        // a3_rsVector
                                        Vector a3_rsVector = ZeroVector(mDim);
                                        noalias(a3_rsVector) = pow(a3,-1)*( aa3_rsVector-a3_rs*a3Vector ) + pow(a3,-2)*( 2*a3Vector*a3_r*a3_s - 
                                                            a3_s*aa3_rVector - a3_r*aa3_sVector );

                                        ///////////////////////////////////////////////////////////////////////////////////
                                        /////////////// second derivative of curvature changes ///////////////////////////
                                        temp(I*mDim+i, J*mDim+j) += (-( MathUtils<double>::Dot3(a_abrVector,a3_sVector)
                                                                         +MathUtils<double>::Dot3(a_absVector,a3_rVector)
                                                                         +MathUtils<double>::Dot3(a_ab[alpha][beta],a3_rsVector) ) )
                                                            *MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]) ; 
                                    }
                                }
                            }
                        }
                    }
                }

                noalias(kLocalVector_rs[gamma][delta])= temp;
            }
        }

        kLocalVector_rs[0][1] *= 2;
        kLocalVector_rs[1][0] *= 2;
    }


    void KirchhoffLoveLargeDeformationShell::CreatingBmatrix(Matrix& BMatrix, const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& LocalStrainVector_r)
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




} // Namespace Kratos
