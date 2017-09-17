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

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)
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
        mDN_DX.resize(integration_points.size());  
                 
        for (unsigned int i = 0; i < integration_points.size(); ++i)
        {
            mInvJ0[i].resize(mDim, mDim, false);
            noalias(mInvJ0[i]) = ZeroMatrix(mDim, mDim);

            mN[i].resize(mNumberOfNodes);
            noalias(mN[i]) = ZeroVector(mNumberOfNodes);

            mDN_De[i].resize(mNumberOfNodes, mDim);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes, mDim);

            mDN_DX[i].resize(mNumberOfNodes, mDim);
            noalias(mDN_DX[i]) = ZeroMatrix(mNumberOfNodes, mDim);
        }
        

        mDetJ0.resize(integration_points.size(), false);        
        // TODO remove the storage for Jacobian to save memory
        noalias(mDetJ0) = ZeroVector(integration_points.size());

        mIntegrationWeight.resize(integration_points.size());


        // calculate the Jacobian
        mJ0.resize(integration_points.size());
        mJ0 = mpIsogeometricGeometry->Jacobian0(mJ0, mThisIntegrationMethod);

     
        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++ PointNumber)
        {       
            mN[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsValues( mN[PointNumber] , integration_points[PointNumber]);
            mDN_De[PointNumber] = mpIsogeometricGeometry->ShapeFunctionsLocalGradients( mDN_De[PointNumber], integration_points[PointNumber]);

            MathUtils<double>::InvertMatrix(mJ0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);

            // compute the gradient w.r.t global coordinates
            noalias(mDN_DX[PointNumber]) = prod(mDN_De[PointNumber], mInvJ0[PointNumber]);

            //getting informations for integration
            mIntegrationWeight[PointNumber] = integration_points[PointNumber].Weight();
        }
   


        mIsInitialized = true;

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


            //////////////////////////////////////////////////////////////
            //1. compute normal forces
            Vector nVector;  // normal forces
            ///// 1.1 membrane strain
            Vector eVector;
            //Vector eVector1;
            /////// i. covariant base vectors
            std::vector<Vector> a;
            std::vector<Vector> A;  // covarianti base vector of undeformed configuration
            CovariantBaseVector( A, mDN_De[PointNumber], mNodalCoordinates);
            CovariantBaseVector(a,mDN_De[PointNumber], mNodalCoordinates, CurrentDisplacement);
            /////// ii. membrane strain
            //computeMembraneStrain( eVector, a, A);
            ///// 1.2 Tangantial material stiffness
            Matrix TanC;
            computeTangentMaterialStiffness(TanC, A);
            ///// 1.3 normal forces
            //computeNormalForces(nVector, TanC,  eVector);

            //2. compute membrane B operator
            //CovariantBaseVector(a, DN_De, mNodalCoordinates, CurrentDisplacement);
            Matrix Bm;
            computeMembraneBMatrix(Bm, mDN_De[PointNumber], a);

            computeStrain(eVector,  Bm,  CurrentDisplacement);
            computeNormalForces(nVector, TanC,  eVector);


            //3. compute bending moment
            Vector moVector; //  bending moment
            ///3.1 compute curvature changes
            Vector kVector;
            //Vector kVector1;
            ////i. normal vector
            Vector a3Vector, A3Vector, aa3Vector, AA3Vector; double a3, A3;
            ReferencedNormalDirector( A3Vector,  AA3Vector, A3, A);
            DeformedNormalDirector(a3Vector, aa3Vector, a3, a);
            ////ii. derivative of covariant base vectors
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > A_ab;
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > a_ab;
            DerivativeCovariantBaseVector( A_ab, D2N_De2, mNodalCoordinates);
            DerivativeCovariantBaseVector( a_ab, D2N_De2, mNodalCoordinates, CurrentDisplacement);     

            //4. compute bending B operator
            Matrix Bb;
            computeBendingBMatrix(Bb, a, a3Vector, aa3Vector, a3, a_ab, mDN_De[PointNumber], D2N_De2);

            ////iii. curvature changes 
            computeCurvatureChange(kVector,A_ab,a_ab,a3Vector,A3Vector );

            //computeStrain(kVector,  Bb,  CurrentDisplacement);
            computeBendingMoments(moVector, TanC, kVector);

            //////////////////////////////////////////////////////////
            ///////// this section for nonlinear computation /////////
            // compute e_rs
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>> eVector_rs;
            SecondDerivativeMembraneStrain_rs( eVector_rs, mDN_De[PointNumber]);

            // compute k_rs
            std::vector<Vector> e; //basis vector
            UnitBaseVectors(e);

            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>> a_abr; // a_{alpha beta,r}
            SecondDerivativeCovariantBaseVector_ur( a_abr, D2N_De2, e);

            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>> kVector_rs;
            SecondDerivativeCurvatureChange_rs( kVector_rs, mDN_De[PointNumber], D2N_De2, a_ab, a, a3Vector, aa3Vector, a3, e);
            /////////////////////////////////////////////////////////
    
            ////////////////////////// compute LHS
            if(CalculateStiffnessMatrixFlag==true)
            {
                AddLinearMembraneStiffness(rLeftHandSideMatrix, TanC, Bm, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                AddLinearBendingStiffness(rLeftHandSideMatrix, TanC, Bb, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);

                AddNonlinearMembraneStiffness(rLeftHandSideMatrix, nVector, eVector_rs,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                AddNonlinearBendingStiffness(rLeftHandSideMatrix, moVector, kVector_rs,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
            }

            ////////////////////////// compute RHS   
            if(CalculateResidualVectorFlag==true)
            {
                /////////// compute f_ext
                AddExternalForces(rRightHandSideVector, mN[PointNumber], mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                /////////////////////////

                // add f_int
                AddMembraneInternalForces(rRightHandSideVector, nVector, Bm, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);  
                AddBendingInternalForces(rRightHandSideVector, moVector, Bb, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                
            }

            
        }
        ///////////////////////////////////////////////////end loop over integration points

        if(CalculateResidualVectorFlag == true)
        {
            for ( unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node )
            {
                if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_X))
                {
                    double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                    for( unsigned int i = 0; i < mNumberOfDof; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim) * temp;
                }
                if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Y))
                {
                    double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                    for( unsigned int i = 0; i < mNumberOfDof; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 1) * temp;
                }
                if((*mpIsogeometricGeometry)[node].IsFixed(DISPLACEMENT_Z))
                {
                    double temp = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                    for( unsigned int i = 0; i < mNumberOfDof; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 2) * temp;
                }
            }
        }
        

        KRATOS_CATCH("")
    }

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution

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

    ///////////////////// add right hand side contribution
    void KirchhoffLoveLargeDeformationShell::AddExternalForces(VectorType& RightHandSideVector, const Vector& N, 
            double DetJ, double Weight)
    {
        const Vector BodyForce = GetProperties()[BODY_FORCE];
        for(unsigned int I=0; I < mNumberOfNodes; ++I)
            for(unsigned int i=0; i<mDim; ++i)
                RightHandSideVector(I*mDim + i) += N(I)*BodyForce(i)*DetJ*Weight;
    }

    void KirchhoffLoveLargeDeformationShell::AddMembraneInternalForces(VectorType& RightHandSideVector
                                , const Vector& nVector, const Matrix& Bm, double DetJ, double Weight)
    {
        noalias(RightHandSideVector) -= prod(trans(Bm), nVector)*DetJ*Weight;
    }

    void KirchhoffLoveLargeDeformationShell::AddBendingInternalForces(VectorType& RightHandSideVector
                                , const Vector& moVector, const Matrix& Bb, double DetJ, double Weight)
    {
        noalias(RightHandSideVector) -= prod(trans(Bb), moVector)*DetJ*Weight;
    }
    

    ////////////////////////////////////////////////////////////////////////////////////////// end 


    ///////////////////// tangent stiffness
    void KirchhoffLoveLargeDeformationShell::computeTangentMaterialStiffness(Matrix& TanC, std::vector<Vector>& A)
    {
        TanC = ZeroMatrix(mStrainSize, mStrainSize);

        Matrix Aab(2,2);
        Aab(0, 0) = inner_prod(A[0], A[0]);
        Aab(0, 1) = inner_prod(A[0], A[1]);
        Aab(1, 0) = inner_prod(A[1], A[0]);
        Aab(1, 1) = inner_prod(A[1], A[1]);
    
        double temp_ab;
        Matrix InvAab(2,2);
        MathUtils<double>::InvertMatrix(Aab, InvAab, temp_ab);
    
        // material parameters
        double E = GetProperties()[YOUNG_MODULUS];
        double NU = GetProperties()[POISSON_RATIO];

    
        double aux = E/(1-pow(NU,2) );
    
        if (TanC.size1() != 3 || TanC.size2() != 3)
            TanC.resize(3, 3, false);
    
            TanC(0,0) = aux*pow(InvAab(0,0), 2);
            TanC(0,1) = aux*( NU*InvAab(0,0)*InvAab(1,1) + (1-NU)*pow(InvAab(0,1), 2) ); 
            TanC(0,2) = aux*InvAab(0,0)*InvAab(0,1);
            TanC(1,0) = TanC(0,1);
            TanC(1,1) = aux*pow(InvAab(1,1), 2); 
            TanC(1,2) = aux*InvAab(1,1)*InvAab(0,1);
            TanC(2,0) = TanC(0,2);
            TanC(2,1) = TanC(1,2);
            TanC(2,2) = aux*0.5*( (1-NU)*InvAab(0,0)*InvAab(1,1) + (1+NU)*pow(InvAab(0,1), 2) ) ;
    }
    

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
    
    void KirchhoffLoveLargeDeformationShell::computeMembraneStrain(Vector& eVector,  std::vector<Vector>& a,  std::vector<Vector>& A)
    {
        eVector.resize(mDim);
        eVector(0) = 0.5*( MathUtils<double>::Dot3(a[0],a[0]) - MathUtils<double>::Dot3(A[0], A[0]) );
        eVector(1) = 0.5*( MathUtils<double>::Dot3(a[1],a[1]) - MathUtils<double>::Dot3(A[1], A[1]) );
        eVector(2) = ( MathUtils<double>::Dot3(a[0],a[1]) - MathUtils<double>::Dot3(A[0], A[1]) );
    }

    void KirchhoffLoveLargeDeformationShell::computeStrain(Vector& StrainVector,  const Matrix& B,  const Matrix& Displacements)
    {
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = 3;
   
        StrainVector = ZeroVector(strain_size);

        for(unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
            for(unsigned int j = 0; j < strain_size; ++j)
                for(unsigned int k = 0; k < dim; ++k)
                    StrainVector[j] += B(j, i*dim + k) * Displacements(i, k);
    }
    
    
    

    void KirchhoffLoveLargeDeformationShell::computeCurvatureChange(Vector& kVector,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
        ,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& a_ab,  Vector a3Vector,  Vector A3Vector )
    {
        kVector.resize(mDim);
        kVector(0) = MathUtils<double>::Dot3(A_ab[0][0], A3Vector) - MathUtils<double>::Dot3(a_ab[0][0], a3Vector) ;
        kVector(1) = MathUtils<double>::Dot3(A_ab[1][1], A3Vector) - MathUtils<double>::Dot3(a_ab[1][1], a3Vector) ;
        kVector(2) = 2*( MathUtils<double>::Dot3(A_ab[0][1], A3Vector) - MathUtils<double>::Dot3(a_ab[0][1], a3Vector) ) ;
    }
    ///////////////////// this section contains all B matrix
    void KirchhoffLoveLargeDeformationShell::computeMembraneBMatrix(Matrix& Bm, const Matrix& DN_De, const std::vector<Vector>& a)
    {
        if(Bm.size1() != mStrainSize)
            Bm.resize(mStrainSize,mNumberOfDof);
        noalias(Bm) = ZeroMatrix(mStrainSize, mNumberOfDof);

        for(unsigned int I=0; I< mNumberOfNodes; ++I)
        {
            for(unsigned int i=0; i< mDim; ++i)
            {
                Bm(0, I*mDim + i) = DN_De(I, 0)*a[0](i);
                Bm(1, I*mDim + i) = DN_De(I, 1)*a[1](i);
                Bm(2, I*mDim + i) =  DN_De(I,0)*a[1](i) + DN_De(I,1)*a[0](i);
            }
        }
    }

    void KirchhoffLoveLargeDeformationShell::computeBendingBMatrix(Matrix& Bb,  std::vector<Vector>& a
        , Vector& a3Vector, Vector& aa3Vector, double a3 ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                 ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2)
    {
        if(Bb.size1() != mStrainSize)
            Bb.resize(mStrainSize, mNumberOfDof);
        noalias(Bb) = ZeroMatrix(mStrainSize, mNumberOfDof);

        std::vector<Vector> ei;
        this->UnitBaseVectors(ei);

        std::vector<Vector> eixa2(3);
        eixa2[0] = MathUtils<double>::CrossProduct(ei[0],a[1]);
        eixa2[1] = MathUtils<double>::CrossProduct(ei[1],a[1]);
        eixa2[2] = MathUtils<double>::CrossProduct(ei[2],a[1]);

        std::vector<Vector> a1xei(3);
        a1xei[0] = MathUtils<double>::CrossProduct(a[0],ei[0]);
        a1xei[1] = MathUtils<double>::CrossProduct(a[0],ei[1]);
        a1xei[2] = MathUtils<double>::CrossProduct(a[0],ei[2]);

        for(unsigned int i=0; i< mNumberOfNodes; ++i)
        {
            for(unsigned int j=0; j< mDim; ++j)
            {


                Vector a3_rVector(mDim), aa3_rVector(mDim) ; // need to declare the size of vector
                double a3_r;

                noalias(aa3_rVector) = ( DN_De(i,0)*eixa2[j] + DN_De(i,1)*a1xei[j] );
                a3_r = MathUtils<double>::Dot3(aa3_rVector, a3Vector);   // not use noalias here
                a3_rVector  = (aa3_rVector*a3 - aa3Vector*a3_r)*pow(a3,-2);
       
                

                Bb(0, i*mDim + j) = - D2N_De2[i](0,0)*a3Vector(j) - MathUtils<double>::Dot3(a_ab[0][0], a3_rVector);
                Bb(1, i*mDim + j) = - D2N_De2[i](1,1)*a3Vector(j) - MathUtils<double>::Dot3(a_ab[1][1], a3_rVector);
                Bb(2, i*mDim + j) = 2*(- D2N_De2[i](0,1)*a3Vector(j) - MathUtils<double>::Dot3(a_ab[0][1], a3_rVector) );


            }
        }
    }    
    //////////////////// second derivatives of membrane strains and curvature changes
    void KirchhoffLoveLargeDeformationShell::SecondDerivativeMembraneStrain_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eVector_rs,
        const Matrix& DN_De)
    {
        if (eVector_rs.size() != 2)
            eVector_rs.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(eVector_rs[i].size() !=2)
                eVector_rs[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (eVector_rs[i][j].size1() != mNumberOfDof)
                {
                    eVector_rs[i][j].resize(mNumberOfDof, mNumberOfDof);
                    eVector_rs[i][j] = ZeroMatrix(mNumberOfDof, mNumberOfDof);
                }
            }
        }

        for(unsigned int a=0; a< 2; ++a)
        {
            for(unsigned int b=0; b<2; ++b)
            {
                for(unsigned int I=0; I<mNumberOfNodes; ++I)
                {
                    for(unsigned int J=0; J<mNumberOfNodes; ++J)
                    {
                        for(unsigned int i=0;i<mDim; ++i)
                        {
                            for(unsigned int j=0; j<mDim; ++j)
                            {
                                eVector_rs[a][b](I*mDim+i, J*mDim+j)= 0.5*( DN_De(I,a)*DN_De(J,b) + DN_De(I,b)*DN_De(J,a) )*this->KroneckerDelta(i,j); 
                            }
                        }
                    }
                }
            }
        }

        eVector_rs[0][1] *= 2;
        eVector_rs[1][0] *= 2;
    }

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeCurvatureChange_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kVector_rs
        ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2
        ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > a_ab
        ,std::vector<Vector>& a,Vector& a3Vector, Vector& aa3Vector, const double a3,  std::vector<Vector>& e)
    {
        if (kVector_rs.size() != 2)
            kVector_rs.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(kVector_rs[i].size() !=2)
                kVector_rs[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (kVector_rs[i][j].size1() != mNumberOfDof)
                {
                    kVector_rs[i][j].resize(mNumberOfDof, mNumberOfDof);
                    kVector_rs[i][j] = ZeroMatrix(mNumberOfDof, mNumberOfDof);
                }
            }
        }



        for(unsigned int alpha=0; alpha< 2; ++alpha)
        {
            for(unsigned int beta=0; beta<2; ++beta)
            {
                for(unsigned int I=0; I<mNumberOfNodes; ++I)
                {
                    for(unsigned int J=0; J<mNumberOfNodes; ++J)
                    {
                        for(unsigned int i=0;i<mDim; ++i)
                        {
                            for(unsigned int j=0; j<mDim; ++j)
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
                                kVector_rs[alpha][beta](I*mDim+i, J*mDim+j)= -( MathUtils<double>::Dot3(a_abrVector,a3_sVector)
                                                                         +MathUtils<double>::Dot3(a_absVector,a3_rVector)
                                                                         +MathUtils<double>::Dot3(a_ab[alpha][beta],a3_rsVector) ) ; 
                            }
                        }
                    }
                }
            }
        }

        kVector_rs[0][1] *= 2;
        kVector_rs[1][0] *= 2;
    }


    ///////////////////// covariant base vectors
    void KirchhoffLoveLargeDeformationShell::CovariantBaseVector(std::vector<Vector>& a, const Matrix& DN_De
                    , const std::vector<Vector> X, const Matrix& u )
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
    //////////////////////////////////////////////////////////////////////// end

    ///////////////////// derivative of covariant base vectors w.r.t curvilinear coordinates
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
    ////////////////////////////////////////////////////////////// end

    ///////////////////// derivative of covariant base vectors w.r.t nodal displacements
    void KirchhoffLoveLargeDeformationShell::DerivativeCovariantBaseVector_u(std::vector<Matrix>& a_ar, const Matrix& DN_De,std::vector<Vector>& e)
    {
        if (a_ar.size() != 2)
            a_ar.resize(2);

        for(unsigned int alpha=0; alpha< 2; ++alpha)
        {
            if (a_ar[alpha].size1() != mDim)
                a_ar[alpha].resize(mDim, mNumberOfDof);
            a_ar[alpha] = ZeroMatrix(mDim, mNumberOfDof);
        }

        // note that a_ar means the derivative of covariant base vectors a_alpha w.r.t nodal displacement u_r
        for(unsigned int alpha = 0; alpha < 2; ++alpha )
        {
            for(unsigned int a=0; a < mNumberOfNodes; ++a)
            {
                for(unsigned int i=0; i < mDim; ++i)
                {
                    // a_ar : a_alpha,r
                    a_ar[alpha](0, a*mDim +i) = DN_De(a,alpha)*e[i](0);
                    a_ar[alpha](1, a*mDim +i) = DN_De(a,alpha)*e[i](1);
                    a_ar[alpha](2, a*mDim +i) = DN_De(a,alpha)*e[i](2);
                }
            }
        }
    }

    void KirchhoffLoveLargeDeformationShell::DerivativeNumeratorNomalDirector(Matrix& aa3_r,  std::vector<Vector>& a, 
                                                                            std::vector<Vector>& e, const Matrix DN_De)
    {
        if(aa3_r.size1() != mStrainSize)
            aa3_r.resize(mStrainSize, mNumberOfDof);

        std::vector<Vector> eixa2 , eixa1;
        eixa1.resize(3);
        eixa2.resize(3);

        for(unsigned int i=0; i<3; i++)
        {
            eixa1[i].resize(mDim);
            eixa2[i].resize(mDim);
        }


        for(unsigned int i=0; i<3; i++)
        {
            eixa1[i] = MathUtils<double>::CrossProduct(e[i], a[0]);
            eixa2[i] = MathUtils<double>::CrossProduct(e[i], a[1]);
        }


        for(unsigned int j=0; j< mStrainSize; j++)
        {
            for(unsigned int I=0; I< mNumberOfNodes; I++)
            {
                for(unsigned int i=0; i< mDim; i++)
                {
                    aa3_r(j, I*mDim + i) = DN_De(I,0)*eixa2[i](j) + DN_De(I,1)*eixa1[i](j);
                }
            }
        }

    }

    void KirchhoffLoveLargeDeformationShell::DerivativeDenominatorNormalDirector(Matrix& a3_r,  Matrix& aa3_r,  Vector& a3Vector )
    {
        a3_r = ZeroMatrix(1, mNumberOfDof);

        for(unsigned int I=0; I < mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i< mDim; i++)
            {
                a3_r(0, I*mDim + i) = aa3_r(0, I*mDim +i)*a3Vector(0) + aa3_r(1, I*mDim +i)*a3Vector(1) + aa3_r(2, I*mDim +i)*a3Vector(2);
            }
        }
    }



    void KirchhoffLoveLargeDeformationShell::SecondDerivativeCovariantBaseVector_ur(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> >& a_abr
        , const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector> e)
    {
        if (a_abr.size() != 2)
            a_abr.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(a_abr[i].size() !=2)
                a_abr[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (a_abr[i][j].size1() != mStrainSize)
                {
                    a_abr[i][j].resize(mStrainSize, mNumberOfDof);
                    a_abr[i][j] = ZeroMatrix(mStrainSize, mNumberOfDof);
                }
            }
        }

        for(unsigned int alpha=0; alpha<2; ++alpha)
        {
            for(unsigned int beta=0; beta<2; ++beta)
            {
                for(unsigned int I=0; I< mNumberOfNodes; ++I)
                {
                    for(unsigned int i=0; i< mDim; ++i)
                    {
                        a_abr[alpha][beta](0, I*mDim + i) += D2N_De2[I](alpha,beta)*e[i](0);
                        a_abr[alpha][beta](1, I*mDim + i) += D2N_De2[I](alpha,beta)*e[i](1);
                        a_abr[alpha][beta](2, I*mDim + i) += D2N_De2[I](alpha,beta)*e[i](2);
                    }
                }
            }
        }
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


} // Namespace Kratos
