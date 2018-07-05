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
#include "custom_elements/K_L_linear_shell.h"
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


namespace Kratos
{
    extern Variable<Vector> STRESSES;
    extern Variable<array_1d<double, 3> > PRESCRIBED_DELTA_DISPLACEMENT;

    //************************************************************************************
    //***** Constructor and Destructor ***************************************************
    //************************************************************************************
    KirchhoffLoveLinearShell::KirchhoffLoveLinearShell
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

    KirchhoffLoveLinearShell::KirchhoffLoveLinearShell(IndexType NewId,
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
    KirchhoffLoveLinearShell::~KirchhoffLoveLinearShell()
    {
    }


    //********************************************************
    //**** Operations ****************************************
    //********************************************************

    Element::Pointer KirchhoffLoveLinearShell::Create(IndexType NewId,
        NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KirchhoffLoveLinearShell(NewId,
                                GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer KirchhoffLoveLinearShell::Create(IndexType NewId,
        GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KirchhoffLoveLinearShell(NewId,
                                                            pGeom, pProperties));
    }

    //************************************************************************************
    //************************************************************************************
    /**
    * returns the used integration method
    */
    KirchhoffLoveLinearShell::IntegrationMethod KirchhoffLoveLinearShell::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    //************************************************************************************
    //************************************************************************************
    /**
    * Setting up the EquationIdVector
    */
    void KirchhoffLoveLinearShell::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        DofsVectorType ElementalDofList;
        GetDofList(ElementalDofList, rCurrentProcessInfo);

        if (rResult.size() != ElementalDofList.size())
            rResult.resize(ElementalDofList.size(), false);

        for(unsigned int i = 0; i < ElementalDofList.size(); ++i)
            rResult[i] = ElementalDofList[i]->EquationId();
    }

//************************************************************************************
//************************************************************************************
/**
* Setting up the DOF list
*/
    void KirchhoffLoveLinearShell::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    void KirchhoffLoveLinearShell::Initialize()
    {
        KRATOS_TRY //EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

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

            mDN_De[i].resize(mNumberOfNodes, 2);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes, 2);

        }
        

        mDetJ0.resize(integration_points.size(), false);        
        // TODO remove the storage for Jacobian to save memory
        noalias(mDetJ0) = ZeroVector(integration_points.size());

        mIntegrationWeight.resize(integration_points.size());

        if (mConstitutiveLawVector.size() != integration_points.size())
            mConstitutiveLawVector.resize(integration_points.size());
        
        for(unsigned int i=0; i< integration_points.size(); ++i)
        {
            if(mConstitutiveLawVector[i].size() !=3)
                mConstitutiveLawVector[i].resize(3);
        }

        for (unsigned int i = 0; i < integration_points.size(); i++)
        {
            for(unsigned int j= 0; j< 3; j++)
            {
                mConstitutiveLawVector[i][j] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[i][j]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
     

            }   
       
        }

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
        //KRATOS_WATCH(mTotalDomainInitialSize)
        //cout<< " print mDomain Size" << endl;
            // initialize the stresses container
            mCurrentStresses.resize(integration_points.size());
            for(unsigned int i = 0; i < integration_points.size(); ++i)
                mCurrentStresses[i].resize(mStrainSize);


        // integration points and weights used for thickness calculation
        mIntegrationPoint1D.resize(3);
        mIntegrationPoint1D(0) = -std::sqrt(3.00 / 5.00) ;
        mIntegrationPoint1D(1) = 0.00;
        mIntegrationPoint1D(2) = std::sqrt(3.00 / 5.00);
               
        mWeight1D.resize(3);
        mWeight1D(0) = mWeight1D(2) = 5.00 / 9.00; mWeight1D(1) = 8.00 / 9.00;
        mDetJ1D = 0.5*mThickness;
   
        mIsPlasticExistence = false;
        mIsInitialized = true;


        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    void KirchhoffLoveLinearShell::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif
    
        const GeometryType::IntegrationPointsArrayType& integration_points =
        mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        for (unsigned int i = 0; i < integration_points.size(); i++)
        {
                for(unsigned int j= 0; j< 3; j++)
                {
                    mConstitutiveLawVector[i][j]->FinalizeSolutionStep( GetProperties(),GetGeometry(),row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),CurrentProcessInfo );
         
                }   
           
        }  
    }

    
    void KirchhoffLoveLinearShell::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
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
    void KirchhoffLoveLinearShell::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
    void KirchhoffLoveLinearShell::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
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
    void KirchhoffLoveLinearShell::CalculateAll( MatrixType& rLeftHandSideMatrix, 
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


        Vector NodalDisplacement = ZeroVector(mNumberOfDof);
        for(unsigned int I=0; I<mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                NodalDisplacement(I*mDim +i) = CurrentDisplacement(I,i);
            }
        }

        if(mIsPlasticExistence == false)
        {
        ///////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////// loop over integration points
        for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
        {

            //////////// get shape function values and their derivatives
            ShapeFunctionsSecondDerivativesType D2N_De2;
            D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    

            // first derivative of displacement
            std::vector<Vector> u_a;
            FirstDerivativeDisplacement_a( u_a,  mDN_De[PointNumber] , CurrentDisplacement);
            
            // second derivative of displacement
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>> u_ab;
            SecondDerivativeDisplacement_ab(u_ab, D2N_De2, CurrentDisplacement);

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

            Matrix Bab(2,2);
            CovariantCurvatureCoefficient(Bab , A_ab, A3Vector);


            //KRATOS_WATCH(Aab)
      
            //KRATOS_WATCH(Bab)
          

            // contravariant base vectors
            std::vector<Vector> AA;
            ContravariantBaseVector( AA, A, Aab);

            /////////////////////////////////////////////
            ///// local Cartesian basis /////////////////
            std::vector<Vector> EE;
            LocalCartesianBasisVector(EE,  A);


            ///// 1.1 membrane strain
            Matrix eTensor(2,2);
            computeMembraneStrain(eTensor, A,  u_a );
            // transform tensor to vector (not necessarily)
            Matrix eeTensor(2,2);
            LocalTransformationOfTensor(eeTensor , eTensor ,  EE, AA);
            Vector eeVector = ZeroVector(3);
            eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
               
            // curvature change
            Matrix kTensor(2,2);
            computeCurvatureChange(kTensor, A, A3, A3Vector, u_a, u_ab, A_ab);

            Matrix kkTensor(2,2);
            LocalTransformationOfTensor(kkTensor , kTensor ,  EE, AA);
            Vector kkVector = ZeroVector(3);
            kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);
            

                
            // first derivative of membrane strain
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > eLocalVector_r;
            FirstDerivativeLocalMembraneStrain_r(eLocalVector_r,  mDN_De[PointNumber], A, EE, AA) ;

            // membrane B matrix
            Matrix BBm(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBm, eLocalVector_r);
            
            // first derivative of curvature change
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > kLocalVector_r;
            FirstDerivativeLocalCurvatureChange_r(kLocalVector_r, A, A3, A3Vector, mDN_De[PointNumber], D2N_De2,A_ab, EE, AA);

            // bending B matrix
            Matrix BBb(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBb, kLocalVector_r);


            // material matrix
            Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
            CalculateElasticMatrix(D0, mE, mNU);


            // stress resultants
            Vector nVector= ZeroVector(mStrainSize);
            Vector moVector= ZeroVector(mStrainSize);
            noalias(nVector) = mThickness*prod(D0, eeVector);
            noalias(moVector) = pow(mThickness,3)/12.0*prod(D0, kkVector);

            Matrix nTensor(2,2) ;
            Matrix moTensor(2,2);
            SD_MathUtils<int>::VectorToTensor(nVector,nTensor);
            SD_MathUtils<int>::VectorToTensor(moVector,moTensor);

            //////////////////////////////////
            // compute the contribution to LHS
            AddStiffnessMatrixComponents(rLeftHandSideMatrix, D0, BBm, BBm, mThickness*detJA, mIntegrationWeight[PointNumber]);
            AddStiffnessMatrixComponents(rLeftHandSideMatrix, D0, BBb, BBb,  pow(mThickness,3.0)/12.0*detJA, mIntegrationWeight[PointNumber]);


            //////////////////////////////////
            // compute the contribution to RHS
            AddInternalForces( rRightHandSideVector, nTensor, eLocalVector_r, detJA, mIntegrationWeight[PointNumber]);
            AddInternalForces( rRightHandSideVector, moTensor, kLocalVector_r,detJA, mIntegrationWeight[PointNumber]);
            AddExternalForces(rRightHandSideVector,  mN[PointNumber],  detJA, mIntegrationWeight[PointNumber] );

                

        }// loop over integration points

        } // end elastic computation
        else // plasticity exits
        {
            for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
            {
                        //////////// get shape function values and their derivatives
                        ShapeFunctionsSecondDerivativesType D2N_De2;
                        D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    
            
                        // first derivative of displacement
                        std::vector<Vector> u_a;
                        FirstDerivativeDisplacement_a( u_a,  mDN_De[PointNumber] , CurrentDisplacement);
                        
                        // second derivative of displacement
                        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>> u_ab;
                        SecondDerivativeDisplacement_ab(u_ab, D2N_De2, CurrentDisplacement);
            
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
            
                        Matrix Bab(2,2);
                        CovariantCurvatureCoefficient(Bab , A_ab, A3Vector);
            
            
                        //KRATOS_WATCH(Aab)
                  
                        //KRATOS_WATCH(Bab)
                      
            
                        // contravariant base vectors
                        std::vector<Vector> AA;
                        ContravariantBaseVector( AA, A, Aab);
            
                        /////////////////////////////////////////////
                        ///// local Cartesian basis /////////////////
                        std::vector<Vector> EE;
                        LocalCartesianBasisVector(EE,  A);
            
            
                        ///// 1.1 membrane strain
                        Matrix eTensor(2,2);
                        computeMembraneStrain(eTensor, A,  u_a );
                        // transform tensor to vector (not necessarily)
                        Matrix eeTensor(2,2);
                        LocalTransformationOfTensor(eeTensor , eTensor ,  EE, AA);
                        Vector eeVector = ZeroVector(3);
                        eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
                           
                        // curvature change
                        Matrix kTensor(2,2);
                        computeCurvatureChange(kTensor, A, A3, A3Vector, u_a, u_ab, A_ab);
            
                        Matrix kkTensor(2,2);
                        LocalTransformationOfTensor(kkTensor , kTensor ,  EE, AA);
                        Vector kkVector = ZeroVector(3);
                        kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);

                                    // first derivative of membrane strain
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > eLocalVector_r;
                FirstDerivativeLocalMembraneStrain_r(eLocalVector_r,  mDN_De[PointNumber], A, EE, AA) ;

                // membrane B matrix
                Matrix BBm(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBm, eLocalVector_r);
            
                // first derivative of curvature change
                boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > kLocalVector_r;
                FirstDerivativeLocalCurvatureChange_r(kLocalVector_r, A, A3, A3Vector, mDN_De[PointNumber], D2N_De2,A_ab, EE, AA);

                // bending B matrix
                Matrix BBb(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBb, kLocalVector_r);


                Vector nVector= ZeroVector(mStrainSize);
                Vector moVector= ZeroVector(mStrainSize);
                Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
                Matrix D1 = ZeroMatrix(mStrainSize,mStrainSize);
                Matrix D2 = ZeroMatrix(mStrainSize,mStrainSize);

                for(unsigned int PointNumber1D = 0; PointNumber1D < mIntegrationPoint1D.size(); PointNumber1D++)
                {

                    double theta3i = 0.5*mThickness*mIntegrationPoint1D(PointNumber1D);

	                //////////////////////////////////////////////////////////////
	                //////// membrane and bending stress resultants //////////////
                    // 1. Total Strain tensor
	                Vector GLStrainVector = ZeroVector(mStrainSize);
		            noalias(GLStrainVector) +=eeVector;
                    noalias(GLStrainVector) += theta3i*kkVector;
	    

                    //////////////////////////////////////////////////////////////
                    ///////////// tangential material tensor TanC ////////////////
                    Matrix TanC = ZeroMatrix(mStrainSize,mStrainSize);
                    Vector PK2StressVector(mStrainSize);

                    mConstitutiveLawVector[PointNumber][PointNumber1D]->CalculateMaterialResponse(
                        GLStrainVector,
                        ZeroMatrix( 1 ),
                        PK2StressVector,
                        TanC,
                        rCurrentProcessInfo,
                        GetProperties(),
                        GetGeometry(),
                        mN[PointNumber],
                        true,
                        true,
                        true
                    );

           



                   
                    ///////////////////////////////////////////////////////////////////////
                    ///////////////// stress resultants //////////////////////////////////
                    // 3. membrane stress vector
                    noalias(nVector) += PK2StressVector*mWeight1D(PointNumber1D)*mDetJ1D;
                
                    // 4. bending stress vector
                    noalias(moVector) += PK2StressVector*theta3i*mWeight1D(PointNumber1D)*mDetJ1D;
                    //KRATOS_WATCH(BendingStressVector)
                
                    ////////////////////////////////////////////
                    ///////// tangent stiffness matrix /////////
                    noalias(D0) += mWeight1D(PointNumber1D)*mDetJ1D*TanC;
                    noalias(D1) += mWeight1D(PointNumber1D)*mDetJ1D*theta3i*TanC;
                    noalias(D2) += mWeight1D(PointNumber1D)*mDetJ1D*pow(theta3i,2)*TanC;      



                } // 1d integration points

                Matrix nTensor(2,2) ;
                Matrix moTensor(2,2);
                SD_MathUtils<int>::VectorToTensor(nVector,nTensor);
                SD_MathUtils<int>::VectorToTensor(moVector,moTensor);

                if(CalculateStiffnessMatrixFlag == true)
                {
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D0, BBm, BBm, detJA, mIntegrationWeight[PointNumber]);
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D1, BBm, BBb, detJA, mIntegrationWeight[PointNumber]);
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D1, BBb, BBm, detJA, mIntegrationWeight[PointNumber]);
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D2, BBb, BBb, detJA, mIntegrationWeight[PointNumber]);
                   
                }
    
                if(CalculateResidualVectorFlag == true)
                {
                      

                    AddInternalForces( rRightHandSideVector, nTensor, eLocalVector_r, detJA, mIntegrationWeight[PointNumber]);
                    AddInternalForces( rRightHandSideVector, moTensor, kLocalVector_r,detJA, mIntegrationWeight[PointNumber]);
                    AddExternalForces(rRightHandSideVector,  mN[PointNumber],  detJA, mIntegrationWeight[PointNumber] );
                }


            }// loop over integration points

        }// plasticity is chosen
      
     
        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif

        KRATOS_CATCH("")
    }

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void KirchhoffLoveLinearShell::AddInternalForces(VectorType& RightHandSideVector, const Matrix& StressResultants,
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight)
    {
       noalias(RightHandSideVector) -= StressResultants(0,0)*StrainVector_r[0][0]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(0,1)*StrainVector_r[0][1]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(1,0)*StrainVector_r[1][0]*DetJ*Weight;
       noalias(RightHandSideVector) -= StressResultants(1,1)*StrainVector_r[1][1]*DetJ*Weight;
    }

    void KirchhoffLoveLinearShell::AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight )
    {
        const Vector BodyForce = GetProperties()[BODY_FORCE];

        for(unsigned int I=0; I < mNumberOfNodes; ++I)
            for(unsigned int i=0; i<mDim; ++i)
                RightHandSideVector(I*mDim + i) += N(I)*BodyForce(i)*DetJ*Weight;
    }

    void KirchhoffLoveLinearShell::AddStiffnessMatrixComponents(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += prod( trans(BlhsMatrix), Matrix(prod(Di, BrhsMatrix)) )*DetJ*Weight;
    }

    ////////////////////////////////////////////////////////////////////////
    /////////////////// membrane strain

    void KirchhoffLoveLinearShell::computeMembraneStrain(Matrix& eTensor,  std::vector<Vector>& A, std::vector<Vector>& u_a )
    {
        eTensor.resize(2,2);
        eTensor = ZeroMatrix(2,2);

        for(int alpha =0 ; alpha<2; alpha++)
        {
            for(int beta=0; beta<2; beta++)
            {
                eTensor(alpha,beta) =  0.5*( MathUtils<double>::Dot3(A[beta],u_a[alpha]) + MathUtils<double>::Dot3(A[alpha], u_a[beta]) )  ;
            }
        }


    }

    ///////////////// curvature change
    void KirchhoffLoveLinearShell::computeCurvatureChange(Matrix& kTensor, std::vector<Vector>& A, double& A3, Vector& A3Vector, std::vector<Vector>& u_a, 
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& u_ab, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab)
    {
        kTensor.resize(2,2);
        kTensor = ZeroMatrix(2,2);

        Vector A2xA3(mDim);
        A2xA3 = MathUtils<double>::CrossProduct(A[1], A3Vector);

        Vector A3xA1(mDim);
        A3xA1 = MathUtils<double>::CrossProduct(A3Vector, A[0]);

        for(int alpha =0 ; alpha<2; alpha++)
        {
            for(int beta=0; beta<2; beta++)
            {
                
                Vector A_abxA2(mDim);
                A_abxA2 = MathUtils<double>::CrossProduct(A_ab[alpha][beta], A[1]);

                Vector A1xA_ab(mDim);
                A1xA_ab = MathUtils<double>::CrossProduct(A[0], A_ab[alpha][beta]);


                kTensor(alpha,beta) =   -MathUtils<double>::Dot3(A3Vector, u_ab[alpha][beta]) + (1/A3)*( MathUtils<double>::Dot3(A_abxA2,u_a[0]) 
                   + MathUtils<double>::Dot3(A1xA_ab,u_a[1]) + MathUtils<double>::Dot3(A_ab[alpha][beta], A3Vector)*
                     ( MathUtils<double>::Dot3( A2xA3,u_a[0]) + MathUtils<double>::Dot3(A3xA1,u_a[1])) ) ;
            }
        }

     
    }


    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
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



    void KirchhoffLoveLinearShell::ReferenceNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A)
    {
        AA3Vector = MathUtils<double>::CrossProduct(A[0], A[1]);
        A3 = MathUtils<double>::Norm3(AA3Vector);
        A3Vector = AA3Vector/A3;
    }

    void KirchhoffLoveLinearShell::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
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

        noalias(AA[0]) = AAab(0,0)*A[0] + AAab(0,1)*A[1] ;

        noalias(AA[1]) = AAab(1,0)*A[0] + AAab(1,1)*A[1] ;
    }


    void KirchhoffLoveLinearShell::DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab, 
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

    void KirchhoffLoveLinearShell::FirstDerivativeDisplacement_a(std::vector<Vector>& u_a, const Matrix& DN_De , const Matrix& u)
    {
        u_a.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
                u_a[i].resize(mDim);
                u_a[i] = ZeroVector(mDim);
        }

        for(unsigned alpha=0; alpha < 2; ++alpha)
        {
            for(unsigned int i=0; i<mDim; ++i)
            {
                for(unsigned int I=0; I<mNumberOfNodes; ++I)
                {
                    // compute a1, a2
                    u_a[alpha](i) += DN_De(I,alpha)*u(I,i) ;
                }
            }
        }    
    }

    void KirchhoffLoveLinearShell::ContinuumCovariantBaseVector(std::vector<Vector>& g, std::vector<Vector>& a
        , std::vector<Vector>& a3Vector_a, Vector& a3Vector, double& theta3)
        {
            if (g.size() != 3)
                g.resize(3);
    
            for(unsigned int i=0; i< 3; ++i)
            {
                if (g[i].size() != 3)
                {
                    g[i].resize(3);
                    g[i] = ZeroVector(3);
                }
            }
    
    
    
    
            for(unsigned int alpha=0; alpha<2; alpha++)
            {
    
                noalias(g[alpha]) = a[alpha] + theta3*a3Vector_a[alpha];
            }
    
            noalias(g[2]) = a3Vector;
    
        }

    void KirchhoffLoveLinearShell::DerivativeReferenceNormalDirector_a(std::vector<Vector>& A3Vector_a, std::vector<Vector>& A
        , Vector& AA3Vector, double& A3, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab)
    {
            if (A3Vector_a.size() != 2)
                A3Vector_a.resize(2);
    
            for(unsigned int i=0; i< 2; ++i)
            {
                A3Vector_a[i].resize(mDim);
                A3Vector_a[i] = ZeroVector(mDim);
                  
            }
    
            for (unsigned int alpha=0; alpha<2; alpha++)
            {    
                // compute aa3_a
                Vector AA3Vector_a = ZeroVector(mDim);
                noalias(AA3Vector_a) = MathUtils<double>::CrossProduct(A_ab[0][alpha], A[1]);
                noalias(AA3Vector_a) += MathUtils<double>::CrossProduct(A[0], A_ab[1][alpha]);
    
                double A3_a = 0.0;
                A3_a = MathUtils<double>::Dot3(AA3Vector, AA3Vector_a)/A3;
    
                // compute derivative normal director w.r.t alpha
                A3Vector_a[alpha] = (A3*AA3Vector_a - AA3Vector*A3_a)*pow(A3,-2.0);
            }
    
    }
    /////////////////////////////////////////////////////////////////////////
    ///// first derivative of membrane strain and necessary components //////
    /////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::FirstDerivativeLocalMembraneStrain_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& eLocalVector_r
        , const Matrix& DN_De,  std::vector<Vector>& A, std::vector<Vector>& EE, std::vector<Vector>& AA)
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
                                temp(I*mDim + i) += ( 0.5*( DN_De(I,alpha)*A[beta](i) + DN_De(I,beta)*A[alpha](i) ) )*
                                                MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                            }
                        }
                    }
                }

                noalias(eLocalVector_r[gamma][delta])= temp;
            }
        }


    


        //eLocalVector_r[0][1] *= 2.0;
        //eLocalVector_r[1][0] *= 2.0;

    }

    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::FirstDerivativeLocalCurvatureChange_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& kLocalVector_r
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


        Vector A2xA3(mDim);
        A2xA3 = MathUtils<double>::CrossProduct(A[1], A3Vector);

        Vector A3xA1(mDim);
        A3xA1 = MathUtils<double>::CrossProduct(A3Vector, A[0]);

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


    

        //kLocalVector_r[0][1] *= 2.0;
        //kLocalVector_r[1][0] *= 2.0;
    }

    void KirchhoffLoveLinearShell::SecondDerivativeDisplacement_ab(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& u_ab
        , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const Matrix& u)
    {
        if (u_ab.size() != 2)
            u_ab.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            if(u_ab[i].size() !=2)
                u_ab[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                if (u_ab[i][j].size() != mDim)
                {
                    u_ab[i][j].resize(mDim);
                    u_ab[i][j] = ZeroVector(mDim);
                }
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
                        u_ab[alpha][beta](i) += D2N_De2[I](alpha,beta)*u(I,i) ;
                    }
                }
            }
        }   
    }

    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
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


    void KirchhoffLoveLinearShell::CovariantCurvatureCoefficient(Matrix& Bab
        , boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab, Vector& A3Vector)
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



    //////////////////////////////////////////////////////////
    ////////////////////// material matrix ///////////////////
    //////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::CalculateElasticMatrix(Matrix& C, const double& E, const double& NU)
    {
        unsigned int strain_size = 3;
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size);

            double c1 = E /( 1.00 - pow(NU,2) );
            double c2 = E*NU /( 1.00 - pow(NU,2) );
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




    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////
    void KirchhoffLoveLinearShell::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A)
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



    void KirchhoffLoveLinearShell::LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector<Vector>& EE, std::vector<Vector>& AA)
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
                        temp += M(alpha,beta)*MathUtils<double>::Dot3(EE[gamma],AA[alpha])*MathUtils<double>::Dot3(AA[beta],EE[delta]);
                    }
                }
                T(gamma,delta)= temp;
            }
        }
    }


    

    void KirchhoffLoveLinearShell::UnitBaseVectors(std::vector<Vector>& e)
    {
        if (e.size() != 3)
            e.resize(3);

        e[0]=e[1]=e[2]=ZeroVector(mDim);
        e[0](0)=e[1](1)=e[2](2)= 1.0;
    }

    void KirchhoffLoveLinearShell::CreatingBmatrix(Matrix& BMatrix, const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& LocalStrainVector_r)
    {
        if(BMatrix.size1() != mStrainSize)
            BMatrix.resize(mStrainSize, mNumberOfDof);
        noalias(BMatrix) = ZeroMatrix(mStrainSize, mNumberOfDof);

        for (int i=0; i< mNumberOfDof; i++)
        {
            BMatrix(0,i) = LocalStrainVector_r[0][0](i);
            BMatrix(1,i) = LocalStrainVector_r[1][1](i);
            BMatrix(2,i) = 2.0*LocalStrainVector_r[0][1](i);
        }

    }

    int KirchhoffLoveLinearShell::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if(mTotalDomainInitialSize < 0.0)
        {
            std::stringstream ss;
            ss << "error on element -> " << this->Id() << ": ";
            ss << "Domain size can not be less than 0. Please check Jacobian. mTotalDomainInitialSize = " << mTotalDomainInitialSize;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }

        return 0;
    }


    void KirchhoffLoveLinearShell::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        GetValueOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    void KirchhoffLoveLinearShell::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {

        // resize as needed
        if ( rValues.size() != mNumberOfIntegrationPoint )
            rValues.resize( mNumberOfIntegrationPoint );

        if(rVariable == STRESSES)
        {
            for ( unsigned int i = 0; i < mNumberOfIntegrationPoint; ++i )
            {
                if ( rValues[i].size() != mStrainSize )
                    rValues[i].resize( mStrainSize );
                noalias( rValues[i] ) = mCurrentStresses[i];
            }
        }
    }


    


} // Namespace Kratos
            /*for(unsigned int PointNumber1D = 0; PointNumber1D < mIntegrationPoint1D.size(); PointNumber1D++)
            {

                double theta3i = 0.5*mThickness*mIntegrationPoint1D(PointNumber1D);

	            //////////////////////////////////////////////////////////////
	            //////// membrane and bending stress resultants //////////////
                // 1. Total Strain tensor
	            Vector GLStrainVector = ZeroVector(mStrainSize);
		        noalias(GLStrainVector) +=eeVector;
                noalias(GLStrainVector) += theta3i*kkVector;
	            //Matrix GLStrainTensor = ZeroMatrix(2,2);
		        //noalias(GLStrainTensor) +=eeTensor;
                //noalias(GLStrainTensor) += theta3i*kkTensor;

                Matrix GLStrainTensor(2, 2);
                IsotropicTensorUtility<2>::StrainVectorToTensor(GLStrainVector, GLStrainTensor);
                //////////////////////////////////////////////////////////////
                ///////////// tangential material tensor TanC ////////////////
                Matrix TanC = ZeroMatrix(mStrainSize,mStrainSize);
                Vector PK2StressVector(mStrainSize);

                //KRATOS_WATCH(GLStrainVector)
                //KRATOS_WATCH(TanC)

                if(mIsPlasticExistence == false)
                {
                    //////// 2.1 Positive part of total strain tensor

                    
                    Vector e;
                    std::vector<Matrix> E;
                    SDParameters decomp_params;
                    decomp_params = IsotropicTensorUtility<2>::SpectralDecomposition(GLStrainTensor, e, E);

                    double TraceStrain = 0.0;
                    for(unsigned int i=0; i< 2; i++)
                        TraceStrain += GLStrainTensor(i,i);

                    ///////// Positive part 
                    Matrix GLStrainTensor_plus(2, 2);
                    MyIsotropicTensorFunction_plus<2> Phi_plus;
                    IsotropicTensorUtility<2>::Value(GLStrainTensor, e, E, decomp_params, Phi_plus, GLStrainTensor_plus);

                    /////// 2.3 calculate Positive part of Cauchy stress tensor
                    Matrix PK2StressTensor_plus = ZeroMatrix(2, 2);
                    if (TraceStrain >= 0)
                        noalias(PK2StressTensor_plus) += ( mLambda - pow(mLambda,2)/(mLambda+2*mMu) )*TraceStrain*IdentityMatrix(2);
                    noalias(PK2StressTensor_plus) += 2*mMu*GLStrainTensor_plus ;


                    Matrix GLStrainTensor_minus(2, 2);
                    MyIsotropicTensorFunction_minus<2> Phi_minus;
                    IsotropicTensorUtility<2>::Value(GLStrainTensor, e, E, decomp_params, Phi_minus, GLStrainTensor_minus);

                    /////// 2.3 calculate negative part of Cauchy stress tensor
                    Matrix PK2StressTensor_minus = ZeroMatrix(2, 2);
                    if (TraceStrain < 0)
                        noalias(PK2StressTensor_minus) += ( mLambda - pow(mLambda,2)/(mLambda+2*mMu) )*TraceStrain*IdentityMatrix(2);
                    noalias(PK2StressTensor_minus) += 2*mMu*GLStrainTensor_minus ;

                
                    Matrix PK2StressTensor = ZeroMatrix(2,2);
                    noalias(PK2StressTensor) += PK2StressTensor_plus;
                    noalias(PK2StressTensor) += PK2StressTensor_minus;
                  

                    // compute stress vector
                    IsotropicTensorUtility<2>::StressTensorToVector(PK2StressTensor, PK2StressVector);

                    





    

                    Matrix Ct_Matrix_plus;
                    IsotropicTensorUtility<2>::Fourth_Order_Tensor Ct_Tensor_plus;
                    IsotropicTensorUtility<2>::InitializeFourthOrderTensor(Ct_Tensor_plus);
                    IsotropicTensorUtility<2>::Derivative(GLStrainTensor, e, E, decomp_params, Phi_plus, Ct_Tensor_plus);
                    IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(Ct_Tensor_plus, Ct_Matrix_plus);

                    IsotropicTensorUtility<2>::Fourth_Order_Tensor IxITensor;
                    IsotropicTensorUtility<2>::InitializeFourthOrderTensor(IxITensor);

                    double alphaI = 1;
                    IsotropicTensorUtility<2>::OuterProductFourthOrderTensor(alphaI, IdentityMatrix(2), IdentityMatrix(2), IxITensor);

                    Matrix IxIMatrix(mStrainSize,mStrainSize);
                    IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(IxITensor, IxIMatrix);


                    if(TraceStrain >= 0)
                    {
                        noalias(TanC) += ( mLambda - pow(mLambda,2)/(mLambda + 2*mMu) )*IxIMatrix;
                    }
                    noalias(TanC) += 2*mMu*Ct_Matrix_plus;

                    Matrix Ct_Matrix_minus;
                    IsotropicTensorUtility<2>::Fourth_Order_Tensor Ct_Tensor_minus;
                    IsotropicTensorUtility<2>::InitializeFourthOrderTensor(Ct_Tensor_minus);
                    IsotropicTensorUtility<2>::Derivative(GLStrainTensor, e, E, decomp_params, Phi_minus, Ct_Tensor_minus);
                    IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(Ct_Tensor_minus, Ct_Matrix_minus);

                    if(TraceStrain < 0)
                    {
                        noalias(TanC) += ( mLambda - pow(mLambda,2)/(mLambda + 2*mMu) )*IxIMatrix;
                    }
                    noalias(TanC) += 2*mMu*Ct_Matrix_minus;
                    

                }
                else
                {
                    Vector GLStrainVector = ZeroVector(3);
                    GLStrainVector = SD_MathUtils<int>::TensorToStrainVector( GLStrainTensor );
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->CalculateMaterialResponse(
                        GLStrainVector,
                        ZeroMatrix( 1 ),
                        PK2StressVector,
                        TanC,
                        rCurrentProcessInfo,
                        GetProperties(),
                        GetGeometry(),
                        mN[PointNumber],
                        true,
                        true,
                        true
                    );

                }



                   
                ///////////////////////////////////////////////////////////////////////
                ///////////////// stress resultants //////////////////////////////////
                // 3. membrane stress vector
                noalias(nVector) += PK2StressVector*mWeight1D(PointNumber1D)*mDetJ1D;
                
                // 4. bending stress vector
                noalias(moVector) += PK2StressVector*theta3i*mWeight1D(PointNumber1D)*mDetJ1D;
                //KRATOS_WATCH(BendingStressVector)
                
                ////////////////////////////////////////////
                ///////// tangent stiffness matrix /////////
                noalias(D0) += mWeight1D(PointNumber1D)*mDetJ1D*TanC;
                noalias(D1) += mWeight1D(PointNumber1D)*mDetJ1D*theta3i*TanC;
                noalias(D2) += mWeight1D(PointNumber1D)*mDetJ1D*pow(theta3i,2)*TanC;      



            } // 1d integration points
            


            // store the stresses
            noalias(mCurrentStresses[PointNumber]) = nVector;
            */

            /*
        for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
        {

            //////////// get shape function values and their derivatives
            ShapeFunctionsSecondDerivativesType D2N_De2;
            D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    

            /////// i. covariant base vectors
            std::vector<Vector> A;  // covarianti base vector of undeformed configuration
            CovariantBaseVector( A, mDN_De[PointNumber], mNodalCoordinates);

            //3. compute bending moment
            ////i. normal vector
            Vector A3Vector, AA3Vector; double A3;
            ReferencedNormalDirector( A3Vector,  AA3Vector, A3, A);

            ////ii. derivative of covariant base vectors
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > A_ab;
            DerivativeCovariantBaseVector( A_ab, D2N_De2, mNodalCoordinates);

            Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
            CalculateElasticMatrix(D0, mE, mNU);

            std::vector<Vector> A3Vector_a;
            DerivativeReferenceNormalDirector_a(A3Vector_a, A , AA3Vector, A3,  A_ab);

            // first derivative of displacement
            std::vector<Vector> u_a;
            FirstDerivativeDisplacement_a( u_a,  mDN_De[PointNumber] , CurrentDisplacement);
            
            // second derivative of displacement
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>> u_ab;
            SecondDerivativeDisplacement_ab(u_ab, D2N_De2, CurrentDisplacement);

            
            /////////////////////////////////////////////
            ///// local transformation /////////////////
            Matrix Aab(2,2);
            CovariantMetricCoefficient( Aab, A);
            //KRATOS_WATCH(Aab)


            std::vector<Vector> AA;
            ContravariantBaseVector( AA, A, Aab);

            std::vector<Vector> EE;
            LocalCartesianBasisVector(EE,  A);

            ///// 1.1 membrane strain
            Matrix eTensor(2,2);
            computeMembraneStrain(eTensor, A,  u_a );
                                           
            // curvature change
            Matrix kTensor(2,2);
            computeCurvatureChange(kTensor, A, A3, A3Vector, u_a, u_ab, A_ab);
                

            
                            // transform strain to local form
            //                Matrix eeTensor = ZeroMatrix(2,2);
            //                LocalTransformationOfTensor(eeTensor , eTensor,  EE, AA);
                            // transform tensor to vector (not necessarily)
            //                Vector eeVector(3);
            //                eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
            
             //               Matrix kkTensor = ZeroMatrix(2,2);
             //               LocalTransformationOfTensor(kkTensor , kTensor,  EE, AA);
                            // transform tensor to vector (not necessarily)
            //                Vector kkVector(3);
            //                kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);

            //                KRATOS_WATCH(kkTensor)
            //                KRATOS_WATCH(eeTensor)

            // first derivative of membrane strain
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > eLocalVector_r;
            FirstDerivativeLocalMembraneStrain_r(eLocalVector_r,  mDN_De[PointNumber], A, EE, AA) ;

            // membrane B matrix
            Matrix BBm(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBm, eLocalVector_r);
                        
            // first derivative of curvature change
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > kLocalVector_r;
            FirstDerivativeLocalCurvatureChange_r(kLocalVector_r, A, A3, A3Vector, u_a, mDN_De[PointNumber], D2N_De2,A_ab, EE, AA);
            
            // bending B matrix
            Matrix BBb(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBb, kLocalVector_r);

            Vector nVector= ZeroVector(mStrainSize);
            Vector moVector= ZeroVector(mStrainSize);

            //noalias(nVector) = mThickness*prod(D0,eeVector) ;
            //noalias(moVector) = pow(mThickness,3)/12.0*prod(D0,kkVector) ;

            for(unsigned int PointNumber1D = 0; PointNumber1D < mIntegrationPoint1D.size(); PointNumber1D++)
            {
                double theta3i = 0.5*mThickness*mIntegrationPoint1D(PointNumber1D);
                
                
                // covarianti base vector of undeformed configuration
                std::vector<Vector> G;
                ContinuumCovariantBaseVector( G,  A, A3Vector_a, A3Vector, theta3i);

                /////////////////////////////////////////////
                ///// local transformation /////////////////
                Matrix Gab(2,2);
                CovariantMetricCoefficient( Gab, G);

             

                Matrix Inv_Gab = ZeroMatrix(2,2); // non sense
                double DetGab ;
                MathUtils<double>::InvertMatrix(Gab, Inv_Gab,  DetGab);
                //KRATOS_WATCH(DetAab)
                //KRATOS_WATCH(sqrt(DetAab))
            
                std::vector<Vector> GG;
                ContravariantBaseVector( GG, G, Gab);

                std::vector<Vector> EE_c;
                LocalCartesianBasisVector(EE_c,  G);

                                                            // transform strain to local form
                                                            Matrix eeTensor = ZeroMatrix(2,2);
                                                            LocalTransformationOfTensor(eeTensor , eTensor,  EE_c, GG);
                                                            // transform tensor to vector (not necessarily)
                                                            Vector eeVector(3);
                                                            eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
                                            
                                                            Matrix kkTensor = ZeroMatrix(2,2);
                                                            LocalTransformationOfTensor(kkTensor , kTensor,  EE_c, GG);
                                                            // transform tensor to vector (not necessarily)
                                                            Vector kkVector(3);
                                                            kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);
                
                                             
                
                
                                                
                                                            noalias(nVector) += mWeight1D(PointNumber1D)*mDetJ1D*prod(D0,eeVector) ;
                                                            noalias(moVector) += mWeight1D(PointNumber1D)*mDetJ1D*pow(theta3i,2)*prod(D0,kkVector) ;

            

                // Green-Largrange strain tensor
                //Matrix GLStrainTensor = ZeroMatrix(2,2);
                //noalias(GLStrainTensor) = eTensor + theta3i*kTensor;

                // transform strain to local form
                //Matrix GLStrainTensor_local = ZeroMatrix(2,2);
                //LocalTransformationOfTensor(GLStrainTensor_local , GLStrainTensor ,  EE_c, GG);

                //Vector GLStrainVector_local = ZeroVector(3);
                //GLStrainVector_local= MathUtils<double>::StrainTensorToVector( GLStrainTensor_local, mStrainSize );

    
                // PK2 stress vector
                //Vector PK2StressVector = ZeroVector(3);
                //noalias(PK2StressVector) = prod(D0, GLStrainVector_local);


                //noalias(nVector) += mWeight1D(PointNumber1D)*mDetJ1D*PK2StressVector;
                //noalias(moVector) += mWeight1D(PointNumber1D)*mDetJ1D*theta3i*PK2StressVector;

                
            } // loop over 1D integration points
            
            if(CalculateStiffnessMatrixFlag == true)
            {
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D0, BBm, BBm, mThickness*mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddStiffnessMatrixComponents(rLeftHandSideMatrix, D0, BBb, BBb, pow(mThickness,3.0)/12.0*mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
            }

            //////////////////////////////////
            // compute the contribution to RHS
            if(CalculateResidualVectorFlag == true)
            {   
                    
                    AddInternalForces( rRightHandSideVector, nVector, eLocalVector_r, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddInternalForces( rRightHandSideVector, moVector, kLocalVector_r, mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                    AddExternalForces(rRightHandSideVector,  mN[PointNumber],mDetJ0[PointNumber], mIntegrationWeight[PointNumber] );

            }



        }// loop over integration points*/
