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


//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_LEVEL3
//#define DEBUG_LEVEL4
//#define DEBUG_LEVEL5
//#define DEBUG_LEVEL6

//#define ENABLE_PROFILING


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

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
            mInvJ0[i].resize(3, 2, false);
            noalias(mInvJ0[i]) = ZeroMatrix(3, 2);

            mN[i].resize(mNumberOfNodes);
            noalias(mN[i]) = ZeroVector(mNumberOfNodes);

            mDN_De[i].resize(mNumberOfNodes, 2);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes, 2);

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


        

        // integration points and weights used for thickness calculation
        mIntegrationPoint1D.resize(3);
        mIntegrationPoint1D(0) = -std::sqrt(3.00 / 5.00) ;
        mIntegrationPoint1D(1) = 0.00;
        mIntegrationPoint1D(2) = std::sqrt(3.00 / 5.00);
               
        mWeight1D.resize(3);
        mWeight1D(0) = mWeight1D(2) = 5.00 / 9.00; mWeight1D(1) = 8.00 / 9.00;
        mDetJ1D = 0.5*mThickness;
   


        mIsInitialized = true;

        mIsIsotropicMaterial = false;

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

    void KirchhoffLoveLargeDeformationShell::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif
    
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            for(unsigned int j=0; j<3; j++)
                mConstitutiveLawVector[i][j]->FinalizeSolutionStep( GetProperties(),GetGeometry(),row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),CurrentProcessInfo );
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif
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


        #ifdef ENABLE_PROFILING
        std::vector<double> timing(10);
        std::fill(timing.begin(), timing.end(), 0.0);
        double start;
        #endif


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
        for(unsigned int node =0; node < mNumberOfNodes ;++node)
            noalias(row(CurrentDisplacement, node)) = (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);
        

        if(mIsIsotropicMaterial == true)
        {
        ///////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////// loop over integration points
        for(unsigned int PointNumber=0; PointNumber < integration_points.size(); ++PointNumber)
        {
            #ifdef ENABLE_PROFILING
            start = OpenMPUtils::GetCurrentTime();
            #endif

            //////////// get shape function values and their derivatives
            ShapeFunctionsSecondDerivativesType D2N_De2;
            D2N_De2 = mpIsogeometricGeometry->ShapeFunctionsSecondDerivatives(D2N_De2, integration_points[PointNumber]);    

            #ifdef ENABLE_PROFILING
            timing[0] += OpenMPUtils::GetCurrentTime() - start;
            start = OpenMPUtils::GetCurrentTime();
            #endif

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

            //for(unsigned int alpha =0; alpha<2; alpha++)
            //{
              
             
             //       KRATOS_WATCH(A[alpha])
             //       KRATOS_WATCH(A3Vector)
             //       KRATOS_WATCH(a[alpha])
             //       KRATOS_WATCH(a3Vector)
              
            //}

            #ifdef ENABLE_PROFILING
            timing[1] += OpenMPUtils::GetCurrentTime() - start;
            start = OpenMPUtils::GetCurrentTime();
            #endif

            ////ii. derivative of covariant base vectors
            std::vector<std::vector<Vector> > u_ab;
            std::vector<std::vector<Vector> >  a_ab;
            std::vector<std::vector<Vector> >  A_ab;
            DerivativeReferenceCovariantBaseVector(A_ab, D2N_De2, mNodalCoordinates);
            SecondDerivativeDisplacement_ab( u_ab,  D2N_De2 , CurrentDisplacement);
            DerivativeDeformedCovariantBaseVector( a_ab , A_ab,  u_ab);

            /*for(unsigned int alpha =0; alpha<2; alpha++)
            {
                for(unsigned int beta=0; beta<2; beta++)
                {
                    KRATOS_WATCH(u_ab[alpha][beta])
                    KRATOS_WATCH(A_ab[alpha][beta])
                    KRATOS_WATCH(a_ab[alpha][beta])
                }
            }*/
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

            
            // aa3_rsVector
            std::vector<std::vector<Vector>> aa3_rsVector;
            SecondDerivativeNonNormalizedDirector_rs( aa3_rsVector , a_ar);
            //SecondDerivativeNonNormalizedDirector_rs(aa3_rsVector,mDN_De[PointNumber], UnitBasisVector);

            // a3_rs
            std::vector<std::vector<double>> a3_rs;
            SecondDerivativeDirectorNorm_rs(a3_rs , aa3_rsVector,  aa3Vector, a3Vector,  a3, aa3_rVector);

            // a3_rsVector
            std::vector<std::vector<Vector>> a3_rsVector;
            SecondDerivativeDirector_rs( a3_rsVector,  aa3_rsVector, a3_rs, a3Vector, a3,  aa3_rVector, a3_r);
            
            
            
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
            Matrix eTensor;
            Matrix kTensor;
            computeMembraneStrain(eTensor, Aab,  aab  );
            computeCurvatureChange(kTensor,  Bab, bab);

            // transform strains to local form
            Matrix eeTensor = ZeroMatrix(2,2);
            Vector eeVector(3);
            LocalTransformationOfTensor(eeTensor , eTensor,  TransformationCoeff);
            eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
                            
            Matrix kkTensor = ZeroMatrix(2,2);
            Vector kkVector(3);
            LocalTransformationOfTensor(kkTensor , kTensor, TransformationCoeff);
            kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);

            // material matrix
            Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
            CalculateElasticMatrix(D0,mE, mNU);

            // stress resultants
            Vector nVector= ZeroVector(mStrainSize);
            noalias(nVector) = mThickness*prod(D0,eeVector) ;

            Vector moVector= ZeroVector(mStrainSize);
            noalias(moVector) = pow(mThickness,3.0)/12.0*prod(D0,kkVector) ;

            // first derivative of membrane strains w.r.t displacement
            std::vector<std::vector<Vector> > eLocalVector_r;
            FirstDerivativeLocalMembraneStrain_r(eLocalVector_r , a_ar, a, TransformationCoeff);

            // membrane B matrix
            Matrix BBm(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBm, eLocalVector_r);
    
            // first derivative of bending strains w.r.t displacement
            std::vector<std::vector<Vector> > kLocalVector_r;
            FirstDerivativeLocalCurvatureChange_r( kLocalVector_r,  a_abr, a_ab, a3Vector , a3_rVector, TransformationCoeff);

            // bending bending B matrix
            Matrix BBb(mStrainSize, mNumberOfDof);
            CreatingBmatrix( BBb, kLocalVector_r);

            #ifdef ENABLE_PROFILING
            timing[2] += OpenMPUtils::GetCurrentTime() - start;
            start = OpenMPUtils::GetCurrentTime();
            #endif

            ///////////////////////////////////////////////
            ///// nonlinear part of stiffness matrix //////
            ///////////////////////////////////////////////
            std::vector<std::vector<Matrix>> eLocalVector_rs;
            SecondDerivativeLocalMembraneStrain_rs(eLocalVector_rs, a_ar, TransformationCoeff);

            #ifdef ENABLE_PROFILING
            timing[3] += OpenMPUtils::GetCurrentTime() - start;
            start = OpenMPUtils::GetCurrentTime();
            #endif

            std::vector<std::vector<Matrix>> kLocalVector_rs;
            SecondDerivativeLocalCurvatureChange_rs(kLocalVector_rs , a_abr, a3_rVector,  a_ab, a3_rsVector, TransformationCoeff);

            #ifdef ENABLE_PROFILING
            timing[4] += OpenMPUtils::GetCurrentTime() - start;
            start = OpenMPUtils::GetCurrentTime();
            #endif

            if(CalculateStiffnessMatrixFlag == true)
            {
                AddLinearStiffnessMatrix(rLeftHandSideMatrix, D0, BBm, BBm, mThickness*detJA , mIntegrationWeight[PointNumber]);
                AddLinearStiffnessMatrix(rLeftHandSideMatrix, D0, BBb, BBb, pow(mThickness,3.0)/12.0*detJA , mIntegrationWeight[PointNumber]);
         
                AddNonlinearStiffnessMatrix(rLeftHandSideMatrix,nVector, eLocalVector_rs, detJA , mIntegrationWeight[PointNumber]);
                AddNonlinearStiffnessMatrix(rLeftHandSideMatrix,moVector, kLocalVector_rs,  detJA , mIntegrationWeight[PointNumber]);
            }

            if(CalculateResidualVectorFlag == true)
            {
                  
                AddInternalForces(rRightHandSideVector , nVector,eLocalVector_r, detJA , mIntegrationWeight[PointNumber]);
                AddInternalForces(rRightHandSideVector , moVector,kLocalVector_r, detJA , mIntegrationWeight[PointNumber]);
        
                ////////// add external forces to RHS
                AddExternalForces(rRightHandSideVector, mN[PointNumber],detJA , mIntegrationWeight[PointNumber]);
            }

            if(CalculateResidualVectorFlag == true)
            {
                // modify the right hand side to account for prescribed displacement
                // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
                // // However, I have to temporarily disable it to keep the consistency.
                for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                {
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                        for( unsigned int i = 0; i < mNumberOfDof; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim) * temp;
                    }
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                        for( unsigned int i = 0; i < mNumberOfDof; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 1) * temp;
                    }
                    if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) )
                    {
                        double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                        for( unsigned int i = 0; i < mNumberOfDof; ++i )
                            rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 2) * temp;
                    }
                }
            }
            



        }// loop over integration points

        } //isotropic material is chosen
        else // plasticity is chosen
        {
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
                std::vector<std::vector<Vector> >  a_ab;
                std::vector<std::vector<Vector> >  A_ab;
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
    
                
                // aa3_rsVector
                std::vector<std::vector<Vector>> aa3_rsVector;
                SecondDerivativeNonNormalizedDirector_rs( aa3_rsVector , a_ar);
    
                // a3_rs
                std::vector<std::vector<double>> a3_rs;
                SecondDerivativeDirectorNorm_rs(a3_rs , aa3_rsVector,  aa3Vector, a3Vector,  a3, aa3_rVector);
    
                // a3_rsVector
                std::vector<std::vector<Vector>> a3_rsVector;
                SecondDerivativeDirector_rs( a3_rsVector,  aa3_rsVector, a3_rs, a3Vector, a3,  aa3_rVector, a3_r);
                
                
                
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
                Matrix eTensor;
                Matrix kTensor;
                computeMembraneStrain(eTensor, Aab,  aab  );
                computeCurvatureChange(kTensor,  Bab, bab);
    
                // transform strains to local form
                Matrix eeTensor = ZeroMatrix(2,2);
                Vector eeVector(3);
                LocalTransformationOfTensor(eeTensor , eTensor,  TransformationCoeff);
                eeVector = SD_MathUtils<int>::TensorToStrainVector( eeTensor);
                                
                Matrix kkTensor = ZeroMatrix(2,2);
                Vector kkVector(3);
                LocalTransformationOfTensor(kkTensor , kTensor, TransformationCoeff);
                kkVector = SD_MathUtils<int>::TensorToStrainVector( kkTensor);
    
                // first derivative of membrane strains w.r.t displacement
                std::vector<std::vector<Vector> > eLocalVector_r;
                FirstDerivativeLocalMembraneStrain_r(eLocalVector_r , a_ar, a, TransformationCoeff);
    
                // membrane B matrix
                Matrix BBm(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBm, eLocalVector_r);
        
                // first derivative of bending strains w.r.t displacement
                std::vector<std::vector<Vector> > kLocalVector_r;
                FirstDerivativeLocalCurvatureChange_r( kLocalVector_r,  a_abr, a_ab, a3Vector , a3_rVector, TransformationCoeff);
    
                // bending bending B matrix
                Matrix BBb(mStrainSize, mNumberOfDof);
                CreatingBmatrix( BBb, kLocalVector_r);
        
                ///////////////////////////////////////////////
                ///// nonlinear part of stiffness matrix //////
                ///////////////////////////////////////////////
                std::vector<std::vector<Matrix>> eLocalVector_rs;
                SecondDerivativeLocalMembraneStrain_rs(eLocalVector_rs, a_ar, TransformationCoeff);
    
                std::vector<std::vector<Matrix>> kLocalVector_rs;
                SecondDerivativeLocalCurvatureChange_rs(kLocalVector_rs , a_abr, a3_rVector,  a_ab, a3_rsVector, TransformationCoeff);
    

                //CalculateElasticMatrix(TanC, mE, mNU);

                Vector nVector= ZeroVector(mStrainSize);
                Vector moVector= ZeroVector(mStrainSize);
                Matrix D0 = ZeroMatrix(mStrainSize,mStrainSize);
                Matrix D1 = ZeroMatrix(mStrainSize,mStrainSize);
                Matrix D2 = ZeroMatrix(mStrainSize,mStrainSize);
                    
                for(unsigned int PointNumber1D = 0; PointNumber1D < mIntegrationPoint1D.size(); PointNumber1D++)
                {
    
                    double theta3i = 0.5*mThickness*mIntegrationPoint1D(PointNumber1D);

                    std::vector<Vector> a3Vector_a;
                    std::vector<Vector> A3Vector_a;

                    //KRATOS_WATCH(a_ab)
                    
                    DerivativeDeformedNormalDirector_a(a3Vector_a, a, aa3Vector, a3, a_ab);
                    DerivativeReferenceNormalDirector_a(A3Vector_a, A , AA3Vector, A3,  A_ab);

                    // covarianti base vector of undeformed configuration
                    std::vector<Vector> g;
                    std::vector<Vector> G;  
                    ContinuumCovariantBaseVector( g,  a, a3Vector_a, a3Vector, theta3i);
                    ContinuumCovariantBaseVector( G,  A, A3Vector_a, A3Vector, theta3i);
             
                    // metric coefficents
                    Matrix Gab(2,2);
                    CovariantMetricCoefficient( Gab, G);
        
                    // contravariant base vectors
                    std::vector<Vector> GG;
                    ContravariantBaseVector( GG, G, Gab);

                    // local Cartesian basis
                    std::vector<Vector> EE_c;
                    LocalCartesianBasisVector(EE_c,  G, G[2]);


                    //cout<< "parameters are set" <<endl;
                    // set covariant basis
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(CURRENT_BASE_VECTOR_1, g[0], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(CURRENT_BASE_VECTOR_2, g[1], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(CURRENT_BASE_VECTOR_3, a3Vector, rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_BASE_VECTOR_1, G[0], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_BASE_VECTOR_2, G[1], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_BASE_VECTOR_3, G[2], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_CONTRA_BASE_VECTOR_1, GG[0], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_CONTRA_BASE_VECTOR_2, GG[1], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(REF_CONTRA_BASE_VECTOR_3, G[2], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(LOCAL_CARTESIAN_VECTOR_1, EE_c[0], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(LOCAL_CARTESIAN_VECTOR_2, EE_c[1], rCurrentProcessInfo);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->SetValue(LOCAL_CARTESIAN_VECTOR_3, EE_c[2], rCurrentProcessInfo);
          

                    
                    Matrix TanC(3,3);
                    Vector PK2StressVector(3);
                    mConstitutiveLawVector[PointNumber][PointNumber1D]->CalculateMaterialResponse(
                        ZeroVector(3),
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

                    //KRATOS_WATCH(TanC)
                    //KRATOS_WATCH(PK2StressVector)

                    ///////////////////////////////////////////////////////////////////////
                    ///////////////// stress resultants //////////////////////////////////
                    // 3. membrane stress vector
                    noalias(nVector) += PK2StressVector*mWeight1D(PointNumber1D)*mDetJ1D;
    
                    // 4. bending stress vector
                    noalias(moVector) += PK2StressVector*theta3i*mWeight1D(PointNumber1D)*mDetJ1D;

                    ////////////////////////////////////////////
                    ///////// tangent stiffness matrix /////////
                    noalias(D0) += mWeight1D(PointNumber1D)*mDetJ1D*TanC;
                    noalias(D1) += mWeight1D(PointNumber1D)*mDetJ1D*theta3i*TanC;
                    noalias(D2) += mWeight1D(PointNumber1D)*mDetJ1D*pow(theta3i,2)*TanC;    

                }

                if(CalculateStiffnessMatrixFlag == true)
                {
                    AddLinearStiffnessMatrix(rLeftHandSideMatrix, D0, BBm, BBm, detJA, mIntegrationWeight[PointNumber]);
                    AddLinearStiffnessMatrix(rLeftHandSideMatrix, D1, BBm, BBb, detJA, mIntegrationWeight[PointNumber]);
                    AddLinearStiffnessMatrix(rLeftHandSideMatrix, D1, BBb, BBm, detJA, mIntegrationWeight[PointNumber]);
                    AddLinearStiffnessMatrix(rLeftHandSideMatrix, D2, BBb, BBb, detJA, mIntegrationWeight[PointNumber]);

                    AddNonlinearStiffnessMatrix(rLeftHandSideMatrix,nVector, eLocalVector_rs, detJA , mIntegrationWeight[PointNumber]);
                    AddNonlinearStiffnessMatrix(rLeftHandSideMatrix,moVector, kLocalVector_rs,  detJA , mIntegrationWeight[PointNumber]);
                }
    
                if(CalculateResidualVectorFlag == true)
                {
                      
                    AddInternalForces(rRightHandSideVector , nVector,eLocalVector_r, detJA , mIntegrationWeight[PointNumber]);
                    AddInternalForces(rRightHandSideVector , moVector,kLocalVector_r, detJA , mIntegrationWeight[PointNumber]);
            
                    ////////// add external forces to RHS
                    AddExternalForces(rRightHandSideVector, mN[PointNumber],detJA , mIntegrationWeight[PointNumber]);
                }

                if(CalculateResidualVectorFlag == true)
                {
                    // modify the right hand side to account for prescribed displacement
                    // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
                    // // However, I have to temporarily disable it to keep the consistency.
                    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                    {
                        if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                        {
                            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                            for( unsigned int i = 0; i < mNumberOfDof; ++i )
                                rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim) * temp;
                        }
                        if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                        {
                            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                            for( unsigned int i = 0; i < mNumberOfDof; ++i )
                                rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 1) * temp;
                        }
                        if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) )
                        {
                            double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                            for( unsigned int i = 0; i < mNumberOfDof; ++i )
                                rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * mDim + 2) * temp;
                        }
                    }
                }



            }// loop over integration points

        }// plasticity is chosen

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry internal data
        mpIsogeometricGeometry->Clean();
        #endif

        #ifdef ENABLE_PROFILING
        double total_time = std::accumulate(timing.begin(), timing.end(), 0.0);
        std::cout << "Timing in percentage:" << std::endl;
        for (std::size_t i = 0; i < timing.size(); ++i)
            std::cout << " " << i << ": " << timing[i]/total_time*100 << " %" << std::endl;
        #endif

        KRATOS_CATCH("")
    }

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void KirchhoffLoveLargeDeformationShell::AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
        std::vector<std::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight)
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


    void KirchhoffLoveLargeDeformationShell::AddLinearStiffnessMatrix(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += prod( trans(BlhsMatrix), Matrix(prod(Di, BrhsMatrix)) )*DetJ*Weight;
    }

    void KirchhoffLoveLargeDeformationShell::AddNonlinearStiffnessMatrix( MatrixType& LeftHandSideMatrix,const Vector& StressResultants
        , const std::vector<std::vector<Matrix>>& StrainVector_rs, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += StressResultants(0)*StrainVector_rs[0][0]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += StressResultants(1)*StrainVector_rs[1][1]*DetJ*Weight;
        noalias(LeftHandSideMatrix) += StressResultants(2)*StrainVector_rs[0][1]*DetJ*Weight;
    }
    ////////////////////////////////////////////////////////
    ///////////////////// strain tensors ///////////////////
    ////////////////////////////////////////////////////////
    void KirchhoffLoveLargeDeformationShell::computeMembraneStrain(Matrix& eTensor,  Matrix& Aab, Matrix& aab  )
    {
        eTensor.resize(2,2);
        eTensor = ZeroMatrix(2,2);
    
        noalias(eTensor) = 0.5*(aab - Aab);
    }
    
    
    void KirchhoffLoveLargeDeformationShell::computeCurvatureChange(Matrix& kTensor, Matrix& Bab, Matrix& bab)
    {
        kTensor.resize(2,2);
        kTensor = ZeroMatrix(2,2);
    
        noalias(kTensor) = Bab - bab;
    
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLargeDeformationShell::DeformedCovariantBaseVector(std::vector<Vector>& a, std::vector<Vector>& A
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
    

    void KirchhoffLoveLargeDeformationShell::ReferenceCovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
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


    void KirchhoffLoveLargeDeformationShell::FirstDerivativeDisplacement_a(std::vector<Vector>& u_a, const Matrix& DN_De , const Matrix& u)
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

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeDisplacement_ab(std::vector<std::vector<Vector> >& u_ab
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

    void KirchhoffLoveLargeDeformationShell::NormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3,  std::vector<Vector>& a)
    {
        aa3Vector = MathUtils<double>::CrossProduct(a[0], a[1]);
        a3 = MathUtils<double>::Norm3(aa3Vector);
        a3Vector = aa3Vector/a3;
    }

  

    void KirchhoffLoveLargeDeformationShell::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
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

    void KirchhoffLoveLargeDeformationShell::DerivativeReferenceCovariantBaseVector(std::vector<std::vector<Vector> >& A_ab, 
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

    void KirchhoffLoveLargeDeformationShell::DerivativeDeformedCovariantBaseVector(std::vector<std::vector<Vector> >& a_ab
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


    void KirchhoffLoveLargeDeformationShell::DerivativeDeformedNormalDirector_a(std::vector<Vector>& a3Vector_a, std::vector<Vector>& a
        , Vector& aa3Vector, double& a3, std::vector<std::vector<Vector> >& a_ab)
    {
        if (a3Vector_a.size() != 2)
            a3Vector_a.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            a3Vector_a[i].resize(mDim);
            a3Vector_a[i] = ZeroVector(mDim);
              
        }

        for (unsigned int alpha=0; alpha<2; alpha++)
        {

            // compute aa3_a
            Vector aa3Vector_a = ZeroVector(mDim);
            noalias(aa3Vector_a) = MathUtils<double>::CrossProduct(a_ab[0][alpha], a[1])+ MathUtils<double>::CrossProduct(a[0], a_ab[1][alpha]);

            double a3_a ;
            a3_a = MathUtils<double>::Dot3(aa3Vector, aa3Vector_a)/a3;

            // compute derivative normal director w.r.t alpha
            a3Vector_a[alpha] = (a3*aa3Vector_a - aa3Vector*a3_a)*pow(a3,-2.0);
        }

    }

    void KirchhoffLoveLargeDeformationShell::DerivativeReferenceNormalDirector_a(std::vector<Vector>& A3Vector_a, std::vector<Vector>& A
        , Vector& AA3Vector, double& A3, std::vector<std::vector<Vector> >& A_ab)
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
 

    void KirchhoffLoveLargeDeformationShell::ContinuumCovariantBaseVector(std::vector<Vector>& g, std::vector<Vector>& a
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

    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void KirchhoffLoveLargeDeformationShell::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
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

    void KirchhoffLoveLargeDeformationShell::CovariantCurvatureCoefficient(Matrix& Bab
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

    //////////////////////////////////////////////////////////
    ////////////////////// material matrix ///////////////////
    //////////////////////////////////////////////////////////

    void KirchhoffLoveLargeDeformationShell::CalculateElasticMatrix(Matrix& C, const double& E, const double& NU)
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

    /////////////////////////////////////////////////////////////////////////
    ///// first derivative of membrane strain and necessary components //////
    /////////////////////////////////////////////////////////////////////////
    // variation of displacement
    void KirchhoffLoveLargeDeformationShell::DerivativeDisplacement_r(std::vector<std::vector<Vector>>& u_r, const Vector& N, std::vector<Vector>& UnitBasisVector)
    {
        u_r.resize(mNumberOfNodes);
        for(unsigned int I=0; I< mNumberOfNodes; I++)
        {
            u_r[I].resize(mDim);

            for(unsigned int i=0; i<mDim; i++)
            {
                u_r[I][i].resize(mDim);
                u_r[I][i] = ZeroVector(mDim);
            }
        }

        for(unsigned int I=0; I < mNumberOfNodes; I++)
        {
            for(unsigned int i=0; i<mDim; i++)
            {
                u_r[I][i] = N(I)*UnitBasisVector[i];
            }
        }


    }

    // variation of base vectors
    void KirchhoffLoveLargeDeformationShell::DerivativeCovariantBaseVector_r(std::vector< std::vector<std::vector<Vector>> >& a_ar, const Matrix& DN_De, std::vector<Vector>& UnitBasisVector)
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

    // second derivative of base vector w.r.t u_r
    void KirchhoffLoveLargeDeformationShell::DerivativeCovariantBaseVector_abr( std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
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

    void KirchhoffLoveLargeDeformationShell::FirstDerivativeLocalMembraneStrain_r(std::vector<std::vector<Vector> >& eLocalVector_r
        , std::vector< std::vector<std::vector<Vector>> >& a_ar,  std::vector<Vector>& a, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {
        eLocalVector_r.resize(2);
        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            eLocalVector_r[alpha].resize(2);
            for(unsigned int beta=0; beta<2; beta++)
            {
                eLocalVector_r[alpha][beta] = ZeroVector(mNumberOfDof);
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
                                // compute a_11
                                temp(I*mDim + i) += 0.5*( MathUtils<double>::Dot3(a_ar[alpha][I][i],a[beta]) + MathUtils<double>::Dot3(a_ar[beta][I][i],a[alpha]) ) *
                                TransformationCoeff[gamma][delta][alpha][beta] ;
                            }
                        }
                    }
                }

                noalias(eLocalVector_r[gamma][delta])= temp;
            }
        }

        eLocalVector_r[0][1] *= 2.0;
        eLocalVector_r[1][0] *= 2.0;
    }

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeLocalMembraneStrain_rs(std::vector<std::vector<Matrix>>& eLocalVector_rs,
        std::vector< std::vector<std::vector<Vector>> >& a_ar, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {
        eLocalVector_rs.resize(2);
        for(unsigned int alpha=0; alpha< 2; alpha++)
        {
            eLocalVector_rs[alpha].resize(2);

            for(unsigned int beta=0; beta<2; ++beta)
            {

                    eLocalVector_rs[alpha][beta].resize(mNumberOfDof, mNumberOfDof);
                    eLocalVector_rs[alpha][beta] = ZeroMatrix(mNumberOfDof, mNumberOfDof);
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
                                        // compute a_11
                                        temp(I*mDim+i, J*mDim+j) +=   0.5*( MathUtils<double>::Dot3(a_ar[alpha][I][i],a_ar[beta][J][j])
                                                    + MathUtils<double>::Dot3(a_ar[alpha][J][j],a_ar[beta][I][i]) )*TransformationCoeff[gamma][delta][alpha][beta] ;
                                    }
                                }
                            }
                        }
                    }
                }

                noalias(eLocalVector_rs[gamma][delta])= temp;
            }
        }

        eLocalVector_rs[0][1] *= 2.0;
        eLocalVector_rs[1][0] *= 2.0;
    }

    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void KirchhoffLoveLargeDeformationShell::DerivativeNonNormalizedDirector_r( std::vector<std::vector<Vector>>& aa3_rVector
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

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeNonNormalizedDirector_rs(
          std::vector< std::vector<Vector> >& aa3_rsVector
        , std::vector< std::vector<std::vector<Vector>> >& a_ar)
    {
        aa3_rsVector.resize(mNumberOfDof);
        for(unsigned int r = 0; r< mNumberOfDof; r++)
        {   
            aa3_rsVector[r].resize(mNumberOfDof);
            for(unsigned int s = 0; s< mNumberOfDof; s++)
            {
                aa3_rsVector[r][s].resize(mDim);
                aa3_rsVector[r][s] = ZeroVector(mDim);
            
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
                        aa3_rsVector[I*mDim + i][J*mDim + j] = MathUtils<double>::CrossProduct(a_ar[0][I][i], a_ar[1][J][j]) + MathUtils<double>::CrossProduct(a_ar[0][J][j], a_ar[1][I][i]);
                    }
                }
            }
        }
        

    }
    

    void KirchhoffLoveLargeDeformationShell::DerivativeDirectorNorm_r(std::vector<std::vector<double>>& a3_r
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

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeDirectorNorm_rs(
          std::vector<std::vector<double>>& a3_rs
        , std::vector<std::vector<Vector>>& aa3_rsVector
        , Vector& aa3Vector, Vector& a3Vector, double& a3,  std::vector<std::vector<Vector>>& aa3_rVector)
    {
        a3_rs.resize(mNumberOfDof);
        for(unsigned int r = 0; r< mNumberOfDof; r++)
        {   
            a3_rs[r].resize(mNumberOfDof);
            for(unsigned int s = 0; s< mNumberOfDof; s++)
            {
                a3_rs[r][s] = 0.0;
        
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
                        a3_rs[I*mDim + i][J*mDim + j] = ( MathUtils<double>::Dot3(aa3_rsVector[I*mDim+i][J*mDim +j], aa3Vector) 
                         + MathUtils<double>::Dot3(aa3_rVector[I][i], aa3_rVector[J][j])
                         - MathUtils<double>::Dot3(aa3_rVector[I][i], a3Vector)*MathUtils<double>::Dot3(aa3_rVector[J][j], a3Vector) )/a3;
                    }
                }
            }
        }
    
    }

    void KirchhoffLoveLargeDeformationShell::DerivativeDirector_r( std::vector<std::vector<Vector>>& a3_rVector
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

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeDirector_rs(
          std::vector<std::vector<Vector>>& a3_rsVector
        , std::vector<std::vector<Vector>>& aa3_rsVector
        , std::vector<std::vector<double>>& a3_rs, Vector& a3Vector, double& a3
        ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r)
    {
        a3_rsVector.resize(mNumberOfDof);
        for(unsigned int r = 0; r< mNumberOfDof; r++)
        {   
            a3_rsVector[r].resize(mNumberOfDof);
            for(unsigned int s = 0; s< mNumberOfDof; s++)
            {
                a3_rsVector[r][s].resize(mDim);
                a3_rsVector[r][s] = ZeroVector(mDim);
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
                        a3_rsVector[I*mDim + i][J*mDim +j] = (aa3_rsVector[I*mDim + i][J*mDim +j] - a3_rs[I*mDim + i][J*mDim +j]*a3Vector)/a3 
                            + (2.0*a3_r[I][i]*a3_r[J][j]*a3Vector - a3_r[I][i]*aa3_rVector[J][j] - a3_r[J][j]*aa3_rVector[I][i])*pow(a3,-2.0);
                    }
                }
            }
        }
    
    }

    void KirchhoffLoveLargeDeformationShell::FirstDerivativeLocalCurvatureChange_r(std::vector<std::vector<Vector> >& kLocalVector_r
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

    void KirchhoffLoveLargeDeformationShell::SecondDerivativeLocalCurvatureChange_rs(std::vector<std::vector<Matrix> >& kLocalVector_rs
        ,std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr, std::vector<std::vector<Vector>>& a3_rVector
        , std::vector<std::vector<Vector> >& a_ab, std::vector<std::vector<Vector>>& a3_rsVector
        , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
    {

        kLocalVector_rs.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            kLocalVector_rs[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                    kLocalVector_rs[i][j].resize(mNumberOfDof, mNumberOfDof);
                    noalias(kLocalVector_rs[i][j]) = ZeroMatrix(mNumberOfDof, mNumberOfDof);
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
                                           + MathUtils<double>::Dot3(a_ab[alpha][beta], a3_rsVector[I*mDim+i][J*mDim +j]) )
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
    void KirchhoffLoveLargeDeformationShell::CreatingBmatrix(Matrix& BMatrix, const std::vector<std::vector<Vector> >& LocalStrainVector_r)
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

    double KirchhoffLoveLargeDeformationShell::KroneckerDelta(int i, int j)
    {
        return (i=j) ? 1.0 : 0.0 ;
    }

    void KirchhoffLoveLargeDeformationShell::LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff)
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

    void KirchhoffLoveLargeDeformationShell::UnitBaseVectors(std::vector<Vector>& e)
    {
            if (e.size() != 3)
                e.resize(3);
    
            e[0]=e[1]=e[2]=ZeroVector(mDim);
            e[0](0)=e[1](1)=e[2](2)= 1.0;
    }
    
    

    void KirchhoffLoveLargeDeformationShell::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A, Vector& A3Vector)
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



    void KirchhoffLoveLargeDeformationShell::LocalTransformationCoefficient(std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff, std::vector<Vector>& EE, std::vector<Vector>& AA)
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



} // Namespace Kratos
