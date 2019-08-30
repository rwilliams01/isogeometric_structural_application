//
//   Project Name:        Kratos
//   Last Modified by:    $Author: DongGiang $
//   Date:                $Date: 25 August 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes
#include <cmath>
#include <iostream>


// External includes

// Project includes
#include "custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell_rev2.h"
#include "isogeometric_structural_application/isogeometric_structural_application.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "utilities/math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"
using namespace std;



#define ENABLE_BEZIER_GEOMETRY

//#define DEBUG1
//#define DEBUG2
//#define DEBUG3




//#define DEBUG_DKGQ

namespace Kratos
{
    //extern Variable<Vector> STRESSES;

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)
    //************************************************************************************
    //***** Constructor and Destructor ***************************************************
    //************************************************************************************
    KinematicLinearKirchhoffLoveIsogeometricShellRev2::KinematicLinearKirchhoffLoveIsogeometricShellRev2
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

    KinematicLinearKirchhoffLoveIsogeometricShellRev2::KinematicLinearKirchhoffLoveIsogeometricShellRev2(IndexType NewId,
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
    KinematicLinearKirchhoffLoveIsogeometricShellRev2::~KinematicLinearKirchhoffLoveIsogeometricShellRev2()
    {
    }


    //********************************************************
    //**** Operations ****************************************
    //********************************************************

    Element::Pointer KinematicLinearKirchhoffLoveIsogeometricShellRev2::Create(IndexType NewId,
        NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KinematicLinearKirchhoffLoveIsogeometricShellRev2(NewId,
                                GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer KinematicLinearKirchhoffLoveIsogeometricShellRev2::Create(IndexType NewId,
        GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new KinematicLinearKirchhoffLoveIsogeometricShellRev2(NewId,
                                                            pGeom, pProperties));
    }

    //************************************************************************************
    //************************************************************************************
    /**
    * returns the used integration method
    */
    KinematicLinearKirchhoffLoveIsogeometricShellRev2::IntegrationMethod KinematicLinearKirchhoffLoveIsogeometricShellRev2::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    //************************************************************************************
    //************************************************************************************
    /**
    * Setting up the EquationIdVector
    */
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
                this->GetValue(NURBS_WEIGHTS),
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

            mInvJ0[i].resize(2, 2, false);
            noalias(mInvJ0[i]) = ZeroMatrix(2, 2);

            mN[i].resize(mNumberOfNodes);
            noalias(mN[i]) = ZeroVector(mNumberOfNodes);

            mDN_De[i].resize(mNumberOfNodes, 2);
            noalias(mDN_De[i]) = ZeroMatrix(mNumberOfNodes, 2);

            mDN_DX[i].resize(mNumberOfNodes, 2);
            noalias(mDN_DX[i]) = ZeroMatrix(mNumberOfNodes, 2);
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

            Matrix JtJ = prod(trans(mJ0[PointNumber]), mJ0[PointNumber]);
            mDetJ0[PointNumber] = sqrt(MathUtils<double>::Det(JtJ));


            Matrix J0_temp (2, 2);
            J0_temp(0,0) = mJ0[PointNumber](0,0);
            J0_temp(0,1) = mJ0[PointNumber](0,1);
            J0_temp(1,0) = mJ0[PointNumber](1,0);
            J0_temp(1,1) = mJ0[PointNumber](1,1);



            MathUtils<double>::InvertMatrix(J0_temp, mInvJ0[PointNumber], DetJ_temp);

            //MathUtils<double>::InvertMatrix(mJ0[PointNumber], mInvJ0[PointNumber], DetJ_temp);

            // compute the gradient w.r.t global coordinates
            noalias(mDN_DX[PointNumber]) = prod(mDN_De[PointNumber], mInvJ0[PointNumber]);



            //getting informations for integration
            mIntegrationWeight[PointNumber] = integration_points[PointNumber].Weight();

        }



        mIsInitialized = true;


        KRATOS_CATCH( "" )
    }




    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
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
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CalculateAll( MatrixType& rLeftHandSideMatrix,
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

            double E = GetProperties()[YOUNG_MODULUS];
            double NU = GetProperties()[POISSON_RATIO];

            //2. compute membrane B operator
            Matrix Bm;
            computeMembraneBMatrix( Bm, mDN_De[PointNumber], A);

            ///// 1.1 membrane strain
            Vector eVector;
            computeStrain(eVector,  Bm,  CurrentDisplacement);

            //4. compute bending B operator
            Matrix Bb;
            computeBendingBMatrix( Bb, A3Vector, A3, A, mDN_De[PointNumber], D2N_De2,A_ab);
            //CalculateBendingBOperator(Bb,A,mDN_De[PointNumber],D2N_De2 );

            ////iii. curvature changes
            Vector kVector;
            computeStrain(kVector,  Bb,  CurrentDisplacement);

            /////////////////////////////////////////////
            ///// local transformation /////////////////
            Matrix Aab(2,2);
            CovariantMetricCoefficient( Aab, A);

            std::vector<Vector> AA;
            ContravariantBaseVector( AA, A, Aab);

            std::vector<Vector> EE;
            LocalCartesianBasisVector(EE,  A, A3Vector);

            Matrix T(3,3);
            TransformationTensor( T,  EE,  AA);

            Vector eeVector(3);
            noalias(eeVector) = prod(T, eVector);

            Vector kkVector(3);
            noalias(kkVector) = prod(T, kVector);

            Matrix BBm(mStrainSize, mNumberOfDof);
            noalias(BBm) = prod(T, Bm);

            Matrix BBb(mStrainSize, mNumberOfDof);
            noalias(BBb) = prod(T, Bb);


            /*KRATOS_WATCH(T)
            KRATOS_WATCH(eeVector)
            KRATOS_WATCH(eVector)
            KRATOS_WATCH(kkVector)
            KRATOS_WATCH(kVector)
            KRATOS_WATCH(BBm)
            KRATOS_WATCH(Bm)
            KRATOS_WATCH(BBb)
            KRATOS_WATCH(Bb)
            */

            Matrix TanC;
            CalculateElasticMatrix(TanC,  mE, mNU);

            //computeTangentMaterialStiffness(TanC, A);

            Vector nVector;  // normal forces
            computeNormalForces(nVector, TanC,  eeVector);
            //computeStrain(kVector,  Bb,  CurrentDisplacement);
            Vector moVector; //  bending moment
            computeBendingMoments(moVector, TanC, kkVector);

            ///////////////////////////////////////////////
            ///////////////////////////////////////////////

            #ifdef DEBUG2
            KRATOS_WATCH(Bm)
            KRATOS_WATCH(Bb)
            char debug2;
            cin>> debug2;
            #endif DEBUG2


            ////////////////////////// compute LHS
            if(CalculateStiffnessMatrixFlag==true)
            {
                AddLinearMembraneStiffness(rLeftHandSideMatrix, TanC, BBm,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                //AddGeometricallyMembraneStiffness(rLeftHandSideMatrix, nVector, DN_De);
                AddLinearBendingStiffness(rLeftHandSideMatrix, TanC, BBb,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
                //Add   GeometricallyBendingStiffness(rLeftHandSideMatrix, moVector, DN_De);
            }


            ////////////////////////// compute RHS
            if(CalculateResidualVectorFlag==true)
            {
                //////////// add membrane internal forces to RHS
                AddInternalForces(rRightHandSideVector, nVector, BBm,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);

                /////////// add bending internal forces to RHS
                AddInternalForces(rRightHandSideVector, moVector, BBb,  mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);

                ////////// add external forces to RHS
                AddExternalForces(rRightHandSideVector, mN[PointNumber], mDetJ0[PointNumber], mIntegrationWeight[PointNumber]);
            }

            #ifdef DEBUG3
            KRATOS_WATCH(rLeftHandSideMatrix)
            KRATOS_WATCH(rRightHandSideVector)
            char debug3;
            cin >> debug3;
            #endif
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
        //KRATOS_WATCH(rRightHandSideVector)
        //KRATOS_WATCH(rLeftHandSideMatrix)


        KRATOS_CATCH("")
    }

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::AddLinearMembraneStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC
                                                , const Matrix& Bm, const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) +=  mThickness*DetJ*Weight*prod( trans(Bm), Matrix(prod(TanC, Bm)) );
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::AddLinearBendingStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC
                                                , const Matrix& Bb,  const double& DetJ, const double& Weight)
    {
        noalias(LeftHandSideMatrix) += pow(mThickness,3)/12*DetJ*Weight*prod( trans(Bb), Matrix(prod(TanC, Bb)) );
    }

    ////////////////////////// add RHS components ////////////////////////
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::AddInternalForces(VectorType& RightHandSideVector
        , const Vector& StressResultants, const Matrix& BMatrix, const double& DetJ, const double& Weight)
    {
        noalias(RightHandSideVector) -= prod(trans(BMatrix), StressResultants)*DetJ*Weight;
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight )
    {
        const Vector BodyForce = GetProperties()[BODY_FORCE];

        for(unsigned int I=0; I < mNumberOfNodes; ++I)
            for(unsigned int i=0; i<mDim; ++i)
                RightHandSideVector(I*mDim + i) += N(I)*BodyForce(i)*DetJ*Weight;
    }


    ////////////////////////////////////////////////////////////////////////
        ///////////////////// normal forces
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeNormalForces(Vector& nVector, const Matrix& C, const Vector& eVector)
    {
        nVector = ZeroVector(mStrainSize);
        noalias(nVector) = prod(C, mThickness*eVector);
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeBendingMoments(Vector& moVector, const Matrix& C, const Vector& kVector)
    {
        moVector = ZeroVector(mStrainSize);
        noalias(moVector) = prod(C, pow(mThickness,3)/12*kVector);
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeStrain(Vector& StrainVector,  const Matrix& B,  const Matrix& Displacements)
    {
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = 3;

        StrainVector = ZeroVector(strain_size);

        for(unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
            for(unsigned int j = 0; j < strain_size; ++j)
                for(unsigned int k = 0; k < dim; ++k)
                    StrainVector[j] += B(j, i*dim + k) * Displacements(i, k);
    }

        ////////// B matrices
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeMembraneBMatrix(Matrix& Bm, const Matrix& DN_De, const std::vector<Vector>& A)
    {
        if(Bm.size1() != mStrainSize)
            Bm.resize(mStrainSize,mNumberOfDof);
        noalias(Bm) = ZeroMatrix(mStrainSize, mNumberOfDof);

        for(unsigned int I=0; I< mNumberOfNodes; ++I)
        {
            for(unsigned int i=0; i< mDim; ++i)
            {
                Bm(0, I*mDim + i) = DN_De(I, 0)*A[0](i);
                Bm(1, I*mDim + i) = DN_De(I, 1)*A[1](i);
                Bm(2, I*mDim + i) = DN_De(I,0)*A[1](i) + DN_De(I,1)*A[0](i);
            }
        }
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeBendingBMatrix(Matrix& Bb, Vector& A3Vector, double& A3,  std::vector<Vector>& A,
        const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2,
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab )
    {
        if(Bb.size1() != mStrainSize)
            Bb.resize(mStrainSize,mNumberOfDof);
        noalias(Bb) = ZeroMatrix(mStrainSize, mNumberOfDof);

        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>> R_abxA2;
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>> A1xR_ab;
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> R_abA3;

        R_abxA2.resize(2);
        A1xR_ab.resize(2);
        R_abA3.resize(2);

        for(unsigned int i=0; i< 2; ++i)
        {
            R_abxA2[i].resize(2);
            A1xR_ab[i].resize(2);
            R_abA3[i].resize(2);

            for(unsigned int j=0; j<2; ++j)
            {
                R_abxA2[i][j].resize(mDim);
                A1xR_ab[i][j].resize(mDim);
                R_abA3[i][j] = 0;
            }
        }

        for(unsigned int i=0; i<2; i++)
        {
            for(unsigned int j=0; j<2; j++)
            {
                R_abxA2[i][j] = MathUtils<double>::CrossProduct(A_ab[i][j], A[1]);
                A1xR_ab[i][j] = MathUtils<double>::CrossProduct(A[0], A_ab[i][j]);
                R_abA3[i][j] = MathUtils<double>::Dot3(A_ab[i][j], A3Vector);
            }
        }

        Vector A2xA3 = MathUtils<double>::CrossProduct(A[1], A3Vector);
        Vector A3xA1 = MathUtils<double>::CrossProduct(A3Vector, A[0]);

        for(unsigned int I=0; I< mNumberOfNodes; ++I)
        {
            for(unsigned int i=0; i< mDim; ++i)
            {
                Bb(0, I*mDim + i) = -A3Vector(i)*D2N_De2[I](0,0) + 1/A3*( R_abxA2[0][0](i)*DN_De(I,0) + A1xR_ab[0][0](i)*DN_De(I,1)
                                            + R_abA3[0][0]*( A2xA3(i)*DN_De(I,0) + A3xA1(i)*DN_De(I,1) ) );
                Bb(1, I*mDim + i) = -A3Vector(i)*D2N_De2[I](1,1) + 1/A3*( R_abxA2[1][1](i)*DN_De(I,0) + A1xR_ab[1][1](i)*DN_De(I,1)
                                            + R_abA3[1][1]*( A2xA3(i)*DN_De(I,0) + A3xA1(i)*DN_De(I,1) ) );
                Bb(2, I*mDim + i) = 2*( -A3Vector(i)*D2N_De2[I](0,1) + 1/A3*( R_abxA2[0][1](i)*DN_De(I,0) + A1xR_ab[0][1](i)*DN_De(I,1)
                                            + R_abA3[0][1]*( A2xA3(i)*DN_De(I,0) + A3xA1(i)*DN_De(I,1) ) ) );
            }
        }
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CalculateBendingBOperator(
        Matrix& Bb,
        std::vector<Vector>& G,
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
         Vector G3(3);
         Vector G1xG2;
         G1xG2=MathUtils<double>::CrossProduct(G[0],G[1]);
         G3 = G1xG2/MathUtils<double>::Norm3(G1xG2);
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

            #ifdef  DEBUG_LEVEL2
            KRATOS_WATCH(G11)
            KRATOS_WATCH(G22)
            KRATOS_WATCH(G12)
            #endif

         //////// cross product of two vectors G1,1 x G2 ; G1 x G1,1; ...
         Vector G11xG2, G22xG2, G12xG2;
         Vector G1xG11, G1xG22, G1xG12 ;
         Vector G2xG3, G3xG3 ;

         G11xG2 =MathUtils<double>::CrossProduct(G11 ,G[1]);
         G22xG2 =MathUtils<double>::CrossProduct(G22 ,G[1]);
         G12xG2 =MathUtils<double>::CrossProduct(G12 ,G[1]);
         G1xG11 =MathUtils<double>::CrossProduct(G[0] ,G11);
         G1xG22 =MathUtils<double>::CrossProduct(G[0] ,G22);
         G1xG12 =MathUtils<double>::CrossProduct(G[0] ,G12);
         G2xG3 =MathUtils<double>::CrossProduct(G[1] ,G3);
         G3xG3 =MathUtils<double>::CrossProduct(G3 ,G3);

         /////// compute Bb
         Bb(strain_size, mat_size);
         for(unsigned int i=0; i< number_of_ctrl_points; i++)
         {
             for(unsigned int j=0; j< dim; j++)
             {
                 Bb(0, i*dim+j) = -D2N_De2[i](0,0)*G3(j) + ( 1/MathUtils<double>::Norm3(G1xG2) )*( DN_De(i,0)*G11xG2(j) + DN_De(i,1)*G1xG11(j)
                                    + MathUtils<double>::Dot3(G3, G11)*( DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG3(j)) );
                 Bb(1, i*dim+j) = -D2N_De2[i](1,1)*G3(j) + ( 1/MathUtils<double>::Norm3(G1xG2) )*( DN_De(i,0)*G22xG2(j) + DN_De(i,1)*G1xG22(j)
                                    + MathUtils<double>::Dot3(G3, G22)*(  DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG3(j) ) );
                 Bb(2, i*dim+j) = 2*( -D2N_De2[i](0,1)*G3(j) + ( 1/MathUtils<double>::Norm3(G1xG2) )*( DN_De(i,0)*G12xG2(j) + DN_De(i,1)*G1xG12(j)
                                    + MathUtils<double>::Dot3(G3, G12)*(  DN_De(i,0)*G2xG3(j) + DN_De(i,1)*G3xG3(j) ) ) );
             }
         }
    }// Bending B operator

    ///////////////////// tangent stiffness
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::computeTangentMaterialStiffness(Matrix& TanC, std::vector<Vector>& A)
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

    /////////////////////////////////////////////////////////////
    //////////// shell analysis utilities //////////////////////
    ///////////////////// covariant base vectors
    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::ReferencedNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A)
    {
        AA3Vector = MathUtils<double>::CrossProduct(A[0], A[1]);
        A3 = MathUtils<double>::Norm3(AA3Vector);
        A3Vector = AA3Vector/A3;
    }

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab,
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A)
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab)
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A, const Vector& A3Vector)
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::TransformationTensor(Matrix& T, std::vector<Vector>& EE, std::vector<Vector> A)
    {
        if (T.size1() != 3)
            T.resize(3, 3);

        T(0,0) = MathUtils<double>::Dot3(EE[0], A[0])*MathUtils<double>::Dot3(A[0], EE[0]);
        T(0,1) = MathUtils<double>::Dot3(EE[0], A[1])*MathUtils<double>::Dot3(A[1], EE[0]);
        T(0,2) = MathUtils<double>::Dot3(EE[0], A[0])*MathUtils<double>::Dot3(A[1], EE[0]);

        T(1,0) = MathUtils<double>::Dot3(EE[1], A[0])*MathUtils<double>::Dot3(A[0], EE[1]);
        T(1,1) = MathUtils<double>::Dot3(EE[1], A[1])*MathUtils<double>::Dot3(A[1], EE[1]);
        T(1,2) = MathUtils<double>::Dot3(EE[1], A[0])*MathUtils<double>::Dot3(A[1], EE[1]);

        T(2,0) = MathUtils<double>::Dot3(EE[0], A[0])*MathUtils<double>::Dot3(A[0], EE[1]);
        T(2,1) = MathUtils<double>::Dot3(EE[0], A[1])*MathUtils<double>::Dot3(A[1], EE[1]);
        T(2,2) = MathUtils<double>::Dot3(EE[0], A[0])*MathUtils<double>::Dot3(A[1], EE[1]);
    }



    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::CalculateElasticMatrix(Matrix& C, const double& E, const double& NU)
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

    void KinematicLinearKirchhoffLoveIsogeometricShellRev2::UnitBaseVectors(std::vector<Vector>& e)
    {
        if (e.size() != 3)
            e.resize(3);

        e[0]=e[1]=e[2]=ZeroVector(mDim);
        e[0](0)=e[1](1)=e[2](2)= 1.0;
    }


} // Namespace Kratos
