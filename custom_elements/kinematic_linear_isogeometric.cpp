/*
 ==============================================================================
 KratosStructuralApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel, Hoang Giang Bui
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 janosch.stascheit@rub.de
 nagel@sd.rub.de
 hgbk2008@gmail.com
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
 *   Date:                $Date: 2013 Sep 10 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/

// System includes
// External includes
#include <boost/timer.hpp>
// Project includes
#include "includes/define.h"
#include "custom_elements/kinematic_linear_isogeometric.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

//#define DEBUG_LEVEL1
// #define DEBUG_LEVEL2

#define ENUMERATE_DOF_LIST_COUPLE
//#define ENUMERATE_DOF_LIST_SEQUENTIAL //just for debugging and compare with geopde

//#define CHECK_JACOBIAN

#define ENABLE_PROFILING

// #define IGNORE_NEGATIVE_JACOBIAN

namespace Kratos
{

KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_DAMPING_ALPHA)
KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_DAMPING_BETA)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)

KinematicLinearIsogeometric::KinematicLinearIsogeometric(IndexType NewId,
        GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
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

/**
 * A simple kinematic linear 3D element for the solution
 * of the momentum balance in structural mechanics.
 * This element is used for students training at the Ruhr University Bochum.
 * Therefore it may includes comments that are obvious for the
 * experienced user.
 */
KinematicLinearIsogeometric::KinematicLinearIsogeometric(IndexType NewId,
        GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
{
    mIsInitialized = false;
//    mpIsogeometricGeometry =
//        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

Element::Pointer KinematicLinearIsogeometric::Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
            new KinematicLinearIsogeometric(NewId,
                    GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer KinematicLinearIsogeometric::Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
            new KinematicLinearIsogeometric(NewId,
                    pGeom, pProperties));
}

KinematicLinearIsogeometric::~KinematicLinearIsogeometric()
{
}

/**
 * Initialization of the element, called at the begin of each simulation.
 * Membervariables and the Material law are initialized here
 */
void KinematicLinearIsogeometric::Initialize()
{
    KRATOS_TRY		//EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )


        if (mIsInitialized)
        {
            //dimension of the problem
            unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize(mpIsogeometricGeometry->size(), dim, false);

            for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
                for (unsigned int i = 0; i < dim; ++i) // hbui edited
                    mInitialDisp(node, i) =
                        (*mpIsogeometricGeometry)[node].GetSolutionStepValue(
                            DISPLACEMENT)[i];

            return;
        }

        ///////////////////////////////////////////////////////////////
        // One time initialization
        ///////////////////////////////////////////////////////////////

        ////////////////////Initialize geometry_data/////////////////////////////
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_1))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_2))
//        KRATOS_WATCH(GetValue(NURBS_KNOTS_3))
//        KRATOS_WATCH(GetValue(NURBS_WEIGHT))
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
             and this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
             and this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
        {
            Vector rowPtr = this->GetValue( EXTRACTION_OPERATOR_CSR_ROWPTR ); // must be 0-base
            Vector colInd = this->GetValue( EXTRACTION_OPERATOR_CSR_COLIND ); // must be 0-base
            Vector values = this->GetValue( EXTRACTION_OPERATOR_CSR_VALUES );
//            int m = mpIsogeometricGeometry->size();
//            int n;
//            unsigned int dim = mpIsogeometricGeometry->NURBS_WorkingSpaceDimension();
//            if(dim == 2)
//                n = (1 + this->GetValue(NURBS_DEGREE_1)) * (1 + this->GetValue(NURBS_DEGREE_2));
//            else if(dim == 3)
//                n = (1 + this->GetValue(NURBS_DEGREE_1)) * (1 + this->GetValue(NURBS_DEGREE_2)) * (1 + this->GetValue(NURBS_DEGREE_3));
//            KRATOS_WATCH(m)
//            KRATOS_WATCH(n)
//            KRATOS_WATCH(rowPtr)
//            KRATOS_WATCH(colInd)
//            KRATOS_WATCH(values)
            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
//            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(m, n, rowPtr, colInd, values);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The extraction operator was not given for element", Id())

        // initialize the geometry
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

        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

        InitializeJacobian();
        ////////////////////End Initialize geometry_data/////////////////////////////

        //Set Up Initial displacement for StressFreeActivation of Elements
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        mInitialDisp.resize(mpIsogeometricGeometry->size(), dim, false);

        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
            for (unsigned int i = 0; i < dim; ++i)
                mInitialDisp(node, i) =
                    (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT)[i];

        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // law was at each integration point
        const GeometryType::IntegrationPointsArrayType& integration_points =
                mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);
        if (mConstitutiveLawVector.size() != integration_points.size())
        {
            mConstitutiveLawVector.resize(integration_points.size());
        }

        InitializeMaterial();

        mIsInitialized = true;

    KRATOS_CATCH( "" )
}

void KinematicLinearIsogeometric::InitializeJacobian1D()
{
    //TODO
}

void KinematicLinearIsogeometric::InitializeJacobian()
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

    // TODO remove the storage for Jacobian to save memory
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
        #ifdef IGNORE_NEGATIVE_JACOBIAN
        mTotalDomainInitialSize += fabs(mDetJ0[PointNumber]) * IntegrationWeight;
        #else
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
        #endif
    }
}

/**
 * Calculate double Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param output Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable, Vector& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    if (Output.size() != mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size())
        Output.resize(
            mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size(),
            false
        );

    for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii)
        Output[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, Output[ii]);
}

/**
 * Calculate Vector Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param output Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable, std::vector<Vector>& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    GetValueOnIntegrationPoints(rVariable, Output, rCurrentProcessInfo);
}

/**
 * Calculate Matrix Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param output Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable, std::vector<Matrix>& Output,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        //Initialize local variables
        Matrix B(strain_size, number_of_ctrl_points * dim);
        Matrix TanC(strain_size, strain_size);
        Vector StrainVector(strain_size);
        Vector StressVector(strain_size);
        Matrix DN_DX(number_of_ctrl_points, dim);
        Matrix CurrentDisp(number_of_ctrl_points, dim);

        const GeometryType::IntegrationPointsArrayType& integration_points =
                mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        //calculating shape function values and local gradients
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        if (Output.size() != integration_points.size())
            Output.resize(integration_points.size());


        //Current displacements
        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
            noalias(row(CurrentDisp, node)) =
                (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

        //Declaration of the integration weight
        //    double Weight;

        //loop over all integration points
        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            noalias(DN_DX) = prod(DN_De[PointNumber], mInvJ0[PointNumber]);
            //Initializing B_Operator at the current integration point
            CalculateBoperator(B, DN_DX);

            //calculate strain
            CalculateStrain(B, CurrentDisp, StrainVector);

            //calculate material response
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector, ZeroMatrix(1), StressVector, TanC,
                rCurrentProcessInfo, GetProperties(), (*mpIsogeometricGeometry),
                row(Ncontainer, PointNumber), true, 0, true
            );

            if (Output[PointNumber].size2() != StrainVector.size())
                Output[PointNumber].resize(1, StrainVector.size(), false);

            if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
            {
                for (unsigned int ii = 0; ii < StrainVector.size(); ++ii)
                    Output[PointNumber](0, ii) = StrainVector[ii];
            }
            else if (rVariable == PK2_STRESS_TENSOR)
            {
                for (unsigned int ii = 0; ii < StrainVector.size(); ++ii)
                    Output[PointNumber](0, ii) = StressVector[ii];
            }
            else if (rVariable == INSITU_STRESS)
            {
                Vector dummy = row(Output[PointNumber], 0);
                row(Output[PointNumber], 0) =
                    mConstitutiveLawVector[PointNumber]->GetValue(
                        INSITU_STRESS, dummy
                    );
            }
        }

    KRATOS_CATCH( "" )
}

/**
 * Initialization of the Material law at each integration point
 */
void KinematicLinearIsogeometric::InitializeMaterial()
{
    KRATOS_TRY

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
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
//				mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
//				mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
//            std::cout << "consitutive law vector " << i << " received clone" << std::endl;
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(),
                (*mpIsogeometricGeometry),
                row(Ncontainer, i)
            );
//            std::cout << "consitutive law vector " << i << " is initialized" << std::endl;
        }

    KRATOS_CATCH( "" )
}

void KinematicLinearIsogeometric::ResetConstitutiveLaw()
{
    KRATOS_TRY

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

    KRATOS_CATCH( "" )
}

/**
 * THIS is the main method here the integration in space (loop over the integration points) is done,
 * the algorithmic tangent and the (inner and outer) load vector is computed
 * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_ctrl_points*dim)*(number_of_ctrl_points*dim)
 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_ctrl_points*dim)
 * @param rCurrentProcessInfo
 * @param CalculateStiffnessMatrixFlag true: algorithmic tangent has to be computed
 * @param CalculateResidualVectorFlag true: load vector has to be computed
 */
void KinematicLinearIsogeometric::CalculateAll(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

        unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        unsigned int mat_size = number_of_ctrl_points * dim;

//        for(int i = 0; i < number_of_ctrl_points; ++i)
//        {
//            KRATOS_WATCH(mpIsogeometricGeometry->GetPoint(i).X());
//            KRATOS_WATCH(mpIsogeometricGeometry->GetPoint(i).Y());
              //Note that X() and Y() are absolute coordinate of a point (which includes displacement)
//        }

        //Initialize local variables
        Matrix B(strain_size, mat_size);
        Matrix TanC(strain_size, strain_size);
        Vector StrainVector(strain_size);
        Vector StressVector(strain_size);
        Matrix DN_DX(number_of_ctrl_points, dim);
        Matrix CurrentDisp(number_of_ctrl_points, dim);

        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            //resize the LHS=StiffnessMatrix if its size is not correct
            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }

        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size, false);

            noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
        }

        const GeometryType::IntegrationPointsArrayType& integration_points =
                mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod);

        //calculating shape function values and local gradients
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        //Current displacements
        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
            noalias(row(CurrentDisp, node)) =
                (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

        //auxiliary terms
        Vector BodyForce;

        #ifdef DEBUG_LEVEL1
        KRATOS_WATCH(integration_points.size())
        KRATOS_WATCH(mpIsogeometricGeometry->size())
        EquationIdVectorType EquationId;
        EquationIdVector(EquationId, rCurrentProcessInfo);
        std::cout << "EquationId: ";
        for(int i = 0; i < EquationId.size(); ++i)
        {
            std::cout << EquationId[i] << " ";
        }
        std::cout << std::endl;
        #endif

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
//        #pragma omp parallel for // currently bug happen on this. Probably because constitutive law is not updated correctly
        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            noalias(DN_DX) = prod(DN_De[PointNumber], mInvJ0[PointNumber]);
            //Initializing B_Operator at the current integration point
            CalculateBoperator(B, DN_DX);

            #ifdef DEBUG_LEVEL1
            KRATOS_WATCH(B)
            #endif

            //calculate strain
            CalculateStrain(B, CurrentDisp, StrainVector);

            //calculate stress
            //             mConstitutiveLawVector[PointNumber]->UpdateMaterial( StrainVector,
            //                     GetProperties(), (*mpIsogeometricGeometry), row(Ncontainer,PointNumber), rCurrentProcessInfo );
            //             mConstitutiveLawVector[PointNumber]->CalculateStress(StrainVector, StressVector);
            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            mConstitutiveLawVector[PointNumber]->SetValue(PARENT_ELEMENT_ID, this->Id(), rCurrentProcessInfo);
            mConstitutiveLawVector[PointNumber]->SetValue(INTEGRATION_POINT_INDEX, PointNumber, rCurrentProcessInfo);
            #endif

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector, ZeroMatrix(1), StressVector, TanC,
                rCurrentProcessInfo, GetProperties(), (*mpIsogeometricGeometry),
                row(Ncontainer, PointNumber), true,
                (int) CalculateStiffnessMatrixFlag, true
            );

            #ifdef DEBUG_LEVEL1
            KRATOS_WATCH(TanC)
            #endif

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight =
                integration_points[PointNumber].Weight();

            if (dim == 2)
                IntToReferenceWeight *= GetProperties()[THICKNESS];

            if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
            {
                //calculate material tangent
                //                 mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(StrainVector, TanC);
                //calculate stiffness matrix
                #ifdef IGNORE_NEGATIVE_JACOBIAN
                noalias(rLeftHandSideMatrix) += prod(trans(B), (IntToReferenceWeight * fabs(mDetJ0[PointNumber])) * Matrix(prod(TanC, B)));
                #else
                noalias(rLeftHandSideMatrix) += prod(trans(B), (IntToReferenceWeight * mDetJ0[PointNumber]) * Matrix(prod(TanC, B)));
                #endif
                //CalculateStiffnesMatrix(rLeftHandSideMatrix, TanC, msB, Weight, mDetJ0[PointNumber]);

                #ifdef DEBUG_LEVEL1
                KRATOS_WATCH(integration_points[PointNumber].Weight())
                KRATOS_WATCH(GetProperties()[THICKNESS])
                #endif
            }

            if (CalculateResidualVectorFlag == true)
            {
                //contribution to external forces
                if(GetProperties().Has(BODY_FORCE))
                {
                    BodyForce = GetProperties()[BODY_FORCE];
                    #ifdef IGNORE_NEGATIVE_JACOBIAN
                    CalculateAndAdd_ExtForceContribution(
                        row(Ncontainer, PointNumber),
                        rCurrentProcessInfo,
                        BodyForce,
                        rRightHandSideVector,
                        IntToReferenceWeight,
                        fabs(mDetJ0[PointNumber])
                    );
                    #else
                    CalculateAndAdd_ExtForceContribution(
                        row(Ncontainer, PointNumber),
                        rCurrentProcessInfo,
                        BodyForce,
                        rRightHandSideVector,
                        IntToReferenceWeight,
                        mDetJ0[PointNumber]
                    );
                    #endif
                }

                //contribution of gravity (if there is)
                if(GetProperties().Has(DENSITY) && GetProperties().Has(GRAVITY))
                {
                    #ifdef IGNORE_NEGATIVE_JACOBIAN
                    AddBodyForcesToRHS(
                        rRightHandSideVector,
                        row(Ncontainer, PointNumber),
                        IntToReferenceWeight,
                        fabs(mDetJ0[PointNumber])
                    );
                    #else
                    AddBodyForcesToRHS(
                        rRightHandSideVector,
                        row(Ncontainer, PointNumber),
                        IntToReferenceWeight,
                        mDetJ0[PointNumber]
                    );
                    #endif
                }

                #ifdef IGNORE_NEGATIVE_JACOBIAN
                AddInternalForcesToRHS(
                    rRightHandSideVector,
                    B,
                    StressVector,
                    IntToReferenceWeight,
                    fabs(mDetJ0[PointNumber])
                );
                #else
                AddInternalForcesToRHS(
                    rRightHandSideVector,
                    B,
                    StressVector,
                    IntToReferenceWeight,
                    mDetJ0[PointNumber]
                );
                #endif
            }
        } //loop over integration points

        // modify the right hand side to account for prescribed displacement
//        if(CalculateStiffnessMatrixFlag && CalculateResidualVectorFlag)
//        {
//            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
//            {
//                if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
//                {
//                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
//                    for( unsigned int i = 0; i < mat_size; ++i )
//                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
//                }
//                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
//                {
//                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
//                    for( unsigned int i = 0; i < mat_size; ++i )
//                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
//                }
//                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) && dim == 3)
//                {
//                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
//                    for( unsigned int i = 0; i < mat_size; ++i )
//                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
//                }
//            }
//        }

        #ifdef DEBUG_LEVEL2
        KRATOS_WATCH(rLeftHandSideMatrix)
        #endif

    KRATOS_CATCH( "" )
}

/**
 * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
 * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_ctrl_points*dim)
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::CalculateRightHandSide(
        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo,
        CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag
    );
}

/**
 * THIS method is called from the scheme at the start of each solution step
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::InitializeSolutionStep(
        ProcessInfo& CurrentProcessInfo)
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
        mConstitutiveLawVector[Point]->InitializeSolutionStep(
            GetProperties(),
            (*mpIsogeometricGeometry),
            row(Ncontainer, Point),
            CurrentProcessInfo
        );
    }
}

void KinematicLinearIsogeometric::InitializeNonLinearIteration(
        ProcessInfo& CurrentProcessInfo)
{
    //reset all resistant forces at node
    for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
    {
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue(REACTION_X) = 0.0;
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue(REACTION_Y) = 0.0;
        (*mpIsogeometricGeometry)[i].GetSolutionStepValue(REACTION_Z) = 0.0;
    }

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
        mConstitutiveLawVector[Point]->InitializeNonLinearIteration(
            GetProperties(),
            (*mpIsogeometricGeometry),
            row(Ncontainer, Point),
            CurrentProcessInfo
        );
    }
}

/**
 * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
 * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = true
 * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_ctrl_points*dim)*(number_of_ctrl_points*dim)
 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_ctrl_points*dim)
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
        CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag
    );
}

void KinematicLinearIsogeometric::FinalizeNonLinearIteration(
        ProcessInfo& CurrentProcessInfo)
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

    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim + 1) / 2;
    unsigned int mat_size = number_of_ctrl_points * dim;
    Matrix B(strain_size, mat_size);
    Matrix CurrentDisp(number_of_ctrl_points, dim);
    Vector StrainVector(strain_size);
    Matrix DN_DX(number_of_ctrl_points, dim);

    //Current displacements
    for (unsigned int node = 0; node < number_of_ctrl_points; ++node)
        noalias(row(CurrentDisp, node)) =
            (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

    for (unsigned int Point = 0; Point < NumberOfIntegrationPoints; ++Point)
    {
        if(mConstitutiveLawVector[Point]->Has(POST_STRAIN_VECTOR))
        {
            //Initializing B_Operator at the current integration point
            noalias( DN_DX ) = prod( DN_De[Point], mInvJ0[Point] );
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            //set the strain vector to the constitutive law
            mConstitutiveLawVector[Point]->SetValue(POST_STRAIN_VECTOR, StrainVector, CurrentProcessInfo);
        }

        mConstitutiveLawVector[Point]->FinalizeNonLinearIteration(
            GetProperties(),
            (*mpIsogeometricGeometry),
            row(Ncontainer, Point),
            CurrentProcessInfo
        );
    }
}

/**
 * THIS method is called from the scheme after each solution step, here the time step
 * start and end point variables can be transferred n --> n+1
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::FinalizeSolutionStep(
        ProcessInfo& CurrentProcessInfo)
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
}

void KinematicLinearIsogeometric::MassMatrix(MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //lumped
//        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();
//        unsigned int NumberOfNodes = mpIsogeometricGeometry->size();
//        unsigned int mat_size = dimension * NumberOfNodes;

//        if (rMassMatrix.size1() != mat_size)
//            rMassMatrix.resize(mat_size, mat_size, false);

//        rMassMatrix = ZeroMatrix(mat_size, mat_size);

//        double TotalMass = 0.0;
//        if (GetValue(USE_DISTRIBUTED_PROPERTIES))
//        {
//            TotalMass = mTotalDomainInitialSize * GetValue(DENSITY);
//        }
//        else
//            TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

//        if (dimension == 2)
//            TotalMass *= GetProperties()[THICKNESS];

//        Vector LumpFact;

//        LumpFact = mpIsogeometricGeometry->LumpingFactors(LumpFact);

//        for (unsigned int i = 0; i < NumberOfNodes; ++i)
//        {
//            double temp = LumpFact[i] * TotalMass;

//            for (unsigned int j = 0; j < dimension; ++j)
//            {
//                unsigned int index = i * dimension + j;
//                rMassMatrix(index, index) = temp;
//            }
//        }

    //consistent mass matrix
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim + 1) / 2;
    unsigned int mat_size = number_of_ctrl_points * dim;

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
void KinematicLinearIsogeometric::DampMatrix(MatrixType& rDampMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_ctrl_points * dim;

        if (rDampMatrix.size1() != mat_size)
            rDampMatrix.resize(mat_size, mat_size, false);

        noalias(rDampMatrix) = ZeroMatrix(mat_size, mat_size);

        Matrix StiffnessMatrix = ZeroMatrix(mat_size, mat_size);

        Vector RHS_Vector = ZeroVector(mat_size);

        //rayleigh damping
        CalculateAll(StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false);

        double alpha = 0.001;
        double beta = 0.001;

        if(GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
            alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        if(GetProperties().Has(RAYLEIGH_DAMPING_BETA))
            beta = GetProperties()[RAYLEIGH_DAMPING_BETA];

        MassMatrix(rDampMatrix, rCurrentProcessInfo);

        noalias(rDampMatrix) = alpha * rDampMatrix + beta * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void KinematicLinearIsogeometric::GetValuesVector(Vector& values, int Step)
{
    const unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int mat_size = number_of_ctrl_points * dim;

    if (values.size() != mat_size)
        values.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
    {
        unsigned int index = i * dim;
        values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(DISPLACEMENT_X, Step);
        values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue( DISPLACEMENT_Y, Step);

        if (dim == 3)
            values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(DISPLACEMENT_Z, Step);
    }
}

//************************************************************************************
//************************************************************************************
void KinematicLinearIsogeometric::GetFirstDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int mat_size = number_of_ctrl_points * dim;

    if (values.size() != mat_size)
        values.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
    {
        unsigned int index = i * dim;
        values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(VELOCITY_X, Step);
        values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(VELOCITY_Y, Step);

        if (dim == 3)
            values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(VELOCITY_Z, Step);
    }
}

//************************************************************************************
//************************************************************************************
void KinematicLinearIsogeometric::GetSecondDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    const unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int mat_size = number_of_ctrl_points * dim;

    if (values.size() != mat_size)
        values.resize(mat_size, false);

    for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
    {
        unsigned int index = i * dim;
        values[index] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(ACCELERATION_X, Step);
        values[index + 1] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(ACCELERATION_Y, Step);

        if (dim == 3)
            values[index + 2] = (*mpIsogeometricGeometry)[i].GetSolutionStepValue(ACCELERATION_Z, Step);
    }
}

/**
 * returns the used integration method
 */
KinematicLinearIsogeometric::IntegrationMethod KinematicLinearIsogeometric::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

/**
 * not used
 */
void KinematicLinearIsogeometric::CalculateAndAddExtForceContribution(
    const Vector& N, const ProcessInfo& CurrentProcessInfo,
    Vector& BodyForce, VectorType& rRightHandSideVector, double weight
)
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

/**
 * Informations for the builder and solver to assemble the global vectors and matrices.
 * Here a Vector containing the EquationIds of the differnt Dofs is created
 * @param rResult Vector of the EquationIds
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
)
{
    DofsVectorType ElementalDofList;
    GetDofList(ElementalDofList, rCurrentProcessInfo);

    if (rResult.size() != ElementalDofList.size())
        rResult.resize(ElementalDofList.size(), false);

    for(unsigned int i = 0; i < ElementalDofList.size(); ++i)
        rResult[i] = ElementalDofList[i]->EquationId();
}

/**
 * Informations for the builder and solver to assemble the global vectors and matrices.
 * Here a Container containing the pointers rto the DOFs of this element is created
 * @param ElementalDofList Container with of the DOFs associated with the nodes
 *                           of this element
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();

    rElementalDofList.resize(0);

    #if defined(ENUMERATE_DOF_LIST_COUPLE)
    for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
    {
        rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_Y));
        if (dim == 3)
            rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_Z));
    }
    #elif defined(ENUMERATE_DOF_LIST_SEQUENTIAL)
    for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
        rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_X));
    for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
        rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_Y));
    if(dim == 3)
    {
        for (unsigned int i = 0; i < mpIsogeometricGeometry->size(); ++i)
            rElementalDofList.push_back((*mpIsogeometricGeometry)[i].pGetDof(DISPLACEMENT_Z));
    }
    #endif
}

/**
 * Adds the Body Forces to the load vector
 * @param R RHS Vector
 * @param N_DISP shape function values at the current integration points
 * @param Weight current integration weight
 * @param detJ current Determinant of the Jacobian
 */
inline void KinematicLinearIsogeometric::AddBodyForcesToRHS(
    Vector& R,
    const Vector& N_DISP,
    double Weight,
    double detJ
)
{
    KRATOS_TRY

        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        double density = 0.0;
        Vector gravity(dim);
        noalias(gravity) = GetProperties()[GRAVITY];
        density = GetProperties()[DENSITY];

        for (unsigned int prim = 0; prim < mpIsogeometricGeometry->size(); ++prim)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                R(prim * dim + i) += N_DISP(prim) * density * gravity(i) * detJ * Weight;
            }
        }

    KRATOS_CATCH( "" )
}

inline void KinematicLinearIsogeometric::CalculateAndAdd_ExtForceContribution(
    const Vector& N, const ProcessInfo& CurrentProcessInfo,
    Vector& BodyForce, VectorType& rRightHandSideVector, double weight,
    double detJ
)
{
    KRATOS_TRY

        unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
        unsigned int dimension = mpIsogeometricGeometry->WorkingSpaceDimension();

        for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
        {
            int index = dimension * i;
            for (unsigned int j = 0; j < dimension; ++j)
                rRightHandSideVector[index + j] += weight * detJ * N[i] * BodyForce[j];
        }

    KRATOS_CATCH( "" )
}

/**
 * Adds the Internal Forces to the load vector
 * @param R RHS Vector
 * @param B_Operator B-Operator at the current integration point
 * @param StressVector current stress vector
 * @param Weight current integration weight
 * @param detJ current Determinant of the Jacobian
 */
//	void KinematicLinearIsogeometric::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
//	{
//		KRATOS_TRY
//		unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
//		unsigned int strain_size = dim * (dim + 1) / 2;
//		for ( unsigned int prim = 0; prim < mpIsogeometricGeometry->size(); ++prim )
//		{
//			for ( unsigned int i = 0; i < dim; ++i )
//			{
//				for ( unsigned int gamma = 0; gamma < strain_size; ++gamma ) // hbui edited
//				{
//					R( prim * dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) *
//							detJ * Weight ); //hbui edited
//				}
//			}
//		}
//		//         noalias(R) -= detJ*Weight* prod(trans(B_Operator),StressVector);
//		KRATOS_CATCH( "" )
//	}
void KinematicLinearIsogeometric::AddInternalForcesToRHS(
    Vector& R,
    const Matrix& B_Operator,
    Vector& StressVector,
    double Weight,
    double detJ
)
{
    KRATOS_TRY

        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        double InternalForces[dim];

        for (unsigned int prim = 0; prim < mpIsogeometricGeometry->size(); ++prim)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                InternalForces[i] = 0.0;
                for (unsigned int gamma = 0; gamma < strain_size; ++gamma)
                {
                    InternalForces[i] += B_Operator(gamma, prim * dim + i) * StressVector(gamma) * detJ * Weight;
                }

                R(prim * dim + i) -= InternalForces[i];

                if (i == 0)
                    (*mpIsogeometricGeometry)[prim].GetSolutionStepValue(REACTION_X) += InternalForces[i];

                if (i == 1)
                    (*mpIsogeometricGeometry)[prim].GetSolutionStepValue(REACTION_Y) += InternalForces[i];

                if (i == 2)
                    (*mpIsogeometricGeometry)[prim].GetSolutionStepValue(REACTION_Z) += InternalForces[i];
            }
        }

    KRATOS_CATCH( "" )
}

/**
 * Adds the Contribution of the current quadrature point to the load vector
 * @param K LHS Matrix
 * @param tan_C 6*6 algorithmic tangent of the materia law (derivation of stresses
 *               regarding strains
 * @param B_Operator B-Operator at the current integration point
 * @param Weight current integration weight
 * @param detJ current Determinant of the Jacobian
 */
void KinematicLinearIsogeometric::CalculateStiffnesMatrix(
    Matrix& K,
    const Matrix& tan_C,
    const Matrix& B_Operator,
    double Weight,
    double detJ
)
{
    KRATOS_TRY

        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        for (unsigned int prim = 0; prim < mpIsogeometricGeometry->size(); ++prim)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int sec = 0; sec < mpIsogeometricGeometry->size(); ++sec)
                {
                    for (unsigned int j = 0; j < dim; ++j)
                    {
                        for (unsigned int alpha = 0; alpha < strain_size; ++alpha)
                            for (unsigned int beta = 0; beta < strain_size; ++beta)
                                K(prim * dim + i, sec * dim + j) += B_Operator(alpha, dim * prim + i)
                                        * tan_C(alpha, beta)
                                        * B_Operator(beta, dim * sec + j)
                                        * detJ * Weight;
                    }
                }
            }
        }

        //         noalias(K)-= prod(trans(B_Operator),(Weight*detJ)*Matrix(prod(C,B_Operator)) );

    KRATOS_CATCH( "" )
}

/**
 * Computes the strain vector
 */
void KinematicLinearIsogeometric::CalculateStrain(
    const Matrix& B,
    const Matrix& Displacements,
    Vector& StrainVector
)
{
    KRATOS_TRY
        unsigned int Dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = Dim * (Dim + 1) / 2;
        noalias(StrainVector) = ZeroVector(strain_size);

        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
        {
            for (unsigned int item = 0; item < strain_size; ++item)
                for (unsigned int dim = 0; dim < Dim; ++dim)
                    StrainVector[item] += B(item, Dim * node + dim) * (Displacements(node, dim) - mInitialDisp(node, dim));
        }

    KRATOS_CATCH( "" )
}

/**
 * Computes the B-Operator at the current quadrature point
 * @param B_Operator current B-operator
 * @param DN_DX shape function values at the current integration point
 */
void KinematicLinearIsogeometric::CalculateBoperator(
    Matrix& B_Operator,
    const Matrix& DN_DX
)
{
    KRATOS_TRY

        const unsigned int number_of_ctrl_points = mpIsogeometricGeometry->PointsNumber();

        unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        if(B_Operator.size1() != strain_size || B_Operator.size2() != number_of_ctrl_points * dim)
            B_Operator.resize(strain_size, number_of_ctrl_points * dim);
        noalias(B_Operator) = ZeroMatrix(strain_size, number_of_ctrl_points * dim);

        if (dim == 2)
        {
            for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
            {
                B_Operator(0, i * 2) = DN_DX(i, 0);
                B_Operator(1, i * 2 + 1) = DN_DX(i, 1);
                B_Operator(2, i * 2) = DN_DX(i, 1);
                B_Operator(2, i * 2 + 1) = DN_DX(i, 0);
            }
        }
        else if (dim == 3)
        {
            for (unsigned int i = 0; i < number_of_ctrl_points; ++i)
            {
                B_Operator(0, i * 3) = DN_DX(i, 0);
                B_Operator(1, i * 3 + 1) = DN_DX(i, 1);
                B_Operator(2, i * 3 + 2) = DN_DX(i, 2);
                B_Operator(3, i * 3) = DN_DX(i, 1);
                B_Operator(3, i * 3 + 1) = DN_DX(i, 0);
                B_Operator(4, i * 3 + 1) = DN_DX(i, 2);
                B_Operator(4, i * 3 + 2) = DN_DX(i, 1);
                B_Operator(5, i * 3) = DN_DX(i, 2);
                B_Operator(5, i * 3 + 2) = DN_DX(i, 0);
            }
        }

    KRATOS_CATCH( "" )
}

void KinematicLinearIsogeometric::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if (rVariable == PK2_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    return;
}

/**
 * Calculate Matrix Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 /**
 * Calculate Vector Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if (rValues.size() != mConstitutiveLawVector.size())
        rValues.resize(mConstitutiveLawVector.size());

    unsigned int dim = mpIsogeometricGeometry->WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim + 1) / 2;
    unsigned int number_of_ctrl_points = mpIsogeometricGeometry->size();
    unsigned int mat_size = number_of_ctrl_points * dim;

    if (rVariable == MATERIAL_PARAMETERS)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
    }

    if (rVariable == INSITU_STRESS || rVariable == PRESTRESS)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            if (rValues[i].size() != strain_size)
                rValues[i].resize(strain_size);
            noalias(rValues[i]) = mConstitutiveLawVector[i]->GetValue(PRESTRESS, rValues[i]);
        }
    }

    if (rVariable == PLASTIC_STRAIN_VECTOR)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            if (rValues[i].size() != strain_size)
                rValues[i].resize(strain_size);
            noalias(rValues[i]) = mConstitutiveLawVector[i]->GetValue(PLASTIC_STRAIN_VECTOR, rValues[i]);
        }
    }

    if (rVariable == STRESSES)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            if (rValues[i].size() != strain_size)
                rValues[i].resize(strain_size);
            noalias(rValues[i]) = mConstitutiveLawVector[i]->GetValue(STRESSES, rValues[i]);
        }
    }

    if (rVariable == RECOVERY_STRESSES)
    {
        /////////////////////////////////////////////////////////////////////////
        //// Calculate recover stresses
        /////////////////////////////////////////////////////////////////////////
        if (GetProperties().Has(STRESS_RECOVERY_TYPE) == true)
        {
            int Type = GetProperties()[STRESS_RECOVERY_TYPE];

            if (Type == 0)
            {
                // no recovery
                GetValueOnIntegrationPoints(STRESSES, rValues, rCurrentProcessInfo);
            }
            else if (Type == 1)
            {
                // new recovery method from Bathe

                int ExpansionLevel = GetProperties()[NEIGHBOUR_EXPANSION_LEVEL];

//    		        BatheRecoverStressUtility StressUtils(ExpansionLevel);
//    		        StressUtils.CalculateImprovedStressOnIntegrationPoints( *this, rValues, rCurrentProcessInfo );
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "The stress recovery type is not supported on element", Id());

        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The stress recovery method is not defined for element", Id());
    }

    if (rVariable == INTERNAL_VARIABLES)
    {
        KRATOS_WATCH("INTERNAL_VARIABLE are disabled for this element");
    }

    if (rVariable == STRAIN)
    {
        // calculate shape function values and local gradients
        Matrix B(strain_size, mat_size);
        Vector StrainVector(strain_size);
        Matrix DN_DX(number_of_ctrl_points, dim);
        Matrix CurrentDisp(number_of_ctrl_points, dim);

        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        // extract current displacements
        for (unsigned int node = 0; node < mpIsogeometricGeometry->size(); ++node)
            noalias(row(CurrentDisp, node)) =
                (*mpIsogeometricGeometry)[node].GetSolutionStepValue(DISPLACEMENT);

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            if (rValues[i].size() != 6)
                rValues[i].resize(6);

            // compute B_Operator at the current integration point
            noalias(DN_DX) = prod(DN_De[i], mInvJ0[i]);
            CalculateBoperator(B, DN_DX);

            // compute the strain at integration point
            CalculateStrain(B, CurrentDisp, StrainVector);
            if(dim == 2)
            {
                rValues[i](0) = StrainVector(0);
                rValues[i](1) = StrainVector(1);
                rValues[i](2) = 0.0; // note: it's only correct for plane strain, TODO: we must make this available for plane stress constitutive law
                rValues[i](3) = StrainVector(2);
                rValues[i](4) = 0.0;
                rValues[i](5) = 0.0;
            }
            else if(dim == 3)
                noalias(rValues[i]) = StrainVector;
        }
    }
}

/**
 * Calculate double Variables at each integration point, used for postprocessing etc.
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector to store the values on the qudrature points, output of the method
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if (rValues.size() != mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size())
        rValues.resize(mpIsogeometricGeometry->IntegrationPoints(mThisIntegrationMethod).size());

    //reading integration points and local gradients
    for (unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point)
    {
        rValues[Point] = mConstitutiveLawVector[Point]->GetValue(rVariable, rValues[Point]);
    }
}

/**
 * Set a Matrix Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector of the values on the quadrature points
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

/**
 * Set a Vector Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector of the values on the quadrature points
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

/**
 * Set a Double Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValue value on the quadrature points
 * @param rCurrentProcessInfo
 */
void KinematicLinearIsogeometric::SetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void KinematicLinearIsogeometric::SetValueOnIntegrationPoints(
        const Kratos::Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CONSTITUTIVE_LAW)
    {
        //calculating shape functions values
        GeometryType::ShapeFunctionsGradientsType DN_De;
        Matrix Ncontainer;

        mpIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Ncontainer,
            DN_De,
            mThisIntegrationMethod
        );

        for (unsigned int i = 0; i < rValues.size(); ++i)
        {
            mConstitutiveLawVector[i] = rValues[i];
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(),
                (*mpIsogeometricGeometry),
                row(Ncontainer, i)
            );
        }
    }
}

int KinematicLinearIsogeometric::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        KRATOS_WATCH("At Check")

        unsigned int dimension = this->mpIsogeometricGeometry->WorkingSpaceDimension();

        //verify valid id
        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Invalid element ->", this->Id());
        }

        //verify valid domain size
        #ifndef IGNORE_NEGATIVE_JACOBIAN
        if (mTotalDomainInitialSize < 0)
        {
            std::stringstream ss;
            ss << "error on element -> " << this->Id() << ": ";
            ss << "Domain size can not be less than 0. Please check Jacobian. mTotalDomainInitialSize = " << mTotalDomainInitialSize;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "");
        }
        #endif

        //verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "constitutive law not provided for property ",
                    this->GetProperties().Id());
        }

        //verify that the constitutive law has the correct dimension
        if (dimension == 2)
        {
            if (this->GetProperties().Has(THICKNESS) == false)
                KRATOS_THROW_ERROR(std::logic_error,
                        "THICKNESS not provided for element", this->Id());
        }

        //check constitutive law
        int ok = 0;
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
        {
            ok = mConstitutiveLawVector[i]->Check(GetProperties(),
                    (*mpIsogeometricGeometry), rCurrentProcessInfo);
            if (ok != 0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the consitutive law", i)
            }

//			if( mConstitutiveLawVector[i]->IsIncremental() )
//				KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//			if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//				KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//			if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//				KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
        }

        //check Jacobian (just for debugging)
        //check Jacobian  should be detected by Area() ot Volume()
//        #ifdef CHECK_JACOBIAN
//        GeometryType::CoordinatesArrayType P;
//
//        P[0] = 0.0;
//        P[1] = 0.0;
//        P[2] = 0.0;
//
//        double J0 = mpIsogeometricGeometry->DeterminantOfJacobian( P );
//
//        if(J0 < 0.0)
//        {
//            KRATOS_THROW_ERROR(std::logic_error, "Negative Jacobian is detected", __FUNCTION__)
//        }
//        #endif

//        KRATOS_WATCH("Check completed")

        return ok;

    KRATOS_CATCH( "" );

}

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
#undef ENUMERATE_DOF_LIST_SEQUENTIAL
#undef ENUMERATE_DOF_LIST_COUPLE
#undef CHECK_JACOBIAN
#undef ENABLE_PROFILING
#undef IGNORE_NEGATIVE_JACOBIAN
