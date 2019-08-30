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
 *   Date:                $Date: 2013 Dec 15 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/

// System includes
// External includes
#include <boost/timer.hpp>
// Project includes
#include "includes/define.h"
#include "custom_elements/kinematic_linear_nurbs.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2

#define ENUMERATE_DOF_LIST_COUPLE
//#define ENUMERATE_DOF_LIST_SEQUENTIAL //just for debugging and compare with geopde

//#define CHECK_JACOBIAN

#define ENABLE_PROFILING

namespace Kratos
{

KinematicLinearNURBS::KinematicLinearNURBS(IndexType NewId,
        GeometryType::Pointer pGeometry) :
        KinematicLinear2(NewId, pGeometry)
{
    KRATOS_WATCH("At KinematicLinearNURBS Default Constructor")
    mIsInitialized = false;
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

/**
 * A simple kinematic linear 3D element for the solution
 * of the momentum balance in structural mechanics.
 * This element is used for students training at the Ruhr University Bochum.
 * Therefore it may includes comments that are obvious for the
 * experienced user.
 * hbui modified to work in 2D
 */
KinematicLinearNURBS::KinematicLinearNURBS(IndexType NewId,
        GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        KinematicLinear2(NewId, pGeometry, pProperties)
{
    KRATOS_WATCH("At KinematicLinearNURBS Constructor")
    mIsInitialized = false;
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

Element::Pointer KinematicLinearNURBS::Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
    KRATOS_WATCH("At KinematicLinearNURBS::Create")
    GeometryType::Pointer pGeometry = GetGeometry().Create(ThisNodes);
    KRATOS_WATCH("At KinematicLinearNURBS::Create created geometry")
    KRATOS_WATCH(*pProperties)
    Element::Pointer pElement = Element::Pointer(
            new KinematicLinearNURBS(NewId, pGeometry, pProperties));
    KRATOS_WATCH("KinematicLinearNURBS::Create completed")
    return pElement;
}

Element::Pointer KinematicLinearNURBS::Create(IndexType NewId,
        NodesArrayType const& ThisNodes) const
{
    return Element::Pointer(
            new KinematicLinearNURBS(NewId,
                    GetGeometry().Create(ThisNodes)));
}

KinematicLinearNURBS::~KinematicLinearNURBS()
{
}

/**
 * Initialization of the element, called at the begin of each simulation.
 * Membervariables and the Material law are initialized here
 */
void KinematicLinearNURBS::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY		//EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

        KRATOS_WATCH("At Initialize")

        if (mIsInitialized)
        {
            //dimension of the problem
            unsigned int dim = GetGeometry().WorkingSpaceDimension();

            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize(GetGeometry().size(), dim, false);

            for (unsigned int node = 0; node < GetGeometry().size(); node++)
                for (unsigned int i = 0; i < dim; i++) // hbui edited
                    mInitialDisp(node, i) =
                            GetGeometry()[node].GetSolutionStepValue(
                                    DISPLACEMENT)[i];

            return;
        }

        ///////////////////////////////////////////////////////////////
        // One time initialisation
        ///////////////////////////////////////////////////////////////

        ////////////////////Initialize geometry_data/////////////////////////////
        #ifdef ENABLE_PROFILING
        double start_compute, end_compute;
        start_compute = OpenMPUtils::GetCurrentTime();
        #endif

        //initialize the nurbs geometry
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

        mThisIntegrationMethod = mpIsogeometricGeometry->GetDefaultIntegrationMethod(); //default method

        InitializeJacobian();

        #ifdef ENABLE_PROFILING
        end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GenerateGeometryData for element " << Id() << " completed: " << end_compute - start_compute << " s" << std::endl;
        #endif
        ////////////////////End Initialize geometry_data/////////////////////////////

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize(GetGeometry().size(), dim, false);

        for (unsigned int node = 0; node < GetGeometry().size(); node++)
            for (unsigned int i = 0; i < dim; ++i)
                mInitialDisp(node, i) =
                        GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT)[i];

        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // law was at each integration point
        const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints(mThisIntegrationMethod);
        if (mConstitutiveLawVector.size() != integration_points.size())
        {
            mConstitutiveLawVector.resize(integration_points.size());
        }

        InitializeMaterial();

        mIsInitialized = true;

    KRATOS_CATCH( "" )
}

void KinematicLinearNURBS::InitializeJacobian1D()
{
    //TODO
}

void KinematicLinearNURBS::InitializeJacobian()
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //number of integration points used, mThisIntegrationMethod refers to the
    //integration method defined in the constructor
    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints(mThisIntegrationMethod);

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
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);

    //calculating the inverse Jacobian
    #pragma omp parallel for
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix(J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);
        //calculating the total area
        #pragma omp atomic
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }
}

int KinematicLinearNURBS::Check(const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "Element found with Id 0 or negative", "");
        }

        if (this->GetGeometry().Length() < 0)
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Length can not be less than 0. Please check Jacobian.", __FUNCTION__);
        }

        if (this->GetGeometry().Area() < 0)
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area can not be less than 0. Please check Jacobian.", __FUNCTION__);
        }

        if (this->GetGeometry().Volume() < 0)
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Volume can not be less than 0. Please check Jacobian.", __FUNCTION__);
        }

        //verify that the constitutive law exists
        if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "constitutive law not provided for property ",
                    this->GetProperties().Id());
        }

        //Verify that the body force is defined
        if (this->GetProperties().Has(BODY_FORCE) == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,
                    "BODY_FORCE not provided for property ",
                    this->GetProperties().Id())
        }

        //verify that the constitutive law has the correct dimension
        if (dimension == 2)
        {
            if (this->GetProperties().Has(THICKNESS) == false)
                KRATOS_THROW_ERROR(std::logic_error,
                        "THICKNESS not provided for element ", this->Id());

            if (this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() != 3)
                KRATOS_THROW_ERROR(std::logic_error,
                        "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ",
                        this->Id());
        }
        else if (dimension == 3)
        {
            if (this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() != 6)
                KRATOS_THROW_ERROR(std::logic_error,
                        "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ",
                        this->Id());
        }

        //check constitutive law
        int ok = 0;
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            ok = mConstitutiveLawVector[i]->Check(GetProperties(),
                    GetGeometry(), rCurrentProcessInfo);
            if (ok != 0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the consitutive law", __FUNCTION__)
            }

//			if( mConstitutiveLawVector[i]->IsIncremental() )
//				KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//			if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//				KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//			if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//				KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
        }

        //check Jacobian
        #ifdef CHECK_JACOBIAN
        GeometryType::CoordinatesArrayType P;

        P[0] = 0.0;
        P[1] = 0.0;
        P[2] = 0.0;

        double J0 = GetGeometry().DeterminantOfJacobian( P );

        if(J0 < 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Negative Jacobian is detected", __FUNCTION__)
        }
        #endif

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

