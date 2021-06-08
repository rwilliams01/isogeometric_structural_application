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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 23 Oct 2014 $
*   Revision:            $Revision: 1.1 $
*
* ***************************************************************************************/


// System includes
// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "custom_conditions/line_force_isogeometric.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "structural_application/structural_application_variables.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "structural_application/structural_application.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

//#define ENABLE_PROFILING

namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
LineForceIsogeometric::LineForceIsogeometric()
{
}

// Constructor
LineForceIsogeometric::LineForceIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());

}

// Constructor
LineForceIsogeometric::LineForceIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGetGeometry());
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineForceIsogeometric::Create( IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties ) const
{
return Condition::Pointer( new LineForceIsogeometric( NewId, pGeom, pProperties ) );
}


//***********************************************************************************
//***********************************************************************************
// Destructor
LineForceIsogeometric::~LineForceIsogeometric()
{
}

//***********************************************************************************
//***********************************************************************************
void LineForceIsogeometric::EquationIdVector( EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    DofsVectorType ConditionalDofList;
    GetDofList(ConditionalDofList, rCurrentProcessInfo);

    if (rResult.size() != ConditionalDofList.size())
    rResult.resize(ConditionalDofList.size(), false);

    for(unsigned int i = 0; i < ConditionalDofList.size(); ++i)
    {
        rResult[i] = ConditionalDofList[i]->EquationId();
    }



    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineForceIsogeometric::GetDofList( DofsVectorType& ElementalDofList,
const ProcessInfo& rCurrentProcessInfo ) const 
{
    const unsigned int dim = 3;

    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************
void LineForceIsogeometric::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        ////////////////////Initialize geometry_data/////////////////////////////
        #ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        #endif

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
            manual_initilization = true;
        }

//        KRATOS_WATCH(manual_initilization)
//        KRATOS_WATCH(ExtractionOperator)

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

        // integration rule
        if(this->Has( INTEGRATION_ORDER ))
        {
            if(this->GetValue(INTEGRATION_ORDER) == 1)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 2)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 4)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "LineForceIsogeometric element does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
        }
        else if(GetProperties().Has( INTEGRATION_ORDER ))
        {
            if(GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "LineForceIsogeometric element does not support for integration points", GetProperties()[INTEGRATION_ORDER])
        }
        else
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1; // default method

        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GenerateGeometryData for condition " << Id() << " completed: " << end_compute - start_compute << " s" << std::endl;
        #endif

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void LineForceIsogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void LineForceIsogeometric::CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    rDampMatrix.resize(0,0,false);
    KRATOS_CATCH("")
}
//***********************************************************************************
//***********************************************************************************
void LineForceIsogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag,
                                bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = 3;
    unsigned int MatSize = number_of_nodes * dim;

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true )
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true )
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS
    }

    //KRATOS_WATCH(typeid(*mpIsogeometricGeometry).name())

    //reading integration points
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

    //loop over integration points
    Vector Load( dim );
    Vector temp( dim );
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        // compute Load vector at integration point
        noalias( Load ) = ZeroVector( dim );
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            noalias( temp ) = ( GetGeometry()[n] ).GetSolutionStepValue(FACE_LOAD);
//            KRATOS_WATCH(temp)

            for ( unsigned int i = 0; i < dim; ++i )
            {
                Load( i ) += temp( i ) * Ncontainer( PointNumber, n );
            }
        }

        // integration weight
        double IntegrationWeight = integration_points[PointNumber].Weight();

        // compute the tangential vector
        Vector t = ZeroVector( dim );
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 0 );
            t[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 0 );
            t[2] += GetGeometry().GetPoint( n ).Z0() * DN_De[PointNumber]( n, 0 );
        }

//        KRATOS_WATCH(t)
        //calculating length
        double dL = sqrt( t[0] * t[0] + t[1] * t[1] + t[2] * t[2] );

//        KRATOS_WATCH(Ncontainer)
//        KRATOS_WATCH(Load)
//        KRATOS_WATCH(IntegrationWeight)
//        KRATOS_WATCH(dL)

        // contribute to RHS vector
        if ( CalculateResidualVectorFlag == true )
        {
            for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
                for ( unsigned int i = 0; i < dim; ++i )
                    rRightHandSideVector( prim*dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dL;
        }
    }
//    KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
