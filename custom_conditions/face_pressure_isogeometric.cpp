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
*   Date:                $Date: 27 Oct 2014 $
*   Revision:            $Revision: 1.0 $
*
* ***************************************************************************************/


// System includes
// External includes
#include <boost/timer.hpp>
// Project includes
#include "includes/define.h"
#include "custom_conditions/face_pressure_isogeometric.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
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
FacePressureIsogeometric::FacePressureIsogeometric()
{
}

// Constructor
FacePressureIsogeometric::FacePressureIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

// Constructor
FacePressureIsogeometric::FacePressureIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer FacePressureIsogeometric::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FacePressureIsogeometric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer FacePressureIsogeometric::Create( IndexType NewId,
                                        GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FacePressureIsogeometric( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
FacePressureIsogeometric::~FacePressureIsogeometric()
{
}

//***********************************************************************************
//***********************************************************************************
void FacePressureIsogeometric::EquationIdVector( EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo )
{
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
}

//***********************************************************************************
//***********************************************************************************
void FacePressureIsogeometric::GetDofList( DofsVectorType& ElementalDofList,
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
void FacePressureIsogeometric::Initialize()
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
             and this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
             and this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
        {
            Vector rowPtr = this->GetValue( EXTRACTION_OPERATOR_CSR_ROWPTR ); // must be 0-base
            Vector colInd = this->GetValue( EXTRACTION_OPERATOR_CSR_COLIND ); // must be 0-base
            Vector values = this->GetValue( EXTRACTION_OPERATOR_CSR_VALUES );
            ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
            manual_initilization = true;
        }
//        else
//            KRATOS_THROW_ERROR(std::logic_error, "The extraction operator was not given for element", Id())
//        KRATOS_WATCH(ExtractionOperator)

        // initialize the geometry
        if(manual_initilization)
        {
            int num_integration_method = 2; // by default compute two integration rules
            if( GetProperties().Has(NUM_IGA_INTEGRATION_METHOD) )
                num_integration_method = GetProperties()[NUM_IGA_INTEGRATION_METHOD];
            mpIsogeometricGeometry->AssignGeometryData(
                this->GetValue(NURBS_KNOTS_1),
                this->GetValue(NURBS_KNOTS_2),
                this->GetValue(NURBS_KNOTS_3),
                this->GetValue(NURBS_WEIGHT),
                ExtractionOperator,
                this->GetValue(NURBS_DEGREE_1),
                this->GetValue(NURBS_DEGREE_2),
                this->GetValue(NURBS_DEGREE_3),
                num_integration_method
            );
        }

        mThisIntegrationMethod =
//            GetGeometry().GetDefaultIntegrationMethod(); //default method
            GeometryData::GI_GAUSS_1;
        
        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GenerateGeometryData for condition " << Id() << " completed: " << end_compute - start_compute << " s" << std::endl;
        #endif
        
    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void FacePressureIsogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void FacePressureIsogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag,
                                bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = mpIsogeometricGeometry->size();
    const unsigned int dim = 3;
    unsigned int mat_size = number_of_nodes * dim;

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true )
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true )
    {
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

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
    Vector t1(3), t2(3), v3(3);
    const double& P = this->GetValue(PRESSURE);
//    KRATOS_WATCH(P)
//    KRATOS_WATCH(Ncontainer)

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        t1 = ZeroVector( 3 );//first tangential vector
        t2 = ZeroVector( 3 );//second tangential vector
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t1[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 0 );
            t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 0 );
            t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_De[PointNumber]( n, 0 );
            t2[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 1 );
            t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 1 );
            t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_De[PointNumber]( n, 1 );
        }

        //calculating normal
        v3[0] = t1[1] * t2[2] - t1[2] * t2[1];

        v3[1] = t1[2] * t2[0] - t1[0] * t2[2];

        v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

        double dA = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] );
        v3 *= 1.0 / dA; //normalize v3
        
        //calculating load vector
        Vector Load = P * v3;

        double IntegrationWeight = integration_points[PointNumber].Weight();

//        KRATOS_WATCH(Load)
//        KRATOS_WATCH(IntegrationWeight)
//        KRATOS_WATCH(dA)

        // RIGHT HAND SIDE VECTOR
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
                for ( unsigned int i = 0; i < dim; ++i )
                    rRightHandSideVector( prim*dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dA;
        }
    }

//    if(CalculateResidualVectorFlag)
//        KRATOS_WATCH(rRightHandSideVector)

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.

#undef ENABLE_PROFILING

