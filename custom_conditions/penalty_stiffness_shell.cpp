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
#include "custom_conditions/penalty_stiffness_shell.h"
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
PenaltyStiffnessShell::PenaltyStiffnessShell()
{
}

// Constructor
PenaltyStiffnessShell::PenaltyStiffnessShell( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

// Constructor
PenaltyStiffnessShell::PenaltyStiffnessShell( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry =
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer PenaltyStiffnessShell::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new PenaltyStiffnessShell( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer PenaltyStiffnessShell::Create( IndexType NewId,
                                        GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new PenaltyStiffnessShell( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
PenaltyStiffnessShell::~PenaltyStiffnessShell()
{
}

//***********************************************************************************
//***********************************************************************************
void PenaltyStiffnessShell::EquationIdVector( EquationIdVectorType& rResult,
                                    const ProcessInfo& rCurrentProcessInfo ) const
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
void PenaltyStiffnessShell::GetDofList( DofsVectorType& ElementalDofList,
                              const ProcessInfo& rCurrentProcessInfo ) const
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
void PenaltyStiffnessShell::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
//        else
//            KRATOS_THROW_ERROR(std::logic_error, "The extraction operator was not given for element", Id())
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
void PenaltyStiffnessShell::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo )
{

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );
}

//***********************************************************************************
//***********************************************************************************
void PenaltyStiffnessShell::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY


    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int mat_size = 6;

    //resizing as needed the LHS
    if ( rLeftHandSideMatrix.size1() != mat_size )
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS

    //resizing as needed the RHS
    if ( rRightHandSideVector.size() != mat_size )
        rRightHandSideVector.resize( mat_size, false );

    rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS

     // Current displacements
     Matrix CurrentDisp(number_of_nodes, 3);
     Vector DispVector(number_of_nodes*3);
     // extract current displacements
     for(unsigned int node = 0; node < number_of_nodes; ++node)
        noalias(row(CurrentDisp, node)) = GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

     for(unsigned int i = 0; i < number_of_nodes; ++i)
        for(unsigned int j = 0; j < 3; ++j)
             DispVector(i*3+j) = CurrentDisp(i,j);

    const double& W = GetProperties()[INITIAL_PENALTY];

    rLeftHandSideMatrix(0,0) = W;
    rLeftHandSideMatrix(0,3) = -W;
    rLeftHandSideMatrix(1,1) = W;
    rLeftHandSideMatrix(1,4) = -W;
    rLeftHandSideMatrix(2,2) = W;
    rLeftHandSideMatrix(2,5) = -W;
    rLeftHandSideMatrix(3,0) = -W;
    rLeftHandSideMatrix(3,3) = W;
    rLeftHandSideMatrix(4,1) = -W;
    rLeftHandSideMatrix(4,4) = W;
    rLeftHandSideMatrix(5,2) = -W;
    rLeftHandSideMatrix(5,5) = W;

    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,DispVector);

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.

#undef ENABLE_PROFILING

