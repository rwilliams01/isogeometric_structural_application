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
*   Date:                $Date: 29 Jan 2015 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/legacy_structural_app_vars.h"
#include "custom_conditions/master_contact_face_3D_isogeometric.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
MasterContactFace3DIsogeometric::MasterContactFace3DIsogeometric( IndexType NewId,
                                                                  GeometryType::Pointer pGeometry )
: Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
    GetValue( IS_CONTACT_SLAVE  )  = 0;
    GetValue( IS_CONTACT_MASTER )  = 1;
    
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************
MasterContactFace3DIsogeometric::MasterContactFace3DIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                                                                  PropertiesType::Pointer pProperties )
: Condition( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

Condition::Pointer MasterContactFace3DIsogeometric::Create( IndexType NewId,
                                                            NodesArrayType const& ThisNodes,
                                                            PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new MasterContactFace3DIsogeometric( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Condition::Pointer MasterContactFace3DIsogeometric::Create( IndexType NewId,
                                                            GeometryType::Pointer pGeometry,
                                                            PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new MasterContactFace3DIsogeometric( NewId, pGeometry, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
void MasterContactFace3DIsogeometric::Initialize()
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
                this->GetValue(NURBS_WEIGHT),
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

        GetValue( LAMBDAS ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( LAMBDAS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2, false );
        noalias( GetValue( LAMBDAS_T ) ) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( GAPS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS ).resize( GetGeometry().IntegrationPoints().size() , false );
        noalias( GetValue( DELTA_LAMBDAS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2, false );
        noalias( GetValue( DELTA_LAMBDAS_T ) ) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( PENALTY ) ) =  ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( PENALTY_T ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( IS_CONTACT_MASTER ) = 1;
        GetValue( IS_CONTACT_SLAVE ) = 0;
    
        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GenerateGeometryData for condition " << Id() << " completed: " << end_compute - start_compute << " s" << std::endl;
        #endif
        
    KRATOS_CATCH("")
}

/**
* returns closest point on current condition element with regard to given point in global coordinates
* @param rResultGlobal a Point in global coordinates being overwritten by the desired information
* @param rResultLocal a Point in global coordinates being overwritten by the desired information
* @param rSlaveContactGlobalPoint the point in global coordinates the closest point on the current condition element is to
* @param rCandidateGlobal the closest node to rSlaveContactGlobalPoint on current
* surface
* be calculated for
* @return true if an orthogonal projection of the given point lies within the boundaries of the current
* condition element
 */
bool MasterContactFace3DIsogeometric::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal,
                                                    GeometryType::CoordinatesArrayType& rResultLocal,
                                                    const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                                                    const GeometryType::CoordinatesArrayType& rCandidateGlobal )
{
    return false;
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void MasterContactFace3DIsogeometric::CalculateRightHandSide( VectorType& rRightHandSideVector,
                                                              ProcessInfo& rCurrentProcessInfo )
{
    unsigned int ndof = GetGeometry().size() * 3;

    if ( rRightHandSideVector.size() != ndof )
        rRightHandSideVector.resize( ndof, false );

    rRightHandSideVector = ZeroVector( ndof );
}

//************************************************************************************
//************************************************************************************

/**
 * calculates this contact element's local contributions
 */
void MasterContactFace3DIsogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                            VectorType& rRightHandSideVector,
                                                            ProcessInfo& rCurrentProcessInfo )
{
    unsigned int ndof = GetGeometry().size() * 3;

    if ( rRightHandSideVector.size() != ndof )
        rRightHandSideVector.resize( ndof, false );

    rRightHandSideVector = ZeroVector( ndof );

    if ( rLeftHandSideMatrix.size1() != ndof )
        rLeftHandSideMatrix( ndof, ndof );

    rLeftHandSideMatrix = ZeroMatrix( ndof, ndof );
}

//************************************************************************************
//************************************************************************************
/**
 * calculates the contact related contributions to the system
 * Does nothing as assembling is to be switched to linking objects
 */
void MasterContactFace3DIsogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo,
                                                    bool CalculateStiffnessMatrixFlag,
                                                    bool CalculateResidualVectorFlag )
{
}

//***********************************************************************
//***********************************************************************
/**
 * System matrix contribution due to contact energy
 * TO BE IMPLEMENTED
 */
void MasterContactFace3DIsogeometric::CalculateAndAddKc( Matrix& K,
                                                         const Vector& N,
                                                         double weight,
                                                         double dA,
                                                         double penalty,
                                                         Vector v )
{
}

//***********************************************************************
//***********************************************************************
/**
 *
 */
void MasterContactFace3DIsogeometric::CalculateAndAdd_PressureForce( Vector& residualvector,
                                                                     const Vector& N,
                                                                     Vector& v3,
                                                                     double pressure,
                                                                     double weight, double dA )
{
}

//************************************************************************************
//************************************************************************************
/**
 * REMOVED: the DOFs are managed by the linking conditions
 */
void MasterContactFace3DIsogeometric::EquationIdVector( EquationIdVectorType& rResult,
                                                        ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
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

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
/**
 * REMOVED: the DOFs are managed by the linking conditions
 */
void MasterContactFace3DIsogeometric::GetDofList( DofsVectorType& ConditionalDofList,
                                                  ProcessInfo& CurrentProcessInfo )
{
    ConditionalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

void MasterContactFace3DIsogeometric::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
                                                                   std::vector<array_1d<double, 3> >& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& sf_gradients = GetGeometry().ShapeFunctionsLocalGradients();

    if ( rVariable == NORMAL )
    {
        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //setting up result matrix
            Matrix T( 2, 3 );
            noalias( T ) = ZeroMatrix( 2, 3 );
            //shape function gradients
            Matrix DN = sf_gradients[PointNumber];
            //calculating tangential vectors

            for ( unsigned int n = 0; n < GetGeometry().PointsNumber(); n++ )
            {
                T( 0, 0 ) += ( GetGeometry().GetPoint( n ).X0()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) ) * DN( n, 0 );
                T( 0, 1 ) += ( GetGeometry().GetPoint( n ).Y0()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) ) * DN( n, 0 );
                T( 0, 2 ) += ( GetGeometry().GetPoint( n ).Z0()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) ) * DN( n, 0 );
                T( 1, 0 ) += ( GetGeometry().GetPoint( n ).X0()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) ) * DN( n, 1 );
                T( 1, 1 ) += ( GetGeometry().GetPoint( n ).Y0()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) ) * DN( n, 1 );
                T( 1, 2 ) += ( GetGeometry().GetPoint( n ).Z()
                               + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) ) * DN( n, 1 );
            }

            Vector Result( 3 );

            //calculating normal vector
            Result[0] = T( 0, 1 ) * T( 1, 2 ) - T( 0, 2 ) * T( 1, 1 );
            Result[1] = T( 0, 2 ) * T( 1, 0 ) - T( 0, 0 ) * T( 1, 2 );
            Result[2] = T( 0, 0 ) * T( 1, 1 ) - T( 0, 1 ) * T( 1, 0 );
            SD_MathUtils<double>::Normalize( Result );

            rValues[PointNumber][0] = Result[0];
            rValues[PointNumber][1] = Result[1];
            rValues[PointNumber][2] = Result[2];
        }
    }
}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int MasterContactFace3DIsogeometric::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}


/// Print information about this object.

void MasterContactFace3DIsogeometric::PrintInfo( std::ostream& rOStream ) const
{
    rOStream << "Condition #" << Id();
}

/// Print object's data.

void MasterContactFace3DIsogeometric::PrintData( std::ostream& rOStream ) const
{
    rOStream << "MasterContactFace3DIsogeometric" << std::endl;
}


} // Namespace Kratos
