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
*   Date:                $Date: 25 Jan 2015 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/legacy_structural_app_vars.h"
#include "custom_conditions/master_contact_face_3D_isogeometric.h"
#include "custom_conditions/slave_contact_face_3D_isogeometric.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

extern Variable<double> CONTACT_PENETRATION;

//************************************************************************************
//************************************************************************************
SlaveContactFace3DIsogeometric::SlaveContactFace3DIsogeometric( IndexType NewId,
                                                                GeometryType::Pointer pGeometry )
: Condition( NewId, pGeometry )
{
    mpIsogeometricGeometry = 
        boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************
SlaveContactFace3DIsogeometric::SlaveContactFace3DIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                                                                PropertiesType::Pointer pProperties )
: Condition( NewId, pGeometry, pProperties )
{
    mpIsogeometricGeometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(pGeometry);
}

Condition::Pointer SlaveContactFace3DIsogeometric::Create( IndexType NewId,
                                                           NodesArrayType const& ThisNodes,
                                                           PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new SlaveContactFace3DIsogeometric( NewId, (*mpIsogeometricGeometry).Create( ThisNodes ), pProperties ) );
}

Condition::Pointer SlaveContactFace3DIsogeometric::Create( IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new SlaveContactFace3DIsogeometric( NewId, pGeometry, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
void SlaveContactFace3DIsogeometric::Initialize()
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

        mpMasterElements = ContactMasterContainerType::Pointer( new ContactMasterContainerType() );
        GetValue( LAMBDAS ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        noalias( GetValue( LAMBDAS ) ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( LAMBDAS_T ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), 2 , false );
        noalias( GetValue( LAMBDAS_T ) ) = ZeroMatrix( (*mpIsogeometricGeometry).IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        noalias( GetValue( GAPS ) ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        noalias( GetValue( DELTA_LAMBDAS ) ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), 2, false );
        noalias( GetValue( DELTA_LAMBDAS_T ) ) = ZeroMatrix( (*mpIsogeometricGeometry).IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size() , false );
        noalias( GetValue( PENALTY ) ) =  ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size() , false );
        noalias( GetValue( PENALTY_T ) ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( IS_CONTACT_SLAVE ) = 1;
        GetValue( IS_CONTACT_MASTER ) = 0;

        GetValue( NORMAL_STRESS ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        GetValue( NORMAL_STRESS ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( TANGENTIAL_STRESS ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        GetValue( TANGENTIAL_STRESS ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
        GetValue( STICK ).resize( (*mpIsogeometricGeometry).IntegrationPoints().size(), false );
        GetValue( STICK ) = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );

        GetValue( NORMAL_CONTACT_STRESS ) = 0.0;
        GetValue( TANGENTIAL_CONTACT_STRESS ) = 0.0;
        GetValue( CONTACT_STICK ) = 0.0;
    
        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GenerateGeometryData for condition " << Id() << " completed: " << end_compute - start_compute << " s" << std::endl;
        #endif
        
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void SlaveContactFace3DIsogeometric::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                                                   Vector& Output,
                                                                   const ProcessInfo& rCurrentProcessInfo )
{
    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = (*mpIsogeometricGeometry).IntegrationPoints();

    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size(), false );

    double result = 0.0;

    double result_friction = 0.0;

    double reference = 0.0;

    if ( rVariable == NORMAL_CONTACT_STRESS || rVariable == TANGENTIAL_CONTACT_STRESS )
    {
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            Matrix TSlave = TangentialVectors_inOrigin( integration_points[PointNumber] );

            Vector vSlaveNonNormalized = ZeroVector( 3 );

            vSlaveNonNormalized[0] = TSlave( 0, 1 ) * TSlave( 1, 2 ) - TSlave( 0, 2 ) * TSlave( 1, 1 );
            vSlaveNonNormalized[1] = TSlave( 0, 2 ) * TSlave( 1, 0 ) - TSlave( 0, 0 ) * TSlave( 1, 2 );
            vSlaveNonNormalized[2] = TSlave( 0, 0 ) * TSlave( 1, 1 ) - TSlave( 0, 1 ) * TSlave( 1, 0 );

            double dASlave = MathUtils<double>::Norm3( vSlaveNonNormalized );

            result += ( this->GetValue( NORMAL_STRESS )[PointNumber] ) * ( (*mpIsogeometricGeometry).IntegrationPoints()[PointNumber].Weight() ) * dASlave;
            result_friction += ( this->GetValue( TANGENTIAL_STRESS )[PointNumber] ) * ( (*mpIsogeometricGeometry).IntegrationPoints()[PointNumber].Weight() ) * dASlave;
            reference += (*mpIsogeometricGeometry).IntegrationPoints()[PointNumber].Weight() * dASlave;
        }
    }

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        if ( rVariable == NORMAL_CONTACT_STRESS )
        {

            Output[PointNumber] = result / reference;

        }

        if ( rVariable == TANGENTIAL_CONTACT_STRESS )
        {
            Output[PointNumber] =  result_friction / reference;
        }

        if ( rVariable == CONTACT_PENETRATION )
        {
            //Output[PointNumber]=  this->GetValue(STICK)[PointNumber] ;
            if( this->GetValue( GAPS )[PointNumber] > 0.0 )
            {
                Output[PointNumber] = this->GetValue( GAPS )[PointNumber];
            }
            else
                Output[PointNumber] =  0.0 ;
        }
    }
}

void SlaveContactFace3DIsogeometric::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                                  std::vector<double>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo )
{
    Vector result = ZeroVector( (*mpIsogeometricGeometry).IntegrationPoints().size() );
    CalculateOnIntegrationPoints( rVariable, result, rCurrentProcessInfo );

    if ( rValues.size() != (*mpIsogeometricGeometry).IntegrationPoints().size() )
        rValues.resize( (*mpIsogeometricGeometry).IntegrationPoints().size() );

    for ( unsigned int i = 0; i < result.size(); i++ )
        rValues[i] = result[i];
}

//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void SlaveContactFace3DIsogeometric::CalculateRightHandSide( VectorType& rRightHandSideVector,
                                                             ProcessInfo& rCurrentProcessInfo )
{
    unsigned int ndof = (*mpIsogeometricGeometry).size() * 3;

    if ( rRightHandSideVector.size() != ndof )
        rRightHandSideVector.resize( ndof, false );

    rRightHandSideVector = ZeroVector( ndof );
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void SlaveContactFace3DIsogeometric::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                           VectorType& rRightHandSideVector,
                                                           ProcessInfo& rCurrentProcessInfo )
{
    unsigned int ndof = (*mpIsogeometricGeometry).size() * 3;

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
 * does nothing as assembling is switched to link condition
 */
void SlaveContactFace3DIsogeometric::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                                   VectorType& rRightHandSideVector,
                                                   ProcessInfo& rCurrentProcessInfo,
                                                   bool CalculateStiffnessMatrixFlag,
                                                   bool CalculateResidualVectorFlag )
{
} // CalculateAll

//***********************************************************************
//***********************************************************************
/**
 * System matrix contribution due to contact energy
 * TODO: implement mixed elementary contribution
 */
void SlaveContactFace3DIsogeometric::CalculateAndAddKc( Matrix& K,
                                                        const Vector& N,
                                                        double weight,
                                                        double dA,
                                                        Vector v )
{
}

//***********************************************************************
//***********************************************************************
/**
 * TO BE TESTED!!!
 */
void SlaveContactFace3DIsogeometric::CalculateAndAdd_PressureForce( Vector& residualvector,
                                                                    const Vector& N,
                                                                    Vector& v3,
                                                                    double pressure,
                                                                    double weight,
                                                                    double DetJ )
{
}

//************************************************************************************
//************************************************************************************
/**
 * REMOVED: the DOFs are managed by the linking conditions
 */
void SlaveContactFace3DIsogeometric::EquationIdVector( EquationIdVectorType& rResult,
                                                       ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    unsigned int number_of_nodes = (*mpIsogeometricGeometry).size();
    unsigned int dim = number_of_nodes * 3;

    if ( rResult.size() != dim )
        rResult.resize( dim );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = (*mpIsogeometricGeometry)[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void SlaveContactFace3DIsogeometric::MasterElementsEquationIdVectors( EquationIdVectorContainerType& rResult,
                                                                      ProcessInfo& rCurrentProcessInfo )
{
}

//************************************************************************************
//************************************************************************************
/**
 * REMOVED: the DOFs are managed by the linking conditions
 */
void SlaveContactFace3DIsogeometric::GetDofList( DofsVectorType& ConditionalDofList,
                                                 ProcessInfo& CurrentProcessInfo )
{
    ConditionalDofList.resize( 0 );

    for ( unsigned int i = 0; i < (*mpIsogeometricGeometry).size(); i++ )
    {
        ConditionalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_Y ) );
        ConditionalDofList.push_back( (*mpIsogeometricGeometry)[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

Matrix SlaveContactFace3DIsogeometric::TangentialVectors_inOrigin( const GeometryType::CoordinatesArrayType& rPoint )
{
    //setting up result matrix
    Matrix T = ZeroMatrix( 2, 3 );

    //shape function gradients
    Matrix DN = ZeroMatrix( (*mpIsogeometricGeometry).PointsNumber(), 2 );
    (*mpIsogeometricGeometry).ShapeFunctionsLocalGradients( DN, rPoint );

    //calculating tangential vectors
    for ( unsigned int n = 0; n < (*mpIsogeometricGeometry).PointsNumber(); ++n )
    {
        T( 0, 0 ) += (*mpIsogeometricGeometry).GetPoint( n ).X0() * DN( n, 0 );
        T( 0, 1 ) += (*mpIsogeometricGeometry).GetPoint( n ).Y0() * DN( n, 0 );
        T( 0, 2 ) += (*mpIsogeometricGeometry).GetPoint( n ).Z0() * DN( n, 0 );
        T( 1, 0 ) += (*mpIsogeometricGeometry).GetPoint( n ).X0() * DN( n, 1 );
        T( 1, 1 ) += (*mpIsogeometricGeometry).GetPoint( n ).Y0() * DN( n, 1 );
        T( 1, 2 ) += (*mpIsogeometricGeometry).GetPoint( n ).Z0() * DN( n, 1 );
    }

    return( T );
}

/**
* This function provides the place to perform checks on the completeness of the input.
* It is designed to be called only once (or anyway, not often) typically at the beginning
* of the calculations, so to verify that nothing is missing from the input
* or that no common error is found.
* @param rCurrentProcessInfo
*/
int SlaveContactFace3DIsogeometric::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

/// Print information about this object.
void SlaveContactFace3DIsogeometric::PrintInfo( std::ostream& rOStream ) const
{
    rOStream << "Condition #" << Id();
}

/// Print object's data.
void SlaveContactFace3DIsogeometric::PrintData( std::ostream& rOStream ) const
{
    rOStream << "SlaveContactFace3DIsogeometric" << std::endl;
}

} // Namespace Kratos
