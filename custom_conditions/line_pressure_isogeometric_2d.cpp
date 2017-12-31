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
*   Date:                $Date: 24 Oct 2014 $
*   Revision:            $Revision: 1.1 $
*
* ***************************************************************************************/


// System includes
// External includes
#include <boost/timer.hpp>
// Project includes
#include "includes/define.h"
#include "line_pressure_isogeometric_2d.h"
#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"


//#define DEBUG_LEVEL1


namespace Kratos
{
//***********************************************************************************
//***********************************************************************************
// -------- //
//  PUBLIC  //
// -------- //

// Constructor
LinePressureIsogeometric2D::LinePressureIsogeometric2D()
{
}

// Constructor
LinePressureIsogeometric2D::LinePressureIsogeometric2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : LineForceIsogeometric( NewId, pGeometry )
{
}

// Constructor
LinePressureIsogeometric2D::LinePressureIsogeometric2D( IndexType NewId, GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties )
    : LineForceIsogeometric( NewId, pGeometry, pProperties )
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LinePressureIsogeometric2D::Create( IndexType NewId,
                                        NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LinePressureIsogeometric2D( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LinePressureIsogeometric2D::Create( IndexType NewId,
                                        GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LinePressureIsogeometric2D( NewId, pGeom, pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
LinePressureIsogeometric2D::~LinePressureIsogeometric2D()
{
}

//***********************************************************************************
//***********************************************************************************
void LinePressureIsogeometric2D::GetDofList( DofsVectorType& ConditionalDofList,
                              ProcessInfo& rCurrentProcessInfo )
{
    ConditionalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
    {
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
    }
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
void LinePressureIsogeometric2D::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag,
                                bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = 2;
    unsigned int MatSize = number_of_nodes * dim;
    //resizing as needed the LHS

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
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
        integration_points
    );
    
    #ifdef DEBUG_LEVEL1
    KRATOS_WATCH(mThisIntegrationMethod)
    KRATOS_WATCH(integration_points.size())
    KRATOS_WATCH(Ncontainer)
    #endif

    const double& P = this->GetValue(PRESSURE);
//    KRATOS_WATCH(P)

    // loop over integration points
    Vector Load( dim );
    Vector t( dim );
    double dL;
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        // compute integration weight
        double IntegrationWeight = integration_points[PointNumber].Weight();

        IntegrationWeight *= GetProperties()[THICKNESS];

        // compute length
        noalias( t ) = ZeroVector( dim );//tangential vector
        for ( unsigned int n = 0; n < GetGeometry().size(); ++n )
        {
            t[0] += GetGeometry().GetPoint( n ).X0() * DN_De[PointNumber]( n, 0 );
            t[1] += GetGeometry().GetPoint( n ).Y0() * DN_De[PointNumber]( n, 0 );
        }
        dL = norm_2(t);
//        KRATOS_WATCH(t)

        //calculating load 
        Load[0] = -P*t[1]/dL;
        Load[1] = P*t[0]/dL;

        // contribute to RIGHT HAND SIDE VECTOR
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
                for ( unsigned int i = 0; i < dim; ++i )
                    rRightHandSideVector( prim * dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dL;
        }
    }

    #ifdef DEBUG_LEVEL1
    KRATOS_WATCH(rRightHandSideVector)
    #endif

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.

#undef DEBUG_LEVEL1

