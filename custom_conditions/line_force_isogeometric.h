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
*   Date:                $Date: 23 Oct 2014 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_ISOGEOMETRIC_APP_LINE_FORCE_ISOGEOMETRIC_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APP_LINE_FORCE_ISOGEOMETRIC_H_INCLUDED


// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"

namespace Kratos
{

extern Variable<array_1d<double, 3> > FACE_LOAD;

/*
 * Implement tht force on line condition in 3D. The force at the integration point is interpolated from the FACE_LOAD value at the control points.
 * Therefore this condition works in both 2D and 3D. In 2D, the thickness is not incorporated. If thickness must be accounted, use LineForceIsogeometric2D.
 */
class LineForceIsogeometric : public Condition
{
public:

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    // Counted pointer of LineForceIsogeometric
    KRATOS_CLASS_POINTER_DEFINITION( LineForceIsogeometric );

    // Constructor void
    LineForceIsogeometric();

    // Constructor using an array of nodes
    LineForceIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    LineForceIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~LineForceIsogeometric();

    // Name Operations

    virtual Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties ) const;

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo );

    virtual void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

protected:

    IntegrationMethod mThisIntegrationMethod;

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    virtual void CalculateAll( MatrixType& rLeftHandSideMatrix,
                   VectorType& rRightHandSideVector,
                   const ProcessInfo& rCurrentProcessInfo,
                   bool CalculateStiffnessMatrixFlag,
                   bool CalculateResidualVectorFlag );

private:
    ///@name Static Member Variables

    /// privat variables

    // privat name Operations


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class LineForceIsogeometric.

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_ISOGEOMETRIC_H_INCLUDED  defined
