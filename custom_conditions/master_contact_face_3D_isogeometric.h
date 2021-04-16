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

#if !defined(KRATOS_MASTER_CONTACT_FACE_3D_ISOGEOMETRIC_H_INCLUDED )
#define  KRATOS_MASTER_CONTACT_FACE_3D_ISOGEOMETRIC_H_INCLUDED



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


/**
 * Contact surface element for 3D contact problems.
 * Defines a facet of a 3D-Element as a contact surface for
 * master contact surfaces
 * adapted from face2D originally written by riccardo.
 */
class MasterContactFace3DIsogeometric : public Condition
{
    typedef Condition BaseType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef PointerVectorSet< EquationIdVectorType, IndexedObject>
    EquationIdVectorContainerType;
    typedef BaseType::MatrixType LHS_ContributionType;
    typedef PointerVectorSet< LHS_ContributionType, IndexedObject> LHS_ContainerType;
    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

public:
    // Counted pointer of MasterContactFace3DIsogeometric
    KRATOS_CLASS_POINTER_DEFINITION( MasterContactFace3DIsogeometric );

    /**
     * Default constructor.
     */
    MasterContactFace3DIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry );
    MasterContactFace3DIsogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties );
    /**
     * Destructor.
     */
    virtual ~MasterContactFace3DIsogeometric() {}

    /**
     * Operations.
     */
    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties ) const;

    Condition::Pointer Create( IndexType NewId,
                               GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties ) const;

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

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
    bool ClosestPoint(
        GeometryType::CoordinatesArrayType& rResultGlobal,
        GeometryType::CoordinatesArrayType& rResultLocal,
        const GeometryType::CoordinatesArrayType& rCandidateGlobal,
        const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint
    );

    /**
     * applies the contact stress from a dedicated slave condition.
     * @param coords the coordinates of the slave condition's partner point on the
     * current master condition in local coordinates
     * @param Stress the value of the current contact stress in current point
     * @param Weight the integration weight in slave element's integration point
     * @param dA the differential area in slave element's integration point
     * @param NormalDirection the normal vector on current master surface
     */
    void AddContactStress( Vector coords, const double Stress, double Weight, double dA,
                           const Vector NormalDirection, double Penalty, double Gap,
                           EquationIdVectorType SlaveEquationId, Vector SlaveNcontainer );

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo );

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo ) const;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo ) const;

    void GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo );

    /**
    * This function provides the place to perform checks on the completeness of the input.
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    /**
     * Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const;

    /**
     * Print object's data.
     */
    virtual void PrintData(std::ostream& rOStream) const;

protected:

private:

    IntegrationMethod mThisIntegrationMethod;
    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    void CalculateAndAddKc( Matrix& K,
                            const Vector& N,
                            double pressure,
                            double weight,
                            double penalty,
                            Vector v
                          );

    void CalculateAndAdd_PressureForce( Vector& residualvector,
                                        const Vector& N,
                                        Vector& v3,
                                        double pressure,
                                        double weight,
                                        double DetJ );

    /**
     * Assignment operator.
     * (DEACTIVATED)
     */
    //MasterContactFace3DIsogeometric& operator=(const MasterContactFace3DIsogeometric& rOther);


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    MasterContactFace3DIsogeometric() {};

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

    /**
     * Copy constructor.
     * (DEACTIVATED)
     */
    //MasterContactFace3DIsogeometric(const MasterContactFace3DIsogeometric& rOther);

}; // Class MasterContactFace3DIsogeometric

}  // namespace Kratos.

#endif // KRATOS_MASTER_CONTACT_FACE_3D_ISOGEOMETRIC_H_INCLUDED  defined
