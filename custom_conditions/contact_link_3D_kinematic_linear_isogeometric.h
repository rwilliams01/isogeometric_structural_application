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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 30 Jan 2015 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_CONTACT_LINK_3D_KINEMATIC_LINEAR_ISOGEOMETRIC_H_INCLUDED )
#define  KRATOS_CONTACT_LINK_3D_KINEMATIC_LINEAR_ISOGEOMETRIC_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "master_contact_face_3D_isogeometric.h"
#include "slave_contact_face_3D_isogeometric.h"


namespace Kratos
{
/**
 * Contact link element for 3D contact problems.
 * This condition links two contact surfaces (one master and
 * one slave element) to a condition which may be assembled
 * at once by the builder.
 * Within this condition all system contribution with regard to
 * the gap function, the contact stress, the normal vectors
 * and their linearizations are calculated.
 */
class ContactLink3D_Kinematic_Linear_Isogeometric : public Condition
{
public:
    // Counted pointer of ContactLink3D_Kinematic_Linear_Isogeometric
    KRATOS_CLASS_POINTER_DEFINITION(ContactLink3D_Kinematic_Linear_Isogeometric);

    /**
     * Default constructor.
     */
    ContactLink3D_Kinematic_Linear_Isogeometric( IndexType NewId, GeometryType::Pointer pGeometry);

    ContactLink3D_Kinematic_Linear_Isogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                                                 PropertiesType::Pointer pProperties );


    ContactLink3D_Kinematic_Linear_Isogeometric( IndexType NewId, GeometryType::Pointer pGeometry,
                                                 PropertiesType::Pointer pProperties,
                                                 Condition::Pointer Master,
                                                 Condition::Pointer Slave,
                                                 Point<3>& MasterContactLocalPoint,
                                                 Point<3>& SlaveContactLocalPoint,
                                                 int SlaveIntegrationPointIndex );
    /**
     * Destructor.
     */
    virtual ~ContactLink3D_Kinematic_Linear_Isogeometric();

    /**
     * Operations.
     */

    Condition::Pointer Create( IndexType NewId,
                               NodesArrayType const& ThisNodes,
                               PropertiesType::Pointer pProperties ) const;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo );

    void DampMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix(
            MatrixType& rDampMatrix,
            const ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo )const;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo ) const;

    /**
     * Turn back information as a string.
     * (DEACTIVATED)
     */
    //std::string Info();

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
    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

    void CalculateAndAdd_RHS( Vector& residualvector,
                              const Vector& NMaster,
                              const Vector& NSlave,
                              const Vector& vMaster,
                              const Matrix& T,
                              const Vector& tangentialStresses,
                              double Gap,
                              double normalStress,
                              double SlaveIntegrationWeight,
                              double dASlave);

    /**
     * This function calculates updates the local and global coordinates
     * of the master contact partner in order to follow the movement of
     * the slave surface along the master surface
     */
    void UpdateMasterLocalPoint( );


    Vector NormalVector( Condition::Pointer Surface,
                         const GeometryType::CoordinatesArrayType& LocalPoint );

    Matrix TangentialVectors( Condition::Pointer Surface,
                              const GeometryType::CoordinatesArrayType& LocalPoint );

    Matrix TangentialVectorsGlobal( Condition::Pointer Surface,
                                    const GeometryType::CoordinatesArrayType& LocalPoint );

    Matrix TangentialVectors_inOrigin( Condition::Pointer Surface,
                                       const GeometryType::CoordinatesArrayType& rPoint );

    Point<3>& GlobalCoordinates(Condition::Pointer Surface, Point<3>& rResult, Point<3> const& LocalCoordinates);


    Vector GetRelativTangentialVelocity(Matrix& T);

    Vector GetRelativVelocity();
    /**
     * Assignment operator.
     * (DEACTIVATED)
     */
    //ContactLink3D_Kinematic_Linear_Isogeometric& operator=(const ContactLink3D_Kinematic_Linear_Isogeometric& rOther);

    /**
     * Copy constructor.
     * (DEACTIVATED)
     */
    //ContactLink3D_Kinematic_Linear_Isogeometric(const ContactLink3D_Kinematic_Linear_Isogeometric& rOther);


    /**
     * private members
     */

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    ContactLink3D_Kinematic_Linear_Isogeometric() {}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    Vector mvMaster;
    Matrix mTMaster;
}; // Class ContactLink3D_Kinematic_Linear_Isogeometric

}  // namespace Kratos.

#endif // KRATOS_CONTACT_LINK_3D_KINEMATIC_LINEAR_ISOGEOMETRIC_H_INCLUDED  defined
