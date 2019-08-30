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
//   Date:                $Date: 2013 Dec 15 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_KINEMATIC_LINEAR_NURBS_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_NURBS_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"
#include "structural_application/custom_elements/kinematic_linear.h"

namespace Kratos
{

///@name Kratos Globals
///@{

extern Variable<Vector> PRESTRESS;
extern Variable<Vector> PLASTIC_STRAIN_VECTOR;
extern Variable<Vector> STRESSES;
extern Variable<double> PRESTRESS_FACTOR;
extern Variable<double> OVERCONSOLIDATION_RATIO;
extern Variable<int> PARENT_ELEMENT_ID;
extern Variable<int> INTEGRATION_POINT_INDEX;
extern Variable<int> NEIGHBOUR_EXPANSION_LEVEL;
extern Variable<Vector> RECOVERY_STRESSES;
extern Variable<int> STRESS_RECOVERY_TYPE;

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 KinematicLinearNURBS is designed to be a general linear strutural element support for both 2D and 3D.

 */

class KinematicLinearNURBS : public KinematicLinear
{

public:
    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    KRATOS_CLASS_POINTER_DEFINITION( KinematicLinearNURBS );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KinematicLinearNURBS( IndexType NewId, GeometryType::Pointer pGeometry );
    KinematicLinearNURBS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~KinematicLinearNURBS();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const;

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes ) const;

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void InitializeJacobian1D();

    void InitializeJacobian();

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
      virtual std::string Info() const
      {
        return "KinematicLinearNURBS";
      }

    /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
        rOStream << "KinematicLinearNURBS";
      }

    /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
        rOStream << GetGeometry();
      }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    KinematicLinearNURBS()
    {}

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //KinematicLinearNURBS& operator=(const KinematicLinearNURBS& rOther);

    /// Copy constructor.
    //KinematicLinearNURBS(const KinematicLinearNURBS& rOther);

    ///@}

}; // Class KinematicLinearNURBS

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
 KinematicLinearNURBS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
 const KinematicLinearNURBS& rThis)
 {
 rThis.PrintInfo(rOStream);
 rOStream << std::endl;
 rThis.PrintData(rOStream);

 return rOStream;
 }*/
///@}
}
  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR2_INCLUDED defined

