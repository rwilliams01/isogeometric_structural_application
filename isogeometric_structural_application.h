//   
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: May 31, 2016 $
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   28/7/2015: create appname_application.h

#if !defined(KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_H_INCLUDED)
#define KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

#include "custom_elements/kinematic_linear_nurbs.h"
#include "custom_elements/kinematic_linear_isogeometric.h"
#include "custom_elements/total_lagrangian_isogeometric.h"
#include "custom_elements/unsaturated_soils_element_2phase_small_strain_isogeometric.h"
#include "custom_conditions/line_force_isogeometric.h"
#include "custom_conditions/line_force_isogeometric_2d.h"
#include "custom_conditions/face_load_isogeometric.h"
#include "custom_conditions/face_pressure_isogeometric.h"
#include "custom_conditions/slave_contact_face_3D_isogeometric.h"
#include "custom_conditions/master_contact_face_3D_isogeometric.h"
#include "isogeometric_application/custom_geometries/geo_1d_nurbs.h"
#include "isogeometric_application/custom_geometries/geo_2d_nurbs.h"
#include "isogeometric_application/custom_geometries/geo_3d_nurbs.h"
#include "isogeometric_application/custom_geometries/geo_1d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier_3.h"
#include "isogeometric_application/custom_geometries/geo_3d_bezier.h"
#include "structural_application/custom_elements/kinematic_linear.h"
#include "structural_application/custom_elements/total_lagrangian.h"
#include "structural_application/custom_elements/unsaturated_soils_element_2phase_small_strain.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KratosIsogeometricStructuralApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{
        
        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosIsogeometricStructuralApplication);

        ///@}
        ///@name Life Cycle
        ///@{ 

        /// Default constructor.
        KratosIsogeometricStructuralApplication();

        /// Destructor.
        virtual ~KratosIsogeometricStructuralApplication(){}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register();

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
            return "Application of Isogeometric Analysis in Structural Mechanics";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << "in KratosIsogeometricStructuralApplication:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
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


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


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


//        const KinematicLinearNURBS mKinematicLinearGeo1dNURBS; //TODO: doesn't work, segmentation fault error
//        const KinematicLinearNURBS mKinematicLinearGeo2dNURBS;
//        const KinematicLinearNURBS mKinematicLinearGeo3dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo1dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo2dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo3dNURBS;
        const KinematicLinearIsogeometric mKinematicLinearGeo1dBezier;
        const KinematicLinearIsogeometric mKinematicLinearGeo2dBezier;
        const KinematicLinearIsogeometric mKinematicLinearGeo3dBezier;
        const TotalLagrangianIsogeometric mTotalLagrangianGeo3dBezier;
        const UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier;

        const LineForceIsogeometric mLineLoadNURBS;
        const LineForceIsogeometric2D mLineLoadNURBS2D;
        const LineForceIsogeometric mLineLoadBezier;
        const LineForceIsogeometric2D mLineLoadBezier2D;

        const FaceLoadIsogeometric mFaceLoadNURBS;
        const FaceLoadIsogeometric mFaceLoadBezier;
        const FacePressureIsogeometric mFacePressureNURBS;
        const FacePressureIsogeometric mFacePressureBezier;

        const MasterContactFace3DIsogeometric mMasterContactFace3DBezier;
        const SlaveContactFace3DIsogeometric mSlaveContactFace3DBezier;

        const KinematicLinear mKinematicLinearBezier2D;
        const KinematicLinear mKinematicLinearBezier3D;
        const KinematicLinear mTotalLagrangianBezier2D;
        const KinematicLinear mTotalLagrangianBezier3D;
        const UnsaturatedSoilsElement_2phase_SmallStrain mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D;



        ///@}
        ///@name Private Operators
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
        KratosIsogeometricStructuralApplication& operator=(KratosIsogeometricStructuralApplication const& rOther);

        /// Copy constructor.
        KratosIsogeometricStructuralApplication(KratosIsogeometricStructuralApplication const& rOther);


        ///@}

    }; // Class KratosIsogeometricStructuralApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_H_INCLUDED defined

