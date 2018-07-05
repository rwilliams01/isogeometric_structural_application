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
#include "custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell.h"
#include "custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell_rev2.h"
#include "custom_elements/K_L_large_deformation_shell.h"
#include "custom_elements/convective_linear_KL_shell.h"
#include "custom_elements/linear_bending_strip.h"
#include "custom_elements/nonlinear_bending_strip.h"
#include "custom_conditions/line_force_isogeometric.h"
#include "custom_conditions/line_force_isogeometric_2d.h"
#include "custom_conditions/line_pressure_isogeometric_2d.h"
#include "custom_conditions/face_load_isogeometric.h"
#include "custom_conditions/face_pressure_isogeometric.h"
#include "custom_conditions/slave_contact_face_3D_isogeometric.h"
#include "custom_conditions/master_contact_face_3D_isogeometric.h"
#include "custom_conditions/penalty_stiffness_shell.h"
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
#include "structural_application/custom_conditions/elastic_constraint.h"
#ifdef PLATE_AND_SHELL_APPLICATION_IS_ON
#include "plate_and_shell_application/custom_elements/kirchhoff_love_linear_shell.h"
#ifdef PLATE_AND_SHELL_APPLICATION_USE_ADOL_C
#include "plate_and_shell_application/custom_elements/kirchhoff_love_shell_ad.h"
#include "plate_and_shell_application/custom_elements/kirchhoff_love_shell_2_ad.h"
#endif
#endif

namespace Kratos
{

KRATOS_DEFINE_VARIABLE( Vector, REF_BASE_VECTOR_1 )
KRATOS_DEFINE_VARIABLE( Vector, REF_BASE_VECTOR_2 )
KRATOS_DEFINE_VARIABLE( Vector, REF_BASE_VECTOR_3 )
KRATOS_DEFINE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_1 )
KRATOS_DEFINE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_2 )
KRATOS_DEFINE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_3 )
KRATOS_DEFINE_VARIABLE( Vector, CURRENT_BASE_VECTOR_1 )
KRATOS_DEFINE_VARIABLE( Vector, CURRENT_BASE_VECTOR_2 )
KRATOS_DEFINE_VARIABLE( Vector, CURRENT_BASE_VECTOR_3 )
KRATOS_DEFINE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_1 )
KRATOS_DEFINE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_2 )
KRATOS_DEFINE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_3 )
    
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

        // element with old interface
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
        const KinematicLinearKirchhoffLoveIsogeometricShell mKinematicLinearKirchhoffLoveIsogeometricShellBezier2D3;
        const KinematicLinearKirchhoffLoveIsogeometricShell mKinematicLinearKirchhoffLoveIsogeometricShellRev2Bezier2D3;

        const LineForceIsogeometric mLineLoadNURBS;
        const LineForceIsogeometric2D mLineLoadNURBS2D;
        const LinePressureIsogeometric2D mLinePressureNURBS2D;
        const LineForceIsogeometric mLineLoadBezier;
        const LineForceIsogeometric2D mLineLoadBezier2D;
        const LinePressureIsogeometric2D mLinePressureBezier2D;

        const FaceLoadIsogeometric mFaceLoadNURBS;
        const FaceLoadIsogeometric mFaceLoadBezier;
        const FacePressureIsogeometric mFacePressureNURBS;
        const FacePressureIsogeometric mFacePressureBezier2D3;

        const MasterContactFace3DIsogeometric mMasterContactFaceBezier2D3;
        const SlaveContactFace3DIsogeometric mSlaveContactFaceBezier2D3;

        // the elements below use new interface, which means it exploits the NURBS geometry in the kernel so the element is not needed to modify to work with NURBS
        const KinematicLinear mKinematicLinearBezier2D;
        const KinematicLinear mKinematicLinearBezier3D;
        const KinematicLinear mTotalLagrangianBezier2D;
        const KinematicLinear mTotalLagrangianBezier3D;
        const UnsaturatedSoilsElement_2phase_SmallStrain mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D;
        #ifdef PLATE_AND_SHELL_APPLICATION_IS_ON
        const KirchhoffLoveLinearShell mKirchhoffLoveLinearShellBezier2D3;
        #ifdef PLATE_AND_SHELL_APPLICATION_USE_ADOL_C
        const KirchhoffLoveShellAD mKirchhoffLoveShellADBezier2D3;
        const KirchhoffLoveShell2AD mKirchhoffLoveShell2ADBezier2D3;
        #endif
        #endif
        const KirchhoffLoveLargeDeformationShell mKirchhoffLoveLargeDeformationShellBezier2D3;

        const ElasticConstraint mElasticFaceConstraintBezier2D3;

        const PenaltyStiffnessShell mPenaltyStiffnessShell3D2N;

        const ConvectiveLinearKirchhoffLoveShell mConvectiveLinearKirchhoffLoveShellBezier2D3;
        const LinearBendingStrip mLinearBendingStripBezier2D3;
        const NonLinearBendingStrip mNonLinearBendingStripBezier2D3;

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

