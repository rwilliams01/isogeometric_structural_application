//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: May 31, 2016 $
//   Revision:            $Revision: 1.0 $
//
// 


// System includes


// External includes


// Project includes
#include "isogeometric_structural_application.h"
#include "geometries/geometry.h"
#include "geometries/line_3d_2.h"
#include "isogeometric_application/custom_geometries/geo_1d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier_3.h"
#include "isogeometric_application/custom_geometries/geo_3d_bezier.h"

namespace Kratos
{

    KRATOS_CREATE_VARIABLE( Vector, REF_BASE_VECTOR_1 )
    KRATOS_CREATE_VARIABLE( Vector, REF_BASE_VECTOR_2 )
    KRATOS_CREATE_VARIABLE( Vector, REF_BASE_VECTOR_3 )
    KRATOS_CREATE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_1 )
    KRATOS_CREATE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_2 )
    KRATOS_CREATE_VARIABLE( Vector, REF_CONTRA_BASE_VECTOR_3 )
    KRATOS_CREATE_VARIABLE( Vector, CURRENT_BASE_VECTOR_1 )
    KRATOS_CREATE_VARIABLE( Vector, CURRENT_BASE_VECTOR_2 )
    KRATOS_CREATE_VARIABLE( Vector, CURRENT_BASE_VECTOR_3 )
    KRATOS_CREATE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_1 )
    KRATOS_CREATE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_2 )
    KRATOS_CREATE_VARIABLE( Vector, LOCAL_CARTESIAN_VECTOR_3 )

    KratosIsogeometricStructuralApplication::KratosIsogeometricStructuralApplication()
    :
      mKinematicLinearGeo1dNURBS( 0, Element::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo2dNURBS( 0, Element::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo3dNURBS( 0, Element::GeometryType::Pointer( new Geo3dNURBS<Node<3> >() ) )
    , mKinematicLinearGeo1dBezier( 0, Element::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mKinematicLinearGeo2dBezier( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mKinematicLinearGeo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mTotalLagrangianGeo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mLineLoadNURBS( 0, Condition::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mLineLoadNURBS2D( 0, Condition::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mLinePressureNURBS2D( 0, Condition::GeometryType::Pointer( new Geo1dNURBS<Node<3> >() ) )
    , mLineLoadBezier( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mLineLoadBezier2D( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mLinePressureBezier2D( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mFaceLoadNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFaceLoadBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mFacePressureNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFacePressureBezier2D3( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mMasterContactFaceBezier2D3( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mSlaveContactFaceBezier2D3( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mElasticFaceConstraintBezier2D3( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mElasticFaceSpringsBezier2D3( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKinematicLinearBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mKinematicLinearBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mTotalLagrangianBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mTotalLagrangianBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mKinematicLinearKirchhoffLoveIsogeometricShellBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKinematicLinearKirchhoffLoveIsogeometricShellRev2Bezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKirchhoffLoveLargeDeformationShellBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    #ifdef PLATE_AND_SHELL_APPLICATION_IS_ON
    , mKirchhoffLoveLinearShellBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    #ifdef PLATE_AND_SHELL_APPLICATION_USE_ADOL_C
    , mKirchhoffLoveShellADBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKirchhoffLoveShell2ADBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    #endif
    #endif
    , mPenaltyStiffnessShell3D2N( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mConvectiveLinearKirchhoffLoveShellBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mLinearBendingStripBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mNonLinearBendingStripBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    {}

    void KratosIsogeometricStructuralApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosIsogeometricStructuralApplication... " << std::endl;

        // register variables
        KRATOS_REGISTER_VARIABLE( REF_BASE_VECTOR_1 )
        KRATOS_REGISTER_VARIABLE( REF_BASE_VECTOR_2 )
        KRATOS_REGISTER_VARIABLE( REF_BASE_VECTOR_3 )
        KRATOS_REGISTER_VARIABLE( REF_CONTRA_BASE_VECTOR_1 )
        KRATOS_REGISTER_VARIABLE( REF_CONTRA_BASE_VECTOR_2 )
        KRATOS_REGISTER_VARIABLE( REF_CONTRA_BASE_VECTOR_3 )
        KRATOS_REGISTER_VARIABLE( CURRENT_BASE_VECTOR_1 )
        KRATOS_REGISTER_VARIABLE( CURRENT_BASE_VECTOR_2 )
        KRATOS_REGISTER_VARIABLE( CURRENT_BASE_VECTOR_3 )
        KRATOS_REGISTER_VARIABLE( LOCAL_CARTESIAN_VECTOR_1 )
        KRATOS_REGISTER_VARIABLE( LOCAL_CARTESIAN_VECTOR_2 )
        KRATOS_REGISTER_VARIABLE( LOCAL_CARTESIAN_VECTOR_3 )

        // register elements
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo1dNURBS", mKinematicLinearGeo1dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo2dNURBS", mKinematicLinearGeo2dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo3dNURBS", mKinematicLinearGeo3dNURBS )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo1dBezier", mKinematicLinearGeo1dBezier )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo2dBezier", mKinematicLinearGeo2dBezier )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearGeo3dBezier", mKinematicLinearGeo3dBezier )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianGeo3dBezier", mTotalLagrangianGeo3dBezier )
        KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier", mUnsaturatedSoilsElement_2phase_SmallStrain_Geo3dBezier )

        KRATOS_REGISTER_ELEMENT( "KinematicLinearBezier2D", mKinematicLinearBezier2D )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearBezier3D", mKinematicLinearBezier3D )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianBezier2D", mTotalLagrangianBezier2D )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianBezier3D", mTotalLagrangianBezier3D )
        KRATOS_REGISTER_ELEMENT( "UnsaturatedSoilsElement_2phase_SmallStrainBezier3D", mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearKirchhoffLoveIsogeometricShellBezier2D3", mKinematicLinearKirchhoffLoveIsogeometricShellBezier2D3 )
        KRATOS_REGISTER_ELEMENT( "KinematicLinearKirchhoffLoveIsogeometricShellRev2Bezier2D3", mKinematicLinearKirchhoffLoveIsogeometricShellRev2Bezier2D3 )
        KRATOS_REGISTER_ELEMENT( "KirchhoffLoveLargeDeformationShellBezier2D3", mKirchhoffLoveLargeDeformationShellBezier2D3 )
        #ifdef PLATE_AND_SHELL_APPLICATION_IS_ON
        KRATOS_REGISTER_ELEMENT( "KirchhoffLoveLinearShellBezier2D3", mKirchhoffLoveLinearShellBezier2D3 )
        #ifdef PLATE_AND_SHELL_APPLICATION_USE_ADOL_C
        KRATOS_REGISTER_ELEMENT( "KirchhoffLoveShellADBezier2D3", mKirchhoffLoveShellADBezier2D3 )
        KRATOS_REGISTER_ELEMENT( "KirchhoffLoveShell2ADBezier2D3", mKirchhoffLoveShell2ADBezier2D3 )
        #endif
        #endif

        // register conditions
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS3D", mLineLoadNURBS )
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS2D", mLineLoadNURBS2D )
        KRATOS_REGISTER_CONDITION( "LinePressureNURBS2D", mLinePressureNURBS2D )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier3D", mLineLoadBezier )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier2D", mLineLoadBezier2D )
        KRATOS_REGISTER_CONDITION( "LinePressureBezier2D", mLinePressureBezier2D )
        KRATOS_REGISTER_CONDITION( "FaceLoadNURBS", mFaceLoadNURBS )
        KRATOS_REGISTER_CONDITION( "FaceLoadBezier", mFaceLoadBezier )
        KRATOS_REGISTER_CONDITION( "FacePressureNURBS", mFacePressureNURBS )
        KRATOS_REGISTER_CONDITION( "FacePressureBezier2D3", mFacePressureBezier2D3 )
        KRATOS_REGISTER_CONDITION( "MasterContactFaceBezier2D3", mMasterContactFaceBezier2D3 )
        KRATOS_REGISTER_CONDITION( "SlaveContactFaceBezier2D3", mSlaveContactFaceBezier2D3 )
        KRATOS_REGISTER_CONDITION( "PenaltyStiffnessShell3D2N", mPenaltyStiffnessShell3D2N )
        KRATOS_REGISTER_CONDITION( "ElasticFaceConstraintBezier2D3", mElasticFaceConstraintBezier2D3 )
        KRATOS_REGISTER_CONDITION( "ElasticFaceSpringsBezier2D3", mElasticFaceSpringsBezier2D3 )
    }

} // namespace Kratos

