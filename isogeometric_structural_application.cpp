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
#include "isogeometric_application/custom_geometries/geo_1d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier_3.h"
#include "isogeometric_application/custom_geometries/geo_3d_bezier.h"

namespace Kratos
{
    
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
    , mLineLoadBezier( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mLineLoadBezier2D( 0, Condition::GeometryType::Pointer( new Geo1dBezier<Node<3> >() ) )
    , mFaceLoadNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFaceLoadBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mFacePressureNURBS( 0, Condition::GeometryType::Pointer( new Geo2dNURBS<Node<3> >() ) )
    , mFacePressureBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mMasterContactFace3DBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mSlaveContactFace3DBezier( 0, Condition::GeometryType::Pointer( new Geo2dBezier3<Node<3> >() ) )
    , mKinematicLinearBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mKinematicLinearBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mTotalLagrangianBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<Node<3> >() ) )
    , mTotalLagrangianBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    , mUnsaturatedSoilsElement_2phase_SmallStrainBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<Node<3> >() ) )
    {}

    void KratosIsogeometricStructuralApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosIsogeometricStructuralApplication... " << std::endl;

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

        // register conditions
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS3D", mLineLoadNURBS )
        KRATOS_REGISTER_CONDITION( "LineLoadNURBS2D", mLineLoadNURBS2D )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier3D", mLineLoadBezier )
        KRATOS_REGISTER_CONDITION( "LineLoadBezier2D", mLineLoadBezier2D )
        KRATOS_REGISTER_CONDITION( "FaceLoadNURBS", mFaceLoadNURBS )
        KRATOS_REGISTER_CONDITION( "FaceLoadBezier", mFaceLoadBezier )
        KRATOS_REGISTER_CONDITION( "FacePressureNURBS", mFacePressureNURBS )
        KRATOS_REGISTER_CONDITION( "FacePressureBezier", mFacePressureBezier )
        KRATOS_REGISTER_CONDITION( "MasterContactFace3DBezier", mMasterContactFace3DBezier )
        KRATOS_REGISTER_CONDITION( "SlaveContactFace3DBezier", mSlaveContactFace3DBezier )
    }

} // namespace Kratos

