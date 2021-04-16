//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Apr 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_KINEMATIC_LINEAR_KIRCHHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_KIRCHHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"

namespace Kratos
{


/**
Implementation of the thin shell using Kirchhoff-Love theory and isogeometric analysis
REF: Nguyen Vinh Phu, igafem
 */
class KinematicLinearKirchhoffLoveIsogeometricShell : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    // Counted pointer of KinematicLinearKirchhoffLoveIsogeometricShell
    KRATOS_CLASS_POINTER_DEFINITION(KinematicLinearKirchhoffLoveIsogeometricShell);

    /**
     * Default constructor.
     */
    KinematicLinearKirchhoffLoveIsogeometricShell( IndexType NewId, GeometryType::Pointer pGeometry );
    KinematicLinearKirchhoffLoveIsogeometricShell( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /**
     * Destructor.
     */
    virtual ~KinematicLinearKirchhoffLoveIsogeometricShell();

    /**
     * Operations.
     */

    virtual Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void InitializeJacobian();

    void ResetConstitutiveLaw();

    void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo);

    void CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo );

    void CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo );

    IntegrationMethod GetIntegrationMethod() const;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList( DofsVectorType& ElementalDofList,
                     const ProcessInfo& CurrentProcessInfo) const;

    void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo );

    void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& Output, const ProcessInfo& rCurrentProcessInfo );

    void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo );

    /**
     * Turn back information as a string.
     * (DEACTIVATED)
     */
    //std::string Info();

    /**
     * Print information about this object.
     * (DEACTIVATED)
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;

    /**
     * Print object's data.
     * (DEACTIVATED)
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:

    friend class Serializer;

    // A private default constructor necessary for serialization
    KinematicLinearKirchhoffLoveIsogeometricShell()
    {
    }


    virtual void save ( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
    }

    virtual void load ( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
    }

    //virtual int Check( const ProcessInfo& rCurrentProcessInfo );


private:

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    IntegrationMethod mThisIntegrationMethod;

    double mTotalDomainInitialSize;

    bool mIsInitialized;

    Matrix mInitialDisp;

    void InitializeMaterial();

    ////////////////////////////////////////////////////////////////////////////
	////////////////////// private subroutines /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag = true,
                       bool CalculateResidualVectorFlag = true,
                       bool MaterialUpdateFlag = true );

    void CalculateStrain(Vector& StrainVector, const Matrix& B, const Matrix& Displacements);

    void CalculateAndAddExtForceContribution(VectorType& rRightHandSideVector, const Vector& N,
					     const double& weight, const double& detJ);

    void AddBodyForcesToRHS(VectorType& rRightHandSideVector, const Vector& N, const double& Weight, const double& DetJ);

    void AddInternalForcesToRHS(Vector& rRightHandSideVector, const Matrix& B,
				const Vector& StressVector, const double& Weight, const double& detJ);

    void CalculateMembraneBOperator(Matrix& Bm, const array_1d<double,3>& G1,const array_1d<double,3>& G2, const Matrix& DN_De);

    void CalculateBendingBOperator (Matrix& Bb, const array_1d<double,3>& G1, const array_1d<double,3>& G2, const Matrix& DN_De,
														 const ShapeFunctionsSecondDerivativesType& D2N_De2 ) ;

    void CalculateHookeanMatrix(Matrix& D, const Vector& G1, const Vector& G2, const double& E, const double& NU);

}; // Class KinematicLinearKirchhoffLoveIsogeometricShell

}  // namespace Kratos.


#endif // KRATOS_KINEMATIC_LINEAR_Kirchhoff_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined
