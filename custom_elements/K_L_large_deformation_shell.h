//
//   Project Name:        Kratos
//   Last Modified by:    $Author: DongGiang $
//   Date:                $Date: 25 August 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_KIRCHHOFF_LOVE_LARGE_DEFORMATION_SHELL_H_INCLUDED )
#define  KRATOS_KIRCHHOFF_LOVE_LARGE_DEFORMATION_SHELL_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"



// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"

namespace Kratos
{

class KirchhoffLoveLargeDeformationShell : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    // Counted pointer of KinematicLinearKirchoffLoveIsogeometricShell
    KRATOS_CLASS_POINTER_DEFINITION(KirchhoffLoveLargeDeformationShell);

    /**
     * Default constructor.
     */
    KirchhoffLoveLargeDeformationShell ( IndexType NewId, GeometryType::Pointer pGeometry );
    KirchhoffLoveLargeDeformationShell ( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /**
     * Destructor.
     */
    virtual ~KirchhoffLoveLargeDeformationShell ();

    /**
     * Operations.
     */

    virtual Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;

    ///////////////////////////////////////////////////////////
    IntegrationMethod GetIntegrationMethod() const;

    void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo);

    void GetDofList( DofsVectorType& ElementalDofList,
                         ProcessInfo& CurrentProcessInfo);

    ////////////////////////////////////////////////////////////
    void Initialize();

    void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );


    ////////////////////////////////////////////////////////////
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 ProcessInfo& rCurrentProcessInfo);

protected:

    friend class Serializer;

    // A private default constructor necessary for serialization
    KirchhoffLoveLargeDeformationShell ()
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

    // flags
    bool mIsInitialized;

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    IntegrationMethod mThisIntegrationMethod;

    int mDim ;
    int mNumberOfNodes;
    int mStrainSize;
    int mNumberOfDof;
    double mThickness;

    std::vector<Vector> mNodalCoordinates ;

    // geometric parameters
    unsigned int mNumberOfIntegrationPoint;
    std::vector<Matrix> mInvJ0;
    std::vector<Vector> mN;
    std::vector<Matrix> mDN_De;
    std::vector<Matrix> mDN_DX;
    Vector mDetJ0;
    GeometryType::JacobiansType mJ0;
    Vector mIntegrationWeight;


    ////////////////////////////////////////////////////////////////////////////
	////////////////////// private subroutines /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag = true,
                       bool CalculateResidualVectorFlag = true,
                       bool MaterialUpdateFlag = true );

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void AddLinearMembraneStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC, const Matrix& Bm, const double& DetJ, const double& Weight);

    void AddLinearBendingStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC, const Matrix& Bb,  const double& DetJ, const double& Weight);

    void AddNonlinearMembraneStiffness( MatrixType& LeftHandSideMatrix,const Vector& nVector
        , const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eVector_rs, const double& DetJ, const double& Weight);

    void AddNonlinearBendingStiffness( MatrixType& LeftHandSideMatrix,const Vector& moVector
        , const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kVector_rs, const double& DetJ, const double& Weight);

    ////////////////////// add right hand side contribution
    void AddExternalForces(VectorType& RightHandSideVector, const Vector& N, double DetJ, double Weight );

    void AddMembraneInternalForces(VectorType& RightHandSideVector, const Vector& nVector, const Matrix& Bm, double DetJ, double Weight);

    void AddBendingInternalForces(VectorType& RightHandSideVector, const Vector& moVector, const Matrix& Bb, double DetJ, double Weight);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////// stress resultants
    void computeNormalForces(Vector& nVector, const Matrix& C, const Vector& eVector);

    void computeBendingMoments(Vector& moVector, const Matrix& C, const Vector& kVector);

    //////////////////// material stiffness
    void computeTangentMaterialStiffness(Matrix& TanC, std::vector<Vector>& A);

    //////////////////// Strain Vectors
    void computeMembraneStrain(Vector& eVector,  std::vector<Vector>& a,  std::vector<Vector>& A);
    void computeStrain(Vector& StrainVector,  const Matrix& B,  const Matrix& Displacements);

    //void computeMembraneStrain(Vector& eVector,  std::vector<Vector>& a,  std::vector<Vector>& A);

    void computeCurvatureChange(Vector& kVector,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
        , boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& a_ab, Vector a3Vector,  Vector A3Vector );
    /////////////////////  B matrix
    void computeMembraneBMatrix(Matrix& Bm, const Matrix& DN_De, const std::vector<Vector>& a);

    void computeBendingBMatrix(Matrix& Bb,  std::vector<Vector>& a
        , Vector& a3Vector, Vector& aa3Vector, double a3 ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                 ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2);

    ///////////////////// second derivatives of membrane strains and curvature changes
    void SecondDerivativeMembraneStrain_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eVector_rs,
        const Matrix& DN_De);

    void SecondDerivativeCurvatureChange_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kVector_rs
        ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2
        ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> > a_ab
        ,std::vector<Vector>& a,Vector& a3Vector, Vector& aa3Vector, const double a3, std::vector<Vector>& e);

    ///////////////////// covariant base vectors
    void CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
        , const std::vector<Vector>& X);

    void CovariantBaseVector(std::vector<Vector>& a, const Matrix& DN_De
        , const std::vector<Vector> X, const Matrix& u );

    void DeformedNormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3, std::vector<Vector>& a);

    void ReferencedNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A);

    ///////////////////// derivative of covariant base vectors w.r.t curvilinear coordinates
    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                                    , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X, const Matrix& u);

    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
                                                        , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X);

    ///////////////////// derivative of covariant base vectors w.r.t nodal displacements
    void DerivativeCovariantBaseVector_u(std::vector<Matrix>& a_ar, const Matrix& DN_De,std::vector<Vector>& e);

    void DerivativeNumeratorNomalDirector(Matrix& aa3_r,  std::vector<Vector>& a,
        std::vector<Vector>& e, const Matrix DN_De);

    void DerivativeDenominatorNormalDirector(Matrix& a3_r,  Matrix& aa3_r,  Vector& a3Vector );


    void SecondDerivativeCovariantBaseVector_ur(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> >& a_abr
                                                , const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector> e);


    double KroneckerDelta(int i, int j);

    void UnitBaseVectors(std::vector<Vector>& e);


}; // Class KinematicLinearKirchoffLoveIsogeometricShell

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_KIRCHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined
