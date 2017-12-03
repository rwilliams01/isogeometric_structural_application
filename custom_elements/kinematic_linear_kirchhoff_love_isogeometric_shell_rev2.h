//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: DongGiang $
//   Date:                $Date: 25 August 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_KINEMATIC_LINEAR_KIRCHHOFF_LOVE_ISOGEOMETRIC_SHELL_REV2_H_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_KIRCHHOFF_LOVE_ISOGEOMETRIC_SHELL_REV2_H_INCLUDED


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

class KinematicLinearKirchhoffLoveIsogeometricShellRev2: public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    // Counted pointer of KinematicLinearKirchhoffLoveIsogeometricShell
    KRATOS_CLASS_POINTER_DEFINITION(KinematicLinearKirchhoffLoveIsogeometricShellRev2);

    /** 
     * Default constructor.
     */
    KinematicLinearKirchhoffLoveIsogeometricShellRev2 ( IndexType NewId, GeometryType::Pointer pGeometry );
    KinematicLinearKirchhoffLoveIsogeometricShellRev2( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /**
     * Destructor.
     */
    virtual ~KinematicLinearKirchhoffLoveIsogeometricShellRev2 ();

    /**
     * Operations.
     */

    virtual Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;


    /////////////////////////////////////////////////////////////////////////////////
    /////////////////// getting method //////////////////////////////////////////////

    IntegrationMethod GetIntegrationMethod() const;
    
    void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );
    
    void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

    /////////////////////////////////////////////////////////////////////////////////
    ///////////////////// starting-ending method ////////////////////////////////////
    void Initialize();
    
    void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////// computing method //////////////////////////////////////
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                               VectorType& rRightHandSideVector, 
                               ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                                 ProcessInfo& rCurrentProcessInfo);

protected:

    friend class Serializer;

    // A private default constructor necessary for serialization
    KinematicLinearKirchhoffLoveIsogeometricShellRev2()
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

    double mE, mNU;

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
                       bool CalculateStiffnessMatrixFlag=true ,
                       bool CalculateResidualVectorFlag =true,
                       bool MaterialUpdateFlag=true );

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void AddLinearMembraneStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC, const Matrix& Bm, const double& DetJ, const double& Weight);
    void AddLinearBendingStiffness(MatrixType& LeftHandSideMatrix, const Matrix& TanC, const Matrix& Bb,  const double& DetJ, const double& Weight);

    ////////////////////// add right hand side contribution
    void AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants, const Matrix& BMatrix, const double& DetJ, const double& Weight);
    void AddExternalForces(VectorType& RightHandSideVector, const Vector& N, const double& DetJ, const double& Weight );
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////// stress resultants
    void computeNormalForces(Vector& nVector, const Matrix& C, const Vector& eVector);

    void computeBendingMoments(Vector& moVector, const Matrix& C, const Vector& kVector); 

    void computeStrain(Vector& StrainVector,  const Matrix& B,  const Matrix& Displacements);

    void computeMembraneBMatrix(Matrix& Bm, const Matrix& DN_De, const std::vector<Vector>& A);

    void computeBendingBMatrix(Matrix& Bb,  Vector& A3Vector, double& A3,  std::vector<Vector>& A, 
        const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2,
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab);



    void CalculateBendingBOperator(
        Matrix& Bb,
        std::vector<Vector>& G,
        const Matrix& DN_De,
        const ShapeFunctionsSecondDerivativesType& D2N_De2 );
    //////////////////// material stiffness
    void computeTangentMaterialStiffness(Matrix& TanC, std::vector<Vector>& A);

    //////////////////////////////////////////////////////////////////
    ////////////// shell analysis utilities /////////////////////////
        ///////////////////// covariant base vectors                  
    void CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De, const std::vector<Vector>& X);

    void ReferencedNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A);

    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
            , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X);

    void CalculateElasticMatrix(Matrix& C,const double& E, const double& NU);


    /////////////////////////////////////////////////////////////
    void LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A,const Vector& A3Vector);

    void TransformationTensor(Matrix& T, std::vector<Vector>& EE, std::vector<Vector> A);


    void CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A);
    
    void ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab);



    void UnitBaseVectors(std::vector<Vector>& e);

}; // Class KinematicLinearKirchhoffLoveIsogeometricShell 

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_Kirchhoff_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined 
