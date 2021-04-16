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
#include "custom_utilities/isotropic_tensor_utility.h"


namespace Kratos
{

class KirchhoffLoveLargeDeformationShell : public Element
{
public:

    template<std::size_t TDimension>
    class MyIsotropicTensorFunction_plus : public SingleIsotropicTensorFunction<TDimension>
    {
        virtual double Value(const double& x, const double& y)
        {
           if(x >= 0)
               return x;
           else
               return 0.0;
        }

        virtual double Value(const double& x, const double& y, const double& z)
        {
           if(x >= 0)
               return x;
           else
               return 0.0;
        }

        virtual void Derivative(const double& x, const double& y, Vector& Y)
        {
           if(x >= 0)
               Y[0] = 1.0;
           else
               Y[0] = 0.0;
           Y[1] = 0.0;
        }

        virtual void Derivative(const double& x, const double& y, const double& z, Vector& Y)
        {
           if(x >= 0)
               Y[0] = 1.0;
           else
               Y[0] = 0.0;
           Y[1] = 0.0;
           Y[2] = 0.0;
        }

    };

    template<std::size_t TDimension>
    class MyIsotropicTensorFunction_minus : public SingleIsotropicTensorFunction<TDimension>
    {
        virtual double Value(const double& x, const double& y)
        {
           if(x < 0)
               return x;
           else
               return 0.0;
        }

        virtual double Value(const double& x, const double& y, const double& z)
        {
           if(x < 0)
               return x;
           else
               return 0.0;
        }

        virtual void Derivative(const double& x, const double& y, Vector& Y)
        {
           if(x < 0)
               Y[0] = 1.0;
           else
               Y[0] = 0.0;
           Y[1] = 0.0;
        }

        virtual void Derivative(const double& x, const double& y, const double& z, Vector& Y)
        {
           if(x < 0)
               Y[0] = 1.0;
           else
               Y[0] = 0.0;
           Y[1] = 0.0;
           Y[2] = 0.0;
        }
    };

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
                               const ProcessInfo& rCurrentProcessInfo) const;

    void GetDofList( DofsVectorType& ElementalDofList,
                         const ProcessInfo& CurrentProcessInfo) const;

    ////////////////////////////////////////////////////////////
    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

    ////////////////////////////////////////////////////////////
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo);

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

    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


private:

    // flags
    bool mIsInitialized;

    bool mIsIsotropicMaterial;

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    IntegrationMethod mThisIntegrationMethod;

    boost::numeric::ublas::vector<boost::numeric::ublas::vector<ConstitutiveLaw::Pointer>> mConstitutiveLawVector;

    int mDim ;
    int mNumberOfNodes;
    int mStrainSize;
    int mNumberOfDof;
    double mThickness;
    double mTotalDomainInitialSize;

    // material parameters
    double mE, mNU;


    std::vector<Vector> mNodalCoordinates ;

    // geometric parameters
    unsigned int mNumberOfIntegrationPoint;
    std::vector<Matrix> mInvJ0;
    std::vector<Vector> mN;
    std::vector<Matrix> mDN_De;

    Vector mDetJ0;
    GeometryType::JacobiansType mJ0;
    Vector mIntegrationWeight;

    Vector mIntegrationPoint1D;
    Vector mWeight1D;
    double mDetJ1D;


    ////////////////////////////////////////////////////////////////////////////
	////////////////////// private subroutines /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag = true,
                       bool CalculateResidualVectorFlag = true,
                       bool MaterialUpdateFlag = true );

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    void AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
        std::vector<std::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight);

    void AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight );

    void AddLinearStiffnessMatrix(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight);

    void AddNonlinearStiffnessMatrix( MatrixType& LeftHandSideMatrix,const Vector& StressResultants
            , const std::vector<std::vector<Matrix>>& StrainVector_rs, const double& DetJ, const double& Weight);
    ///////////////////// add left hand side contribution
    void computeMembraneStrain(Matrix& eTensor,  Matrix& Aab, Matrix& aab  );

    void computeCurvatureChange(Matrix& kTensor, Matrix& Bab, Matrix& bab);

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////
    void DeformedCovariantBaseVector(std::vector<Vector>& a, std::vector<Vector>& A
                    , std::vector<Vector>& u_a );

    void ReferenceCovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
                                                                            , const std::vector<Vector>& X);

    void FirstDerivativeDisplacement_a(std::vector<Vector>& u_a, const Matrix& DN_De , const Matrix& u);

    void SecondDerivativeDisplacement_ab(std::vector<std::vector<Vector> >& u_ab
        , const ShapeFunctionsSecondDerivativesType& D2N_De2 , const Matrix& u);

    void NormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3,  std::vector<Vector>& a);

    void ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab);

    void DerivativeReferenceCovariantBaseVector(std::vector<std::vector<Vector> >& A_ab,
        const ShapeFunctionsSecondDerivativesType& D2N_De2, const std::vector<Vector>& X);

    void DerivativeDeformedCovariantBaseVector(std::vector<std::vector<Vector> >& a_ab
        , std::vector<std::vector<Vector> >& A_ab
        , std::vector<std::vector<Vector> >& u_ab);


    void DerivativeDeformedNormalDirector_a(std::vector<Vector>& a3Vector_a, std::vector<Vector>& a
            , Vector& aa3Vector, double& a3, std::vector<std::vector<Vector> >& a_ab);

    void DerivativeReferenceNormalDirector_a(std::vector<Vector>& A3Vector_a, std::vector<Vector>& A
                , Vector& AA3Vector, double& A3, std::vector<std::vector<Vector> >& A_ab);

    void ContinuumCovariantBaseVector(std::vector<Vector>& g, std::vector<Vector>& a
                    , std::vector<Vector>& a3Vector_a, Vector& a3Vector, double& theta3);
    //////////////////////////////////////////////////////////
    ////////////////////// material matrix ///////////////////
    //////////////////////////////////////////////////////////
    void CalculateElasticMatrix(Matrix& C, const double& E, const double& NU);

    /////////////////////////////////////////////////////////////////////////
    ///// first derivative of membrane strain and necessary components //////
    /////////////////////////////////////////////////////////////////////////
    void DerivativeCovariantBaseVector_r(std::vector< std::vector<std::vector<Vector>> >& a_ar, const Matrix& DN_De, std::vector<Vector>& UnitBasisVector);

    void DerivativeDisplacement_r(std::vector<std::vector<Vector>>& u_r, const Vector& N, std::vector<Vector>& UnitBasisVector);

    void DerivativeCovariantBaseVector_abr( std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
            ,const ShapeFunctionsSecondDerivativesType& D2N_De2, std::vector<Vector>& UnitBasisVector);

    void FirstDerivativeLocalMembraneStrain_r(std::vector<std::vector<Vector> >& eLocalVector_r
        , std::vector< std::vector<std::vector<Vector>> >& a_ar,  std::vector<Vector>& a, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);

    void SecondDerivativeLocalMembraneStrain_rs(std::vector<std::vector<Matrix>>& eLocalVector_rs,
            std::vector< std::vector<std::vector<Vector>> >& a_ar, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);
    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void DerivativeDirector_r( std::vector<std::vector<Vector>>& a3_rVector
        ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r, Vector& a3Vector, double& a3);

    void DerivativeDirectorNorm_r(std::vector<std::vector<double>>& a3_r
            , std::vector<std::vector<Vector>>& aa3_rVector, Vector& a3Vector);

    void DerivativeNonNormalizedDirector_r( std::vector<std::vector<Vector>>& aa3_rVector
                , std::vector< std::vector<std::vector<Vector>> >& a_ar, std::vector<Vector>& a);

    void SecondDerivativeDirectorNorm_rs(
                    std::vector<std::vector<double>>& a3_rs
                  , std::vector<std::vector<Vector>>& aa3_rsVector
                  , Vector& aa3Vector, Vector& a3Vector, double& a3,  std::vector<std::vector<Vector>>& aa3_rVector);

    void SecondDerivativeNonNormalizedDirector_rs(
                    std::vector< std::vector<Vector> >& aa3_rsVector
                 , std::vector< std::vector<std::vector<Vector>> >& a_ar);

    void SecondDerivativeDirector_rs(
                std::vector<std::vector<Vector>>& a3_rsVector
              , std::vector<std::vector<Vector>>& aa3_rsVector
              , std::vector<std::vector<double>>& a3_rs, Vector& a3Vector, double& a3
              ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r);

    void FirstDerivativeLocalCurvatureChange_r(std::vector<std::vector<Vector> >& kLocalVector_r
        , std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
        , std::vector<std::vector<Vector> >& a_ab
        , Vector& a3Vector , std::vector<std::vector<Vector>>& a3_rVector
        , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);

    void SecondDerivativeLocalCurvatureChange_rs(std::vector<std::vector<Matrix> >& kLocalVector_rs
            ,std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr, std::vector<std::vector<Vector>>& a3_rVector
            , std::vector<std::vector<Vector> >& a_ab, std::vector<std::vector<Vector>>& a3_rsVector
            , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);
    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A);

    void CovariantCurvatureCoefficient(Matrix& Bab, std::vector<std::vector<Vector> >& A_ab, Vector& A3Vector);

    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////
    void CreatingBmatrix(Matrix& BMatrix
        , const std::vector<std::vector<Vector> >& LocalStrainVector_r);

    double KroneckerDelta(int i, int j);

    void LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);

    void UnitBaseVectors(std::vector<Vector>& e);

    void LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A, Vector& A3Vector);

    void LocalTransformationCoefficient(std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff, std::vector<Vector>& EE, std::vector<Vector>& AA);




    //int Check( const ProcessInfo& rCurrentProcessInfo );

}; // Class KinematicLinearKirchoffLoveIsogeometricShell

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_KIRCHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined
