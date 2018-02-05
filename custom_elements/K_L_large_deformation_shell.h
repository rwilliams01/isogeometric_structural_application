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
#include "phase_field_application/custom_utilities/isotropic_tensor_utility.h"


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

    bool mIsIsotropicMaterial;


    
    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    IntegrationMethod mThisIntegrationMethod;

    int mDim ;
    int mNumberOfNodes;
    int mStrainSize;
    int mNumberOfDof;
    double mThickness;

    // material parameters
    double mE, mNU , mLambda, mMu ,mKappa ;
    

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
    void AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight);

    void AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight );
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////// stress resultants
    void computeNormalForces(Vector& nVector, const Matrix& C, const Vector& eVector);

    
    void computeBendingMoments(Vector& moVector, const Matrix& C, const Vector& kVector); 


    //////////////////// Strain Vectors
    void computeMembraneStrain(Matrix& eTensor,  std::vector<Vector>& a,  std::vector<Vector>& A, std::vector<Vector>& EE, std::vector<Vector>& AA );

    //void computeMembraneStrain(Vector& eVector,  std::vector<Vector>& a,  std::vector<Vector>& A);

    void computeCurvatureChange(Matrix& kTensor,  boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
        , boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& a_ab, Vector& a3Vector,  Vector& A3Vector, std::vector<Vector>& EE, std::vector<Vector>& AA );
    /////////////////////  B matrix

    ///////////////////// first derivarives of strains
    void FirstDerivativeLocalMembraneStrain_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& eVector_r
            , const Matrix& DN_De,  std::vector<Vector>& a, std::vector<Vector>& EE, std::vector<Vector>& AA) ;


    void FirstDerivativeLocalCurvatureChange_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& kLocalVector_r
                , std::vector<Vector>& a
                , Vector& a3Vector, Vector& aa3Vector, double& a3 ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                , const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2, std::vector<Vector>& EE, std::vector<Vector>& AA);

    ///////////////////// second derivatives of membrane strains and curvature changes

    void SecondDerivativeLocalMembraneStrain_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& eLocalVector_rs,
                const Matrix& DN_De, std::vector<Vector>& EE, std::vector<Vector>& AA);

    void SecondDerivativeLocalCurvatureChange_rs(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix>>& kLocalVector_rs
                    ,const Matrix& DN_De, const ShapeFunctionsSecondDerivativesType& D2N_De2
                    ,boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                    ,std::vector<Vector>& a,Vector& a3Vector, Vector& aa3Vector, const double& a3
                    , std::vector<Vector>& EE, std::vector<Vector>& AA);

    ///////////////////// covariant base vectors                  
    void CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De
        , const std::vector<Vector>& X);

    void CovariantBaseVector(std::vector<Vector>& a, const Matrix& DN_De
        , const std::vector<Vector>& X, const Matrix& u );

    void DeformedNormalDirector(Vector& a3Vector, Vector& aa3Vector, double& a3, std::vector<Vector>& a);

    void ReferencedNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A);
 
    ///////////////////// derivative of covariant base vectors w.r.t curvilinear coordinates

    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& a_ab
                                    , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X, const Matrix& u);

    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
                                                        , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X);


    // material matrix
    void CalculateElasticMatrix(Matrix& C,const double& E, const double& NU);


    double KroneckerDelta(int i, int j);

    void UnitBaseVectors(std::vector<Vector>& e);


    //////////////////////////////////////////////////////////////////////////////
    ///////////// transformation operator ///////////////////////////////////////
    void LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A);
 
    void CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A);
    
    void ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab);

    void CreatingBmatrix(Matrix& BMatrix, const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& LocalStrainVector_r);

  
}; // Class KinematicLinearKirchoffLoveIsogeometricShell 

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_KIRCHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined 
