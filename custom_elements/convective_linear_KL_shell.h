//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: DongGiang $
//   Date:                $Date: 25 August 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_CONVECTIVE_LINEAR_KL_SHELL_H_INCLUDED )
#define  KRATOS_CONVECTIVE_LINEAR_KL_SHELL_H_INCLUDED


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
    // variable imported from structural_application
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)
class ConvectiveLinearKirchhoffLoveShell: public Element
{
public:
 

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    // Counted pointer of KinematicLinearKirchoffLoveIsogeometricShell
    KRATOS_CLASS_POINTER_DEFINITION(ConvectiveLinearKirchhoffLoveShell);

    /** 
     * Default constructor.
     */
    ConvectiveLinearKirchhoffLoveShell ( IndexType NewId, GeometryType::Pointer pGeometry );
    ConvectiveLinearKirchhoffLoveShell( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /**
     * Destructor.
     */
    virtual ~ConvectiveLinearKirchhoffLoveShell ();

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
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////// computing method //////////////////////////////////////
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                               VectorType& rRightHandSideVector, 
                               ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                                 ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

   


protected:

    friend class Serializer;

    // A private default constructor necessary for serialization
    ConvectiveLinearKirchhoffLoveShell()
    {
    }


    virtual void save ( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
        rSerializer.save( "mCurrentStresses", mCurrentStresses );
    }

    virtual void load ( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
        rSerializer.load( "mCurrentStresses", mCurrentStresses );
    }

    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


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
    double mTotalDomainInitialSize;

    // material parameters
    double mE, mNU , mLambda, mMu /*,mKappa*/ ;

    std::vector<Vector> mNodalCoordinates ;
    std::vector<Vector> mCurrentStresses;

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

    bool mIsPlasticExistence ;
    boost::numeric::ublas::vector<boost::numeric::ublas::vector<ConstitutiveLaw::Pointer>> mConstitutiveLawVector;


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
    void AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight);

    void AddExternalForces(VectorType& RightHandSideVector, const Vector& N
        , const double& DetJ, const double& Weight );

    void AddStiffnessMatrixComponents(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight);
  
    //////////////////// local strain vectors
    void computeMembraneStrain(Matrix& eTensor,  std::vector<Vector>& A, std::vector<Vector>& u_a);
    
    void computeCurvatureChange(Matrix& kTensor, std::vector<Vector>& A, double& A3, Vector& A3Vector, std::vector<Vector>& u_a, 
        boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& u_ab, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab);
    /////////////////////  B matrix

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// base vectors and their derivatives///////////
    /////////////////////////////////////////////////////////////////////////               
    void CovariantBaseVector(std::vector<Vector>& A, const Matrix& DN_De, const std::vector<Vector>& X);

    void ReferenceNormalDirector(Vector& A3Vector, Vector& AA3Vector, double& A3,  std::vector<Vector>& A);

    void ContravariantBaseVector(std::vector<Vector>& AA, std::vector<Vector>& A, Matrix& Aab);

    void DerivativeCovariantBaseVector(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab
            , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const std::vector<Vector>& X);

    void FirstDerivativeDisplacement_a(std::vector<Vector>& u_a, const Matrix& DN_De , const Matrix& u);

    void SecondDerivativeDisplacement_ab(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& u_ab
        , const ShapeFunctionsSecondDerivativesType& D2N_De2,  const Matrix& u);

    void ContinuumCovariantBaseVector(std::vector<Vector>& g, std::vector<Vector>& a
        , std::vector<Vector>& a3Vector_a, Vector& a3Vector, double& theta3);

    void DerivativeReferenceNormalDirector_a(std::vector<Vector>& A3Vector_a, std::vector<Vector>& A
            , Vector& AA3Vector, double& A3, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab);

    /////////////////////////////////////////////////////////////////////////
    ///// first derivative of membrane strain and necessary components //////
    /////////////////////////////////////////////////////////////////////////
    void FirstDerivativeLocalMembraneStrain_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& eLocalVector_r
        , const Matrix& DN_De,  std::vector<Vector>& A) ;

    ///////////////////////////////////////////////////////////////////////////
    ///// first derivative of curvature changes and necessary components //////
    ///////////////////////////////////////////////////////////////////////////
    void FirstDerivativeLocalCurvatureChange_r(boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& kLocalVector_r
        , std::vector<Vector>& A, double& A3, Vector& A3Vector, const Matrix& DN_De
        , const ShapeFunctionsSecondDerivativesType& D2N_De2, boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector>>& A_ab);

    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A);

    void CovariantCurvatureCoefficient(Matrix& Bab
        , boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& A_ab, Vector& A3Vector);

    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////

    void UnitBaseVectors(std::vector<Vector>& e);

    void LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A);
    
    void LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector<Vector>& EE, std::vector<Vector>& AA);

    void CreatingBmatrix(Matrix& BMatrix, const boost::numeric::ublas::vector<boost::numeric::ublas::vector<Vector> >& LocalStrainVector_r);

    //////////////////////////////////////////////////////////
    ////////////////////// material matrix ///////////////////
    //////////////////////////////////////////////////////////
    void CalculateElasticMatrix(Matrix& C, const double& E, const double& NU,  Matrix& Aab);


}; // Class KinematicLinearKirchoffLoveIsogeometricShell 

}  // namespace Kratos.

#endif // KRATOS_KINEMATIC_LINEAR_KIRCHOFF_LOVE_ISOGEOMETRIC_SHELL_H_INCLUDED defined 
