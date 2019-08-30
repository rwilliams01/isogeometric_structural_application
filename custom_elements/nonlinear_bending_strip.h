/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author: G.D.Huynh $
*   Date:                $Date: 25.042018 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


#if !defined(KRATOS_NONLINEAR_BENDING_STRIP_H_INCLUDED )
#define  KRATOS_NONLINEAR_BENDING_STRIP_H_INCLUDED


// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"
#include "custom_utilities/isotropic_tensor_utility.h"

namespace Kratos
{

/*
 * Remarks: this works with 3d geometry
 */
class NonLinearBendingStrip : public Element
{
public:
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;


    // Counted pointer of FacePressureIsogeometric
    KRATOS_CLASS_POINTER_DEFINITION( NonLinearBendingStrip);

    // Constructor void
    NonLinearBendingStrip();

    // Constructor using an array of nodes
    NonLinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry );

    // Constructor using an array of nodes with properties
    NonLinearBendingStrip( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    // Destructor
    virtual ~NonLinearBendingStrip();

    // Name Operations

    virtual Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties ) const;


    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo );

    virtual void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo );

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo );

protected:



private:
    ///@name Static Member Variables

    /// privat variables
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
    double mE, mNU , mLambda, mMu ;

    std::vector<Vector> mNodalCoordinates ;

    // geometric parameters
    unsigned int mNumberOfIntegrationPoint;
    std::vector<Matrix> mInvJ0;
    std::vector<Vector> mN;
    std::vector<Matrix> mDN_De;
    Vector mDetJ0;
    GeometryType::JacobiansType mJ0;
    Vector mIntegrationWeight;

    // privat name Operations
    void CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag );

    ///////////////////////////////////////// all components of residual vectors and stiffness matrices ///////////
    ///////////////////// add left hand side contribution
    void AddInternalForces(VectorType& RightHandSideVector, const Vector& StressResultants,
        std::vector<std::vector<Vector> >& StrainVector_r, const double& DetJ, const double& Weight);


    void AddLinearStiffnessMatrix(MatrixType& LeftHandSideMatrix, const Matrix& Di, const Matrix& BlhsMatrix
        , const Matrix& BrhsMatrix, const double& DetJ, const double& Weight);

    void AddNonlinearStiffnessMatrix( MatrixType& LeftHandSideMatrix,const Vector& StressResultants
            , const std::vector<std::vector<Matrix>>& StrainVector_rs, const double& DetJ, const double& Weight);

    void computeCurvatureChange(Matrix& kTensor, Matrix& Bab, Matrix& bab);

    void DerivativeDirector_r( std::vector<std::vector<Vector>>& a3_rVector
        ,  std::vector<std::vector<Vector>>& aa3_rVector, std::vector<std::vector<double>>& a3_r, Vector& a3Vector, double& a3);

    void DerivativeDirectorNorm_r(std::vector<std::vector<double>>& a3_r
            , std::vector<std::vector<Vector>>& aa3_rVector, Vector& a3Vector);

    void DerivativeNonNormalizedDirector_r( std::vector<std::vector<Vector>>& aa3_rVector
                , std::vector< std::vector<std::vector<Vector>> >& a_ar, std::vector<Vector>& a);

    void SecondDerivativeDirectorNorm_rs(
                std::vector<std::vector<std::vector<std::vector<double>>> >& a3_rs
              , std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
              , Vector& aa3Vector, Vector& a3Vector, double& a3,  std::vector<std::vector<Vector>>& aa3_rVector);

    void SecondDerivativeNonNormalizedDirector_rs(
                std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
              , std::vector< std::vector<std::vector<Vector>> >& a_ar);

    void DerivativeCovariantBaseVector_abr( std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
                ,const ShapeFunctionsSecondDerivativesType& D2N_De2, std::vector<Vector>& UnitBasisVector);

    void SecondDerivativeDirector_rs(
          std::vector<std::vector<std::vector<std::vector<Vector>>> >& a3_rsVector
        , std::vector<std::vector<std::vector<std::vector<Vector>>> >& aa3_rsVector
        , std::vector<std::vector<std::vector<std::vector<double>>> >& a3_rs
        , Vector& a3Vector, double& a3,  std::vector<std::vector<Vector>>& aa3_rVector
        , std::vector<std::vector<double>>& a3_r);

    void FirstDerivativeLocalCurvatureChange_r(std::vector<std::vector<Vector> >& kLocalVector_r
        , std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr
        , std::vector<std::vector<Vector> >& a_ab
        , Vector& a3Vector , std::vector<std::vector<Vector>>& a3_rVector
        , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);

    void SecondDerivativeLocalCurvatureChange_rs(std::vector<std::vector<Matrix> >& kLocalVector_rs
            ,std::vector< std::vector<std::vector<std::vector<Vector>>> >& a_abr, std::vector<std::vector<Vector>>& a3_rVector
            , std::vector<std::vector<Vector> >& a_ab, std::vector<std::vector<std::vector<std::vector<Vector>>> >& a3_rsVector
            , std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);
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

    void DerivativeCovariantBaseVector_r(std::vector< std::vector<std::vector<Vector>> >& a_ar, const Matrix& DN_De, std::vector<Vector>& UnitBasisVector);
    ///////////////////////////////////////////////////////////////////
    ////////// shell fundamental properties ///////////////////////////
    ///////////////////////////////////////////////////////////////////
    void CovariantMetricCoefficient(Matrix& Aab, std::vector<Vector>& A);
    void CovariantCurvatureCoefficient(Matrix& Bab, std::vector<std::vector<Vector> >& A_ab, Vector& A3Vector);


    //////////////////////////////////////////////////////////////
    //////////// addtional utilities /////////////////////////////
    //////////////////////////////////////////////////////////////

    void UnitBaseVectors(std::vector<Vector>& e);

    void LocalCartesianBasisVector(std::vector<Vector>& EE, std::vector<Vector>& A, Vector& A3Vector);

    void CreatingBmatrix(Matrix& BMatrix, const std::vector<std::vector<Vector> >& LocalStrainVector_r);

    void LocalTransformationOfTensor(Matrix& T, Matrix& M, std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff);

    void LocalTransformationCoefficient(std::vector< std::vector<std::vector<std::vector<double>>> >& TransformationCoeff, std::vector<Vector>& EE, std::vector<Vector>& AA);
    /////////////// material matrix
    void CalculateConstitutiveMatrix(Matrix& D);
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }

}; // class FacePressureIsogeometric.

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_ISOGEOMETRIC_H_INCLUDED  defined
