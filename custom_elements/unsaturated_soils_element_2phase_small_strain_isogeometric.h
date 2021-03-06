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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Oct 2014 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_UNSATURATED_SOILS_ELEMENT_2PHASE_SMALL_STRAIN_ISOGEOMETRIC_INCLUDED )
#define  KRATOS_UNSATURATED_SOILS_ELEMENT_2PHASE_SMALL_STRAIN_ISOGEOMETRIC_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"

namespace Kratos
{

///@name Kratos Globals
///@{

extern Variable<Vector> PRESTRESS;
extern Variable<Vector> PLASTIC_STRAIN_VECTOR;
extern Variable<double> PLASTICITY_INDICATOR;
extern Variable<double> EXCESS_PORE_WATER_PRESSURE;
extern Variable<double> WATER_PRESSURE;
extern Variable<double> WATER_PRESSURE_DT;
extern Variable<double> WATER_PRESSURE_EINS;
extern Variable<double> WATER_PRESSURE_NULL;
extern Variable<double> REFERENCE_WATER_PRESSURE;
extern Variable<Vector> STRESSES;
extern Variable<double> PRESTRESS_FACTOR;
extern Variable<double> OVERCONSOLIDATION_RATIO;
extern Variable<int> PARENT_ELEMENT_ID;
extern Variable<int> INTEGRATION_POINT_INDEX;
extern Variable<int> NEIGHBOUR_EXPANSION_LEVEL;
extern Variable<Vector> RECOVERY_STRESSES;
extern Variable<int> STRESS_RECOVERY_TYPE;
extern Variable<Vector> COORDINATES;
extern Variable<Vector> FLUID_FLOWS;

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */

class UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric : public Element
{

public:
    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw ConstitutiveLawType;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    /// Counted pointer of UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric

    KRATOS_CLASS_POINTER_DEFINITION( UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric( IndexType NewId, GeometryType::Pointer pGeometry );
    UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric();


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    IntegrationMethod GetIntegrationMethod() const;

    Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void InitializeJacobian();

    virtual void ResetConstitutiveLaw();

    virtual void EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const;

    virtual void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const;

    virtual void InitializeSolutionStep( ProcessInfo& CurrentProcessInfo );

    virtual void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo );

    virtual void DampMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo );

    virtual void MassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo);

    virtual void FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    virtual void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

    //************************************************************************************
    virtual void CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValuesVector( Vector& values, int Step );

    virtual void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable, std::vector<ConstitutiveLaw::Pointer>& rValues, const ProcessInfo& rCurrentProcessInfo );

//    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

//    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);


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
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    bool mIsInitialized;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    IntegrationMethod mThisIntegrationMethod;

    double mTotalDomainInitialSize;

    unsigned int mNodesPressMin;
    unsigned int mNodesPressMax;
    unsigned int mNodesDispMin;
    unsigned int mNodesDispMax;
    Matrix mInitialDisp;
    std::vector<double> mReferencePressures;
    bool mIsStabilised;

    std::vector< Matrix > mInvJ0;
    Vector mDetJ0;

    ///@}
    ///@name Private Operators
    ///@{
    /** K += weight*Btrans*D*B */
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
    );

    void InitializeVariables();

    void InitializeMaterial();

    void CalculateAndAddExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& mResidualVector,
        double weight
    );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //ASSEMBLE STIFFNESS-MATRICES AND FORCE-VECTORS

    void AssembleStiffness( Matrix& HelpLeftHandSideMatrix,
                            Matrix& K_UU, Matrix& K_UW,
                            Matrix& K_WU, Matrix& K_WW );

    void AssembleDamping( Matrix& HelpDampingMatrix,
                          Matrix& D_UU, Matrix& D_UW,
                          Matrix&  D_WU, Matrix&  D_WW );

    void AssembleRHS( Vector& rRightHandSideVector, Vector& R_U, Vector& R_W );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //DO TIME INTEGRATION STIFFNESS AND FORCES

    void AssembleTimeSpaceRHSFromSubVectors( VectorType& rRightHandSideVector,
                                             const Vector& R_U, const Vector& R_W );

    void AssembleTimeSpaceStiffnessFromDampSubMatrices( MatrixType& rLeftHandSideMatrix,
                                                        const Matrix& D_UU, const Matrix& D_UW,
                                                        const Matrix& D_WU, const Matrix& D_WW );

    void AssembleTimeSpaceStiffnessFromStiffSubMatrices( MatrixType& rLeftHandSideMatrix,
                                                         const Matrix& K_UU, const Matrix& K_UW,
                                                         const Matrix& K_WU, const Matrix& K_WW );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE FORCEVECTORS DISPLACEMENT

    void AddBodyForcesToRHSVectorU( Vector& R, Vector& N_DISP, double density, double Weight, double detJ );

    void AddInternalForcesToRHSU( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE STIFFNESS MATRICES DISPLACEMENT

    void CalculateStiffnesMatrixUU( Matrix& K, const
                                    Matrix& C, const Matrix& B_Operator, const Matrix& DN_DX_DISP, Vector& N_DISP,
                                    double density, double capillaryPressure, double Weight, double detJ );

    void CalculateStiffnesMatrixUW( Matrix& Help_K_UW, Matrix& tanC_W,
                                    const Matrix& DN_DX_DISP, Vector& N_DISP, Vector& N_PRESS,
                                    double capillaryPressure, double Weight, double DetJ );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE FORCEVECTORS WATER

    void AddInternalForcesToRHSW( Vector& Help_R_U, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double  DetJ );

    void AddInternalForcesToRHSWs( Vector& Help_R_W, Vector& N_PRESS, Vector& N_PRESS_averaged, double capillaryPressure_Dt, double averageCapillaryPressure_Dt, double Weight, double DetJ );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE STIFFNESS MATRICES WATER

    void CalculateStiffnesMatrixWU( Matrix& Help_K_WU, const Matrix& DN_DX, const Matrix& DN_DX_PRESS, Vector& N, double capillaryPressure, double Weight, double DetJ );

    void CalculateStiffnesMatrixWW( Matrix& Help_K_WW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ );

    void CalculateStiffnesMatrixWWs( Matrix& Help_K_WW, Vector& N_PRESS, Vector& N_PRESS_averaged, double Weight, double DetJ );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE DAMPING MATRICES WATER

    void CalculateDampingMatrixWU( Matrix& Help_D_WU, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ );

    void CalculateDampingMatrixWW( Matrix& Help_D_WW, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ );

    void CalculateDampingMatrixWWs( Matrix& Help_D_WW, const Matrix& DN_DX_DISP, Vector& N_PRESS, Vector& N_PRESS_averaged, double Weight, double DetJ );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //PRIMARY VARIABLES AND THEIR DERIVATIVES

    Matrix CalculateDisplacementGradient( const Matrix& DN_DX_DISP );

    void GetDerivativeDPressuresDt( const Vector& N_PRESS, double& capillaryPressure_Dt, double& waterPressure_Dt );

    void GetPressures( const Vector& N, double& capillaryPressure, double& waterPressure );

    double GetDerivativeDCapillaryPressureDt( const Vector& N_PRESS );

    Vector GetGradientWater( const Matrix& DN_DX_PRESS );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //POROSITY AND ITS DERIVATIVES

    double GetPorosity( const Matrix& DN_DX_DISP );

    double GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP );
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //AVERAGED DENSITY
    Vector GetGravity();
    double GetDivU( const Matrix& DN_DX_DISP );
    double GetDerivativeDDivUDt( const Matrix& DN_DX_DISP );
    double GetAveragedDensity( double capillaryPressure, double porosity );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //SATURATION AND ITS DERIVATIVES

    double GetSaturation( double capillaryPressure );

    double GetDerivativeDSaturationDpc( double capillaryPressure );

    double GetSecondDerivativeD2SaturationDpc2( double capillaryPressure );

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //WATER FLOW AND ITS DERIVATIVES

    Vector GetFlowWater( const Matrix& DN_DX, const Matrix& DN_DX_DISP, double capillaryPressure );

    Vector GetDerivativeDWaterFlowDpw( const Matrix& DN_DX, const Matrix& DN_DX_DISP, double capillaryPressure );

    double GetDerivativeDWaterFlowDGradpw( const Matrix& DN_DX_DISP, double capillaryPressure );

    //************************************************************************************
    //************************************************************************************
    //STRESSES, STRAINS AND CONSTITUTIVE MODELL (UNSATURATED CASE)
    //************************************************************************************
    //************************************************************************************

    void CalculateEffectiveStress( Vector& StressVector, Matrix& tanC_W, const double waterPressure );

    void CalculateStressAndTangentialStiffnessUnsaturatedSoils( Vector& StressTensor,
                                                                Matrix& tanC_U,
                                                                Matrix& tanC_W,
                                                                Vector& StrainVector,
                                                                double waterPressure,
                                                                int PointNumber,
                                                                const ProcessInfo& rCurrentProcessInfo );
    //************************************************************************************
    //************************************************************************************
    //STRESSES, STRAINS AND CONSTITUTIVE MODELL
    //************************************************************************************
    //************************************************************************************

    Matrix GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, int PointNumber );

    //************************************************************************************
    //************************************************************************************
    //NONLINEAR CONTRIBUTION OF VOLUMECHANGE
    //************************************************************************************
    //************************************************************************************
//    double Determinant_DeformationTensor(const Matrix& DN_DX_DISP);
    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX );

    void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector );

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric() {}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric& operator=(const UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric& rOther);

    /// Copy constructor.
    //UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric(const UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric& rOther);


    ///@}

}; // Class UnsaturatedSoilsElement_2phase_SmallStrain_Isogeometric

}  // namespace Kratos.

#endif // KRATOS_UNSATURATED_SOILS_ELEMENT_INCLUDED defined


