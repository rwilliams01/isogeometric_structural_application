/*
    see thermal_chem_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Feb 2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_THERMAL_CHEM_APP_COMPUTE_ERROR_ESTIMATE_ELASTICITY_PROCESS_ISOGEOMETRIC_H_INCLUDED )
#define  KRATOS_THERMAL_CHEM_APP_COMPUTE_ERROR_ESTIMATE_ELASTICITY_PROCESS_ISOGEOMETRIC_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "structural_application/custom_utilities/sd_math_utils.h"



namespace Kratos
{

///@name Kratos Globals
///@{

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

class ComputeErrorEstimateElasticityProcessIsogeometric : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeErrorEstimateElasticityProcessIsogeometric
    KRATOS_CLASS_POINTER_DEFINITION(ComputeErrorEstimateElasticityProcessIsogeometric);

    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef typename SparseSpaceType::VectorType SystemVectorType;
    typedef typename LocalSpaceType::VectorType VectorType;
    typedef Element::GeometryType GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeErrorEstimateElasticityProcessIsogeometric(ElementsContainerType& pElements, ProcessInfo& rProcessInfo)
    : Process(), mpElements(pElements), mrProcessInfo(rProcessInfo), mC0_coeff(1.0)
    {}

    ComputeErrorEstimateElasticityProcessIsogeometric(ElementsContainerType& pElements, ProcessInfo& rProcessInfo, const double& C0)
    : Process(), mpElements(pElements), mrProcessInfo(rProcessInfo), mC0_coeff(C0)
    {}

    /// Destructor.
    virtual ~ComputeErrorEstimateElasticityProcessIsogeometric()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        KRATOS_TRY

        Vector RHS_Contribution;
        std::map<std::size_t, double> ErrorEstimate;
//        std::map<std::size_t, double> Weights;
        std::map<std::size_t, std::map<std::size_t, double> > ElementSizes;
        std::map<std::size_t, std::size_t> DofLevel;
        Element::EquationIdVectorType EquationId;
        GeometryData::IntegrationMethod ThisIntegrationMethod;
        Variable<int>& HIERARCHICAL_LEVEL_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("HIERARCHICAL_LEVEL"));
        Variable<double>& NURBS_WEIGHT_var = static_cast<Variable<double>&>(KratosComponents<VariableData>::Get("NURBS_WEIGHT"));
        for (typename ElementsContainerType::ptr_iterator it = mpElements.ptr_begin(); it != mpElements.ptr_end(); ++it)
        {
//            KRATOS_WATCH((*it)->Id())
//            KRATOS_WATCH((*it)->GetValue(HIERARCHICAL_LEVEL_var))

//            if ((*it)->GetValue(HIERARCHICAL_LEVEL_var) != 2)
//                continue;

            /* Compute the estimator */
            ThisIntegrationMethod = (*it)->GetIntegrationMethod();

            GeometryType& rGeometry = (*it)->GetGeometry();

            const unsigned int dim = rGeometry.WorkingSpaceDimension();

            const unsigned int number_of_nodes = rGeometry.size();

            if (RHS_Contribution.size() != number_of_nodes)
                RHS_Contribution.resize(number_of_nodes, false);
            noalias(RHS_Contribution) = ZeroVector(number_of_nodes);

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            rGeometry.Initialize(ThisIntegrationMethod);
            #endif

            //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( ThisIntegrationMethod );
            const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeometry.ShapeFunctionsLocalGradients( ThisIntegrationMethod );
            const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

            //initializing the Jacobian in the reference configuration
            Matrix DeltaPosition(rGeometry.size(), 3);
            for ( unsigned int node = 0; node < rGeometry.size(); ++node )
                noalias( row( DeltaPosition, node ) ) = rGeometry[node].Coordinates() - rGeometry[node].GetInitialPosition();

            GeometryType::JacobiansType J0;
            J0 = rGeometry.Jacobian( J0, ThisIntegrationMethod, DeltaPosition );

            double temp, DetJ0;
            Matrix InvJ0(dim, dim);
            Matrix DN_DX(number_of_nodes, dim);
            std::vector<Vector> D2N_DX2;
            double domain_size = 0.0;
            for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
                noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ0 );

                SD_MathUtils<double>::CalculateSecondDerivatives(rGeometry, D2N_DX2, J0[PointNumber], DN_DX, integration_points[PointNumber]);

                // calculating weights for integration on the reference configuration
                double IntToReferenceWeight = integration_points[PointNumber].Weight();

                // calculate displacement and gradient
                temp = 0.0;
                for (unsigned int node = 0; node < number_of_nodes; ++node)
                    temp += Ncontainer(PointNumber, node) * rGeometry[node].GetSolutionStepValue(TEMPERATURE);

                // modify integration weight in case of 2D
                if ( dim == 2 ) IntToReferenceWeight *= ((*it)->GetProperties()[THICKNESS]);

                // compute domain size
                domain_size += DetJ0 * IntToReferenceWeight;

//KRATOS_WATCH(heat_source)
////KRATOS_WATCH(IntToReferenceWeight)
////KRATOS_WATCH(DetJ0)
//KRATOS_WATCH(DetJ0*IntToReferenceWeight)
//KRATOS_WATCH(row(Ncontainer, PointNumber))
//KRATOS_WATCH(number_of_nodes)

                // TODO add body force
                double aux = 0.0;
                for (unsigned int prim = 0; prim < number_of_nodes; ++prim)
                {
                    const array_1d<double, 3>& disp = rGeometry[prim].GetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int d = 0; d < dim; ++d)
                        aux += D2N_DX2[prim](d) * disp[d];
                }

                for (unsigned int prim = 0; prim < number_of_nodes; ++prim)
                {
                    double tmp = Ncontainer(PointNumber, prim) * std::pow(aux, 2) * DetJ0 * IntToReferenceWeight;
                    RHS_Contribution(prim) += tmp;

                    if (tmp < 0.0)
                    {
                        std::cout << "At element " << (*it)->Id() << ", point " << PointNumber << ", node index " << prim << ":" << std::endl;
                        KRATOS_WATCH(Ncontainer)
                        KRATOS_WATCH(Ncontainer(PointNumber, prim))
                        KRATOS_WATCH(DetJ0)
                        KRATOS_WATCH(IntToReferenceWeight)
                    }
                }
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            rGeometry.Clean();
            #endif

            double element_size = std::pow(domain_size, 1.0/dim) * std::sqrt(dim);
            ElementSizes[(*it)->GetValue(HIERARCHICAL_LEVEL_var)][(*it)->Id()] = element_size;

//            KRATOS_WATCH(RHS_Contribution)

            for (unsigned int prim = 0; prim < number_of_nodes; ++prim)
            {
                std::size_t node_id = rGeometry[prim].Id();

                std::map<std::size_t, double>::iterator ite = ErrorEstimate.find(node_id);

                if (ite == ErrorEstimate.end())
                    ErrorEstimate[node_id] = RHS_Contribution[prim];
                else
                    ite->second += RHS_Contribution[prim];

//                Weights[node_id] = rGeometry[prim].GetValue(NURBS_WEIGHT_var);

                std::map<std::size_t, std::size_t>::iterator it2 = DofLevel.find(node_id);
                if (it2 == DofLevel.end())
                {
                    DofLevel[node_id] = rGeometry[prim].GetValue(HIERARCHICAL_LEVEL_var);
                }
                else
                {
                    if (it2->second < rGeometry[prim].GetValue(HIERARCHICAL_LEVEL_var))
                        it2->second = rGeometry[prim].GetValue(HIERARCHICAL_LEVEL_var);
                }
            }
        }

//        std::cout << "ErrorEstimate at 262:" << std::endl;
//        for (std::map<std::size_t, double>::iterator it = ErrorEstimate.begin(); it != ErrorEstimate.end(); ++it)
//        {
//            std::cout << it->first << ": " << it->second << std::endl;
//        }

        std::map<std::size_t, double> MaxElementSizePerLevel;
        for (std::map<std::size_t, std::map<std::size_t, double> >::iterator it = ElementSizes.begin();
            it != ElementSizes.end(); ++it)
        {
            std::size_t level = it->first;
            MaxElementSizePerLevel[level] = -1.0e99;
            for (std::map<std::size_t, double>::iterator it2 = it->second.begin();
                it2 != it->second.end(); ++it2)
            {
                if (it2->second > MaxElementSizePerLevel[level])
                    MaxElementSizePerLevel[level] = it2->second;
            }
        }

        for (std::map<std::size_t, double>::iterator it = ErrorEstimate.begin(); it != ErrorEstimate.end(); ++it)
        {
//            std::cout << "dof " << it->first << " level: " << DofLevel[it->first] << std::endl;
//            std::cout << "dof " << it->first << " level element size: " << MaxElementSizePerLevel[DofLevel[it->first]] << std::endl;
            it->second = mC0_coeff * MaxElementSizePerLevel[DofLevel[it->first]] * std::sqrt(it->second);
        }

        double err_norm = 0.0;
//        KRATOS_WATCH(ErrorEstimate.size())
        for (std::map<std::size_t, double>::iterator it = ErrorEstimate.begin(); it != ErrorEstimate.end(); ++it)
        {
//            std::cout << it->first << ": " << it->second
////                      << ", w: " << Weights[it->first]
////                      << ", v: " << it->second/Weights[it->first]
//                      << std::endl;
            err_norm += std::pow(it->second, 2);
        }
//        KRATOS_WATCH(err_norm)
        err_norm = std::sqrt(err_norm);
        std::cout << "Error estimator: " << err_norm << std::endl;

        for (typename ElementsContainerType::ptr_iterator it = mpElements.ptr_begin(); it != mpElements.ptr_end(); ++it)
        {
            Element::GeometryType& rGeometry = (*it)->GetGeometry();

            for (std::size_t i = 0; i < rGeometry.size(); ++i)
            {
                rGeometry[i].GetSolutionStepValue(ERROR_RATIO) = ErrorEstimate[rGeometry[i].Id()];
            }
        }

        KRATOS_CATCH("")
    }

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
    virtual std::string Info() const
    {
        return "ComputeErrorEstimateElasticityProcessIsogeometric";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeErrorEstimateElasticityProcessIsogeometric";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ElementsContainerType& mpElements;
    ProcessInfo& mrProcessInfo;
    double mC0_coeff;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeErrorEstimateElasticityProcessIsogeometric& operator=(ComputeErrorEstimateElasticityProcessIsogeometric const& rOther);

    /// Copy constructor.
    //ComputeErrorEstimateElasticityProcessIsogeometric(ComputeErrorEstimateElasticityProcessIsogeometric const& rOther);


    ///@}

}; // Class ComputeErrorEstimateElasticityProcessIsogeometric

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeErrorEstimateElasticityProcessIsogeometric& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeErrorEstimateElasticityProcessIsogeometric& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_THERMAL_CHEM_APP_COMPUTE_ERROR_ESTIMATE_ELASTICITY_PROCESS_ISOGEOMETRIC_H_INCLUDED defined 


