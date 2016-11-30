//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_POSTERIORI_ESTIMATOR_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_POSTERIORI_ESTIMATOR_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes 
#include <omp.h>
#include "boost/progress.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "isogeometric_application/custom_geometries/isogeometric_geometry.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication

///@{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
extern Variable<bool> HAS_STRAIN_AT_NODE;
extern Variable<bool> HAS_STRESSES_AT_NODE;


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
/*** Detail class definition.
 */
class IsogeometricPosterioriEstimator
{
public:
    ///@name Type Definitions
    ///@{
    typedef typename ModelPart::NodesContainerType NodesArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsArrayType;

    typedef typename ModelPart::ConditionsContainerType ConditionsArrayType;

    typedef typename Element::GeometryType GeometryType;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SerialSparseSpaceType;

    typedef UblasSpace<double, Matrix, Vector> SerialDenseSpaceType;

    typedef LinearSolver<SerialSparseSpaceType, SerialDenseSpaceType> LinearSolverType;

    /// Pointer definition of IsogeometricPosterioriEstimator
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricPosterioriEstimator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricPosterioriEstimator()
    {
    }

    /// Destructor.
    virtual ~IsogeometricPosterioriEstimator()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Simple posteriori error estimator based on nodal stress & strain
    /// Reference: Matlab Implementation of the Finite Element Method in Elasticity, Alberty et al
    double ComputeSimplePosterioriError(ModelPart& r_model_part);
    void ComputeSimplePosterioriErrorOnNodes(const Variable<double>& rThisVariable,
                                             ModelPart& r_model_part,
                                             LinearSolverType::Pointer pSolver);

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
        std::stringstream buffer;
        buffer << "A collection of posteriori estimators for isogeometric method";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Calculate vector variable (Stress, Strain) at a local point
    Vector& CalculateOnPoint(const Variable<Vector>& rVariable,
                             Vector& rResult,
                             Element::Pointer& pElement,
                             const CoordinatesArrayType& rCoordinates)
    {
        Vector N;
        pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);
        
        for(unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
        {
            Vector& NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
            
            if(i == 0)
            {
                if(rResult.size() != NodalValues.size())
                    rResult.resize(NodalValues.size());
                noalias(rResult) = N( i ) * NodalValues;
            }
            else
            {
                noalias(rResult) += N( i ) * NodalValues;
            }
        }
        return rResult;
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while ( i != endit && (*i) != candidate)
        {
            ++i;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    void ConstructMatrixStructure (
        SerialSparseSpaceType::MatrixType& A,
        ElementsArrayType& rElements,
        std::map<unsigned int, unsigned int> MapNodeIdToVec,
        ProcessInfo& CurrentProcessInfo
    )
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        Element::EquationIdVectorType ids;
        for(typename ElementsArrayType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; ++i_element)
        {
            ids.resize((i_element)->GetGeometry().size());
            for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                ids[i] = MapNodeIdToVec[(i_element)->GetGeometry()[i].Id()];

            for(std::size_t i = 0 ; i < ids.size() ; ++i)
            {
                if(ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for(std::size_t j = 0 ; j < ids.size() ; ++j)
                    {
                        if(ids[j] < equation_size)
                            AddUnique(row_indices, ids[j]);
                    }
                }
            }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; ++i)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
#ifndef _OPENMP
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
        for( int k=0; k < number_of_threads; ++k )
        {
            #pragma omp parallel
            if( omp_get_thread_num() == k )
            {
                for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; i++ )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(unsigned int i = 1; i < number_of_threads; ++i)
            partitions[i] = partitions[i-1] + partition_size ;
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IsogeometricPosterioriEstimator& operator=(IsogeometricPosterioriEstimator const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricPosterioriEstimator(IsogeometricPosterioriEstimator const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricPosterioriEstimator

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, IsogeometricPosterioriEstimator& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const IsogeometricPosterioriEstimator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
#undef DEBUG_MULTISOLVE
#undef DEBUG_GENERATE_MESH
#undef ENABLE_PROFILING

#endif
