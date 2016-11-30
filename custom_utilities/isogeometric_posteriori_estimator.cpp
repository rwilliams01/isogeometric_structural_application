//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 May 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#include <cmath>
#include "includes/deprecated_variables.h"
#include "isogeometric_posteriori_estimator.h"
#include "utilities/math_utils.h"

// #define DEBUG_ESTIMATOR

namespace Kratos
{

double IsogeometricPosterioriEstimator::ComputeSimplePosterioriError(ModelPart& r_model_part)
{
    #ifdef DEBUG_ESTIMATOR
    std::ofstream logfile;
    logfile.open("estimator.log");

    // export the nodal stress
//    Vector DummyStressVector(6);
//    Vector DummyStrainVector(6);
//    for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
//    {
//        noalias(DummyStressVector) = (*it).GetSolutionStepValue(STRESSES);
//        logfile << "Node " << (*it).Id() << " stresses: " << DummyStressVector << std::endl;
//        noalias(DummyStrainVector) = (*it).GetSolutionStepValue(STRAIN);
//        logfile << "Node " << (*it).Id() << " strain: " << DummyStrainVector << std::endl;
//    }
    #endif

    // check if model_part has nodal stresses and nodal strain
    if(!r_model_part.GetProcessInfo()[HAS_STRESSES_AT_NODE])
        KRATOS_THROW_ERROR(std::logic_error, "Nodal stress is not yet available", "")

    if(!r_model_part.GetProcessInfo()[HAS_STRAIN_AT_NODE])
        KRATOS_THROW_ERROR(std::logic_error, "Nodal strain is not yet available", "")

    double posteriori_error_nom = 0.0;
    double posteriori_error_denom = 0.0;

    // get the element array from model_part
    ElementsArrayType& ElementsArray = r_model_part.Elements();

    //create a partition of the elements
    int number_of_threads = omp_get_max_threads();
    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

    KRATOS_WATCH( number_of_threads )
    KRATOS_WATCH( element_partition )

    #pragma omp parallel for
    for(int k = 0; k < number_of_threads; ++k)
    {
        #ifdef DEBUG_ESTIMATOR
        std::ofstream logfile_thread;
        std::stringstream ss;
        ss << "estimator_thread_" << k << ".log";
        logfile_thread.open(ss.str().c_str());
        #endif

        typename ElementsArrayType::ptr_iterator it_begin = ElementsArray.ptr_begin() + element_partition[k];
        typename ElementsArrayType::ptr_iterator it_end = ElementsArray.ptr_begin() + element_partition[k + 1];

        Vector StrainVector;
        Vector StressVector;
        double dV, DetJ, Aux1, Aux2;
        double posteriori_error_nom_thread = 0.0;
        double posteriori_error_denom_thread = 0.0;
        for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
        {
            if(!(*it)->GetValue(IS_INACTIVE))
            {
                #ifdef DEBUG_ESTIMATOR
                logfile_thread << "Element " << (*it)->Id() << ":" << std::endl;
                #endif
            
                // get the integration_points of the element
                const GeometryType::IntegrationPointsArrayType& integration_points =
                (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                // get the stress at integration_points
                std::vector<Vector> IntegrationPointStresses;
                (*it)->GetValueOnIntegrationPoints(STRESSES, IntegrationPointStresses, r_model_part.GetProcessInfo());

                // get the strain at integration_points
                std::vector<Vector> IntegrationPointStrains;
                (*it)->GetValueOnIntegrationPoints(STRAIN, IntegrationPointStrains, r_model_part.GetProcessInfo());

                #ifdef DEBUG_ESTIMATOR
                for(unsigned int i = 0; i < IntegrationPointStresses.size(); ++i)
                    logfile_thread << "stress " << i << ":" << IntegrationPointStresses[i] << std::endl;
                for(unsigned int i = 0; i < IntegrationPointStrains.size(); ++i)
                    logfile_thread << "strain " << i << ":" << IntegrationPointStrains[i] << std::endl;
                #endif

                // compute the  Jacobian
                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                
                #ifdef DEBUG_ESTIMATOR
                logfile_thread << "J:" << J << std::endl;
                #endif
                
                // iterate through integration_points
                for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
                {
                    // compute the nodal stress at integration point
                    StressVector = CalculateOnPoint(STRESSES, StressVector, (*it), integration_points[PointNumber]);

                    // compute the nodal strain at integration point
                    StrainVector = CalculateOnPoint(STRAIN, StrainVector, (*it), integration_points[PointNumber]);

                    // contribute to the posteriori error
                    DetJ = MathUtils<double>::Det(J[PointNumber]);
                    dV = DetJ * integration_points[PointNumber].Weight();
                    Aux1 = fabs(inner_prod(StressVector - IntegrationPointStresses[PointNumber], StrainVector - IntegrationPointStrains[PointNumber])) * dV;
                    Aux2 = fabs(inner_prod(StressVector, StrainVector)) * dV;

                    posteriori_error_nom_thread += Aux1;
                    posteriori_error_denom_thread += Aux2;
                }
            }
        }
        #ifdef DEBUG_ESTIMATOR
        std::cout << "thread " << k << ": " << posteriori_error_nom_thread << "," << posteriori_error_denom_thread << std::endl;
        #endif
        #pragma omp critical
        {
            posteriori_error_nom += posteriori_error_nom_thread;
            posteriori_error_denom += posteriori_error_denom_thread;
        }
        #ifdef DEBUG_ESTIMATOR
        logfile_thread.close();
        #endif
    }

    #ifdef DEBUG_ESTIMATOR
    logfile.close();
    #endif

    return sqrt(posteriori_error_nom / posteriori_error_denom);
}

void IsogeometricPosterioriEstimator::ComputeSimplePosterioriErrorOnNodes(
    const Variable<double>& rThisVariable,
    ModelPart& r_model_part,
    LinearSolverType::Pointer pSolver)
{
    // check if model_part has nodal stresses and nodal strain
    if(!r_model_part.GetProcessInfo()[HAS_STRESSES_AT_NODE])
        KRATOS_THROW_ERROR(std::logic_error, "Nodal stress is not yet available", "")

    if(!r_model_part.GetProcessInfo()[HAS_STRAIN_AT_NODE])
        KRATOS_THROW_ERROR(std::logic_error, "Nodal strain is not yet available", "")

    // get the element array from model_part
    ElementsArrayType& ElementsArray = r_model_part.Elements();

    //create a partition of the elements
    int number_of_threads = omp_get_max_threads();
    vector<unsigned int> element_partition;
    CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

    KRATOS_WATCH( number_of_threads )
    KRATOS_WATCH( element_partition )

    // create and initialize matrix and vectors for transferring the error to node
    int NumberOfNodes = r_model_part.NumberOfNodes();
    SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
    noalias(M)= ZeroMatrix(NumberOfNodes, NumberOfNodes);

    SerialSparseSpaceType::VectorType g(NumberOfNodes);
    noalias(g)= ZeroVector(NumberOfNodes);

    SerialSparseSpaceType::VectorType b(NumberOfNodes);
    noalias(b)= ZeroVector(NumberOfNodes);

    // create a map from node Id to matrix/vector row
    std::map<unsigned int, unsigned int> MapNodeIdToVec;
    unsigned int cnt = 0;
    for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
        MapNodeIdToVec[it->Id()] = cnt++;

    // create the structure for M a priori
    ConstructMatrixStructure(M, ElementsArray, MapNodeIdToVec, r_model_part.GetProcessInfo());

    //create the array of lock for matrix/vector assembly
    std::vector< omp_lock_t > lock_array(NumberOfNodes);
    for(unsigned int i = 0; i < NumberOfNodes; ++i)
        omp_init_lock(&lock_array[i]);

    #pragma omp parallel for
    for(int k = 0; k < number_of_threads; ++k)
    {
        typename ElementsArrayType::ptr_iterator it_begin = ElementsArray.ptr_begin() + element_partition[k];
        typename ElementsArrayType::ptr_iterator it_end = ElementsArray.ptr_begin() + element_partition[k + 1];
        
        Vector StrainVector;
        Vector StressVector;
        double dV, DetJ, Aux;
        unsigned int row, col;
        for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
        {
            if(!(*it)->GetValue(IS_INACTIVE))
            {
                // get the integration_points of the element
                const GeometryType::IntegrationPointsArrayType& integration_points =
                (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                // get the stress at integration_points
                std::vector<Vector> IntegrationPointStresses;
                (*it)->GetValueOnIntegrationPoints(STRESSES, IntegrationPointStresses, r_model_part.GetProcessInfo());

                // get the strain at integration_points
                std::vector<Vector> IntegrationPointStrains;
                (*it)->GetValueOnIntegrationPoints(STRAIN, IntegrationPointStrains, r_model_part.GetProcessInfo());

                // compute the  Jacobian
                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                
                // compute the shape function values
                GeometryType::ShapeFunctionsGradientsType DN_De;
                Matrix Ncontainer;
                IsogeometricGeometryType& rIsogeometricGeometry = dynamic_cast<IsogeometricGeometryType&>((*it)->GetGeometry());
                rIsogeometricGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                    Ncontainer,
                    DN_De,
                    (*it)->GetIntegrationMethod()
                );

                // iterate through integration_points
                for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
                {
                    // compute the nodal stress at integration point
                    StressVector = CalculateOnPoint(STRESSES, StressVector, (*it), integration_points[PointNumber]);

                    // compute the nodal strain at integration point
                    StrainVector = CalculateOnPoint(STRAIN, StrainVector, (*it), integration_points[PointNumber]);

                    // contribute to the posteriori error
                    DetJ = MathUtils<double>::Det(J[PointNumber]);
                    dV = DetJ * integration_points[PointNumber].Weight();
                    Aux = fabs(inner_prod(StressVector - IntegrationPointStresses[PointNumber], StrainVector - IntegrationPointStrains[PointNumber])) * dV;

                    for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                    {
                        row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];
                        omp_set_lock(&lock_array[row]);
                        b(row) += Aux * Ncontainer(PointNumber, prim) * dV;
                        for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                        {
                            col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                            M(row, col) += Ncontainer(PointNumber, prim) * Ncontainer(PointNumber, sec) * dV;
                        }
                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
            else
            {
                // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                {
                    row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];
                    omp_set_lock(&lock_array[row]);
                    for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                    {
                        col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                        if(col == row)
                            M(row, col) += 1.0;
                    }
                    omp_unset_lock(&lock_array[row]);
                }
            }
        }
    }

    for(unsigned int i = 0; i < NumberOfNodes; ++i)
        omp_destroy_lock(&lock_array[i]);

    // solver the system
    pSolver->Solve(M, g, b);

    // transfer the solution to the nodal variables
    for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
    {
        unsigned int row = MapNodeIdToVec[it->Id()];
        it->GetSolutionStepValue(rThisVariable) = g(row);
    }
    std::cout << "Compute posteriori error for " << rThisVariable.Name() << " completed" << std::endl;
}

}

#undef DEBUG_ESTIMATOR

