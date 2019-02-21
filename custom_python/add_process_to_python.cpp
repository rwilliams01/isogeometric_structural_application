// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Feb 2019 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_process_to_python.h"
#include "custom_processes/compute_error_estimate_elasticity_process_isogeometric.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void IsogeometricStructuralApplication_AddProcessToPython()
{

    class_<ComputeErrorEstimateElasticityProcessIsogeometric, ComputeErrorEstimateElasticityProcessIsogeometric::Pointer, bases<Process>, boost::noncopyable>
    ("ComputeErrorEstimateElasticityProcessIsogeometric", init<ModelPart::ElementsContainerType&, ProcessInfo&>())
    .def(init<ModelPart::ElementsContainerType&, ProcessInfo&, const double&>())
    ;

}

}  // namespace Python.

}  // namespace Kratos.

