/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 9, 2013 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/isogeometric_posteriori_estimator.h"
#include "custom_utilities/isogeometric_contact_utility.h"

namespace Kratos
{

namespace Python
{

void IsogeometricStructuralApplication_AddCustomUtilitiesToPython()
{
    using namespace boost::python;
    
    class_<IsogeometricContactUtility, IsogeometricContactUtility::Pointer, boost::noncopyable>
    ("IsogeometricContactUtility", init<int>())
//    .def("CalculateOnIntegrationPoints", &IsogeometricContactUtility::CalculateOnIntegrationPoints)
//    .def("ComputeCOG", &IsogeometricContactUtility::ComputeCOG)
    .def("SetUpContactConditions", &IsogeometricContactUtility::SetUpContactConditions)
    .def("Update", &IsogeometricContactUtility::Update)
    .def("IsConverged", &IsogeometricContactUtility::IsConverged)
    .def("Clean", &IsogeometricContactUtility::Clean)
    ;

    class_<IsogeometricPosterioriEstimator, IsogeometricPosterioriEstimator::Pointer, boost::noncopyable>
    ("IsogeometricPosterioriEstimator", init<>())
    .def("ComputeSimplePosterioriError", &IsogeometricPosterioriEstimator::ComputeSimplePosterioriError)
    .def("ComputeSimplePosterioriErrorOnNodes", &IsogeometricPosterioriEstimator::ComputeSimplePosterioriErrorOnNodes)
    ;

}

}  // namespace Python.

} // Namespace Kratos

