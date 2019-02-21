// see thermal_chem_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Feb 2019 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_ADD_PROCESS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_ADD_PROCESS_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  IsogeometricStructuralApplication_AddProcessToPython();

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_ADD_PROCESS_TO_PYTHON_H_INCLUDED  defined
