//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: May 31, 2016 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "isogeometric_structural_application.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosIsogeometricStructuralApplication)
    {

        class_<KratosIsogeometricStructuralApplication, KratosIsogeometricStructuralApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosIsogeometricStructuralApplication");

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
