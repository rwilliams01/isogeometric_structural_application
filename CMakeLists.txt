set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIsogeometricStructuralApplication ****")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR}/applications/isogeometric_application )
include_directories( ${CMAKE_SOURCE_DIR}/applications/structural_application )
#######################################################################

if(ADOL_C_DIR AND (PLATE_AND_SHELL_APPLICATION MATCHES ON))
    add_definitions(-DPLATE_AND_SHELL_APPLICATION_IS_ON)
    add_definitions(-DPLATE_AND_SHELL_APPLICATION_USE_ADOL_C)
endif()

######### add compilation units here
set( KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/isogeometric_structural_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/isogeometric_structural_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/isogeometric_posteriori_estimator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/kinematic_linear_isogeometric.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/kinematic_linear_nurbs.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_2phase_small_strain_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/kinematic_linear_kirchhoff_love_isogeometric_shell_rev2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/K_L_large_deformation_shell.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/K_L_linear_shell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/convective_linear_KL_shell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/linear_bending_strip.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nonlinear_bending_strip.cpp
    )

set( KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_SOURCES ${KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_force_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_force_isogeometric_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_pressure_isogeometric_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_load_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_pressure_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/penalty_stiffness_shell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/master_contact_face_3D_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/slave_contact_face_3D_isogeometric.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_link_3D_kinematic_linear_isogeometric.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosIsogeometricStructuralApplication SHARED ${KRATOS_ISOGEOMETRIC_STRUCTURAL_APPLICATION_SOURCES})
target_link_libraries(KratosIsogeometricStructuralApplication KratosCore)
target_link_libraries(KratosIsogeometricStructuralApplication KratosIsogeometricApplication)
if(PLATE_AND_SHELL_APPLICATION MATCHES ON)
    target_link_libraries(KratosIsogeometricStructuralApplication KratosPlateAndShellApplication)
endif()
target_link_libraries(KratosIsogeometricStructuralApplication KratosStructuralApplication)
set_target_properties(KratosIsogeometricStructuralApplication PROPERTIES PREFIX "")
install(TARGETS KratosIsogeometricStructuralApplication DESTINATION libs )
#######################################################################

## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosIsogeometricStructuralApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################

## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################

## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/IsogeometricStructuralApplication.py" DESTINATION KratosMultiphysics )
#######################################################################
