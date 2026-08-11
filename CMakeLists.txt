cmake_minimum_required(VERSION 3.20)

project(Microphysics
    VERSION 1.0.0
    DESCRIPTION "Building networks in Microphysics with CMake"
    LANGUAGES CXX C)

#----------------------------------------------------------------------------------------------------------------------
# General settings and options
#----------------------------------------------------------------------------------------------------------------------

find_package (Python COMPONENTS Interpreter Development)
find_package(Threads)

include(CTest)
include(CMakePackageConfigHelpers)

#----------------------------------------------------------------------------------------------------------------------
# Network registry
#----------------------------------------------------------------------------------------------------------------------
#
# This file registers only networks self-contained within Microphysics. Networks with paths outside
# Microphysics are registered by the external code instead -- see src/networks/NetworkRegistry.cmake
# for Quokka's.
#
# Fields (unset ones default to empty/FALSE):
#   EOSDIR                    Directory name under EOS/.
#   EOSPARAMFILE              Path to the EOS's _parameters file.
#   PARAMFILE                 The network's unit-test _parameters file. Only read when UNIT_TEST_GATE
#                             is ON, and no problem in this family ships a unit-test _parameters file.
#   NETWORKPARAMFILE          The network's own _parameters file.
#   HAS_NET_FILE              TRUE if the network has a general_null-style .net file for
#                             write_network.py/get_naux.py to process. When FALSE, network_properties.H
#                             comes from the problem's own CMakeLists.txt instead.
#   NETWORKFILE               Path to the .net file (only used if HAS_NET_FILE).
#   NETWORKDIR                Directory holding it (only used if HAS_NET_FILE).
#   UNIT_TEST_GATE            NAME of the BUILD_UNIT_TEST_* cache variable gating whether the unit-test
#                             paramfiles join the write_probin.py call.
#   IS_JAFF                   TRUE if jaff-generated, so write_jaff_network_vars.py runs to derive
#                             NSPEC/SPECIES_ENUM/SPEC_NAMES/CHEM_BANDS.
#   USES_INTEGRATOR_DIRS      TRUE unless the network has no burn_t-based integrator of its own
#                             (gamma_law).
#   EXTRA_SOURCES             Extra .cpp sources.
#   NETWORK_INCLUDE_DIRS      Network-specific include dirs. Order and contents differ per network.

set(_registry_one_value_fields
      EOSDIR EOSPARAMFILE PARAMFILE
      NETWORKPARAMFILE HAS_NET_FILE NETWORKFILE NETWORKDIR
      UNIT_TEST_GATE IS_JAFF USES_INTEGRATOR_DIRS CACHE INTERNAL "")
set(_registry_multi_value_fields EXTRA_SOURCES NETWORK_INCLUDE_DIRS CACHE INTERNAL "")
set(_registry_fields ${_registry_one_value_fields} ${_registry_multi_value_fields} CACHE INTERNAL "")

macro(register_microphysics_network network_name)
  cmake_parse_arguments(_reg "" "${_registry_one_value_fields}" "${_registry_multi_value_fields}" ${ARGN})
  if (DEFINED _reg_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "register_microphysics_network(${network_name}): unrecognized arguments: ${_reg_UNPARSED_ARGUMENTS}")
  endif()

  foreach(_field ${_registry_fields})
    set(${network_name}_${_field} "${_reg_${_field}}" CACHE INTERNAL "")
  endforeach()

  if (NOT network_name IN_LIST KNOWN_NETWORK_NAMES)
    set(KNOWN_NETWORK_NAMES ${KNOWN_NETWORK_NAMES} ${network_name} CACHE INTERNAL "")
  endif()
endmacro()

set(KNOWN_NETWORK_NAMES "" CACHE INTERNAL "")

register_microphysics_network(gamma_law
  EOSDIR            "gamma_law"
  EOSPARAMFILE      "${CMAKE_CURRENT_LIST_DIR}/EOS/gamma_law/_parameters"
  NETWORKPARAMFILE  "${CMAKE_CURRENT_LIST_DIR}/networks/_parameters"
  HAS_NET_FILE      TRUE
  NETWORKFILE       "${CMAKE_CURRENT_LIST_DIR}/networks/general_null/gammalaw.net"
  NETWORKDIR        "${CMAKE_CURRENT_LIST_DIR}/networks/general_null/"
  EXTRA_SOURCES     "${CMAKE_CURRENT_LIST_DIR}/interfaces/eos_data.cpp"
                    "${CMAKE_CURRENT_LIST_DIR}/interfaces/network_initialization.cpp"
  USES_INTEGRATOR_DIRS FALSE
  NETWORK_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/networks"
                        "${CMAKE_CURRENT_LIST_DIR}/networks/general_null"
)

register_microphysics_network(primordial_chem
  EOSDIR            "primordial_chem"
  EOSPARAMFILE      "${CMAKE_CURRENT_LIST_DIR}/EOS/primordial_chem/_parameters"
  PARAMFILE         "${CMAKE_CURRENT_LIST_DIR}/unit_test/burn_cell_primordial_chem/_parameters"
  NETWORKPARAMFILE  "${CMAKE_CURRENT_LIST_DIR}/networks/primordial_chem/_parameters"
  HAS_NET_FILE      TRUE
  NETWORKFILE       "${CMAKE_CURRENT_LIST_DIR}/networks/primordial_chem/pynucastro.net"
  NETWORKDIR        "${CMAKE_CURRENT_LIST_DIR}/networks/primordial_chem/"
  EXTRA_SOURCES     "${CMAKE_CURRENT_LIST_DIR}/interfaces/eos_data.cpp"
                    "${CMAKE_CURRENT_LIST_DIR}/interfaces/network_initialization.cpp"
                    "${CMAKE_CURRENT_LIST_DIR}/networks/primordial_chem/actual_network_data.cpp"
  UNIT_TEST_GATE    "BUILD_UNIT_TEST_PC"
  USES_INTEGRATOR_DIRS TRUE
  NETWORK_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/networks/primordial_chem"
                        "${CMAKE_CURRENT_LIST_DIR}/networks"
)

# Unlike primordial_chem, actual_network_data.cpp is needed here because
# that's where the Semenov opacity table is read in.
register_microphysics_network(metal_chem
  EOSDIR            "metal_chem"
  EOSPARAMFILE      "${CMAKE_CURRENT_LIST_DIR}/EOS/metal_chem/_parameters"
  PARAMFILE         "${CMAKE_CURRENT_LIST_DIR}/unit_test/burn_cell_metal_chem/_parameters"
  NETWORKPARAMFILE  "${CMAKE_CURRENT_LIST_DIR}/networks/metal_chem/_parameters"
  HAS_NET_FILE      TRUE
  NETWORKFILE       "${CMAKE_CURRENT_LIST_DIR}/networks/metal_chem/pynucastro.net"
  NETWORKDIR        "${CMAKE_CURRENT_LIST_DIR}/networks/metal_chem/"
  EXTRA_SOURCES     "${CMAKE_CURRENT_LIST_DIR}/interfaces/eos_data.cpp"
                    "${CMAKE_CURRENT_LIST_DIR}/interfaces/network_initialization.cpp"
                    "${CMAKE_CURRENT_LIST_DIR}/networks/metal_chem/actual_network_data.cpp"
  UNIT_TEST_GATE    "BUILD_UNIT_TEST_MC"
  USES_INTEGRATOR_DIRS TRUE
  NETWORK_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/networks/metal_chem"
                        "${CMAKE_CURRENT_LIST_DIR}/networks"
)

# Builds Microphysics with a given network/EOS for one target in an external code.
#
#   setup_target_for_microphysics_compilation(network_name output_dir
#                                             [INTEGRATOR <name>]
#                                             [<FIELD> <value> ...])
#
# network_name is looked up in the registry; INTEGRATOR picks the integration subdirectory,
# defaulting to Rosenbrock. Any FIELD may be passed directly to override the registry, or to
# describe a network that isn't registered at all.
function(setup_target_for_microphysics_compilation network_name output_dir)

  cmake_parse_arguments(PARSE_ARGV 2 _arg
                        ""
                        "INTEGRATOR;${_registry_one_value_fields}"
                        "${_registry_multi_value_fields}")
  if (DEFINED _arg_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "setup_target_for_microphysics_compilation(${network_name}): unrecognized arguments: ${_arg_UNPARSED_ARGUMENTS}")
  endif()

  set(integrator_name "Rosenbrock")
  if (DEFINED _arg_INTEGRATOR)
    set(integrator_name "${_arg_INTEGRATOR}")
  endif()

  set(_any_field_given FALSE)
  foreach(_field ${_registry_fields})
    if (DEFINED _arg_${_field})
      set(_f_${_field} "${_arg_${_field}}")
      set(_any_field_given TRUE)
    else()
      set(_f_${_field} "${${network_name}_${_field}}")
    endif()
  endforeach()

  if (NOT network_name IN_LIST KNOWN_NETWORK_NAMES AND NOT _any_field_given)
    message(FATAL_ERROR "Given network_name '${network_name}' is not registered (via register_microphysics_network()) "
                         "and no fields were supplied explicitly. Known networks: ${KNOWN_NETWORK_NAMES}")
  endif()

  set(networkparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/_parameters")
  set(integrator_dir "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/${integrator_name}")
  if(NOT EXISTS "${integrator_dir}/actual_integrator.H")
    message(FATAL_ERROR "Given integrator_name ${integrator_name} does not provide ${integrator_dir}/actual_integrator.H")
  endif()
  set(integratorparamfile "${integrator_dir}/_parameters")
  set(integrationparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/_parameters")
  set(unittestparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/unit_test/_parameters")

  set(networkpropfile "${output_dir}/network_properties.H")
  set(networkheadertemplatefile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/network_header.template")

  #DO NOT change the order of the directories below! (NETWORK_INCLUDE_DIRS' own internal
  #order/contents are set per-network in the registry above)
  set(integrator_dirs "")
  if (_f_USES_INTEGRATOR_DIRS)
    set(integrator_dirs ${integrator_dir} ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/utils
                         ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration)
  endif()
  set(_dirs ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/gcem/include
            ${integrator_dirs} ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces
            ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/${_f_EOSDIR}
            ${_f_NETWORK_INCLUDE_DIRS}
            ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/constants)
  # <EOSDIR>_dirs: kept for networks in the same family that already read it by
  # that shared name (see the EOSDIR field comment above). <network_name>_dirs:
  # always exported too, so a caller can use this target's exact network_name
  # without needing to know which EOSDIR family it belongs to.
  set(${_f_EOSDIR}_dirs ${_dirs} PARENT_SCOPE)
  set(${network_name}_dirs ${_dirs} PARENT_SCOPE)

  set(probin_paramfiles "${_f_EOSPARAMFILE} ${_f_NETWORKPARAMFILE} ${networkparamfile} ${integratorparamfile} ${integrationparamfile}")
  if (NOT "${_f_UNIT_TEST_GATE}" STREQUAL "")
    if (${_f_UNIT_TEST_GATE})
      set(probin_paramfiles "${_f_PARAMFILE} ${probin_paramfiles} ${unittestparamfile}")
    endif()
  endif()

  execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/build_scripts/write_probin.py"
                  --pa "${probin_paramfiles}" --use_namespace WORKING_DIRECTORY ${output_dir}/)

  # Same <EOSDIR>_sources / <network_name>_sources duality as <EOSDIR>_dirs above.
  set(_sources ${_f_EXTRA_SOURCES} ${output_dir}/extern_parameters.cpp)
  set(${_f_EOSDIR}_sources ${_sources} PARENT_SCOPE)
  set(${network_name}_sources ${_sources} PARENT_SCOPE)

  if (_f_HAS_NET_FILE)
    execute_process(COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null" python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/get_naux.py" --net "${_f_NETWORKDIR}" --microphysics_path "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/" WORKING_DIRECTORY ${output_dir}/)
    execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${_f_NETWORKFILE}" WORKING_DIRECTORY ${output_dir}/)
  endif()

  if (_f_IS_JAFF)
    set(jaff_network_vars_file "${output_dir}/jaff_network_vars.cmake")
    execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/build_scripts/write_jaff_network_vars.py"
                    --network-parameters "${_f_NETWORKPARAMFILE}"
                    --jaff-toml "${CMAKE_SOURCE_DIR}/src/networks/${network_name}/jaffgen.toml"
                    --output "${jaff_network_vars_file}"
                    WORKING_DIRECTORY ${output_dir}/
                    RESULT_VARIABLE jaff_network_vars_result)
    if (NOT jaff_network_vars_result EQUAL 0)
      message(FATAL_ERROR "write_jaff_network_vars.py failed for network ${network_name}")
    endif()
    set(jaff_network_vars_file "${jaff_network_vars_file}" PARENT_SCOPE)
  endif()

  if (NOT _f_HAS_NET_FILE)
    set(NEXTRASPEC 0)
    set(PROPERTIES "")

    if (_f_IS_JAFF)
      include(${jaff_network_vars_file})
    endif()

    include_directories(
      BEFORE ${_dirs} "${output_dir}"
      "includes/extern_parameters.H" "includes/network_properties.H")

    configure_file(
        ${CMAKE_SOURCE_DIR}/src/networks/network_header.template
        ${networkpropfile}
        @ONLY
    )
  endif()

endfunction()


#set a cache variable that controls whether
#we want to fetch and link amrex or not
set(BUILD_AMReX false CACHE BOOL "Do you want to build and link amrex? (true/false)")
message("Building and linking amrex -- ${BUILD_AMReX}")

if(BUILD_AMReX)
   #fetching amrex
   set (AMReX_TAG "development")
   include(FetchContent)
   FetchContent_Declare(
     amrex
     GIT_REPOSITORY https://github.com/AMReX-Codes/amrex
     GIT_TAG "${AMReX_TAG}"
   )

   # CMake 3.14+
   FetchContent_MakeAvailable(amrex)

   include_directories(${amrex_INCLUDE_DIRS_RET})
   link_libraries(AMReX::amrex)
   message("Built and fetched amrex in Microphysics")
else()
   message("Not building and fetching amrex in Microphysics")
endif()

#set a cache variable that controls whether
#we want to build the unit test or not
set(BUILD_UNIT_TEST_PC false CACHE BOOL "Do you want to build the primordial chem unit test? (true/false)")

set(BUILD_UNIT_TEST_MC false CACHE BOOL "Do you want to build the metal chem unit test? (true/false)")

add_compile_options(-Werror -Wall -Wextra)
add_compile_definitions(NET_LOOP_UNROLL_LEN=1)

#setting sourcefiles and directories needed to make the test here
#so that they are accessible to codes using
#Microphysics as a submodule
#include parent_scope for sourcefiles and directories
#when built as a submodule in an external code
#otherwise these are not visible to the external code
#this will generate a warning but it will build successfully
#do not need unit_test paramfiles when unit_test is not built

if(BUILD_UNIT_TEST_PC AND BUILD_UNIT_TEST_MC)
  message(FATAL_ERROR "Cannot build both primordial chem and metal chem tests at the same time!")
endif()

if(BUILD_UNIT_TEST_PC)
  #Build primordial chem unit test
  setup_target_for_microphysics_compilation("primordial_chem" ${CMAKE_BINARY_DIR})
  include_directories(${primordial_chem_dirs})

  #adding unit_tests as subdirectories
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_test)

elseif(BUILD_UNIT_TEST_MC)
  #Build metal chem unit test
  setup_target_for_microphysics_compilation("metal_chem" ${CMAKE_BINARY_DIR})
  include_directories(${metal_chem_dirs})

  #adding unit_tests as subdirectories
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_test)

  message("Not building primordial chem unit test")
endif()

#Sample command: cmake .. -DBUILD_UNIT_TEST_MC=True -DBUILD_AMReX=True -DAMReX_SPACEDIM=1
