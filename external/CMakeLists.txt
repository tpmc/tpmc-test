include(ExternalProject)

set(python_tpmc_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/python-tpmc")
set(python_tpmc_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/python-tpmc/install")

ExternalProject_Add(python-tpmc
  PREFIX ${python_tpmc_PREFIX}
  GIT_REPOSITORY "https://users.dune-project.org/repositories/projects/dune-mc-git.git"
  GIT_TAG "p/christi/python-tpmc"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Add_Step(python-tpmc build_and_install
  COMMAND python setup.py install --home ${python_tpmc_INSTALL_DIR}
  DEPENDEES update
  WORKING_DIRECTORY ${python_tpmc_PREFIX}/src/python-tpmc)

set(TPMC_INCLUDE_DIR "${python_tpmc_INSTALL_DIR}/include/python" PARENT_SCOPE)
include_directories(${TPMC_INCLUDE_DIR})
set(TPMC_LIBRARIES "-L${python_tpmc_INSTALL_DIR}/lib/python/tpmc/lib -ltpmc_tables" PARENT_SCOPE)
