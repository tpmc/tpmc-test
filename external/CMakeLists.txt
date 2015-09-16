include(ExternalProject)

set(tpmc_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/tpmc")
set(tpmc_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/tpmc/install")

ExternalProject_Add(tpmc
  PREFIX ${tpmc_PREFIX}
  GIT_REPOSITORY "https://www.github.com/tpmc/tpmc.git"
  GIT_TAG "master"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Add_Step(tpmc build_and_install
  COMMAND python setup.py install --home ${tpmc_INSTALL_DIR}
  DEPENDEES update
  WORKING_DIRECTORY ${tpmc_PREFIX}/src/tpmc)

set(TPMC_INCLUDE_DIR "${tpmc_INSTALL_DIR}/lib/python/tpmc/tpmc/include" PARENT_SCOPE)
include_directories(${TPMC_INCLUDE_DIR})
set(TPMC_LIBRARIES "-L${tpmc_INSTALL_DIR}/lib/python/tpmc/lib -ltpmc_tables" PARENT_SCOPE)
