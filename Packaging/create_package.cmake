# Copyright (C) 2005-2017 Centre National d'Etudes Spatiales (CNES)
#
# This file is part of Orfeo Toolbox
#
#     https://www.orfeo-toolbox.org/
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
set(pkg_target_dir   "/opt/maja/${PKG_MAJA_VERSION_STRING}")
if(WITH_VALIDATION)
  set(pkg_target_dir "${PKG_STAGE_DIR}")
endif()
macro(create_package)
  install(CODE "
    message(STATUS \"Creating ${PACKAGE_OUTPUT_FILE}\") 
      execute_process(
        COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/makeself/makeself.sh
        \"--quiet\"
        \"--target\"
        \"${pkg_target_dir}\"
        \"${CMAKE_INSTALL_PREFIX}/${PKG_STAGE_DIR}\"
        \"${PKG_STAGE_DIR}.run\"
        \"${PACKAGE_LONG_NAME} ${PKG_MAJA_VERSION_STRING}\"
        \"./pkgsetup\"
        WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\"
     )
  ")
endmacro()

