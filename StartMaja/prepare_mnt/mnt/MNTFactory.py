#
# Copyright (C) 2020 Centre National d'Etudes Spatiales (CNES)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Author:         Peter KETTIG <peter.kettig@cnes.fr>, Pierre LASSALLE <pierre.lassalle@cnes.fr>
Project:        StartMaja, CNES
Created on:     Tue Sep 11 15:31:00 2018
"""


class MNTFactory:
    """
    Create a given DEM in Maja format
    """
    def __init__(self, site, platform_id, mission_field, mnt_resolutions, **kwargs):
        self.mnt_type = kwargs.get("mnt_type", "srtm")
        self.site = site
        self.plaform_id = platform_id
        self.mission_field = mission_field
        self.mnt_resolutions = mnt_resolutions
        self.coarse_res = kwargs.get("coarse_res", (240, -240))
        if type(self.coarse_res) in [float, int, str]:
            self.coarse_res = (int(self.coarse_res), -1 * int(self.coarse_res))
        if type(self.coarse_res) != tuple:
            raise TypeError("Unknown coarse_res type %s: %s" % (type(self.coarse_res), self.coarse_res))
        self.kwargs = kwargs

    def factory(self):
        """
        Checks the given mnt_type and returns the class accordingly
        :return: A derived MNTBase object, None if mnt_type unknown.
        """
        from StartMaja.prepare_mnt.mnt.SRTM import SRTM
        # TODO Add more options here: ALOS, TDX, EuDEM...
        if self.mnt_type == "srtm":
            # SRTM is distributed in 90m.
            # Thus, all initial calculation has to be done at this resolution:
            self.site.res_x, self.site.res_y = 90, 90
            return SRTM(site=self.site,
                        **self.kwargs).to_maja_format(platform_id=self.plaform_id,
                                                      mission_field=self.mission_field,
                                                      mnt_resolutions=self.mnt_resolutions,
                                                      coarse_res=self.coarse_res)
        return None