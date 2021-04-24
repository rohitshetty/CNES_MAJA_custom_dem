from StartMaja.prepare_mnt.mnt.MNTBase import MNT
import os
import logging
import math
from StartMaja.Common import FileSystem, ImageTools



class NorDEM(MNT):
  """
    Base class to get norway DEM for a given site 
  """

  def __init__(self, site, **kwargs):

    """
    Initialise an norway DEM.

    :param site: The :class:`prepare_mnt.mnt.SiteInfo` struct containing the basic information.
    :param kwargs: Forwarded parameters to :class:`prepare_mnt.mnt.MNTBase`
    """

    import math 
    super(NorDEM, self).__init__(site, **kwargs)
    if not self.dem_version:
        self.dem_version = 4001




  def prepare_mnt(self):
    """ 
      We expect that the norway DEM file would be in /nordem and work from that assumption there, 
    """
    nordem_full_res = os.path.join(self.wdir, "no10m.tif")
    src_data="/nordem/no10m.tif"

    ImageTools.gdal_warp(src_data, dst=nordem_full_res,
                             r="cubic",
                             te=self.site.te_str,
                             t_srs=self.site.epsg_str,
                             tr=self.site.tr_str,
                             dstnodata=0,
                             srcnodata=0,
                             multi=True)

    return nordem_full_res

if __name__ == "__main__":
    pass
