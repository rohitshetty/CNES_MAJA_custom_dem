import csv
import os
import re
import subprocess
import sys

import fiona as fiona
import rasterio
from rasterio.windows import Window
import rasterio.features
import rasterio.warp
import rasterio.mask
from lxml import etree
import shutil
import gdal
from osgeo import ogr, osr


def create_roi(main_image, vector_file, window):
    print('Crop window', window)

    print('Crop input TIF:', main_image)
    with rasterio.open(main_image) as src:
        kwargs = src.meta.copy()
        gt = rasterio.windows.transform(window, src.transform)
        kwargs.update({
            'height': window.height,
            'width': window.width,
            'transform': gt})
        file_out = os.path.join(os.path.dirname(vector_file), 'temp.tif')
        with rasterio.open(file_out, 'w', **kwargs) as dst:
            dst.write(src.read(window=window))

        subprocess.call(['gdaltindex', vector_file, file_out])
        #os.remove(file_out)
        return gt


def update_landsat8(descriptor_file, output_file, gt, crop_size):
    outputs = []
    with open(descriptor_file, newline='') as dfile:
        reader = csv.reader(dfile, delimiter='=')
        for line in reader:
            if line[0].strip() == 'END':
                continue
            key, value = line[0].strip(), line[1].strip()
            if re.match("CORNER_UL_PROJECTION_X_PRODUCT", key):
                outputs.append([line[0], gt[2]])
            elif re.match("CORNER_UL_PROJECTION_Y_PRODUCT", key):
                outputs.append([line[0], gt[5]])
            elif re.match("REFLECTIVE_LINES", key):
                outputs.append([line[0], crop_size])
            elif re.match("REFLECTIVE_SAMPLES", key):
                outputs.append([line[0], crop_size])
            elif re.match("THERMAL_LINES", key):
                outputs.append([line[0], crop_size])
            elif re.match("THERMAL_SAMPLES", key):
                outputs.append([line[0], crop_size])
            else:
                outputs.append(line)

    with open(output_file, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='=', escapechar='', quotechar='', quoting=csv.QUOTE_NONE)
        for out in outputs:
            writer.writerow(out)


def update_s2muscate(descriptor_file, output_file, gt, crop_size):
    print("Fixing S2 muscate xml", descriptor_file, '...')
    root = etree.fromstring(open(descriptor_file).read().encode())
    xpath = "./Geoposition_Informations/Geopositioning/Group_Geopositioning_List"
    node_gps = root.findall(xpath, namespaces=root.nsmap)[0]
    points = node_gps.findall("./Group_Geopositioning", namespaces=root.nsmap)
    for point in points:
        if point.attrib['group_id'] == 'R1':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size)
            node_ncols.text = str(crop_size)
        if point.attrib['group_id'] == 'R2':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size // 2)
            node_ncols.text = str(crop_size // 2)
        if point.attrib['group_id'] == 'R3':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str((crop_size // 2) // 3)
            node_ncols.text = str((crop_size // 2) // 3)
    with open(output_file, 'w') as fp:
        fp.write(etree.tostring(root, pretty_print=True).decode())


def update_venusmuscate(descriptor_file, output_file, gt, crop_size):
    print("Fixing Venus muscate xml", descriptor_file, '...')
    root = etree.fromstring(open(descriptor_file).read().encode())
    xpath = "./Geoposition_Informations/Geopositioning/Group_Geopositioning_List"
    node_gps = root.findall(xpath, namespaces=root.nsmap)[0]
    points = node_gps.findall("./Group_Geopositioning", namespaces=root.nsmap)
    for point in points:
        if point.attrib['group_id'] == 'XS':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size)
            node_ncols.text = str(crop_size)
    with open(output_file, 'w') as fp:
        fp.write(etree.tostring(root, pretty_print=True).decode())


def update_venus(descriptor_file, output_file, gt, crop_size):
    print("Fixing Venus xml", descriptor_file, '...')
    out_f = os.path.join(os.path.split(os.path.split(output_file)[0])[0], os.path.join(os.path.basename(output_file)))
    root = etree.fromstring(open(descriptor_file).read().encode())
    xpath = "./Variable_Header/Specific_Product_Header/Geo_Referencing_Information/Product_Coverage/Cartographic"
    node_carto = root.findall(xpath, namespaces=root.nsmap)[0]
    node_ulx = node_carto.findall("./Upper_Left_Corner/X", namespaces=root.nsmap)[0]
    node_uly = node_carto.findall("./Upper_Left_Corner/Y", namespaces=root.nsmap)[0]
    ulx = gt[2]
    uly = gt[5]
    node_ulx.text = str(ulx)
    node_uly.text = str(ulx)

    node_ulw = node_carto.findall("./Width", namespaces=root.nsmap)[0]
    node_ulh = node_carto.findall("./Height", namespaces=root.nsmap)[0]
    ulw = crop_size
    ulh = crop_size
    node_ulw.text = str(ulw)
    node_ulh.text = str(ulh)

    xpath_size = "./Variable_Header/Specific_Product_Header/Image_Information"
    node_carto_size = root.findall(xpath_size, namespaces=root.nsmap)[0]
    node_lines = node_carto_size.findall("./Size/Lines", namespaces=root.nsmap)[0]
    node_columns = node_carto_size.findall("./Size/Columns", namespaces=root.nsmap)[0]
    node_lines.text = str(crop_size)
    node_columns.text = str(crop_size)

    with open(out_f, 'w') as fp:
        ##print(etree.tostring(root, pretty_print=True))
        fp.write(etree.tostring(root, pretty_print=True).decode())


def update_s2(descriptor_file, output_file, gt, crop_size):
    print("Fixing S2 muscate xml", descriptor_file, '...')
    root = etree.fromstring(open(descriptor_file).read().encode())
    xpath = "./n1:Geometric_Info"
    node_gps = root.findall(xpath, namespaces=root.nsmap)[0]
    points = node_gps.findall("./Size", namespaces=root.nsmap)
    for point in points:
        if point.attrib['resolution'] == '10':
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size)
            node_ncols.text = str(crop_size)
        if point.attrib['resolution'] == '20':
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size // 2)
            node_ncols.text = str(crop_size // 2)
        if point.attrib['resolution'] == '60':
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str((crop_size // 2) // 3)
            node_ncols.text = str((crop_size // 2) // 3)


def update_l8muscate(descriptor_file, output_file, gt, crop_size):
    print("Fixing L8 muscate xml", descriptor_file, '...')
    root = etree.fromstring(open(descriptor_file).read().encode())
    xpath = "./Geoposition_Informations/Geopositioning/Group_Geopositioning_List"
    node_gps = root.findall(xpath, namespaces=root.nsmap)[0]
    points = node_gps.findall("./Group_Geopositioning", namespaces=root.nsmap)
    for point in points:
        if point.attrib['group_id'] == 'XS' or point.attrib['group_id'] == 'TH':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size // 2)
            node_ncols.text = str(crop_size // 2)
        if point.attrib['group_id'] == 'PAN':
            node_ulx = point.findall('./ULX')[0]
            node_uly = point.findall('./ULY')[0]
            ulx = gt[2]
            uly = gt[5]
            node_ulx.text = str(ulx)
            node_uly.text = str(uly)
            node_nrows = point.findall('./NROWS')[0]
            node_ncols = point.findall('./NCOLS')[0]
            node_nrows.text = str(crop_size)
            node_ncols.text = str(crop_size)
    with open(output_file, 'w') as fp:
        fp.write(etree.tostring(root, pretty_print=True).decode())


def crop_images(out_product_path, product_path, vector_file):
    #out_product_path = os.path.join(output_path, os.path.basename(product_path))
    #os.makedirs(out_product_path, exist_ok=True)
    list_cropped_img = []
    for pth in os.listdir(product_path):
        if os.path.isdir(os.path.join(product_path, pth)):
            for f in os.listdir(os.path.join(product_path, pth)):
                file_in = os.path.join(os.path.join(product_path, pth), f)
                if os.path.splitext(f)[1].lower() == '.tif':
                    print("input_file1 : ", file_in)
                    # subprocess.call(['gdal_warp', '-cutline', shp_file, '-crop_to_cutline' file_in, file_out ])
                    with fiona.open(vector_file, "r") as vector_ds:
                        shapes = [feature["geometry"] for feature in vector_ds]
                    file_out = os.path.join(os.path.join(out_product_path, pth), f)
                    if os.path.exists(file_out):
                        os.remove(file_out)
                    with rasterio.open(file_in) as src:
                        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                        out_meta = src.meta
                        out_meta.update({"driver": "GTiff",
                                         "height": out_image.shape[1],
                                         "width": out_image.shape[2],
                                         "transform": out_transform})
                    with rasterio.open(file_out, "w", **out_meta) as dest:
                        dest.write(out_image)
                    list_cropped_img.append(file_out)
        else:
            file_in = os.path.join(product_path, pth)
            if os.path.splitext(pth)[1].lower() == '.tif':
                print("vector file : ", vector_file)
                # Check EPSG
                d_file_in = gdal.Open(file_in)
                proj_file_in = osr.SpatialReference(wkt=d_file_in.GetProjection())
                epsg_file_in = proj_file_in.GetAttrValue('AUTHORITY', 1)

                d_file_in = None
                vds =ogr.Open(vector_file)
                spatialRef = vds.GetLayer().GetNextFeature().GetGeometryRef().GetSpatialReference()
                vproj = osr.SpatialReference(wkt=spatialRef.ExportToPrettyWkt())
                epsg_vector_file = vproj.GetAttrValue('AUTHORITY', 1)
                vds = None
                #d_vector_file = fiona.open(vector_file)
                #spatialRef = d_vector_file.crs
                #epsg_vector_file = spatialRef.get("init").split(":")[1]
                print(epsg_file_in)
                print(epsg_vector_file)

                if epsg_file_in != epsg_vector_file:
                    out_shp = os.path.join(vector_file, epsg_file_in+".shp")
                    print(out_shp)
                    transform_shapefile(int(epsg_vector_file), int(epsg_file_in), os.path.join(vector_file, "roi.shp"), out_shp)
                    vector_file = out_shp



                
                print("input_file 2 : ", file_in)
                # subprocess.call(['gdal_warp', '-cutline', shp_file, '-crop_to_cutline' file_in, file_out ])
                print("new vector file: ", vector_file)
                with fiona.open(vector_file, "r") as vector_ds:
                    shapes = [feature["geometry"] for feature in vector_ds]
                file_out = os.path.join(out_product_path, pth)
                if os.path.exists(file_out):
                    os.remove(file_out)
                with rasterio.open(file_in) as src:
                    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                    out_meta = src.meta
                    out_meta.update({"driver": "GTiff",
                                     "height": out_image.shape[1],
                                     "width": out_image.shape[2],
                                     "transform": out_transform})
                    with rasterio.open(file_out, "w", **out_meta) as dest:
                        dest.write(out_image)
                list_cropped_img.append(file_out)
    return list_cropped_img


def crop_hdr_files(img_list):
    for img in img_list:
        if os.path.exists(img.split('.DBL.TIF')[0] + '.HDR'):
            hdr_file = img.split('.DBL.TIF')[0] + '.HDR'
            print(hdr_file)
            root = etree.fromstring(open(hdr_file).read().encode())
            xpath = "./Variable_Header/Specific_Product_Header/Annex_Information"
            if not root.findall(xpath, namespaces=root.nsmap):
                xpath = "./Variable_Header/Specific_Product_Header/Image_Information"
            node = root.findall(xpath, namespaces=root.nsmap)[0]
            node_lines = node.findall("./Size/Lines", namespaces=root.nsmap)[0]
            node_columns = node.findall("./Size/Columns", namespaces=root.nsmap)[0]
            node_lines.text = str(crop_size)
            node_columns.text = str(crop_size)
            with open(hdr_file, 'w') as fp:
                fp.write(etree.tostring(root, pretty_print=True).decode())


def crop_product(output_path, product_path, vector_file, crop_size):
    print("copy '{}' to '{}'".format(product_path, output_path))
    shutil.copytree(product_path, os.path.join(output_path, os.path.basename(product_path)))
    print("crop_images in '{}'".format(output_path))
    out_product_path = os.path.join(output_path, os.path.basename(product_path))
    os.makedirs(out_product_path, exist_ok=True)
    img_list = crop_images(out_product_path, product_path, vector_file)
    # Venus case
    crop_hdr_files(img_list)


def main(context_dir, product_type, tu_input_dir=None, output_dir=None, offset=60, crop_size=480):


    pats = {
        'LANDSAT8': ['LC08_*', 'LC08_*.*_B1.TIF', 'LC08_*.*MTL.txt', update_landsat8],
        'S2_MUSCATE': ['SENTINEL2B_*.*_V1.0',
                       'SENTINEL2B_*.*R1.tif',
                       'SENTINEL2B_*.*_MTD_*.*.xml$', update_s2muscate],
        'L8_MUSCATE': ['LANDSAT8-*.*_V1.0',
                       'LANDSAT8-*.*XS.tif',
                       'LANDSAT8-*.*_MTD_*.*.xml$', update_l8muscate],
        'VENUS_MUSCATE': ['VENUS-*.*_V1.0',
                          'VENUS-*.*XS.tif',
                          'VENUS-*.*_MTD_*.*.xml$', update_venusmuscate],
        'VENUS': ['VE_VM01_*.*.DBL.DIR',
                  'VE_VM01_*.*TIF',
                  'VE_VM01_VSC_L1VALD_*.*.HDR', update_venus],
        'S2': ['S2A_OPER*.*.SAFE',
               'S2A_OPER_*.*.jp2',
               'S2A_OPER_MTD_*.*.xml$', update_s2],

    }
    window = Window(offset, offset, crop_size, crop_size)
    vector_file = os.path.join(output_dir, 'roi.gml')
    product_dir_pattern = pats[product_type][0]
    main_image_pattern = pats[product_type][1]
    mtd_regex = pats[product_type][2]
    update_descriptor_func = pats[product_type][3]
    products = get_matching_dirs(context_dir, product_dir_pattern)
    images = get_matching_files(products[0], main_image_pattern)

    main_image = images[0]
    gt = create_roi(main_image, vector_file, window)
    ulx = gt[2]
    uly = gt[5]

    for p in os.listdir(context_dir):
        abs_path = os.path.join(context_dir, p)
        if os.path.isdir(abs_path):
            if '_GIP_' in p:
                shutil.copytree(abs_path, os.path.join(output_dir, p))
            elif re.match(product_dir_pattern, p):
                product_path = os.path.join(context_dir, p)
                crop_product(output_dir, product_path, vector_file, crop_size)
                descriptor_file = get_matching_files(product_path, mtd_regex)[0]
                output_file = os.path.join(output_dir,
                                           os.path.basename(product_path),
                                           os.path.basename(descriptor_file))
                update_descriptor_func(descriptor_file, output_file, gt, crop_size)

            elif '_AUX_REFDE2_' in p:
                out_product_path = os.path.join(output_dir, os.path.basename(abs_path))
                os.makedirs(out_product_path, exist_ok=True)
                crop_images(out_product_path, abs_path, vector_file)
            else:
                print("unknown directory", abs_path)
        else:
            if p.endswith('.HDR') and '_AUX_REFDE2_' in p:
                print("Fixing", abs_path, '...')
                root = etree.fromstring(open(abs_path).read().encode())
                xpath = "./Variable_Header/Specific_Product_Header/DEM_Information/Cartographic"
                node_carto = root.findall(xpath, namespaces=root.nsmap)[0]
                node_ulx = node_carto.findall("./Upper_Left_Corner/X", namespaces=root.nsmap)[0]
                node_uly = node_carto.findall("./Upper_Left_Corner/Y", namespaces=root.nsmap)[0]
                node_ulx.text = str(ulx)
                node_uly.text = str(ulx)
                node_lines = node_carto.findall("./Size/Lines", namespaces=root.nsmap)[0]
                node_columns = node_carto.findall("./Size/Columns", namespaces=root.nsmap)[0]
                node_lines.text = str(crop_size)
                node_columns.text = str(crop_size)
                with open(os.path.join(output_dir, p), 'w') as fp:
                    ##print(etree.tostring(root, pretty_print=True))
                    fp.write(etree.tostring(root, pretty_print=True).decode())
            else:
                shutil.copy(abs_path, output_dir)

    print('Cropped context in:', output_dir)
    # file_in = os.path.join(tu_input_dir, 'apTuL1Reader_LANDSAT8_MUSCATE_L2TOAImageOutput.tif')
    if os.path.exists(tu_input_dir):

        tu_output_dir = tu_input_dir.replace(os.sep + "TU" + os.sep, os.sep + "TU.cropped" + os.sep)

        print("copy '{}' to '{}'".format(tu_input_dir, os.path.dirname(tu_output_dir)))
        shutil.copytree(tu_input_dir, tu_output_dir)
        # subprocess.call(['gdalwarp', '-cutline', vector_file, '-crop_to_cutline', file_in, file_out])
        crop_images(tu_output_dir, tu_input_dir, vector_file)
        print('Cropped TU  into: ', tu_output_dir)

    shutil.rmtree(vector_file)
    #os.remove(os.path.splitext(vector_file)[0] + '.xsd')

def get_matching_files(context_dir, pattern):
    results = []

    for subdir, dirs, files in os.walk(context_dir):
        for file in files:
            if re.match(pattern, file):
                results.append(os.path.join(subdir, file))

    if pattern == "VE_VM01_VSC_L1VALD_*.*.HDR":
        context_dir_split = os.path.split(context_dir)[0]
        for s in os.listdir(context_dir_split):
            if re.match(pattern, s):
                results.append(os.path.join(context_dir_split, s))

    if not results:
        raise ValueError("Cannot find file matching {} inside {}".format(
            pattern, context_dir))
    return results


def get_matching_dirs(context_dir, pattern):
    results = []
    for s in os.listdir(context_dir):
        if os.path.isdir(os.path.join(context_dir, s)):
            if re.match(pattern, s):
                results.append(os.path.join(context_dir, s))
    if not results:
        raise ValueError("Cannot find anything matching '{}' in '{}'".format(
            pattern, context_dir))
    return results


def transform_shapefile(in_spatialRef, out_spatial_Ref, in_shp, out_shp):
    print("Reproject the layer : ", in_shp)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # input SpatialReference
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(in_spatialRef)
    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(out_spatial_Ref)
    # CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    # input layer
    inDataSet = driver.Open(in_shp)
    inLayer = inDataSet.GetLayer()
    # utput layer
    if os.path.exists(out_shp):
        driver.DeleteDataSource(out_shp)
    outDataSet = driver.CreateDataSource(out_shp)
    outLayer = outDataSet.CreateLayer("basemap_"+str(out_spatial_Ref), geom_type=ogr.wkbMultiPolygon)
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    # output feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        geom = inFeature.GetGeometryRef()
        geom.Transform(coordTrans)
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        outLayer.CreateFeature(outFeature)
        outFeature = None
        inFeature = inLayer.GetNextFeature()
    #Export projection
    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromEPSG(int(out_spatial_Ref))
    spatialRef.MorphToESRI()
    file = open(os.path.splitext(out_shp)[0]+".prj", 'w')
    file.write(spatialRef.ExportToWkt())
    file.close()
    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None


if __name__ == '__main__':
    context_dir = sys.argv[1]
    product_type = sys.argv[2]
    tu_input_dir = sys.argv[3] if len(sys.argv) > 3 else None
    output_dir = None
    if len(sys.argv) > 4:
        output_dir = sys.argv[4]
    else:
        output_dir = context_dir + ".cropped"
 
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    os.makedirs(output_dir, exist_ok=True)

    offset = 10 if len(sys.argv) < 6 else int(sys.argv[5])
    crop_size = 480 if len(sys.argv) < 7 else int(sys.argv[6])
    print("using crop offset", offset)
    print("using crop size", crop_size)
    main(context_dir, product_type, tu_input_dir, output_dir, offset, crop_size)