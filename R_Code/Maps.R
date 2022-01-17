# Maps #
source ("R_Code/Packages.R")
shape_path <- "/Users/quimbayojp/Dropbox/Manuscripts/2_Manuscripts_Review/Ms Global cleaners/Shapefiles/"
ocean_shapefile <- paste(shape_path, "ne_10m_ocean/ne_10m_ocean.shp", sep="")
land_shapefile  <- paste(shape_path, "ne_10m_land/ne_10m_land.shp", sep="")

layer_ocean <- ogrListLayers(ocean_shapefile)
ogrInfo(ocean_shapefile, layer=layer_ocean)
ocean_poly <- readOGR(ocean_shapefile, layer = layer_ocean)

layer_land <- ogrListLayers(land_shapefile)
ogrInfo(land_shapefile, layer=layer_land)
land_poly <- readOGR(land_shapefile, layer = layer_land)

# Robinson CRS Projection
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
ocean_poly_proj <- spTransform(ocean_poly, CRSobj = robin_crs)
land_poly_proj <- spTransform(land_poly, CRSobj = robin_crs)
