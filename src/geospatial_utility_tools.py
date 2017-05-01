
# a set of functions to calculate the area of a WGS84 pixel analytically
def calculate_area_of_ellipsoidal_slice(lat,a=6378137., b=6356752.3142):

    #convert lat to radians
    f = lat/360.*2*np.pi

    # calculating surface area from equator to latitude
    e = np.sqrt(1 - (b/a)**2)]
    zm = 1 - e*np.sin(f)
    zp = 1 + e*np.sin(f)
    area = np.pi * b**2 * (np.log(zp/zm) / (2*e) + np.sin(f) / (zp*zm))
    return area

def calculate_WGS84_pixel_area(lat1,lat2,long1,long2):

    # distance between east and west boundaries = fraction of whole circle
    q = np.abs(long1-long2)/360.

    # get surface areas from equator to each latitude
    a1 = calculate_area_of_ellipsoidal_slice(lat1)
    a2 = calculate_area_of_ellipsoidal_slice(lat2)

    # difference of these areas multiplied by fraction of 360 degrees taken up by
    # segment bounded by longitude values gives area
    area = q * (np.max([a1,a2]) - np.min(a1,a2))

    return area
