from math import pow, sin, cos, tan, pi, floor, sqrt, isnan, atan2
from fractions import Fraction
"""
usage:  import MapUtils as MU

        v_lat  =  53.16893
        v_long = -1.76163

        c = MU.ConvertLatLongToXY(v_lat, v_long)
		    print(c.OS())
        # 'SK 16032 63561' Arbor Low

author: dave@bankside-computing.co.uk
thanks: Chris Veness @ http://www.movable-type.co.uk/scripts/latlong.html
        I ported the 'maths' from a lot of this to python from his javascript.

tests:     

    def setUp(self):
        self.c = MU.ConvertLatLongToXY(53.16893, -1.76163)  # Arbor Low

    def test_ConvertLatLongToXY_E(self):
        self.assertEqual(self.c.E(), 416032.3538031285)

    def test_ConvertLatLongToXY_N(self):
        self.assertEqual(self.c.N(), 363561.4399606167)

    def test_ConvertLatLongToXY_OS(self):
        self.assertEqual(self.c.OS(), 'SK 16032 63561')


"""

def northingFromOSGrid(p_OSGridRef):
    """
    Takes an OS BNG reference and returns the northing from it.
    northingFromOSGrid('SK 33997 88370') => 433997
    """

    TwoLetters = p_OSGridRef[0:2]

    try:
        numbers = {
            'SV': '0:0', 'SW': '1:0', 'SX': '2:0', 'SY': '3:0', 'SZ': '4:0', 'TV': '5:0', 'TW': '6:0',
            'SQ': '0:1', 'SR': '1:1', 'SS': '2:1', 'ST': '3:1', 'SU': '4:1', 'TQ': '5:1', 'TR': '6:1',
            'SL': '0:2', 'SM': '1:2', 'SN': '2:2', 'SO': '3:2', 'SP': '4:2', 'TL': '5:2', 'TM': '6:2',
            'SF': '0:3', 'SG': '1:3', 'SH': '2:3', 'SJ': '3:3', 'SK': '4:3', 'TF': '5:3', 'TG': '6:3',
            'SA': '0:4', 'SB': '1:4', 'SC': '2:4', 'SD': '3:4', 'SE': '4:4', 'TA': '5:4', 'TB': '6:4',
            'NV': '0:5', 'NW': '1:5', 'NX': '2:5', 'NY': '3:5', 'NZ': '4:5', 'OV': '5:5', 'OW': '6:5',
            'NQ': '0:6', 'NR': '1:6', 'NS': '2:6', 'NT': '3:6', 'NU': '4:6', 'OQ': '5:6', 'OR': '6:6',
            'NL': '0:7', 'NM': '1:7', 'NN': '2:7', 'NO': '3:7', 'NP': '4:7', 'OL': '5:7', 'OM': '6:7',
            'NF': '0:8', 'NG': '1:8', 'NH': '2:8', 'NJ': '3:8', 'NK': '4:8', 'OF': '5:8', 'OG': '6:8',
            'NA': '0:9', 'NB': '1:9', 'NC': '2:9', 'ND': '3:9', 'NE': '4:9', 'OA': '5:9', 'OB': '6:9',
            'HV': '0:10', 'HW': '1:10', 'HX': '2:10', 'HY': '3:10', 'HZ': '4:10', 'JV': '5:10', 'JW': '6:10',
            'HQ': '0:11', 'HR': '1:11', 'HS': '2:11', 'HT': '3:11', 'HU': '4:11', 'JQ': '5:11', 'JR': '6:11',
            'HL': '0:12', 'HM': '1:12', 'HN': '2:12', 'HO': '3:12', 'HP': '4:12', 'JL': '5:12', 'JM': '6:12'
        }[TwoLetters]

        return numbers[0:1] + p_OSGridRef[3:8]


    except KeyError:
        return 'XX'


def eastingFromOSGrid(p_OSGridRef):
    """
    Takes an OS BNG reference and returns the easting from it.
    eastingFromOSGrid('SK 33997 88370') => 388370
    """

    TwoLetters = p_OSGridRef[0:2]

    try:
        numbers = {
            'SV': '0:0', 'SW': '1:0', 'SX': '2:0', 'SY': '3:0', 'SZ': '4:0', 'TV': '5:0', 'TW': '6:0',
            'SQ': '0:1', 'SR': '1:1', 'SS': '2:1', 'ST': '3:1', 'SU': '4:1', 'TQ': '5:1', 'TR': '6:1',
            'SL': '0:2', 'SM': '1:2', 'SN': '2:2', 'SO': '3:2', 'SP': '4:2', 'TL': '5:2', 'TM': '6:2',
            'SF': '0:3', 'SG': '1:3', 'SH': '2:3', 'SJ': '3:3', 'SK': '4:3', 'TF': '5:3', 'TG': '6:3',
            'SA': '0:4', 'SB': '1:4', 'SC': '2:4', 'SD': '3:4', 'SE': '4:4', 'TA': '5:4', 'TB': '6:4',
            'NV': '0:5', 'NW': '1:5', 'NX': '2:5', 'NY': '3:5', 'NZ': '4:5', 'OV': '5:5', 'OW': '6:5',
            'NQ': '0:6', 'NR': '1:6', 'NS': '2:6', 'NT': '3:6', 'NU': '4:6', 'OQ': '5:6', 'OR': '6:6',
            'NL': '0:7', 'NM': '1:7', 'NN': '2:7', 'NO': '3:7', 'NP': '4:7', 'OL': '5:7', 'OM': '6:7',
            'NF': '0:8', 'NG': '1:8', 'NH': '2:8', 'NJ': '3:8', 'NK': '4:8', 'OF': '5:8', 'OG': '6:8',
            'NA': '0:9', 'NB': '1:9', 'NC': '2:9', 'ND': '3:9', 'NE': '4:9', 'OA': '5:9', 'OB': '6:9',
            'HV': '0:10', 'HW': '1:10', 'HX': '2:10', 'HY': '3:10', 'HZ': '4:10', 'JV': '5:10', 'JW': '6:10',
            'HQ': '0:11', 'HR': '1:11', 'HS': '2:11', 'HT': '3:11', 'HU': '4:11', 'JQ': '5:11', 'JR': '6:11',
            'HL': '0:12', 'HM': '1:12', 'HN': '2:12', 'HO': '3:12', 'HP': '4:12', 'JL': '5:12', 'JM': '6:12'
        }[TwoLetters]

        return numbers[2:3] + p_OSGridRef[9:15]

    except KeyError:
        return 'XX'


def whichTwoLetterCode(p_easting, p_northing):
    """ Determines the OS Two Letter Code from the easting and northing.
    e.g. whichTwoLetterCode(420700,364166)
    'SK'

    Args:
      p_easting:  The easting as an int or float
      p_northing: The northing as an int or float

    Returns:
      The OS Two Letter Code
    """

    e = floor(p_easting  / 100000)
    n = floor(p_northing / 100000)
    # print("e={}, n={}".format(e, n))
    both = str(e) + ":" + str(n)

    try:
        return {
            '0:0': 'SV', '1:0': 'SW', '2:0': 'SX', '3:0': 'SY', '4:0': 'SZ', '5:0': 'TV', '6:0': 'TW',
            '0:1': 'SQ', '1:1': 'SR', '2:1': 'SS', '3:1': 'ST', '4:1': 'SU', '5:1': 'TQ', '6:1': 'TR',
            '0:2': 'SL', '1:2': 'SM', '2:2': 'SN', '3:2': 'SO', '4:2': 'SP', '5:2': 'TL', '6:2': 'TM',
            '0:3': 'SF', '1:3': 'SG', '2:3': 'SH', '3:3': 'SJ', '4:3': 'SK', '5:3': 'TF', '6:3': 'TG',
            '0:4': 'SA', '1:4': 'SB', '2:4': 'SC', '3:4': 'SD', '4:4': 'SE', '5:4': 'TA', '6:4': 'TB',
            '0:5': 'NV', '1:5': 'NW', '2:5': 'NX', '3:5': 'NY', '4:5': 'NZ', '5:5': 'OV', '6:5': 'OW',
            '0:6': 'NQ', '1:6': 'NR', '2:6': 'NS', '3:6': 'NT', '4:6': 'NU', '5:6': 'OQ', '6:6': 'OR',
            '0:7': 'NL', '1:7': 'NM', '2:7': 'NN', '3:7': 'NO', '4:7': 'NP', '5:7': 'OL', '6:7': 'OM',
            '0:8': 'NF', '1:8': 'NG', '2:8': 'NH', '3:8': 'NJ', '4:8': 'NK', '5:8': 'OF', '6:8': 'OG',
            '0:9': 'NA', '1:9': 'NB', '2:9': 'NC', '3:9': 'ND', '4:9': 'NE', '5:9': 'OA', '6:9': 'OB',
            '0:10': 'HV', '1:10': 'HW', '2:10': 'HX', '3:10': 'HY', '4:10': 'HZ', '5:10': 'JV', '6:10': 'JW',
            '0:11': 'HQ', '1:11': 'HR', '2:11': 'HS', '3:11': 'HT', '4:11': 'HU', '5:11': 'JQ', '6:11': 'JR',
            '0:12': 'HL', '1:12': 'HM', '2:12': 'HN', '3:12': 'HO', '4:12': 'HP', '5:12': 'JL', '6:12': 'JM',
        }[both]
    except KeyError:
        return 'XX'


def OS_Ref(p_easting, p_northing):
    """ Determines the OS National Grid Reference from the easting and northing.
        It uses a 'Helmert transformation' and returns the WGS84 point.
    """
    return whichTwoLetterCode(p_easting, p_northing) + " " + str(p_easting)[1:6] + " " + str(p_northing)[1:6]


class Vector3D():
    def __init__(self, p_x, p_y, p_z):
        self.x = p_x
        self.y = p_y
        self.z = p_z


class Point():
    def __init__(self, p_latitude, p_longitude, p_datum='WGS84'):
        self.latitude = p_latitude
        self.longitude = p_longitude
        self.datum = p_datum


class ConvertLatLongToXY:
    """
    Accepts a Latitude and Longtitude as 'p_lat' and 'p_lon'

    Methods:
     E()  : Returns the easting
     N()  : Returns the northing
     OS() : Returns the OS Great Britain National Grid reference
    """

    a  = 6377563.396    # Airy 1830
    b  = 6356256.909    # Airy 1830

    # a  = 6378137.000    # GRS80 aka WGS84
    # b  = 6356752.3141   # GRS80 aka WGS84

    F0 = 0.9996012717  # scale factor on central meridian Airy 1830
    N0 = -100000       # northing of true origin
    E0 = 400000        # easting of true origin
    latTrueOrigin = 49 * pi / 180  # Converted to radians  # lat 49° N
    lonTrueOrigin = -2 * pi / 180  # Converted to radians  # long 2° W

    def __init__(self, p_lat, p_lon):

        self._point = Point(p_lat, p_lon)

        self._latitude = self._point.latitude
        self._longitude = self._point.longitude

        # First thing is to convert this lat and lon to use old datum!
        # var oldCartesian = oldLatLon.toCartesian();
        # var newCartesian = oldCartesian.applyTransform();
        # var newLatLon = newCartesian.toLatLonE(toDatum)

        ConvertLatLongToXY.toOSGB36(self, p_old_lat=self._point.latitude,
                                          p_old_lon=self._point.longitude)

        self.lat = ConvertLatLongToXY.to_radians(self, self._point.latitude)    # Converted to radians
        self.lon = ConvertLatLongToXY.to_radians(self, self._point.longitude)   # Converted to radians

    @staticmethod
    def to_radians(self, value):
        return value * pi / 180

    @staticmethod
    def from_radians(self, value):
        return value * 180 / pi

    def toOSGB36(self, p_old_lat, p_old_lon):

        self.v3D = ConvertLatLongToXY.toCartesian(self)
        self.v3D = ConvertLatLongToXY.apply_transform(self)
        self._point = ConvertLatLongToXY.toLatLonE(self)
        return

    def apply_transform(self):
        x1 = self.v3D.x
        y1 = self.v3D.y
        z1 = self.v3D.z

        tx = -446.448
        ty = 125.157
        tz = -542.060

        rx = ConvertLatLongToXY.to_radians(self, -0.1502 / 3600) # normalise seconds to radians
        ry = ConvertLatLongToXY.to_radians(self, -0.2470 / 3600) # normalise seconds to radians
        rz = ConvertLatLongToXY.to_radians(self, -0.8421 / 3600) # normalise seconds to radians

        s1 = 20.4894 / 1e6 + 1 # normalise ppm to(s + 1)

        x2 = tx + x1 * s1 - y1 * rz + z1 * ry
        y2 = ty + x1 * rz + y1 * s1 - z1 * rx
        z2 = tz - x1 * ry + y1 * rx + z1 * s1

        return Vector3D(x2, y2, z2)

    def toLatLonE(self):
        x = float(self.v3D.x)
        y = float(self.v3D.y)
        z = float(self.v3D.z)

        a = 6377563.396
        b = 6356256.909
        f = 1 / 299.3249646

        e2 = 2 * f - f * f       # 1st eccentricity squared = (a²-b²) / a²) / a
        E2 = e2 / (1 - e2)       # 2nd eccentricity squared = (a²-b²) / b²) / b
        p = sqrt((pow(x, 2)) + (pow(y, 2))) # distance from minor axis
        R = sqrt(p * p + z * z)  # polar radius

        # parametric latitude (Bowring eqn 17, replacing tan = z·a / p·b)
        tanB = (b * z) / (a * p) * (1 + e2 * b / R)
        sinB = tanB / sqrt(1 + tanB * tanB)
        cosB = sinB / tanB

        # geodetic latitude(Bowring eqn 18: tanphi = z + E²bsinnB / p - e²cos)o)
        if isnan(cosB):
            phi2 = 0
        else:
            phi2 = atan2(z + E2 * b * sinB * sinB * sinB, p - e2 * a * cosB * cosB * cosB)

        # longitude
        lambda_long = atan2(y, x)

        # height above ellipsoid (Bowring eqn 7) [not currently used]
        # sinphi = sin(phi)
        # cosphi = cos(phi)

        # v1 = a / sqrt(1 - e2 * sinphi * sinphi) # length ophi the normal terminated by the minor axis
        # h = p * cosphi + z * sinphi - (a * a / v1)
        return Point(ConvertLatLongToXY.from_radians(self, phi2), ConvertLatLongToXY.from_radians(self, lambda_long))


    def toCartesian(self):
        h = 0
        a = 6378137
        f = Fraction(1 / 298.257223563)

        phi = ConvertLatLongToXY.to_radians(self, self._point.latitude)
        lamdba = ConvertLatLongToXY.to_radians(self, self._point.longitude)

        sinphi = sin(phi)
        cosphi = cos(phi)
        sinlamdba = sin(lamdba)
        coslamdba = cos(lamdba)

        eSq = (2 * f) - (f * f)

        v = a / sqrt(1 - eSq * sinphi * sinphi)

        _x = (v + h) * cosphi * coslamdba
        _y = (v + h) * cosphi * sinlamdba
        _z = (v * (1 - eSq) + h) * sinphi

        return Vector3D(_x, _y, _z)

    def degrees_to_decimal(self, degrees, minutes, seconds):
        d_secs = seconds / 60
        d_mins = minutes + d_secs
        d_degrees = degrees + (d_mins / 60)
        return (d_degrees)

    def n(self):
        return (ConvertLatLongToXY.a - ConvertLatLongToXY.b) / (ConvertLatLongToXY.a + ConvertLatLongToXY.b)

    def e2(self):
        return (pow(self.a, 2) - pow(self.b, 2)) / pow(self.a, 2)

    def v(self):
        return (self.a * self.F0 * (pow((1.0 - self.e2() * pow(sin(self.lat), 2)), -0.5)))

    def p(self):
        return (self.a * self.F0 * (1 - self.e2()) * pow((1 - self.e2() * (pow(sin(self.lat), 2))), -1.5))

    def n2(self):
        return ((self.v() / self.p()) - 1)

    def M(self):
        # two_thirds = Fraction(2, 3)
        FirstBit = (1 + self.n() + (Fraction(5 / 4) * pow(self.n(), 2)) + (Fraction(5 / 4) * pow(self.n(), 3))) * (
        self.lat - self.latTrueOrigin)

        SecondBit = (((3 * self.n()) + (3 * pow(self.n(), 2)) + Fraction(21 / 8) * pow(self.n(), 3)) * sin(
            (self.lat - self.latTrueOrigin)) * cos((self.lat + self.latTrueOrigin)))
        FirstLine = FirstBit - SecondBit

        ThirdBit = (((Fraction(15 / 8) * pow(self.n(), 2)) +
                     (Fraction(15 / 8) * pow(self.n(), 3))) *
                    sin(2 * (self.lat - self.latTrueOrigin)) * cos(2 * (self.lat + self.latTrueOrigin)))
        FourthBit = (Fraction(35 / 24) * pow(self.n(), 3) * sin(3 * (self.lat - self.latTrueOrigin)) * cos(
            (3 * (self.lat + self.latTrueOrigin))))
        SecondLine = ThirdBit - FourthBit

        FifthBit = self.b * self.F0

        All = FifthBit * (FirstLine + SecondLine)
        return All

    def I(self):
        return self.M() + self.N0

    def II(self):
        return (self.v() / 2) * sin(self.lat) * cos(self.lat)

    def III(self):
        return ((self.v() / 24) * sin(self.lat) * pow(cos(self.lat), 3) * (5 - pow(tan(self.lat), 2) + 9 * self.n2()))

    def IIIA(self):
        return (self.v() / 720) * sin(self.lat) * pow(cos(self.lat), 5) * (
        61 - 58 * pow(tan(self.lat), 2) + pow(tan(self.lat), 4))

    def IV(self):
        return self.v() * cos(self.lat)

    def V(self):
        return (self.v() / 6) * pow(cos(self.lat), 3) * (self.v() / self.p() - pow(tan(self.lat), 2))

    def VI(self):
        return (self.v() / 120) * pow(cos(self.lat), 5) * (5 - (18 * pow(tan(self.lat), 2)) +
                                                           (pow(tan(self.lat), 4) + (14 * self.n2()) - (
                                                           58 * pow(tan(self.lat), 2) * self.n2())))

    def N(self):
        return self.I() + (self.II() * (pow((self.lon - self.lonTrueOrigin), 2))) + (
        self.III() * (pow((self.lon - self.lonTrueOrigin), 4))) + (
               self.IIIA() * (pow((self.lon - self.lonTrueOrigin), 6)))

    def E(self):
        return self.E0 + self.IV() * (self.lon - self.lonTrueOrigin) + self.V() * pow((self.lon - self.lonTrueOrigin),
                                                                                      3) + self.VI() * pow(
            (self.lon - self.lonTrueOrigin), 5)

    def OS(self):
        return whichTwoLetterCode(self.E(), self.N()) + " " + str(self.E())[1:6] + " " + str(self.N())[1:6]
        # return whichTwoLetterCode( floor(self.E()  / 100000), floor(self.N()  / 100000)) + " " + str(self.E())[1:6] + " " + str(self.N())[1:6]

    def __call__(self):
        return self.OS()
