/**
 * M_PNRW: Celestial Transformation Engine Class
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 * Logic: M = P * N * R * W
 */


export class M_PNRW {
    /* These are the standard parameters for the WGS84 geodetic reference system, used in GPS and most modern geodetic calculations */

    // The semi-major axis of the WGS84 ellipsoid, defined as 6,378,137 meters
    static WGS84_A= 6378137

    // The flattening factor of the WGS84 ellipsoid, defined as 1/298.257223563
    static WGS84_F= 1 / 298.257223563

    // The square of the first eccentricity of the WGS84 ellipsoid
    static WGS84_E_SQ= (2 * M_PNRW.WGS84_F) - (M_PNRW.WGS84_F ** 2)

    constructor() {
        // --- CONSTANTS ---
        this.rad_deg= Math.PI / 180
        this.mas_to_rad= (Math.PI / 180) / 3600000

        // Principal Nutation Terms (IAU 2000B simplified)
        this.nutation_terms= [
            { l: 0, lp: 0, f: 0, d: 0, om: 1, psi: -17206, eps: 9205 },
            { l: 0, lp: 0, f: 2, d: -2, om: 2, psi: -1317, eps: 573 },
            { l: 0, lp: 0, f: 2, d: 0, om: 2, psi: -227, eps: 97 },
            { l: 0, lp: 1, f: 0, d: 0, om: 0, psi: 147, eps: 7 }
        ]
    }

    // --- CONVERSION & UTILITY METHODS ---

    /**
     * Calculates the Julian Date from a JavaScript Date object.
     * @param {Date} date - The date to convert.
     * @returns {number} The Julian Date.
     */
    get_julian_date(date) {
        if (!(date instanceof Date)) throw new Error('date must be a Date object')
        if (isNaN(date.getTime())) throw new Error('Invalid date')

        var year= date.getUTCFullYear()
        var month= date.getUTCMonth() + 1
        var day= date.getUTCDate()
        var hours= date.getUTCHours() + (date.getUTCMinutes() / 60) + (date.getUTCSeconds() / 3600)
        
        var jd= (367 * year) - 
                 Math.floor(7 * (year + Math.floor((month + 9) / 12)) / 4) + 
                 Math.floor(275 * month / 9) + 
                 day + 1721013.5 + (hours / 24)

        return jd
    }

    /**
     * Converts ECEF coordinates to a 3D vector.
     * @param {number} lat - Latitude in degrees
     * @param {number} lon - Longitude in degrees
     * @param {number} height - Height above ellipsoid in meters
     * @param {number} a - Semi-major axis of the ellipsoid (meters). Defaults to WGS84.
     * @param {number} e_sq - Square of the first eccentricity. Defaults to WGS84.
     * @returns {number[]} [x, y, z] in meters (ECEF)
     */
    to_ecef(lat, lon, height, a = M_PNRW.WGS84_A, e_sq = M_PNRW.WGS84_E_SQ) {
        var rad_lat= lat * this.rad_deg
        var rad_lon= lon * this.rad_deg
        var n_phi= a / Math.sqrt(1 - e_sq * Math.pow(Math.sin(rad_lat), 2))
        
        var x= (n_phi + height) * Math.cos(rad_lat) * Math.cos(rad_lon)
        var y= (n_phi + height) * Math.cos(rad_lat) * Math.sin(rad_lon)
        var z= (n_phi * (1 - e_sq) + height) * Math.sin(rad_lat)

        return [x, y, z]
    }

    /**
     * Calculates the distance between two 3D vectors. This is needed since in space, 
     * the distance between two points is not the same as the distance on the surface of the Earth.
     * @param {number[]} v1 - The first vector [x, y, z].
     * @param {number[]} v2 - The second vector [x, y, z].
     * @returns {number} The distance between the two vectors.
     */
    get_distance(v1, v2) {
        var dx= v1[0] - v2[0]
        var dy= v1[1] - v2[1]
        var dz= v1[2] - v2[2]

        return Math.sqrt((dx * dx) + (dy * dy) + (dz * dz))
    }

    // --- MATRIX MATH ---

    /**
     * Multiplies two 3x3 matrices.
     * @param {number[][]} a - The first matrix.
     * @param {number[][]} b - The second matrix.
     * @returns {number[][]} The resulting matrix.
     */
    multiply_matrices(a, b) {
        var res= [[0,0,0], [0,0,0], [0,0,0]]

        for (var i= 0; i < 3; i++) {
            for (var j= 0; j < 3; j++) {
                for (var k= 0; k < 3; k++) {
                    res[i][j] += a[i][k] * b[k][j]
                }
            }
        }

        return res
    }

    /**
     * Multiplies a 3x3 matrix by a 3D vector.
     * @param {number[][]} m - The matrix.
     * @param {number[]} v - The vector.
     * @returns {number[]} The resulting vector.
     */
    multiply_matrix_vector(m, v) {
        return [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        ]
    }

    // --- PNRW COMPONENT GENERATORS ---

    /**
     * Generates a rotation matrix for rotation around the y-axis.
     * @param {number} angle - The angle of rotation in radians.
     * @returns {number[][]} The rotation matrix.
     */
    rot_y(angle) {
        var c= Math.cos(angle), s= Math.sin(angle)
        
        return [[c, 0, -s], [0, 1, 0], [s, 0, c]]
    }

    /**
     * Generates a rotation matrix for rotation around the z-axis.
     * @param {number} angle - The angle of rotation in radians.
     * @returns {number[][]} The rotation matrix.
     */
    rot_z(angle) {
        var c= Math.cos(angle), s= Math.sin(angle)

        return [[c, s, 0], [-s, c, 0], [0, 0, 1]]
    }

    /**
     * Generates the Precession matrix (P) using IAU 1976 precession angles.
     * @param {number} t - Time in Julian centuries from J2000.0.
     * @returns {number[][]} The precession matrix.
     */
    get_p_matrix(t) {
        // Proper IAU 1976 Precession angles
        var zeta= (2306.2181 * t) * (this.rad_deg / 3600)
        var z= (2306.2181 * t) * (this.rad_deg / 3600)
        var theta= (2004.3109 * t) * (this.rad_deg / 3600)
        
        // P = Rz(-z) * Ry(theta) * Rz(-zeta)
        var rz_zeta= this.rot_z(-zeta)
        var ry_theta= this.rot_y(theta)
        var rz_z= this.rot_z(-z)
        
        return this.multiply_matrices(rz_z, this.multiply_matrices(ry_theta, rz_zeta))
    }

    /**
     * Generates the Nutation matrix (N) using IAU 2000B nutation terms.
     * @param {number} t - Time in Julian centuries from J2000.0.
     * @returns {number[][]} The nutation matrix.
     */
    get_n_matrix(t) {
        // The 5 fundamental Delaunay arguments
        var l= (134.963 + 477198.868 * t) * this.rad_deg
        var lp= (357.529 + 35999.050 * t) * this.rad_deg
        var f= (93.272 + 483202.018 * t) * this.rad_deg
        var d= (297.850 + 445267.111 * t) * this.rad_deg
        var om= (125.045 - 1934.136 * t) * this.rad_deg

        var d_psi= 0, d_eps= 0
        var eps= (23.43929 - 0.013 * t) * this.rad_deg

        for (var i= 0; i < this.nutation_terms.length; i++) {
            var term= this.nutation_terms[i]
            var arg= (term.l * l) + (term.lp * lp) + (term.f * f) + (term.d * d) + (term.om * om)
            d_psi += term.psi * Math.sin(arg)
            d_eps += term.eps * Math.cos(arg)
        }

        var dp= d_psi * this.mas_to_rad, de= d_eps * this.mas_to_rad

        return [
            [1, -dp * Math.cos(eps), -dp * Math.sin(eps)],
            [dp * Math.cos(eps), 1, -de],
            [dp * Math.sin(eps), de, 1]
        ]
    }

    /**
     * Generates the Rotation matrix (R) using Greenwich Mean Sidereal Time (GMST).
     * @param {number} jd - The Julian Date.
     * @returns {number[][]} The rotation matrix.
     */
    get_r_matrix(jd) {
        // Correct GMST calculation treating time dynamically
        var days_since_j2000= jd - 2451545.0
        var gmst_hours= (18.697374558 + 24.06570982441908 * days_since_j2000) % 24
        if (gmst_hours < 0) gmst_hours += 24
        
        var theta= gmst_hours * 15 * this.rad_deg

        return [
            [Math.cos(theta), Math.sin(theta), 0],
            [-Math.sin(theta), Math.cos(theta), 0],
            [0, 0, 1]
        ]
    }

    /**
     * Generates the Wobble matrix (W) using polar motion coordinates.
     * @param {number} xp - Polar motion in the x-direction (radians).
     * @param {number} yp - Polar motion in the y-direction (radians).
     * @returns {number[][]} The wobble matrix.
     */
    get_w_matrix(xp= 0, yp= 0) {
        // Polar motion — pass xp, yp in radians from IERS bulletin for full accuracy
        return [
            [1,   0,  xp],
            [0,   1, -yp],
            [-xp, yp,  1]
        ]
    }

    // --- MAIN ENGINE ---

    /**
     * Calculates the Universal Address for a given location and time.
     * @param {number} lat - Latitude in degrees (WGS84 ellipsoid)
     * @param {number} lon - Longitude in degrees (WGS84 ellipsoid)
     * @param {number} height - Height above ellipsoid in meters
     * @param {Date} date - The date and time for the transformation
     * @param {number[]} sun_vector - The position vector of the Sun or Barycenter in ICRS [x, y, z] in km
     * @returns {number[]} The Universal Space Vector [x, y, z] in meters
     */
    get_universal_address(lat, lon, height, date, sun_vector) {
        if (typeof lat !== 'number' || lat < -90 || lat > 90) {
            throw new Error('Latitude must be a number between -90 and 90')
        }
        if (typeof lon !== 'number' || lon < -180 || lon > 180) {
            throw new Error('Longitude must be a number between -180 and 180')
        }
        if (typeof height !== 'number' || height < 0) {
            throw new Error('Height must be a number greater than or equal to 0')
        }
        if (!(date instanceof Date)) {
            throw new Error('Date must be a Date object')
        }
        if (!Array.isArray(sun_vector) || sun_vector.length !== 3) {
            throw new Error('Sun vector must be an array of 3 numbers')
        }

        var jd= this.get_julian_date(date)
        // t = Julian centuries from J2000.0
        var t= (jd - 2451545.0) / 36525
        // local_v = ECEF coordinates of the location
        var local_v= this.to_ecef(lat, lon, height)
        
        var p= this.get_p_matrix(t)
        var n= this.get_n_matrix(t)
        var r= this.get_r_matrix(jd)
        var w= this.get_w_matrix()
        // m_pn = P * N
        var m_pn= this.multiply_matrices(p, n)
        // m_pnr = P * N * R
        var m_pnr= this.multiply_matrices(m_pn, r)
        // master_m = P * N * R * W
        var master_m= this.multiply_matrices(m_pnr, w)
        // space_v = master_m * local_v
        var space_v= this.multiply_matrix_vector(master_m, local_v)
        // universal_address = space_v + sun_vector
        return [
            space_v[0] + sun_vector[0],
            space_v[1] + sun_vector[1],
            space_v[2] + sun_vector[2]
        ]
    }
}