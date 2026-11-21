//! M_PNRW: Celestial Transformation Engine
//!
//! This program is free software: you can redistribute it and/or modify
//! it under the terms of the GNU General Public License as published by
//! the Free Software Foundation, either version 3 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU General Public License for more details.
//!
//! You should have received a copy of the GNU General Public License
//! along with this program.  If not, see <https://www.gnu.org/licenses/>.
//!
//! Logic: `M = P * N * R * W`

use std::error::Error;
use std::f64::consts::PI;
use std::fmt;

/// A 3D vector `[x, y, z]`.
pub type Vec3 = [f64; 3];

/// A 3x3 matrix in row-major order.
pub type Mat3 = [[f64; 3]; 3];

/// A principal nutation term (IAU 2000B simplified).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NutationTerm {
    pub l: i32,
    pub lp: i32,
    pub f: i32,
    pub d: i32,
    pub om: i32,
    pub psi: i32,
    pub eps: i32,
}

/// A UTC calendar date and time used for Julian Date calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UtcDateTime {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
}

impl UtcDateTime {
    /// Creates and validates a UTC date and time.
    pub fn new(
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
    ) -> Result<Self, MpnrwError> {
        let dt = Self {
            year,
            month,
            day,
            hour,
            minute,
            second,
        };
        dt.validate()?;
        Ok(dt)
    }

    /// Converts a Unix timestamp (seconds since `1970-01-01T00:00:00Z`) to a `UtcDateTime`.
    pub fn from_unix_timestamp(timestamp_secs: i64) -> Self {
        let days = timestamp_secs.div_euclid(86_400);
        let secs_of_day = timestamp_secs.rem_euclid(86_400) as u32;

        // Howard Hinnant's civil-from-days algorithm
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let doe = (z - era * 146_097) as u32;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let y = (yoe as i64) + era * 400;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let d = doy - (153 * mp + 2) / 5 + 1;
        let m = if mp < 10 { mp + 3 } else { mp - 9 };
        let year = (y + if m <= 2 { 1 } else { 0 }) as i32;

        let hour = (secs_of_day / 3600) as u8;
        let minute = ((secs_of_day % 3600) / 60) as u8;
        let second = (secs_of_day % 60) as u8;

        Self {
            year,
            month: m as u8,
            day: d as u8,
            hour,
            minute,
            second,
        }
    }

    /// Validates date and time field ranges.
    pub fn validate(&self) -> Result<(), MpnrwError> {
        if !(1..=12).contains(&self.month) {
            return Err(MpnrwError::InvalidDate("month must be between 1 and 12"));
        }

        let is_leap = (self.year % 4 == 0 && self.year % 100 != 0) || (self.year % 400 == 0);
        let max_day = match self.month {
            1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
            4 | 6 | 9 | 11 => 30,
            2 if is_leap => 29,
            2 => 28,
            _ => unreachable!(),
        };

        if self.day < 1 || self.day > max_day {
            return Err(MpnrwError::InvalidDate("day is out of range for month"));
        }
        if self.hour > 23 {
            return Err(MpnrwError::InvalidDate("hour must be between 0 and 23"));
        }
        if self.minute > 59 {
            return Err(MpnrwError::InvalidDate("minute must be between 0 and 59"));
        }
        if self.second > 59 {
            return Err(MpnrwError::InvalidDate("second must be between 0 and 59"));
        }

        Ok(())
    }
}

/// Unpacked 64-bit signed integer fields from a 40-byte "Sealed Envelope" buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DepackedVector {
    pub lat: i64,
    pub lon: i64,
    pub height: i64,
    pub date: i64,
    pub sun: i64,
}

/// Errors returned by `M_PNRW` operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MpnrwError {
    InvalidLatitude,
    InvalidLongitude,
    InvalidHeight,
    InvalidDate(&'static str),
    InvalidSunVector,
    BufferTooShort { expected: usize, actual: usize },
}

impl fmt::Display for MpnrwError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLatitude => {
                write!(f, "Latitude must be a finite number between -90 and 90")
            }
            Self::InvalidLongitude => {
                write!(f, "Longitude must be a finite number between -180 and 180")
            }
            Self::InvalidHeight => {
                write!(f, "Height must be a finite number greater than or equal to 0")
            }
            Self::InvalidDate(reason) => write!(f, "Invalid date: {reason}"),
            Self::InvalidSunVector => write!(f, "Sun vector components must be finite numbers"),
            Self::BufferTooShort { expected, actual } => {
                write!(f, "Buffer too short: expected at least {expected} bytes, got {actual}")
            }
        }
    }
}

impl Error for MpnrwError {}

/// Idiomatic Rust alias for [`M_PNRW`].
pub type MPnrw = M_PNRW;

/// `M_PNRW`: Celestial Transformation Engine.
///
/// Transforms WGS84 geodetic coordinates into a celestial space vector via
/// `M = P * N * R * W` (Precession, Nutation, Rotation, Wobble).
#[allow(non_camel_case_types)]
#[derive(Debug, Clone)]
pub struct M_PNRW {
    /// Degrees-to-radians conversion factor (`PI / 180`).
    pub rad_deg: f64,
    /// Milliarcseconds-to-radians conversion factor (`(PI / 180) / 3_600_000`).
    pub mas_to_rad: f64,
    /// Principal Nutation Terms (IAU 2000B simplified).
    pub nutation_terms: [NutationTerm; 4],
}

impl Default for M_PNRW {
    fn default() -> Self {
        Self::new()
    }
}

impl M_PNRW {
    /// The semi-major axis of the WGS84 ellipsoid, defined as 6,378,137 meters.
    pub const WGS84_A: f64 = 6_378_137.0;

    /// The flattening factor of the WGS84 ellipsoid, defined as `1 / 298.257223563`.
    pub const WGS84_F: f64 = 1.0 / 298.257223563;

    /// The square of the first eccentricity of the WGS84 ellipsoid.
    pub const WGS84_E_SQ: f64 = (2.0 * Self::WGS84_F) - (Self::WGS84_F * Self::WGS84_F);

    /// Creates a new `M_PNRW` transformation engine with default constants and IAU 2000B terms.
    pub const fn new() -> Self {
        Self {
            rad_deg: PI / 180.0,
            mas_to_rad: (PI / 180.0) / 3_600_000.0,
            nutation_terms: [
                NutationTerm {
                    l: 0,
                    lp: 0,
                    f: 0,
                    d: 0,
                    om: 1,
                    psi: -17206,
                    eps: 9205,
                },
                NutationTerm {
                    l: 0,
                    lp: 0,
                    f: 2,
                    d: -2,
                    om: 2,
                    psi: -1317,
                    eps: 573,
                },
                NutationTerm {
                    l: 0,
                    lp: 0,
                    f: 2,
                    d: 0,
                    om: 2,
                    psi: -227,
                    eps: 97,
                },
                NutationTerm {
                    l: 0,
                    lp: 1,
                    f: 0,
                    d: 0,
                    om: 0,
                    psi: 147,
                    eps: 7,
                },
            ],
        }
    }

    // --- CONVERSION & UTILITY METHODS ---

    /// Converts the 40-byte "Sealed Envelope" (5 big-endian signed 64-bit integers) back into values.
    pub fn depack_vector(&self, uint8_buffer: &[u8]) -> Result<DepackedVector, MpnrwError> {
        if uint8_buffer.len() < 40 {
            return Err(MpnrwError::BufferTooShort {
                expected: 40,
                actual: uint8_buffer.len(),
            });
        }

        let lat = i64::from_be_bytes(uint8_buffer[0..8].try_into().unwrap());
        let lon = i64::from_be_bytes(uint8_buffer[8..16].try_into().unwrap());
        let height = i64::from_be_bytes(uint8_buffer[16..24].try_into().unwrap());
        let date = i64::from_be_bytes(uint8_buffer[24..32].try_into().unwrap());
        let sun = i64::from_be_bytes(uint8_buffer[32..40].try_into().unwrap());

        Ok(DepackedVector {
            lat,
            lon,
            height,
            date,
            sun,
        })
    }

    /// Calculates the Julian Date from a [`UtcDateTime`].
    pub fn get_julian_date(&self, date: UtcDateTime) -> Result<f64, MpnrwError> {
        date.validate()?;

        let year = date.year as i64;
        let month = date.month as i64;
        let day = date.day as f64;
        let hours =
            (date.hour as f64) + (date.minute as f64 / 60.0) + (date.second as f64 / 3600.0);

        let term1 = (7 * (year + (month + 9).div_euclid(12))).div_euclid(4);
        let term2 = (275 * month).div_euclid(9);

        let jd = (367 * year - term1 + term2) as f64 + day + 1_721_013.5 + (hours / 24.0);

        Ok(jd)
    }

    /// Converts geodetic coordinates on the WGS84 ellipsoid to an ECEF 3D vector `[x, y, z]` in meters.
    pub fn to_ecef(&self, lat: f64, lon: f64, height: f64) -> Vec3 {
        self.to_ecef_with_ellipsoid(lat, lon, height, Self::WGS84_A, Self::WGS84_E_SQ)
    }

    /// Converts geodetic coordinates on a custom ellipsoid to an ECEF 3D vector `[x, y, z]` in meters.
    pub fn to_ecef_with_ellipsoid(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
        a: f64,
        e_sq: f64,
    ) -> Vec3 {
        let rad_lat = lat * self.rad_deg;
        let rad_lon = lon * self.rad_deg;
        let sin_lat = rad_lat.sin();
        let cos_lat = rad_lat.cos();
        let n_phi = a / (1.0 - e_sq * sin_lat.powi(2)).sqrt();

        let x = (n_phi + height) * cos_lat * rad_lon.cos();
        let y = (n_phi + height) * cos_lat * rad_lon.sin();
        let z = (n_phi * (1.0 - e_sq) + height) * sin_lat;

        [x, y, z]
    }

    /// Calculates the Euclidean distance between two 3D vectors.
    pub fn get_distance(&self, v1: Vec3, v2: Vec3) -> f64 {
        let dx = v1[0] - v2[0];
        let dy = v1[1] - v2[1];
        let dz = v1[2] - v2[2];

        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Calculates the position of the Sun in ICRF equatorial coordinates for a given Julian Date.
    /// Returns `[x, y, z]` in kilometers (matching `src.mjs/M_PNRW.mjs`).
    pub fn get_sun_vector(&self, jd: f64) -> Vec3 {
        // t = Julian centuries from J2000.0
        let t = (jd - 2_451_545.0) / 36_525.0;
        // l = mean longitude of the Sun
        let l = (280.460 + 36_000.771 * t) % 360.0;
        // g = mean anomaly of the Sun
        let g = (357.528 + 35_999.050 * t) * self.rad_deg;
        // lambda = ecliptic longitude of the Sun
        let lambda = (l + 1.915 * g.sin() + 0.020 * (2.0 * g).sin()) * self.rad_deg;
        // r_au = distance from the Sun to the Earth in astronomical units
        let r_au = 1.00014 - 0.01671 * g.cos() - 0.00014 * (2.0 * g).cos();
        // eps = obliquity of the ecliptic
        let eps = (23.439 - 0.013 * t) * self.rad_deg;
        // au_km = astronomical unit in kilometers
        let au_km = 149_597_870.7;

        // Rotate ecliptic -> ICRF (equatorial) so frame matches get_universal_address
        [
            r_au * lambda.cos() * au_km,
            r_au * lambda.sin() * eps.cos() * au_km,
            r_au * lambda.sin() * eps.sin() * au_km,
        ]
    }

    // --- MATRIX MATH ---

    /// Multiplies two 3x3 matrices.
    pub fn multiply_matrices(&self, a: Mat3, b: Mat3) -> Mat3 {
        let mut res = [[0.0; 3]; 3];

        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    res[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        res
    }

    /// Multiplies a 3x3 matrix by a 3D vector.
    pub fn multiply_matrix_vector(&self, m: Mat3, v: Vec3) -> Vec3 {
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    }

    // --- PNRW COMPONENT GENERATORS ---

    /// Generates a rotation matrix around the y-axis for `angle` in radians.
    pub fn rot_y(&self, angle: f64) -> Mat3 {
        let (s, c) = angle.sin_cos();

        [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]]
    }

    /// Generates a rotation matrix around the z-axis for `angle` in radians.
    pub fn rot_z(&self, angle: f64) -> Mat3 {
        let (s, c) = angle.sin_cos();

        [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
    }

    /// Generates the Precession matrix (P) using IAU 1976 precession angles.
    /// `t` is time in Julian centuries from J2000.0.
    pub fn get_p_matrix(&self, t: f64) -> Mat3 {
        let arcsec_to_rad = self.rad_deg / 3600.0;
        let zeta = (2306.2181 * t) * arcsec_to_rad;
        let z = (2306.2181 * t) * arcsec_to_rad;
        let theta = (2004.3109 * t) * arcsec_to_rad;

        // P = Rz(-z) * Ry(theta) * Rz(-zeta)
        let rz_zeta = self.rot_z(-zeta);
        let ry_theta = self.rot_y(theta);
        let rz_z = self.rot_z(-z);

        self.multiply_matrices(rz_z, self.multiply_matrices(ry_theta, rz_zeta))
    }

    /// Generates the Nutation matrix (N) using simplified IAU 2000B nutation terms.
    /// `t` is time in Julian centuries from J2000.0.
    pub fn get_n_matrix(&self, t: f64) -> Mat3 {
        // The 5 fundamental Delaunay arguments
        let l = (134.963 + 477_198.868 * t) * self.rad_deg;
        let lp = (357.529 + 35_999.050 * t) * self.rad_deg;
        let f = (93.272 + 483_202.018 * t) * self.rad_deg;
        let d = (297.850 + 445_267.111 * t) * self.rad_deg;
        let om = (125.045 - 1_934.136 * t) * self.rad_deg;

        let mut d_psi = 0.0;
        let mut d_eps = 0.0;
        let eps = (23.43929 - 0.013 * t) * self.rad_deg;

        for term in &self.nutation_terms {
            let arg = (term.l as f64 * l)
                + (term.lp as f64 * lp)
                + (term.f as f64 * f)
                + (term.d as f64 * d)
                + (term.om as f64 * om);
            d_psi += (term.psi as f64) * arg.sin();
            d_eps += (term.eps as f64) * arg.cos();
        }

        let dp = d_psi * self.mas_to_rad;
        let de = d_eps * self.mas_to_rad;
        let (sin_eps, cos_eps) = eps.sin_cos();

        [
            [1.0, -dp * cos_eps, -dp * sin_eps],
            [dp * cos_eps, 1.0, -de],
            [dp * sin_eps, de, 1.0],
        ]
    }

    /// Generates the Rotation matrix (R) using Greenwich Mean Sidereal Time (GMST).
    pub fn get_r_matrix(&self, jd: f64) -> Mat3 {
        let days_since_j2000 = jd - 2_451_545.0;
        let mut gmst_hours = (18.697374558 + 24.06570982441908 * days_since_j2000) % 24.0;
        if gmst_hours < 0.0 {
            gmst_hours += 24.0;
        }

        let theta = gmst_hours * 15.0 * self.rad_deg;
        let (s, c) = theta.sin_cos();

        [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
    }

    /// Generates the Wobble matrix (W) using polar motion coordinates (`xp`, `yp` in radians).
    pub fn get_w_matrix(&self, xp: f64, yp: f64) -> Mat3 {
        [[1.0, 0.0, xp], [0.0, 1.0, -yp], [-xp, yp, 1.0]]
    }

    // --- MAIN ENGINE ---

    /// Calculates the Universal Address `[x, y, z]` for a given location, UTC time, and sun vector.
    ///
    /// - `lat`: Latitude in degrees (`-90.0..=90.0`, WGS84)
    /// - `lon`: Longitude in degrees (`-180.0..=180.0`, WGS84)
    /// - `height`: Height above ellipsoid in meters (`>= 0.0`)
    /// - `date`: UTC date and time for the transformation
    /// - `sun_vector`: Position vector of the Sun or Barycenter in ICRS `[x, y, z]`
    pub fn get_universal_address(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
        date: UtcDateTime,
        sun_vector: Vec3,
    ) -> Result<Vec3, MpnrwError> {
        if !lat.is_finite() || !(-90.0..=90.0).contains(&lat) {
            return Err(MpnrwError::InvalidLatitude);
        }
        if !lon.is_finite() || !(-180.0..=180.0).contains(&lon) {
            return Err(MpnrwError::InvalidLongitude);
        }
        if !height.is_finite() || height < 0.0 {
            return Err(MpnrwError::InvalidHeight);
        }
        if !sun_vector.iter().all(|v| v.is_finite()) {
            return Err(MpnrwError::InvalidSunVector);
        }

        let jd = self.get_julian_date(date)?;
        // t = Julian centuries from J2000.0
        let t = (jd - 2_451_545.0) / 36_525.0;
        // local_v = ECEF coordinates of the location
        let local_v = self.to_ecef(lat, lon, height);

        let p = self.get_p_matrix(t);
        let n = self.get_n_matrix(t);
        let r = self.get_r_matrix(jd);
        let w = self.get_w_matrix(0.0, 0.0);

        // m_pn = P * N
        let m_pn = self.multiply_matrices(p, n);
        // m_pnr = P * N * R
        let m_pnr = self.multiply_matrices(m_pn, r);
        // master_m = P * N * R * W
        let master_m = self.multiply_matrices(m_pnr, w);
        // space_v = master_m * local_v
        let space_v = self.multiply_matrix_vector(master_m, local_v);

        // universal_address = space_v + sun_vector
        Ok([
            space_v[0] + sun_vector[0],
            space_v[1] + sun_vector[1],
            space_v[2] + sun_vector[2],
        ])
    }
}
