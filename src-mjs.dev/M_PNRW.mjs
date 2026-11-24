// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Matthew & Contributors — Free Software (GNU GPLv3)

/**
 * M_PNRW: Celestial Transformation Engine (IEEE-754 float64 reference)
 *
 * Logic: M = P * N * R * W
 *
 * Exposes Dual Surface API: stateless topic-prefixed functional primitives
 * (`mpnrw_*`) and thin `M_PNRW` class wrapper.
 */

// 1p
// (none)

// 2p
// (none)

// 3p
// (none - zero external npm dependencies in resident core)

function invariant(cond, msg) {
    if (!cond) {
        throw new Error(msg)
    }
}

const WGS84_A= 6378137
const WGS84_F= 1 / 298.257223563
const WGS84_E_SQ= (2 * WGS84_F) - (WGS84_F ** 2)

export const MPNRW_CONSTANTS= Object.freeze({
    WGS84_A: WGS84_A,
    WGS84_F: WGS84_F,
    WGS84_E_SQ: WGS84_E_SQ
})

export const MPNRW_NUTATION_TERMS= Object.freeze([
    Object.freeze({ l: 0, lp: 0, f: 0, d: 0, om: 1, psi: -17206, eps: 9205 }),
    Object.freeze({ l: 0, lp: 0, f: 2, d: -2, om: 2, psi: -1317, eps: 573 }),
    Object.freeze({ l: 0, lp: 0, f: 2, d: 0, om: 2, psi: -227, eps: 97 }),
    Object.freeze({ l: 0, lp: 1, f: 0, d: 0, om: 0, psi: 147, eps: 7 })
])

export function mpnrw_create_ctx() {
    return Object.freeze({
        rad_deg: Math.PI / 180,
        mas_to_rad: (Math.PI / 180) / 3600000,
        nutation_terms: MPNRW_NUTATION_TERMS
    })
}

export function mpnrw_depack_vector(uint8_buffer) {
    invariant(uint8_buffer instanceof Uint8Array && uint8_buffer.byteLength >= 40, 'mpnrw_depack_vector: 40-byte Uint8Array required')
    var view= new DataView(uint8_buffer.buffer, uint8_buffer.byteOffset, uint8_buffer.byteLength)
    return {
        lat_bi: view.getBigInt64(0),
        lon_bi: view.getBigInt64(8),
        height_bi: view.getBigInt64(16),
        date_bi: view.getBigInt64(24),
        sun_bi: view.getBigInt64(32)
    }
}

export function mpnrw_julian_date(date) {
    invariant(date instanceof Date, 'date must be a Date object')
    invariant(!isNaN(date.getTime()), 'Invalid date')

    return 2440587.5 + date.getTime() / 86400000
}

export function mpnrw_to_ecef(ctx, lat, lon, height, a= MPNRW_CONSTANTS.WGS84_A, e_sq= MPNRW_CONSTANTS.WGS84_E_SQ) {
    var rad_lat= lat * ctx.rad_deg
    var rad_lon= lon * ctx.rad_deg
    var n_phi= a / Math.sqrt(1 - e_sq * Math.pow(Math.sin(rad_lat), 2))

    var x= (n_phi + height) * Math.cos(rad_lat) * Math.cos(rad_lon)
    var y= (n_phi + height) * Math.cos(rad_lat) * Math.sin(rad_lon)
    var z= (n_phi * (1 - e_sq) + height) * Math.sin(rad_lat)

    return [x, y, z]
}

export function mpnrw_distance(v1, v2) {
    var dx= v1[0] - v2[0]
    var dy= v1[1] - v2[1]
    var dz= v1[2] - v2[2]
    return Math.sqrt((dx * dx) + (dy * dy) + (dz * dz))
}

export function mpnrw_sun_vector(ctx, jd) {
    var t= (jd - 2451545.0) / 36525
    var l= (280.460 + 36000.771 * t) % 360
    var g= (357.528 + 35999.050 * t) * ctx.rad_deg
    var lambda= (l + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * ctx.rad_deg
    var r_au= 1.00014 - 0.01671 * Math.cos(g) - 0.00014 * Math.cos(2 * g)
    var eps= (23.439 - 0.013 * t) * ctx.rad_deg
    var au_km= 149597870.7

    return [
        r_au * Math.cos(lambda) * au_km,
        r_au * Math.sin(lambda) * Math.cos(eps) * au_km,
        r_au * Math.sin(lambda) * Math.sin(eps) * au_km
    ]
}

export function mpnrw_multiply_matrices(a, b) {
    var res= [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for (let i= 0; i < 3; i= i + 1) {
        for (let j= 0; j < 3; j= j + 1) {
            var sum= 0
            for (let k= 0; k < 3; k= k + 1) {
                sum= sum + (a[i][k] * b[k][j])
            }
            res[i][j]= sum
        }
    }
    return res
}

export function mpnrw_multiply_matrix_vector(m, v) {
    return [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
    ]
}

export function mpnrw_transpose_matrix(m) {
    return [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]]
    ]
}

export function mpnrw_invert_matrix(m) {
    var c00= (m[1][1] * m[2][2]) - (m[1][2] * m[2][1])
    var c01= (m[1][2] * m[2][0]) - (m[1][0] * m[2][2])
    var c02= (m[1][0] * m[2][1]) - (m[1][1] * m[2][0])
    var det= (m[0][0] * c00) + (m[0][1] * c01) + (m[0][2] * c02)
    invariant(det !== 0, 'mpnrw_invert_matrix: singular matrix')

    return [
        [c00 / det, ((m[0][2] * m[2][1]) - (m[0][1] * m[2][2])) / det, ((m[0][1] * m[1][2]) - (m[0][2] * m[1][1])) / det],
        [c01 / det, ((m[0][0] * m[2][2]) - (m[0][2] * m[2][0])) / det, ((m[0][2] * m[1][0]) - (m[0][0] * m[1][2])) / det],
        [c02 / det, ((m[0][1] * m[2][0]) - (m[0][0] * m[2][1])) / det, ((m[0][0] * m[1][1]) - (m[0][1] * m[1][0])) / det]
    ]
}

export function mpnrw_from_ecef(ctx, x, y, z, a= MPNRW_CONSTANTS.WGS84_A, e_sq= MPNRW_CONSTANTS.WGS84_E_SQ) {
    var lon_rad= Math.atan2(y, x)
    var p= Math.hypot(x, y)

    if (p === 0) {
        var lat_deg= (z >= 0)
            ? 90
            : -90
        var b= a * Math.sqrt(1 - e_sq)
        return { lat: lat_deg, lon: 0, height: Math.abs(z) - b }
    }

    var lat_rad= Math.atan2(z, p * (1 - e_sq))
    for (let i= 0; i < 6; i= i + 1) {
        var sin_lat= Math.sin(lat_rad)
        var n_phi= a / Math.sqrt(1 - e_sq * sin_lat * sin_lat)
        lat_rad= Math.atan2(z + e_sq * n_phi * sin_lat, p)
    }

    var sin_f= Math.sin(lat_rad)
    var cos_f= Math.cos(lat_rad)
    var n_phi_f= a / Math.sqrt(1 - e_sq * sin_f * sin_f)
    var height= (Math.abs(cos_f) >= Math.abs(sin_f))
        ? (p / cos_f) - n_phi_f
        : (z / sin_f) - (n_phi_f * (1 - e_sq))

    return {
        lat: lat_rad / ctx.rad_deg,
        lon: lon_rad / ctx.rad_deg,
        height: height
    }
}

export function mpnrw_rot_y(angle) {
    var c= Math.cos(angle)
    var s= Math.sin(angle)
    return [[c, 0, -s], [0, 1, 0], [s, 0, c]]
}

export function mpnrw_rot_z(angle) {
    var c= Math.cos(angle)
    var s= Math.sin(angle)
    return [[c, s, 0], [-s, c, 0], [0, 0, 1]]
}

export function mpnrw_p_matrix(ctx, t) {
    var zeta= (2306.2181 * t) * (ctx.rad_deg / 3600)
    var z= (2306.2181 * t) * (ctx.rad_deg / 3600)
    var theta= (2004.3109 * t) * (ctx.rad_deg / 3600)

    var rz_zeta= mpnrw_rot_z(zeta)
    var ry_theta= mpnrw_rot_y(-theta)
    var rz_z= mpnrw_rot_z(z)

    return mpnrw_multiply_matrices(rz_z, mpnrw_multiply_matrices(ry_theta, rz_zeta))
}

export function mpnrw_n_matrix(ctx, t) {
    var l= (134.963 + 477198.868 * t) * ctx.rad_deg
    var lp= (357.529 + 35999.050 * t) * ctx.rad_deg
    var f= (93.272 + 483202.018 * t) * ctx.rad_deg
    var d= (297.850 + 445267.111 * t) * ctx.rad_deg
    var om= (125.045 - 1934.136 * t) * ctx.rad_deg

    var d_psi= 0
    var d_eps= 0
    var eps= (23.43929 - 0.013 * t) * ctx.rad_deg

    for (let i= 0; i < ctx.nutation_terms.length; i= i + 1) {
        var term= ctx.nutation_terms[i]
        var arg= (term.l * l) + (term.lp * lp) + (term.f * f) + (term.d * d) + (term.om * om)
        d_psi= d_psi + term.psi * Math.sin(arg)
        d_eps= d_eps + term.eps * Math.cos(arg)
    }

    var dp= d_psi * ctx.mas_to_rad
    var de= d_eps * ctx.mas_to_rad

    return [
        [1, 0, dp * Math.sin(eps)],
        [0, 1, de],
        [-dp * Math.sin(eps), -de, 1]
    ]
}

export function mpnrw_r_matrix(ctx, jd) {
    var days_since_j2000= jd - 2451545.0
    var gmst_hours= (18.697374558 + 24.06570982441908 * days_since_j2000) % 24
    if (gmst_hours < 0) gmst_hours= gmst_hours + 24

    var theta= gmst_hours * 15 * ctx.rad_deg
    return [
        [Math.cos(theta), -Math.sin(theta), 0],
        [Math.sin(theta), Math.cos(theta), 0],
        [0, 0, 1]
    ]
}

export function mpnrw_w_matrix(xp= 0, yp= 0) {
    return [
        [1, 0, -xp],
        [0, 1, yp],
        [xp, -yp, 1]
    ]
}

export function mpnrw_heliocentric_earth_m(ctx, jd_tt) {
    var t= (jd_tt - 2451545) / 36525
    var sun= mpnrw_sun_vector(ctx, jd_tt)
    return mpnrw_multiply_matrix_vector(mpnrw_p_matrix(ctx, t), sun.map(v => v * -1000))
}

export function mpnrw_universal_address(ctx, lat, lon, height, date, origin_to_earth_m, orientation= {}) {
    invariant(typeof lat === 'number' && lat >= -90 && lat <= 90, 'Latitude must be a number between -90 and 90')
    invariant(typeof lon === 'number' && lon >= -180 && lon <= 180, 'Longitude must be a number between -180 and 180')
    invariant(Number.isFinite(height) && height >= 0, 'Height must be a number greater than or equal to 0')
    invariant(date instanceof Date, 'Date must be a Date object')
    invariant(Array.isArray(origin_to_earth_m) && origin_to_earth_m.length === 3 && Array.from(origin_to_earth_m).every(Number.isFinite), 'Origin offset must be an array of 3 numbers')

    var { dut1= 0, tt_minus_utc= 0, xp= 0, yp= 0 }= orientation
    invariant([dut1, tt_minus_utc, xp, yp].every(Number.isFinite), 'Invalid Earth orientation')
    var jd= mpnrw_julian_date(date)
    var jd_tt= jd + tt_minus_utc / 86400
    var t= (jd_tt - 2451545) / 36525
    var local_v= mpnrw_to_ecef(ctx, lat, lon, height)

    var p= mpnrw_p_matrix(ctx, t)
    var n= mpnrw_n_matrix(ctx, t)
    var r= mpnrw_r_matrix(ctx, jd + dut1 / 86400)
    var w= mpnrw_w_matrix(xp, yp)

    var m_pn= mpnrw_multiply_matrices(p, n)
    var m_pnr= mpnrw_multiply_matrices(m_pn, r)
    var master_m= mpnrw_multiply_matrices(m_pnr, w)
    var space_v= mpnrw_multiply_matrix_vector(master_m, local_v)

    return [
        space_v[0] + origin_to_earth_m[0],
        space_v[1] + origin_to_earth_m[1],
        space_v[2] + origin_to_earth_m[2]
    ]
}

export function mpnrw_decode_universal_address(ctx, address, date, origin_to_earth_m= [0, 0, 0], orientation= {}) {
    invariant(Array.isArray(address) && address.length === 3 && Array.from(address).every(Number.isFinite), 'Address must be an array of 3 numbers')
    invariant(date instanceof Date, 'Date must be a Date object')
    invariant(Array.isArray(origin_to_earth_m) && origin_to_earth_m.length === 3 && Array.from(origin_to_earth_m).every(Number.isFinite), 'Origin offset must be an array of 3 numbers')

    var { dut1= 0, tt_minus_utc= 0, xp= 0, yp= 0 }= orientation
    invariant([dut1, tt_minus_utc, xp, yp].every(Number.isFinite), 'Invalid Earth orientation')
    var jd= mpnrw_julian_date(date)
    var jd_tt= jd + tt_minus_utc / 86400
    var t= (jd_tt - 2451545) / 36525

    var p= mpnrw_p_matrix(ctx, t)
    var n= mpnrw_n_matrix(ctx, t)
    var r= mpnrw_r_matrix(ctx, jd + dut1 / 86400)
    var w= mpnrw_w_matrix(xp, yp)

    var m_pn= mpnrw_multiply_matrices(p, n)
    var m_pnr= mpnrw_multiply_matrices(m_pn, r)
    var master_m= mpnrw_multiply_matrices(m_pnr, w)
    var inv_m= mpnrw_invert_matrix(master_m)

    var space_v= [
        address[0] - origin_to_earth_m[0],
        address[1] - origin_to_earth_m[1],
        address[2] - origin_to_earth_m[2]
    ]
    var local_v= mpnrw_multiply_matrix_vector(inv_m, space_v)
    return mpnrw_from_ecef(ctx, local_v[0], local_v[1], local_v[2])
}

// --- THIN OOP CLASS WRAPPER ---

export class M_PNRW {
    static WGS84_A= MPNRW_CONSTANTS.WGS84_A
    static WGS84_F= MPNRW_CONSTANTS.WGS84_F
    static WGS84_E_SQ= MPNRW_CONSTANTS.WGS84_E_SQ

    #ctx

    constructor() {
        this.#ctx= mpnrw_create_ctx()
        this.rad_deg= this.#ctx.rad_deg
        this.mas_to_rad= this.#ctx.mas_to_rad
        this.nutation_terms= this.#ctx.nutation_terms
    }

    depack_vector(uint8_buffer) {
        return mpnrw_depack_vector(uint8_buffer)
    }

    get_julian_date(date) {
        return mpnrw_julian_date(date)
    }

    to_ecef(lat, lon, height, a= MPNRW_CONSTANTS.WGS84_A, e_sq= MPNRW_CONSTANTS.WGS84_E_SQ) {
        return mpnrw_to_ecef(this.#ctx, lat, lon, height, a, e_sq)
    }

    from_ecef(x, y, z, a= MPNRW_CONSTANTS.WGS84_A, e_sq= MPNRW_CONSTANTS.WGS84_E_SQ) {
        return mpnrw_from_ecef(this.#ctx, x, y, z, a, e_sq)
    }

    transpose_matrix(m) {
        return mpnrw_transpose_matrix(m)
    }

    get_distance(v1, v2) {
        return mpnrw_distance(v1, v2)
    }

    get_heliocentric_earth_m(jd_tt) {
        return mpnrw_heliocentric_earth_m(this.#ctx, jd_tt)
    }

    get_sun_vector(jd) {
        return mpnrw_sun_vector(this.#ctx, jd)
    }

    multiply_matrices(a, b) {
        return mpnrw_multiply_matrices(a, b)
    }

    multiply_matrix_vector(m, v) {
        return mpnrw_multiply_matrix_vector(m, v)
    }

    rot_y(angle) {
        return mpnrw_rot_y(angle)
    }

    rot_z(angle) {
        return mpnrw_rot_z(angle)
    }

    get_p_matrix(t) {
        return mpnrw_p_matrix(this.#ctx, t)
    }

    get_n_matrix(t) {
        return mpnrw_n_matrix(this.#ctx, t)
    }

    get_r_matrix(jd) {
        return mpnrw_r_matrix(this.#ctx, jd)
    }

    get_w_matrix(xp= 0, yp= 0) {
        return mpnrw_w_matrix(xp, yp)
    }

    get_universal_address(lat, lon, height, date, origin_to_earth_m, orientation= {}) {
        return mpnrw_universal_address(this.#ctx, lat, lon, height, date, origin_to_earth_m, orientation)
    }

    decode_universal_address(address, date, origin_to_earth_m= [0, 0, 0], orientation= {}) {
        return mpnrw_decode_universal_address(this.#ctx, address, date, origin_to_earth_m, orientation)
    }
}
