// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Matthew & Contributors — Free Software (GNU GPLv3)

/**
 * M_PNRW_BigInt: Deterministic Fixed-Point Celestial Transformation Engine
 *
 * Logic: M = P * N * R * W
 *
 * Arithmetic: Q.80 signed binary fixed-point using native BigInt.
 * Every operation (+, -, *, /, sqrt, sin, cos) is pure integer arithmetic
 * with zero IEEE-754 floating-point or platform libm dependencies, ensuring
 * bit-for-bit identical results across all CPU architectures and JS engines.
 *
 * Exposes Dual Surface API: stateless topic-prefixed functional primitives
 * (`mpnrw_fp_*`, `mpnrw_bigint_*`) and thin `M_PNRW_BigInt` class wrapper.
 */

// 1p
// (none)

// 2p
// (none)

// 3p
// (none - zero external npm dependencies in resident core)

const SCALE_BITS= 80n
const SCALE= 1n << SCALE_BITS
const HALF_SCALE= 1n << (SCALE_BITS - 1n)

// Internal guard bits for trigonometric Taylor evaluation
const TRIG_GUARD_BITS= 16n
const TRIG_SHIFT= SCALE_BITS + TRIG_GUARD_BITS // 96n
const TRIG_HALF_GUARD= 1n << (TRIG_GUARD_BITS - 1n)

// Pi at 160 fractional bits (for high-precision range reduction of large angles)
const PI_160= 0x3243f6a8885a308d313198a2e03707344a4093822n
const TWO_PI_160= PI_160 << 1n

// Pi constants at Q.80
const PI_80= (PI_160 + (1n << 79n)) >> 80n
const HALF_PI_80= (PI_160 + (1n << 80n)) >> 81n

const MS_PER_DAY= 86_400_000n
// Julian Date at Unix epoch (1970-01-01T00:00:00.000Z) = 2440587.5
const JD_UNIX_EPOCH_80= (2440587n << SCALE_BITS) + HALF_SCALE
// J2000.0 = 2451545.0
const JD_J2000_80= 2451545n << SCALE_BITS

function invariant(cond, msg) {
    if (!cond) {
        throw new Error(msg)
    }
}

function parse_decimal_to_q80(str) {
    var trimmed= str.trim()
    var match= /^([+-])?(?:(\d+)(?:\.(\d*))?|\.(\d+))(?:[eE]([+-]?\d+))?$/.exec(trimmed)
    invariant(match !== null, `Invalid decimal string: "${str}"`)

    var neg= match[1] === '-'
    var int_part= match[2] ?? ''
    var frac_part= match[3] ?? match[4] ?? ''
    var exp= match[5] ? parseInt(match[5], 10) : 0

    var digits= BigInt((int_part + frac_part) || '0')
    if (digits === 0n) return 0n

    var pow10= exp - frac_part.length
    var num= digits << SCALE_BITS
    var den= 1n

    if (pow10 >= 0) {
        num= num * (10n ** BigInt(pow10))
    } else {
        den= 10n ** BigInt(-pow10)
    }

    var q= (num + (den >> 1n)) / den
    return (neg)
        ? -q
        : q
}

// --- FROZEN CONSTANTS (Q.80) ---

const WGS84_A_80= 6378137n << SCALE_BITS
const WGS84_F_80= 0xdbbad4279f09443c5en

export const MPNRW_BIGINT_CONSTANTS= Object.freeze({
    SCALE_BITS: SCALE_BITS,
    SCALE: SCALE,
    PI: PI_80,
    HALF_PI: HALF_PI_80,
    JD_J2000: JD_J2000_80,
    WGS84_A: WGS84_A_80,
    WGS84_F: WGS84_F_80,
    WGS84_E_SQ: (2n * WGS84_F_80) - mpnrw_fp_mul(WGS84_F_80, WGS84_F_80)
})

export const MPNRW_NUTATION_TERMS_BIGINT= Object.freeze([
    Object.freeze({ l: 0n, lp: 0n, f: 0n, d: 0n, om: 1n, psi: -17206n << SCALE_BITS, eps: 9205n << SCALE_BITS }),
    Object.freeze({ l: 0n, lp: 0n, f: 2n, d: -2n, om: 2n, psi: -1317n << SCALE_BITS, eps: 573n << SCALE_BITS }),
    Object.freeze({ l: 0n, lp: 0n, f: 2n, d: 0n, om: 2n, psi: -227n << SCALE_BITS, eps: 97n << SCALE_BITS }),
    Object.freeze({ l: 0n, lp: 1n, f: 0n, d: 0n, om: 0n, psi: 147n << SCALE_BITS, eps: 7n << SCALE_BITS })
])

// --- FIXED-POINT CONVERSION PRIMITIVES ---

export function mpnrw_fp_from_int(n) {
    return BigInt(n) << SCALE_BITS
}

export function mpnrw_fp_from_decimal(str) {
    invariant(typeof str === 'string', 'mpnrw_fp_from_decimal: string required')
    return parse_decimal_to_q80(str)
}

export function mpnrw_fp_from_number(n) {
    invariant(typeof n === 'number' && Number.isFinite(n), 'mpnrw_fp_from_number: finite number required')
    return parse_decimal_to_q80(n.toString())
}

export function mpnrw_fp_from_rational(num, den) {
    var a= BigInt(num)
    var b= BigInt(den)
    invariant(b !== 0n, 'mpnrw_fp_from_rational: division by zero')
    var neg= (a < 0n) !== (b < 0n)
    var abs_a= (a < 0n)
        ? -a
        : a
    var abs_b= (b < 0n)
        ? -b
        : b
    var q= ((abs_a << SCALE_BITS) + (abs_b >> 1n)) / abs_b
    return (neg)
        ? -q
        : q
}

export function mpnrw_fp_to_decimal(fp, decimals= 12) {
    invariant(typeof fp === 'bigint', 'mpnrw_fp_to_decimal: bigint required')
    invariant(Number.isInteger(decimals) && decimals >= 0, 'mpnrw_fp_to_decimal: non-negative integer decimals required')

    var neg= fp < 0n
    var abs= (neg)
        ? -fp
        : fp
    var factor= 10n ** BigInt(decimals)
    var scaled= ((abs * factor) + HALF_SCALE) >> SCALE_BITS

    var int_part= scaled / factor
    if (decimals === 0) {
        var zero_sign= (neg && int_part !== 0n)
            ? '-'
            : ''
        return zero_sign + int_part.toString()
    }

    var frac_part= (scaled % factor).toString().padStart(decimals, '0')
    var sign= (neg && scaled !== 0n)
        ? '-'
        : ''
    return `${sign}${int_part.toString()}.${frac_part}`
}

export function mpnrw_fp_to_number(fp) {
    return Number(mpnrw_fp_to_decimal(fp, 17))
}

// --- FIXED-POINT ARITHMETIC PRIMITIVES (Q.80) ---

export function mpnrw_fp_add(a, b) {
    return a + b
}

export function mpnrw_fp_sub(a, b) {
    return a - b
}

export function mpnrw_fp_mul(a, b) {
    var prod= a * b
    return (prod >= 0n)
        ? (prod + HALF_SCALE) >> SCALE_BITS
        : -((-prod + HALF_SCALE) >> SCALE_BITS)
}

export function mpnrw_fp_div(a, b) {
    invariant(b !== 0n, 'mpnrw_fp_div: division by zero')
    var neg= (a < 0n) !== (b < 0n)
    var abs_a= (a < 0n)
        ? -a
        : a
    var abs_b= (b < 0n)
        ? -b
        : b
    var q= ((abs_a << SCALE_BITS) + (abs_b >> 1n)) / abs_b
    return (neg)
        ? -q
        : q
}

export function mpnrw_fp_mod(a, b) {
    invariant(b !== 0n, 'mpnrw_fp_mod: modulo by zero')
    return a % b
}

export function mpnrw_fp_sqrt(a) {
    invariant(a >= 0n, 'mpnrw_fp_sqrt: negative input')
    if (a === 0n) return 0n

    var n= a << SCALE_BITS
    var bit_len= n.toString(2).length
    var y= 1n << BigInt((bit_len + 1) >> 1)

    while (true) {
        var next= (y + n / y) >> 1n
        if (next >= y) break
        y= next
    }
    return y
}

/**
 * Computes [sin(angle), cos(angle)] for a Q.80 angle in radians.
 * Uses 160-bit Pi range reduction and 96-bit internal Taylor evaluation.
 */
export function mpnrw_fp_sin_cos(angle) {
    var a160= angle << SCALE_BITS
    var half_two_pi= TWO_PI_160 >> 1n
    var k= (a160 >= 0n)
        ? (a160 + half_two_pi) / TWO_PI_160
        : -((-a160 + half_two_pi) / TWO_PI_160)
    var r= (a160 - k * TWO_PI_160) >> SCALE_BITS

    var half_quarter= HALF_PI_80 >> 1n
    var q= (r >= 0n)
        ? (r + half_quarter) / HALF_PI_80
        : -((-r + half_quarter) / HALF_PI_80)
    var u= r - q * HALF_PI_80

    var u96= u << TRIG_GUARD_BITS
    var u2_96= (u96 * u96) >> TRIG_SHIFT

    var sin_sum= u96
    var sin_term= u96
    for (let n= 1n; n <= 12n; n= n + 1n) {
        var s_denom= (2n * n) * (2n * n + 1n)
        sin_term= -((sin_term * u2_96) >> TRIG_SHIFT) / s_denom
        if (sin_term === 0n) break
        sin_sum= sin_sum + sin_term
    }

    var cos_sum= 1n << TRIG_SHIFT
    var cos_term= 1n << TRIG_SHIFT
    for (let n= 1n; n <= 12n; n= n + 1n) {
        var c_denom= (2n * n - 1n) * (2n * n)
        cos_term= -((cos_term * u2_96) >> TRIG_SHIFT) / c_denom
        if (cos_term === 0n) break
        cos_sum= cos_sum + cos_term
    }

    var s= (sin_sum >= 0n)
        ? (sin_sum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS
        : -((-sin_sum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS)
    var c= (cos_sum >= 0n)
        ? (cos_sum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS
        : -((-cos_sum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS)

    var quad= ((q % 4n) + 4n) % 4n
    if (quad === 0n) return [s, c]
    if (quad === 1n) return [c, -s]
    if (quad === 2n) return [-s, -c]
    return [-c, s]
}

export function mpnrw_fp_sin(angle) {
    return mpnrw_fp_sin_cos(angle)[0]
}

export function mpnrw_fp_cos(angle) {
    return mpnrw_fp_sin_cos(angle)[1]
}

// --- STATELESS ENGINE PRIMITIVES ---

export function mpnrw_bigint_create_ctx() {
    var rad_deg= mpnrw_fp_div(PI_80, 180n << SCALE_BITS)
    var mas_to_rad= mpnrw_fp_div(rad_deg, 3_600_000n << SCALE_BITS)

    return Object.freeze({
        rad_deg: rad_deg,
        mas_to_rad: mas_to_rad,
        nutation_terms: MPNRW_NUTATION_TERMS_BIGINT,
        c: Object.freeze({
            days_per_century: 36525n << SCALE_BITS,
            hours_per_day: 24n << SCALE_BITS,
            deg_per_hour: 15n << SCALE_BITS,
            arcsec_to_rad: mpnrw_fp_div(rad_deg, 3600n << SCALE_BITS),

            p_zeta: parse_decimal_to_q80('2306.2181'),
            p_theta: parse_decimal_to_q80('2004.3109'),

            n_l0: parse_decimal_to_q80('134.963'),
            n_l1: parse_decimal_to_q80('477198.868'),
            n_lp0: parse_decimal_to_q80('357.529'),
            n_lp1: parse_decimal_to_q80('35999.050'),
            n_f0: parse_decimal_to_q80('93.272'),
            n_f1: parse_decimal_to_q80('483202.018'),
            n_d0: parse_decimal_to_q80('297.850'),
            n_d1: parse_decimal_to_q80('445267.111'),
            n_om0: parse_decimal_to_q80('125.045'),
            n_om1: parse_decimal_to_q80('-1934.136'),
            n_eps0: parse_decimal_to_q80('23.43929'),
            n_eps1: parse_decimal_to_q80('0.013'),

            gmst0: parse_decimal_to_q80('18.697374558'),
            gmst1: parse_decimal_to_q80('24.06570982441908'),

            s_l0: parse_decimal_to_q80('280.460'),
            s_l1: parse_decimal_to_q80('36000.771'),
            s_g0: parse_decimal_to_q80('357.528'),
            s_g1: parse_decimal_to_q80('35999.050'),
            s_lam1: parse_decimal_to_q80('1.915'),
            s_lam2: parse_decimal_to_q80('0.020'),
            s_r0: parse_decimal_to_q80('1.00014'),
            s_r1: parse_decimal_to_q80('0.01671'),
            s_r2: parse_decimal_to_q80('0.00014'),
            s_eps0: parse_decimal_to_q80('23.439'),
            s_eps1: parse_decimal_to_q80('0.013'),
            au_km: parse_decimal_to_q80('149597870.7'),
            deg_360: 360n << SCALE_BITS
        })
    })
}

export function mpnrw_bigint_depack_vector(uint8_buffer) {
    invariant(uint8_buffer instanceof Uint8Array && uint8_buffer.byteLength >= 40, 'mpnrw_bigint_depack_vector: 40-byte Uint8Array required')
    var view= new DataView(uint8_buffer.buffer, uint8_buffer.byteOffset, uint8_buffer.byteLength)
    return {
        lat_bi: view.getBigInt64(0),
        lon_bi: view.getBigInt64(8),
        height_bi: view.getBigInt64(16),
        date_bi: view.getBigInt64(24),
        sun_bi: view.getBigInt64(32)
    }
}

export function mpnrw_bigint_julian_date(date) {
    var unix_ms= 0n
    if (date instanceof Date) {
        var ms= date.getTime()
        invariant(!Number.isNaN(ms), 'mpnrw_bigint_julian_date: invalid Date')
        unix_ms= BigInt(ms)
    } else if (typeof date === 'bigint') {
        unix_ms= date
    } else {
        invariant(false, 'mpnrw_bigint_julian_date: Date or bigint UTC ms required')
    }

    var days_fp= mpnrw_fp_from_rational(unix_ms, MS_PER_DAY)
    return JD_UNIX_EPOCH_80 + days_fp
}

export function mpnrw_bigint_to_ecef(ctx, lat, lon, height, a= MPNRW_BIGINT_CONSTANTS.WGS84_A, e_sq= MPNRW_BIGINT_CONSTANTS.WGS84_E_SQ) {
    var rad_lat= mpnrw_fp_mul(lat, ctx.rad_deg)
    var rad_lon= mpnrw_fp_mul(lon, ctx.rad_deg)

    var [sin_lat, cos_lat]= mpnrw_fp_sin_cos(rad_lat)
    var [sin_lon, cos_lon]= mpnrw_fp_sin_cos(rad_lon)

    var sin2= mpnrw_fp_mul(sin_lat, sin_lat)
    var n_phi= mpnrw_fp_div(a, mpnrw_fp_sqrt(SCALE - mpnrw_fp_mul(e_sq, sin2)))
    var n_plus_h= n_phi + height

    var x= mpnrw_fp_mul(mpnrw_fp_mul(n_plus_h, cos_lat), cos_lon)
    var y= mpnrw_fp_mul(mpnrw_fp_mul(n_plus_h, cos_lat), sin_lon)
    var z= mpnrw_fp_mul(mpnrw_fp_mul(n_phi, SCALE - e_sq) + height, sin_lat)

    return [x, y, z]
}

export function mpnrw_bigint_distance(v1, v2) {
    var dx= v1[0] - v2[0]
    var dy= v1[1] - v2[1]
    var dz= v1[2] - v2[2]
    return mpnrw_fp_sqrt(mpnrw_fp_mul(dx, dx) + mpnrw_fp_mul(dy, dy) + mpnrw_fp_mul(dz, dz))
}

export function mpnrw_bigint_sun_vector(ctx, jd) {
    var c= ctx.c
    var t= mpnrw_fp_div(jd - JD_J2000_80, c.days_per_century)
    var l= (c.s_l0 + mpnrw_fp_mul(c.s_l1, t)) % c.deg_360
    var g= mpnrw_fp_mul(c.s_g0 + mpnrw_fp_mul(c.s_g1, t), ctx.rad_deg)

    var [sin_g, cos_g]= mpnrw_fp_sin_cos(g)
    var [sin_2g, cos_2g]= mpnrw_fp_sin_cos(g << 1n)

    var lambda= mpnrw_fp_mul(l + mpnrw_fp_mul(c.s_lam1, sin_g) + mpnrw_fp_mul(c.s_lam2, sin_2g), ctx.rad_deg)
    var r_au= c.s_r0 - mpnrw_fp_mul(c.s_r1, cos_g) - mpnrw_fp_mul(c.s_r2, cos_2g)
    var eps= mpnrw_fp_mul(c.s_eps0 - mpnrw_fp_mul(c.s_eps1, t), ctx.rad_deg)

    var [sin_lam, cos_lam]= mpnrw_fp_sin_cos(lambda)
    var [sin_eps, cos_eps]= mpnrw_fp_sin_cos(eps)
    var r_km= mpnrw_fp_mul(r_au, c.au_km)

    return [
        mpnrw_fp_mul(r_km, cos_lam),
        mpnrw_fp_mul(mpnrw_fp_mul(r_km, sin_lam), cos_eps),
        mpnrw_fp_mul(mpnrw_fp_mul(r_km, sin_lam), sin_eps)
    ]
}

export function mpnrw_bigint_multiply_matrices(a, b) {
    var res= [
        [0n, 0n, 0n],
        [0n, 0n, 0n],
        [0n, 0n, 0n]
    ]
    for (let i= 0; i < 3; i= i + 1) {
        for (let j= 0; j < 3; j= j + 1) {
            res[i][j]= mpnrw_fp_mul(a[i][0], b[0][j]) + mpnrw_fp_mul(a[i][1], b[1][j]) + mpnrw_fp_mul(a[i][2], b[2][j])
        }
    }
    return res
}

export function mpnrw_bigint_multiply_matrix_vector(m, v) {
    return [
        mpnrw_fp_mul(m[0][0], v[0]) + mpnrw_fp_mul(m[0][1], v[1]) + mpnrw_fp_mul(m[0][2], v[2]),
        mpnrw_fp_mul(m[1][0], v[0]) + mpnrw_fp_mul(m[1][1], v[1]) + mpnrw_fp_mul(m[1][2], v[2]),
        mpnrw_fp_mul(m[2][0], v[0]) + mpnrw_fp_mul(m[2][1], v[1]) + mpnrw_fp_mul(m[2][2], v[2])
    ]
}

export function mpnrw_bigint_rot_y(angle) {
    var [s, c]= mpnrw_fp_sin_cos(angle)
    return [
        [c, 0n, -s],
        [0n, SCALE, 0n],
        [s, 0n, c]
    ]
}

export function mpnrw_bigint_rot_z(angle) {
    var [s, c]= mpnrw_fp_sin_cos(angle)
    return [
        [c, s, 0n],
        [-s, c, 0n],
        [0n, 0n, SCALE]
    ]
}

export function mpnrw_bigint_p_matrix(ctx, t) {
    var c= ctx.c
    var zeta= mpnrw_fp_mul(mpnrw_fp_mul(c.p_zeta, t), c.arcsec_to_rad)
    var z= zeta
    var theta= mpnrw_fp_mul(mpnrw_fp_mul(c.p_theta, t), c.arcsec_to_rad)

    var rz_zeta= mpnrw_bigint_rot_z(zeta)
    var ry_theta= mpnrw_bigint_rot_y(-theta)
    var rz_z= mpnrw_bigint_rot_z(z)

    return mpnrw_bigint_multiply_matrices(rz_z, mpnrw_bigint_multiply_matrices(ry_theta, rz_zeta))
}

export function mpnrw_bigint_n_matrix(ctx, t) {
    var c= ctx.c
    var l= mpnrw_fp_mul(c.n_l0 + mpnrw_fp_mul(c.n_l1, t), ctx.rad_deg)
    var lp= mpnrw_fp_mul(c.n_lp0 + mpnrw_fp_mul(c.n_lp1, t), ctx.rad_deg)
    var f= mpnrw_fp_mul(c.n_f0 + mpnrw_fp_mul(c.n_f1, t), ctx.rad_deg)
    var d= mpnrw_fp_mul(c.n_d0 + mpnrw_fp_mul(c.n_d1, t), ctx.rad_deg)
    var om= mpnrw_fp_mul(c.n_om0 + mpnrw_fp_mul(c.n_om1, t), ctx.rad_deg)

    var d_psi= 0n
    var d_eps= 0n
    var eps= mpnrw_fp_mul(c.n_eps0 - mpnrw_fp_mul(c.n_eps1, t), ctx.rad_deg)

    for (let i= 0; i < ctx.nutation_terms.length; i= i + 1) {
        var term= ctx.nutation_terms[i]
        var arg= (term.l * l) + (term.lp * lp) + (term.f * f) + (term.d * d) + (term.om * om)
        var [s_arg, c_arg]= mpnrw_fp_sin_cos(arg)
        d_psi= d_psi + mpnrw_fp_mul(term.psi, s_arg)
        d_eps= d_eps + mpnrw_fp_mul(term.eps, c_arg)
    }

    var dp= mpnrw_fp_mul(d_psi, ctx.mas_to_rad)
    var de= mpnrw_fp_mul(d_eps, ctx.mas_to_rad)
    var [sin_eps]= mpnrw_fp_sin_cos(eps)
    var dp_sin= mpnrw_fp_mul(dp, sin_eps)

    return [
        [SCALE, 0n, dp_sin],
        [0n, SCALE, de],
        [-dp_sin, -de, SCALE]
    ]
}

export function mpnrw_bigint_r_matrix(ctx, jd) {
    var c= ctx.c
    var days_since_j2000= jd - JD_J2000_80
    var gmst_hours= (c.gmst0 + mpnrw_fp_mul(c.gmst1, days_since_j2000)) % c.hours_per_day
    if (gmst_hours < 0n) gmst_hours= gmst_hours + c.hours_per_day

    var theta= mpnrw_fp_mul(mpnrw_fp_mul(gmst_hours, c.deg_per_hour), ctx.rad_deg)
    var [s, cos_t]= mpnrw_fp_sin_cos(theta)

    return [
        [cos_t, -s, 0n],
        [s, cos_t, 0n],
        [0n, 0n, SCALE]
    ]
}

export function mpnrw_bigint_w_matrix(xp= 0n, yp= 0n) {
    return [
        [SCALE, 0n, -xp],
        [0n, SCALE, yp],
        [xp, -yp, SCALE]
    ]
}

export function mpnrw_bigint_heliocentric_earth_m(ctx, jd_tt) {
    var t= mpnrw_fp_div(jd_tt - JD_J2000_80, ctx.c.days_per_century)
    var sun= mpnrw_bigint_sun_vector(ctx, jd_tt)
    return mpnrw_bigint_multiply_matrix_vector(mpnrw_bigint_p_matrix(ctx, t), sun.map(v => v * -1000n))
}

export function mpnrw_bigint_universal_address(ctx, lat, lon, height, date, origin_to_earth_m, orientation= {}) {
    invariant(typeof lat === 'bigint' && lat >= (-90n << SCALE_BITS) && lat <= (90n << SCALE_BITS), 'Latitude must be a Q.80 bigint between -90 and 90 degrees')
    invariant(typeof lon === 'bigint' && lon >= (-180n << SCALE_BITS) && lon <= (180n << SCALE_BITS), 'Longitude must be a Q.80 bigint between -180 and 180 degrees')
    invariant(typeof height === 'bigint' && height >= 0n, 'Height must be a non-negative Q.80 bigint')
    invariant((date instanceof Date) || (typeof date === 'bigint'), 'Date must be a Date object or bigint UTC milliseconds')
    invariant(Array.isArray(origin_to_earth_m) && origin_to_earth_m.length === 3 && Array.from(origin_to_earth_m).every((v) => typeof v === 'bigint'), 'Origin offset must be an array of 3 Q.80 bigints')

    var { dut1= 0n, tt_minus_utc= 0n, xp= 0n, yp= 0n }= orientation
    invariant([dut1, tt_minus_utc, xp, yp].every(v => typeof v === 'bigint'), 'Invalid Earth orientation')
    var jd= mpnrw_bigint_julian_date(date)
    var jd_tt= jd + tt_minus_utc / 86400n
    var t= mpnrw_fp_div(jd_tt - JD_J2000_80, ctx.c.days_per_century)
    var local_v= mpnrw_bigint_to_ecef(ctx, lat, lon, height)

    var p= mpnrw_bigint_p_matrix(ctx, t)
    var n= mpnrw_bigint_n_matrix(ctx, t)
    var r= mpnrw_bigint_r_matrix(ctx, jd + dut1 / 86400n)
    var w= mpnrw_bigint_w_matrix(xp, yp)

    var m_pn= mpnrw_bigint_multiply_matrices(p, n)
    var m_pnr= mpnrw_bigint_multiply_matrices(m_pn, r)
    var master_m= mpnrw_bigint_multiply_matrices(m_pnr, w)
    var space_v= mpnrw_bigint_multiply_matrix_vector(master_m, local_v)

    return [
        space_v[0] + origin_to_earth_m[0],
        space_v[1] + origin_to_earth_m[1],
        space_v[2] + origin_to_earth_m[2]
    ]
}

// --- THIN OOP CLASS WRAPPER ---

export class M_PNRW_BigInt {
    static SCALE_BITS= SCALE_BITS
    static SCALE= SCALE
    static PI= PI_80
    static WGS84_A= MPNRW_BIGINT_CONSTANTS.WGS84_A
    static WGS84_F= MPNRW_BIGINT_CONSTANTS.WGS84_F
    static WGS84_E_SQ= MPNRW_BIGINT_CONSTANTS.WGS84_E_SQ

    static fromInt= mpnrw_fp_from_int
    static fromDecimal= mpnrw_fp_from_decimal
    static fromNumber= mpnrw_fp_from_number
    static fromRational= mpnrw_fp_from_rational
    static toDecimal= mpnrw_fp_to_decimal
    static toNumber= mpnrw_fp_to_number

    static fpAdd= mpnrw_fp_add
    static fpSub= mpnrw_fp_sub
    static fpMul= mpnrw_fp_mul
    static fpDiv= mpnrw_fp_div
    static fpMod= mpnrw_fp_mod
    static fpSqrt= mpnrw_fp_sqrt
    static fpSinCos= mpnrw_fp_sin_cos
    static fpSin= mpnrw_fp_sin
    static fpCos= mpnrw_fp_cos

    #ctx

    constructor() {
        this.#ctx= mpnrw_bigint_create_ctx()
        this.rad_deg= this.#ctx.rad_deg
        this.mas_to_rad= this.#ctx.mas_to_rad
        this.nutation_terms= this.#ctx.nutation_terms
    }

    depack_vector(uint8_buffer) {
        return mpnrw_bigint_depack_vector(uint8_buffer)
    }

    get_julian_date(date) {
        return mpnrw_bigint_julian_date(date)
    }

    to_ecef(lat, lon, height, a= M_PNRW_BigInt.WGS84_A, e_sq= M_PNRW_BigInt.WGS84_E_SQ) {
        return mpnrw_bigint_to_ecef(this.#ctx, lat, lon, height, a, e_sq)
    }

    get_distance(v1, v2) {
        return mpnrw_bigint_distance(v1, v2)
    }

    get_heliocentric_earth_m(jd_tt) {
        return mpnrw_bigint_heliocentric_earth_m(this.#ctx, jd_tt)
    }

    get_sun_vector(jd) {
        return mpnrw_bigint_sun_vector(this.#ctx, jd)
    }

    multiply_matrices(a, b) {
        return mpnrw_bigint_multiply_matrices(a, b)
    }

    multiply_matrix_vector(m, v) {
        return mpnrw_bigint_multiply_matrix_vector(m, v)
    }

    rot_y(angle) {
        return mpnrw_bigint_rot_y(angle)
    }

    rot_z(angle) {
        return mpnrw_bigint_rot_z(angle)
    }

    get_p_matrix(t) {
        return mpnrw_bigint_p_matrix(this.#ctx, t)
    }

    get_n_matrix(t) {
        return mpnrw_bigint_n_matrix(this.#ctx, t)
    }

    get_r_matrix(jd) {
        return mpnrw_bigint_r_matrix(this.#ctx, jd)
    }

    get_w_matrix(xp= 0n, yp= 0n) {
        return mpnrw_bigint_w_matrix(xp, yp)
    }

    get_universal_address(lat, lon, height, date, origin_to_earth_m, orientation= {}) {
        return mpnrw_bigint_universal_address(this.#ctx, lat, lon, height, date, origin_to_earth_m, orientation)
    }
}
