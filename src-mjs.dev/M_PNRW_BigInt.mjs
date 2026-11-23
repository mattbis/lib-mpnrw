/**
 * M_PNRW_BigInt: Deterministic Fixed-Point Celestial Transformation Engine
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
 *
 * Arithmetic: Q.80 signed binary fixed-point using native BigInt.
 * Every operation (+, -, *, /, sqrt, sin, cos) is pure integer arithmetic
 * with zero IEEE-754 floating-point or platform libm dependencies, ensuring
 * bit-for-bit identical results across all CPU architectures and JS engines.
 */

const SCALE_BITS = 80n
const SCALE = 1n << SCALE_BITS
const HALF_SCALE = 1n << (SCALE_BITS - 1n)

// Internal guard bits for trigonometric Taylor evaluation
const TRIG_GUARD_BITS = 16n
const TRIG_SHIFT = SCALE_BITS + TRIG_GUARD_BITS // 96n
const TRIG_HALF_GUARD = 1n << (TRIG_GUARD_BITS - 1n)

// Pi at 160 fractional bits (for high-precision range reduction of large angles)
const PI_160 = 0x3243f6a8885a308d313198a2e03707344a4093822n
const TWO_PI_160 = PI_160 << 1n

// Pi constants at Q.80
const PI_80 = (PI_160 + (1n << 79n)) >> 80n
const HALF_PI_80 = (PI_160 + (1n << 80n)) >> 81n

const MS_PER_DAY = 86_400_000n
// Julian Date at Unix epoch (1970-01-01T00:00:00.000Z) = 2440587.5
const JD_UNIX_EPOCH_80 = (2440587n << SCALE_BITS) + HALF_SCALE
// J2000.0 = 2451545.0
const JD_J2000_80 = 2451545n << SCALE_BITS

/**
 * Parses a decimal string (including optional sign and scientific notation)
 * into a Q.80 fixed-point BigInt, rounded half-away-from-zero.
 * @param {string} str
 * @returns {bigint}
 */
function parseDecimalToQ80(str) {
    const trimmed = str.trim()
    const match = /^([+-])?(?:(\d+)(?:\.(\d*))?|\.(\d+))(?:[eE]([+-]?\d+))?$/.exec(trimmed)
    if (!match) {
        throw new Error(`Invalid decimal string: "${str}"`)
    }

    const neg = match[1] === '-'
    const intPart = match[2] ?? ''
    const fracPart = match[3] ?? match[4] ?? ''
    const exp = match[5] ? parseInt(match[5], 10) : 0

    const digits = BigInt((intPart + fracPart) || '0')
    if (digits === 0n) return 0n

    const pow10 = exp - fracPart.length
    let num = digits << SCALE_BITS
    let den = 1n

    if (pow10 >= 0) {
        num *= 10n ** BigInt(pow10)
    } else {
        den = 10n ** BigInt(-pow10)
    }

    const q = (num + (den >> 1n)) / den
    return neg ? -q : q
}

export class M_PNRW_BigInt {
    static SCALE_BITS = SCALE_BITS
    static SCALE = SCALE
    static PI = PI_80

    // --- CONVERSION HELPERS ---

    /**
     * Converts an integer (number or bigint) to Q.80 fixed-point.
     * @param {number|bigint} n
     * @returns {bigint}
     */
    static fromInt(n) {
        return BigInt(n) << SCALE_BITS
    }

    /**
     * Converts an exact decimal string (e.g. "51.477928") to Q.80 fixed-point.
     * @param {string} str
     * @returns {bigint}
     */
    static fromDecimal(str) {
        if (typeof str !== 'string') {
            throw new TypeError('fromDecimal expects a string')
        }
        return parseDecimalToQ80(str)
    }

    /**
     * Converts a finite JS number to Q.80 fixed-point via its decimal representation.
     * @param {number} n
     * @returns {bigint}
     */
    static fromNumber(n) {
        if (typeof n !== 'number' || !Number.isFinite(n)) {
            throw new RangeError('fromNumber expects a finite number')
        }
        return parseDecimalToQ80(n.toString())
    }

    /**
     * Converts a rational num / den into Q.80 fixed-point.
     * @param {bigint|number} num
     * @param {bigint|number} den
     * @returns {bigint}
     */
    static fromRational(num, den) {
        const a = BigInt(num)
        const b = BigInt(den)
        if (b === 0n) throw new RangeError('Division by zero')
        const neg = (a < 0n) !== (b < 0n)
        const absA = a < 0n ? -a : a
        const absB = b < 0n ? -b : b
        const q = ((absA << SCALE_BITS) + (absB >> 1n)) / absB
        return neg ? -q : q
    }

    /**
     * Converts a Q.80 fixed-point BigInt back to a JS number (lossy for display/interop).
     * @param {bigint} fp
     * @returns {number}
     */
    static toNumber(fp) {
        return Number(M_PNRW_BigInt.toDecimal(fp, 17))
    }

    /**
     * Formats a Q.80 fixed-point BigInt as a decimal string without float64 rounding.
     * @param {bigint} fp
     * @param {number} [decimals=12]
     * @returns {string}
     */
    static toDecimal(fp, decimals = 12) {
        if (typeof fp !== 'bigint') throw new TypeError('toDecimal expects a bigint')
        if (!Number.isInteger(decimals) || decimals < 0) {
            throw new RangeError('decimals must be a non-negative integer')
        }

        const neg = fp < 0n
        const abs = neg ? -fp : fp
        const factor = 10n ** BigInt(decimals)
        const scaled = ((abs * factor) + HALF_SCALE) >> SCALE_BITS

        const intPart = scaled / factor
        if (decimals === 0) {
            return (neg && intPart !== 0n ? '-' : '') + intPart.toString()
        }

        const fracPart = (scaled % factor).toString().padStart(decimals, '0')
        const sign = neg && scaled !== 0n ? '-' : ''
        return `${sign}${intPart.toString()}.${fracPart}`
    }

    // --- FIXED-POINT ARITHMETIC PRIMITIVES (Q.80) ---

    static fpAdd(a, b) {
        return a + b
    }

    static fpSub(a, b) {
        return a - b
    }

    static fpMul(a, b) {
        const prod = a * b
        return prod >= 0n
            ? (prod + HALF_SCALE) >> SCALE_BITS
            : -((-prod + HALF_SCALE) >> SCALE_BITS)
    }

    static fpDiv(a, b) {
        if (b === 0n) throw new RangeError('Division by zero')
        const neg = (a < 0n) !== (b < 0n)
        const absA = a < 0n ? -a : a
        const absB = b < 0n ? -b : b
        const q = ((absA << SCALE_BITS) + (absB >> 1n)) / absB
        return neg ? -q : q
    }

    static fpMod(a, b) {
        if (b === 0n) throw new RangeError('Modulo by zero')
        return a % b
    }

    static fpSqrt(a) {
        if (a < 0n) throw new RangeError('Square root of negative fixed-point number')
        if (a === 0n) return 0n

        const n = a << SCALE_BITS
        const bitLen = n.toString(2).length
        let y = 1n << BigInt((bitLen + 1) >> 1)

        while (true) {
            const next = (y + n / y) >> 1n
            if (next >= y) break
            y = next
        }
        return y
    }

    /**
     * Computes [sin(angle), cos(angle)] for a Q.80 angle in radians.
     * Uses 160-bit Pi range reduction and 96-bit internal Taylor evaluation
     * on [-pi/4, pi/4] for full Q.80 determinism and precision.
     * @param {bigint} angle - Angle in Q.80 radians.
     * @returns {[bigint, bigint]} [sin, cos] in Q.80.
     */
    static fpSinCos(angle) {
        // 1. Reduce modulo 2*pi using 160-bit precision
        const a160 = angle << SCALE_BITS
        const halfTwoPi = TWO_PI_160 >> 1n
        const k = a160 >= 0n
            ? (a160 + halfTwoPi) / TWO_PI_160
            : -((-a160 + halfTwoPi) / TWO_PI_160)
        const r = (a160 - k * TWO_PI_160) >> SCALE_BITS // r in [-pi, pi] in Q.80

        // 2. Quadrant reduction into [-pi/4, pi/4]
        const halfQuarter = HALF_PI_80 >> 1n
        const q = r >= 0n
            ? (r + halfQuarter) / HALF_PI_80
            : -((-r + halfQuarter) / HALF_PI_80)
        const u = r - q * HALF_PI_80 // |u| <= pi/4 in Q.80

        // 3. Evaluate sin(u) and cos(u) with 16 guard bits (Q.96)
        const u96 = u << TRIG_GUARD_BITS
        const u2_96 = (u96 * u96) >> TRIG_SHIFT

        // sin(u) = u - u^3/3! + u^5/5! - ... (up to u^25/25!)
        let sinSum = u96
        let sinTerm = u96
        for (let n = 1n; n <= 12n; n++) {
            const denom = (2n * n) * (2n * n + 1n)
            sinTerm = -((sinTerm * u2_96) >> TRIG_SHIFT) / denom
            if (sinTerm === 0n) break
            sinSum += sinTerm
        }

        // cos(u) = 1 - u^2/2! + u^4/4! - ... (up to u^24/24!)
        let cosSum = 1n << TRIG_SHIFT
        let cosTerm = 1n << TRIG_SHIFT
        for (let n = 1n; n <= 12n; n++) {
            const denom = (2n * n - 1n) * (2n * n)
            cosTerm = -((cosTerm * u2_96) >> TRIG_SHIFT) / denom
            if (cosTerm === 0n) break
            cosSum += cosTerm
        }

        const s = sinSum >= 0n
            ? (sinSum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS
            : -((-sinSum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS)
        const c = cosSum >= 0n
            ? (cosSum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS
            : -((-cosSum + TRIG_HALF_GUARD) >> TRIG_GUARD_BITS)

        const quad = ((q % 4n) + 4n) % 4n
        if (quad === 0n) return [s, c]
        if (quad === 1n) return [c, -s]
        if (quad === 2n) return [-s, -c]
        return [-c, s]
    }

    static fpSin(angle) {
        return M_PNRW_BigInt.fpSinCos(angle)[0]
    }

    static fpCos(angle) {
        return M_PNRW_BigInt.fpSinCos(angle)[1]
    }

    // --- WGS84 CONSTANTS (Q.80) ---

    static WGS84_A = 6378137n << SCALE_BITS
    static WGS84_F = 0xdbbad4279f09443c5en
    static WGS84_E_SQ = (2n * M_PNRW_BigInt.WGS84_F) - M_PNRW_BigInt.fpMul(M_PNRW_BigInt.WGS84_F, M_PNRW_BigInt.WGS84_F)

    constructor() {
        // --- CONSTANTS (Q.80) ---
        this.rad_deg = M_PNRW_BigInt.fpDiv(PI_80, 180n << SCALE_BITS)
        this.mas_to_rad = M_PNRW_BigInt.fpDiv(this.rad_deg, 3_600_000n << SCALE_BITS)

        // Pre-parsed decimal constants used across P, N, R, Sun calculations
        this._c = {
            daysPerCentury: 36525n << SCALE_BITS,
            hoursPerDay: 24n << SCALE_BITS,
            degPerHour: 15n << SCALE_BITS,
            arcsecToRad: M_PNRW_BigInt.fpDiv(this.rad_deg, 3600n << SCALE_BITS),

            // Precession (IAU 1976)
            pZeta: parseDecimalToQ80('2306.2181'),
            pTheta: parseDecimalToQ80('2004.3109'),

            // Nutation Delaunay arguments (deg, deg/century)
            nL0: parseDecimalToQ80('134.963'),
            nL1: parseDecimalToQ80('477198.868'),
            nLp0: parseDecimalToQ80('357.529'),
            nLp1: parseDecimalToQ80('35999.050'),
            nF0: parseDecimalToQ80('93.272'),
            nF1: parseDecimalToQ80('483202.018'),
            nD0: parseDecimalToQ80('297.850'),
            nD1: parseDecimalToQ80('445267.111'),
            nOm0: parseDecimalToQ80('125.045'),
            nOm1: parseDecimalToQ80('-1934.136'),
            nEps0: parseDecimalToQ80('23.43929'),
            nEps1: parseDecimalToQ80('0.013'),

            // GMST
            gmst0: parseDecimalToQ80('18.697374558'),
            gmst1: parseDecimalToQ80('24.06570982441908'),

            // Sun vector
            sL0: parseDecimalToQ80('280.460'),
            sL1: parseDecimalToQ80('36000.771'),
            sG0: parseDecimalToQ80('357.528'),
            sG1: parseDecimalToQ80('35999.050'),
            sLam1: parseDecimalToQ80('1.915'),
            sLam2: parseDecimalToQ80('0.020'),
            sR0: parseDecimalToQ80('1.00014'),
            sR1: parseDecimalToQ80('0.01671'),
            sR2: parseDecimalToQ80('0.00014'),
            sEps0: parseDecimalToQ80('23.439'),
            sEps1: parseDecimalToQ80('0.013'),
            auKm: parseDecimalToQ80('149597870.7'),
            deg360: 360n << SCALE_BITS
        }

        // Principal Nutation Terms (IAU 2000B simplified)
        this.nutation_terms = [
            { l: 0n, lp: 0n, f: 0n, d: 0n, om: 1n, psi: -17206n << SCALE_BITS, eps: 9205n << SCALE_BITS },
            { l: 0n, lp: 0n, f: 2n, d: -2n, om: 2n, psi: -1317n << SCALE_BITS, eps: 573n << SCALE_BITS },
            { l: 0n, lp: 0n, f: 2n, d: 0n, om: 2n, psi: -227n << SCALE_BITS, eps: 97n << SCALE_BITS },
            { l: 0n, lp: 1n, f: 0n, d: 0n, om: 0n, psi: 147n << SCALE_BITS, eps: 7n << SCALE_BITS }
        ]
    }

    // --- CONVERSION & UTILITY METHODS ---

    /**
     * Unpacks the 5 BigInt64 fields from a 40-byte "Sealed Envelope" buffer.
     * @param {Uint8Array} uint8_buffer
     * @returns {{ lat_bi: bigint, lon_bi: bigint, height_bi: bigint, date_bi: bigint, sun_bi: bigint }}
     */
    depack_vector(uint8_buffer) {
        const view = new DataView(uint8_buffer.buffer, uint8_buffer.byteOffset, uint8_buffer.byteLength)
        return {
            lat_bi: view.getBigInt64(0),
            lon_bi: view.getBigInt64(8),
            height_bi: view.getBigInt64(16),
            date_bi: view.getBigInt64(24),
            sun_bi: view.getBigInt64(32)
        }
    }

    /**
     * Calculates the Julian Date in Q.80 fixed-point from a Date object or
     * UTC Unix timestamp in milliseconds (bigint).
     * @param {Date|bigint} date
     * @returns {bigint} Julian Date in Q.80.
     */
    get_julian_date(date) {
        let unixMs
        if (date instanceof Date) {
            const ms = date.getTime()
            if (Number.isNaN(ms)) throw new Error('Invalid date')
            unixMs = BigInt(ms)
        } else if (typeof date === 'bigint') {
            unixMs = date
        } else {
            throw new Error('date must be a Date object or bigint UTC milliseconds')
        }

        // Exact: 2440587.5 + unixMs / 86_400_000
        const daysFp = M_PNRW_BigInt.fromRational(unixMs, MS_PER_DAY)
        return JD_UNIX_EPOCH_80 + daysFp
    }

    /**
     * Converts WGS84 geodetic coordinates to ECEF [x, y, z] in Q.80 meters.
     * @param {bigint} lat - Latitude in Q.80 degrees.
     * @param {bigint} lon - Longitude in Q.80 degrees.
     * @param {bigint} height - Height above ellipsoid in Q.80 meters.
     * @param {bigint} [a=M_PNRW_BigInt.WGS84_A]
     * @param {bigint} [e_sq=M_PNRW_BigInt.WGS84_E_SQ]
     * @returns {bigint[]} [x, y, z] in Q.80 meters.
     */
    to_ecef(lat, lon, height, a = M_PNRW_BigInt.WGS84_A, e_sq = M_PNRW_BigInt.WGS84_E_SQ) {
        const { fpMul, fpDiv, fpSqrt, fpSinCos } = M_PNRW_BigInt
        const rad_lat = fpMul(lat, this.rad_deg)
        const rad_lon = fpMul(lon, this.rad_deg)

        const [sin_lat, cos_lat] = fpSinCos(rad_lat)
        const [sin_lon, cos_lon] = fpSinCos(rad_lon)

        const sin2 = fpMul(sin_lat, sin_lat)
        const n_phi = fpDiv(a, fpSqrt(SCALE - fpMul(e_sq, sin2)))
        const n_plus_h = n_phi + height

        const x = fpMul(fpMul(n_plus_h, cos_lat), cos_lon)
        const y = fpMul(fpMul(n_plus_h, cos_lat), sin_lon)
        const z = fpMul(fpMul(n_phi, SCALE - e_sq) + height, sin_lat)

        return [x, y, z]
    }

    /**
     * Calculates the Euclidean distance between two Q.80 3D vectors.
     * @param {bigint[]} v1
     * @param {bigint[]} v2
     * @returns {bigint}
     */
    get_distance(v1, v2) {
        const { fpMul, fpSqrt } = M_PNRW_BigInt
        const dx = v1[0] - v2[0]
        const dy = v1[1] - v2[1]
        const dz = v1[2] - v2[2]
        return fpSqrt(fpMul(dx, dx) + fpMul(dy, dy) + fpMul(dz, dz))
    }

    /**
     * Calculates the position of the Sun in ICRF coordinates for a given Q.80 Julian Date.
     * @param {bigint} jd - Julian Date in Q.80.
     * @returns {bigint[]} [x, y, z] in Q.80 kilometers.
     */
    get_sun_vector(jd) {
        const { fpMul, fpDiv, fpSinCos } = M_PNRW_BigInt
        const c = this._c

        const t = fpDiv(jd - JD_J2000_80, c.daysPerCentury)
        const l = (c.sL0 + fpMul(c.sL1, t)) % c.deg360
        const g = fpMul(c.sG0 + fpMul(c.sG1, t), this.rad_deg)

        const [sin_g, cos_g] = fpSinCos(g)
        const [sin_2g, cos_2g] = fpSinCos(g << 1n)

        const lambda = fpMul(l + fpMul(c.sLam1, sin_g) + fpMul(c.sLam2, sin_2g), this.rad_deg)
        const r_au = c.sR0 - fpMul(c.sR1, cos_g) - fpMul(c.sR2, cos_2g)
        const eps = fpMul(c.sEps0 - fpMul(c.sEps1, t), this.rad_deg)

        const [sin_lam, cos_lam] = fpSinCos(lambda)
        const [sin_eps, cos_eps] = fpSinCos(eps)
        const r_km = fpMul(r_au, c.auKm)

        return [
            fpMul(r_km, cos_lam),
            fpMul(fpMul(r_km, sin_lam), cos_eps),
            fpMul(fpMul(r_km, sin_lam), sin_eps)
        ]
    }

    // --- MATRIX MATH ---

    multiply_matrices(a, b) {
        const { fpMul } = M_PNRW_BigInt
        const res = [
            [0n, 0n, 0n],
            [0n, 0n, 0n],
            [0n, 0n, 0n]
        ]
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                res[i][j] = fpMul(a[i][0], b[0][j]) + fpMul(a[i][1], b[1][j]) + fpMul(a[i][2], b[2][j])
            }
        }
        return res
    }

    multiply_matrix_vector(m, v) {
        const { fpMul } = M_PNRW_BigInt
        return [
            fpMul(m[0][0], v[0]) + fpMul(m[0][1], v[1]) + fpMul(m[0][2], v[2]),
            fpMul(m[1][0], v[0]) + fpMul(m[1][1], v[1]) + fpMul(m[1][2], v[2]),
            fpMul(m[2][0], v[0]) + fpMul(m[2][1], v[1]) + fpMul(m[2][2], v[2])
        ]
    }

    // --- PNRW COMPONENT GENERATORS ---

    rot_y(angle) {
        const [s, c] = M_PNRW_BigInt.fpSinCos(angle)
        return [
            [c, 0n, -s],
            [0n, SCALE, 0n],
            [s, 0n, c]
        ]
    }

    rot_z(angle) {
        const [s, c] = M_PNRW_BigInt.fpSinCos(angle)
        return [
            [c, s, 0n],
            [-s, c, 0n],
            [0n, 0n, SCALE]
        ]
    }

    get_p_matrix(t) {
        const { fpMul } = M_PNRW_BigInt
        const c = this._c

        const zeta = fpMul(fpMul(c.pZeta, t), c.arcsecToRad)
        const z = zeta
        const theta = fpMul(fpMul(c.pTheta, t), c.arcsecToRad)

        const rz_zeta = this.rot_z(-zeta)
        const ry_theta = this.rot_y(theta)
        const rz_z = this.rot_z(-z)

        return this.multiply_matrices(rz_z, this.multiply_matrices(ry_theta, rz_zeta))
    }

    get_n_matrix(t) {
        const { fpMul, fpSinCos } = M_PNRW_BigInt
        const c = this._c

        const l = fpMul(c.nL0 + fpMul(c.nL1, t), this.rad_deg)
        const lp = fpMul(c.nLp0 + fpMul(c.nLp1, t), this.rad_deg)
        const f = fpMul(c.nF0 + fpMul(c.nF1, t), this.rad_deg)
        const d = fpMul(c.nD0 + fpMul(c.nD1, t), this.rad_deg)
        const om = fpMul(c.nOm0 + fpMul(c.nOm1, t), this.rad_deg)

        let d_psi = 0n
        let d_eps = 0n
        const eps = fpMul(c.nEps0 - fpMul(c.nEps1, t), this.rad_deg)

        for (let i = 0; i < this.nutation_terms.length; i++) {
            const term = this.nutation_terms[i]
            const arg = (term.l * l) + (term.lp * lp) + (term.f * f) + (term.d * d) + (term.om * om)
            const [s_arg, c_arg] = fpSinCos(arg)
            d_psi += fpMul(term.psi, s_arg)
            d_eps += fpMul(term.eps, c_arg)
        }

        const dp = fpMul(d_psi, this.mas_to_rad)
        const de = fpMul(d_eps, this.mas_to_rad)
        const [sin_eps, cos_eps] = fpSinCos(eps)
        const dp_cos = fpMul(dp, cos_eps)
        const dp_sin = fpMul(dp, sin_eps)

        return [
            [SCALE, -dp_cos, -dp_sin],
            [dp_cos, SCALE, -de],
            [dp_sin, de, SCALE]
        ]
    }

    get_r_matrix(jd) {
        const { fpMul, fpSinCos } = M_PNRW_BigInt
        const c = this._c

        const days_since_j2000 = jd - JD_J2000_80
        let gmst_hours = (c.gmst0 + fpMul(c.gmst1, days_since_j2000)) % c.hoursPerDay
        if (gmst_hours < 0n) gmst_hours += c.hoursPerDay

        const theta = fpMul(fpMul(gmst_hours, c.degPerHour), this.rad_deg)
        const [s, cos_t] = fpSinCos(theta)

        return [
            [cos_t, s, 0n],
            [-s, cos_t, 0n],
            [0n, 0n, SCALE]
        ]
    }

    get_w_matrix(xp = 0n, yp = 0n) {
        return [
            [SCALE, 0n, xp],
            [0n, SCALE, -yp],
            [-xp, yp, SCALE]
        ]
    }

    // --- MAIN ENGINE ---

    /**
     * Calculates the Universal Address in Q.80 fixed-point.
     * @param {bigint} lat - Latitude in Q.80 degrees [-90, 90].
     * @param {bigint} lon - Longitude in Q.80 degrees [-180, 180].
     * @param {bigint} height - Height above ellipsoid in Q.80 meters (>= 0).
     * @param {Date|bigint} date - Date object or bigint UTC milliseconds.
     * @param {bigint[]} sun_vector - 3-element array of Q.80 BigInts.
     * @returns {bigint[]} Universal Space Vector [x, y, z] in Q.80.
     */
    get_universal_address(lat, lon, height, date, sun_vector) {
        if (typeof lat !== 'bigint' || lat < (-90n << SCALE_BITS) || lat > (90n << SCALE_BITS)) {
            throw new Error('Latitude must be a Q.80 bigint between -90 and 90 degrees')
        }
        if (typeof lon !== 'bigint' || lon < (-180n << SCALE_BITS) || lon > (180n << SCALE_BITS)) {
            throw new Error('Longitude must be a Q.80 bigint between -180 and 180 degrees')
        }
        if (typeof height !== 'bigint' || height < 0n) {
            throw new Error('Height must be a non-negative Q.80 bigint')
        }
        if (!(date instanceof Date) && typeof date !== 'bigint') {
            throw new Error('Date must be a Date object or bigint UTC milliseconds')
        }
        if (!Array.isArray(sun_vector) || sun_vector.length !== 3 || sun_vector.some(v => typeof v !== 'bigint')) {
            throw new Error('Sun vector must be an array of 3 Q.80 bigints')
        }

        const jd = this.get_julian_date(date)
        const t = M_PNRW_BigInt.fpDiv(jd - JD_J2000_80, this._c.daysPerCentury)
        const local_v = this.to_ecef(lat, lon, height)

        const p = this.get_p_matrix(t)
        const n = this.get_n_matrix(t)
        const r = this.get_r_matrix(jd)
        const w = this.get_w_matrix()

        const m_pn = this.multiply_matrices(p, n)
        const m_pnr = this.multiply_matrices(m_pn, r)
        const master_m = this.multiply_matrices(m_pnr, w)
        const space_v = this.multiply_matrix_vector(master_m, local_v)

        return [
            space_v[0] + sun_vector[0],
            space_v[1] + sun_vector[1],
            space_v[2] + sun_vector[2]
        ]
    }
}
