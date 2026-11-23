/**
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
 */

import { M_PNRW } from '../src-mjs.dev/M_PNRW.mjs'
import { M_PNRW_BigInt } from '../src-mjs.dev/M_PNRW_BigInt.mjs'

let passed = 0
let failed = 0

function check(name, cond, detail = '') {
    if (cond) {
        console.log(`  [PASS] ${name}`)
        passed++
    } else {
        console.error(`  [FAIL] ${name}${detail ? ' — ' + detail : ''}`)
        failed++
    }
}

console.log('--- M_PNRW_BigInt Engine Test ---')

const biEngine = new M_PNRW_BigInt()
const fpEngine = new M_PNRW()
const ONE = M_PNRW_BigInt.SCALE

// 1. Fixed-point arithmetic & trig primitives
console.log('\n[Test 1] Fixed-point primitives & trig identity')
{
    const a = M_PNRW_BigInt.fromDecimal('3.5')
    const b = M_PNRW_BigInt.fromDecimal('-2.25')
    check('fpMul(3.5, -2.25) === -7.875', M_PNRW_BigInt.fpMul(a, b) === M_PNRW_BigInt.fromDecimal('-7.875'))
    check('fpDiv(-7.875, 3.5) === -2.25', M_PNRW_BigInt.fpDiv(M_PNRW_BigInt.fromDecimal('-7.875'), a) === b)
    check('fpSqrt(9) === 3', M_PNRW_BigInt.fpSqrt(M_PNRW_BigInt.fromInt(9)) === M_PNRW_BigInt.fromInt(3))

    const angles = ['0', '0.5', '1.5707963267948966', '3.141592653589793', '-2.3456789', '100000.123456']
    let maxTrigErr = 0n
    for (const s of angles) {
        const ang = M_PNRW_BigInt.fromDecimal(s)
        const [sn, cs] = M_PNRW_BigInt.fpSinCos(ang)
        const sumSq = M_PNRW_BigInt.fpMul(sn, sn) + M_PNRW_BigInt.fpMul(cs, cs)
        const err = sumSq > ONE ? sumSq - ONE : ONE - sumSq
        if (err > maxTrigErr) maxTrigErr = err
    }
    // Tolerance: 16 ULPs of Q.80 (~1.3e-23)
    check('sin^2 + cos^2 == 1 within 16 ULPs of Q.80', maxTrigErr <= 16n, `maxErr=${maxTrigErr}`)
}

// 2. Exact Julian Date at J2000.0 and Unix epoch
console.log('\n[Test 2] Exact Julian Date')
{
    const jdJ2000 = biEngine.get_julian_date(new Date('2000-01-01T12:00:00.000Z'))
    check('J2000.0 is exact 2451545.0', jdJ2000 === (2451545n << M_PNRW_BigInt.SCALE_BITS))

    const jdEpoch = biEngine.get_julian_date(0n)
    check('Unix epoch is exact 2440587.5', jdEpoch === M_PNRW_BigInt.fromDecimal('2440587.5'))
}

// 3. Precession matrix orthogonality at Q.80
console.log('\n[Test 3] Precession matrix orthogonality (P * P^T = I)')
{
    const t = M_PNRW_BigInt.fromDecimal('0.2625') // ~2026
    const p = biEngine.get_p_matrix(t)
    let maxOrthoErr = 0n
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            const dot =
                M_PNRW_BigInt.fpMul(p[i][0], p[j][0]) +
                M_PNRW_BigInt.fpMul(p[i][1], p[j][1]) +
                M_PNRW_BigInt.fpMul(p[i][2], p[j][2])
            const expected = i === j ? ONE : 0n
            const err = dot > expected ? dot - expected : expected - dot
            if (err > maxOrthoErr) maxOrthoErr = err
        }
    }
    // Within 64 ULPs of Q.80 (~5e-23)
    check('P * P^T = I within 64 ULPs of Q.80', maxOrthoErr <= 64n, `maxErr=${maxOrthoErr}`)
}

// 4. Parity against f64 M_PNRW across multiple locations & dates
console.log('\n[Test 4] Parity with f64 M_PNRW (sub-millimeter agreement)')
{
    const cases = [
        { name: 'Greenwich J2000', lat: 51.477928, lon: 0, h: 46, date: new Date('2000-01-01T12:00:00Z') },
        { name: 'Stonehenge Equinox 2026', lat: 51.1789, lon: -1.8262, h: 103, date: new Date('2026-03-20T12:00:00Z') },
        { name: 'Sydney 2050', lat: -33.8688, lon: 151.2093, h: 19, date: new Date('2050-07-15T04:30:00Z') }
    ]

    for (const c of cases) {
        const fpVec = fpEngine.get_universal_address(c.lat, c.lon, c.h, c.date, [0, 0, 0])
        const biVec = biEngine.get_universal_address(
            M_PNRW_BigInt.fromNumber(c.lat),
            M_PNRW_BigInt.fromNumber(c.lon),
            M_PNRW_BigInt.fromNumber(c.h),
            c.date,
            [0n, 0n, 0n]
        )

        const diffs = biVec.map((v, idx) => Math.abs(M_PNRW_BigInt.toNumber(v) - fpVec[idx]))
        const maxDiff = Math.max(...diffs)
        // f64 Julian Date (~2.46e6 days) has ~4.5e-10 d quantization (~18 um Earth rotation)
        check(`${c.name} agrees within 1e-4 m`, maxDiff < 1e-4, `maxDiff=${maxDiff}`)
    }
}

// 5. High-fidelity large-scale offset preservation (1 AU in meters + Earth ECEF)
console.log('\n[Test 5] Sub-millimeter fidelity at 1 AU scale (no f64 cancellation)')
{
    const date = new Date('2026-03-20T12:00:00Z')
    const oneAuMeters = M_PNRW_BigInt.fromDecimal('149597870700.0')
    const sunOffset = [oneAuMeters, 0n, 0n]

    const addrZero = biEngine.get_universal_address(
        M_PNRW_BigInt.fromDecimal('51.5073'),
        M_PNRW_BigInt.fromDecimal('-0.1277'),
        M_PNRW_BigInt.fromDecimal('15'),
        date,
        [0n, 0n, 0n]
    )
    const addrAu = biEngine.get_universal_address(
        M_PNRW_BigInt.fromDecimal('51.5073'),
        M_PNRW_BigInt.fromDecimal('-0.1277'),
        M_PNRW_BigInt.fromDecimal('15'),
        date,
        sunOffset
    )

    // Subtracting the 1 AU offset recovers the exact same BigInt bit-for-bit!
    const recoveredX = addrAu[0] - oneAuMeters
    check('1 AU offset subtraction is bit-exact', recoveredX === addrZero[0])
}

// 6. Sun vector parity & deterministic hex snapshot
console.log('\n[Test 6] Sun vector parity & deterministic hex output')
{
    const date = new Date('2026-03-20T12:00:00Z')
    const fpSun = fpEngine.get_sun_vector(fpEngine.get_julian_date(date))
    const biSun = biEngine.get_sun_vector(biEngine.get_julian_date(date))
    const maxSunDiff = Math.max(...biSun.map((v, i) => Math.abs(M_PNRW_BigInt.toNumber(v) - fpSun[i])))
    check('Sun vector agrees with f64 within 1e-3 km', maxSunDiff < 1e-3, `maxDiff=${maxSunDiff}`)

    // Verify deterministic bit-exactness across two fresh instances
    const biEngine2 = new M_PNRW_BigInt()
    const addr1 = biEngine.get_universal_address(
        M_PNRW_BigInt.fromDecimal('51.477928'),
        0n,
        M_PNRW_BigInt.fromDecimal('46'),
        new Date('2000-01-01T12:00:00Z'),
        [0n, 0n, 0n]
    )
    const addr2 = biEngine2.get_universal_address(
        M_PNRW_BigInt.fromDecimal('51.477928'),
        0n,
        M_PNRW_BigInt.fromDecimal('46'),
        new Date('2000-01-01T12:00:00Z'),
        [0n, 0n, 0n]
    )
    check('Two instances produce bit-identical BigInts', addr1.every((v, i) => v === addr2[i]))
}

console.log('\n--- Summary ---')
console.log(`Passed: ${passed} / ${passed + failed}`)
if (failed > 0) process.exit(1)
