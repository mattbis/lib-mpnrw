// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Matthew & Contributors — Free Software (GNU GPLv3)

// 1p
import test from 'node:test'
import assert from 'node:assert/strict'

// 2p
import { M_PNRW } from '../src-mjs.dev/M_PNRW.mjs'
import {
    M_PNRW_BigInt,
    mpnrw_fp_mul,
    mpnrw_fp_div,
    mpnrw_fp_sqrt,
    mpnrw_fp_sin_cos,
    mpnrw_fp_from_int,
    mpnrw_fp_from_decimal,
    mpnrw_fp_from_number,
    mpnrw_fp_to_number
} from '../src-mjs.dev/M_PNRW_BigInt.mjs'

// 3p
// (none)

const ONE= M_PNRW_BigInt.SCALE

test('M_PNRW_BigInt fixed-point primitives and trig identity', () => {
    var a= mpnrw_fp_from_decimal('3.5')
    var b= mpnrw_fp_from_decimal('-2.25')
    assert.equal(mpnrw_fp_mul(a, b), mpnrw_fp_from_decimal('-7.875'))
    assert.equal(mpnrw_fp_div(mpnrw_fp_from_decimal('-7.875'), a), b)
    assert.equal(mpnrw_fp_sqrt(mpnrw_fp_from_int(9)), mpnrw_fp_from_int(3))

    var angles= ['0', '0.5', '1.5707963267948966', '3.141592653589793', '-2.3456789', '100000.123456']
    var max_err= 0n
    for (var s of angles) {
        var ang= mpnrw_fp_from_decimal(s)
        var [sn, cs]= mpnrw_fp_sin_cos(ang)
        var sum_sq= mpnrw_fp_mul(sn, sn) + mpnrw_fp_mul(cs, cs)
        var err= (sum_sq > ONE)
            ? sum_sq - ONE
            : ONE - sum_sq
        if (err > max_err) max_err= err
    }
    assert.ok(max_err <= 16n, `Expected max trig err <= 16 ULPs, got ${max_err}`)
})

test('M_PNRW_BigInt exact Julian Date', () => {
    var engine= new M_PNRW_BigInt()
    var jd_j2000= engine.get_julian_date(new Date('2000-01-01T12:00:00.000Z'))
    assert.equal(jd_j2000, 2451545n << M_PNRW_BigInt.SCALE_BITS)

    var jd_epoch= engine.get_julian_date(0n)
    assert.equal(jd_epoch, mpnrw_fp_from_decimal('2440587.5'))
})

test('M_PNRW_BigInt precession matrix orthogonality', () => {
    var engine= new M_PNRW_BigInt()
    var t= mpnrw_fp_from_decimal('0.2625')
    var p= engine.get_p_matrix(t)
    var max_err= 0n

    for (let i= 0; i < 3; i= i + 1) {
        for (let j= 0; j < 3; j= j + 1) {
            var dot= mpnrw_fp_mul(p[i][0], p[j][0]) +
                mpnrw_fp_mul(p[i][1], p[j][1]) +
                mpnrw_fp_mul(p[i][2], p[j][2])
            var expected= (i === j)
                ? ONE
                : 0n
            var err= (dot > expected)
                ? dot - expected
                : expected - dot
            if (err > max_err) max_err= err
        }
    }
    assert.ok(max_err <= 64n, `Expected orthogonality err <= 64 ULPs, got ${max_err}`)
})

test('M_PNRW_BigInt parity with float64 M_PNRW', () => {
    var bi_engine= new M_PNRW_BigInt()
    var fp_engine= new M_PNRW()
    var cases= [
        { name: 'Greenwich J2000', lat: 51.477928, lon: 0, h: 46, date: new Date('2000-01-01T12:00:00Z') },
        { name: 'Stonehenge Equinox 2026', lat: 51.1789, lon: -1.8262, h: 103, date: new Date('2026-03-20T12:00:00Z') },
        { name: 'Sydney 2050', lat: -33.8688, lon: 151.2093, h: 19, date: new Date('2050-07-15T04:30:00Z') }
    ]

    for (var c of cases) {
        var fp_vec= fp_engine.get_universal_address(c.lat, c.lon, c.h, c.date, [0, 0, 0])
        var bi_vec= bi_engine.get_universal_address(
            mpnrw_fp_from_number(c.lat),
            mpnrw_fp_from_number(c.lon),
            mpnrw_fp_from_number(c.h),
            c.date,
            [0n, 0n, 0n]
        )
        var max_diff= Math.max(...bi_vec.map((v, idx) => Math.abs(mpnrw_fp_to_number(v) - fp_vec[idx])))
        assert.ok(max_diff < 1e-4, `${c.name}: maxDiff=${max_diff}`)
    }
})

test('M_PNRW_BigInt preserves sub-millimeter precision at 1 AU offset', () => {
    var engine= new M_PNRW_BigInt()
    var date= new Date('2026-03-20T12:00:00Z')
    var one_au_meters= mpnrw_fp_from_decimal('149597870700.0')

    var addr_zero= engine.get_universal_address(
        mpnrw_fp_from_decimal('51.5073'),
        mpnrw_fp_from_decimal('-0.1277'),
        mpnrw_fp_from_decimal('15'),
        date,
        [0n, 0n, 0n]
    )
    var addr_au= engine.get_universal_address(
        mpnrw_fp_from_decimal('51.5073'),
        mpnrw_fp_from_decimal('-0.1277'),
        mpnrw_fp_from_decimal('15'),
        date,
        [one_au_meters, 0n, 0n]
    )
    assert.equal(addr_au[0] - one_au_meters, addr_zero[0])
})

test('M_PNRW_BigInt sun vector parity and bit-exact determinism', () => {
    var bi_engine1= new M_PNRW_BigInt()
    var bi_engine2= new M_PNRW_BigInt()
    var fp_engine= new M_PNRW()
    var date= new Date('2026-03-20T12:00:00Z')

    var fp_sun= fp_engine.get_sun_vector(fp_engine.get_julian_date(date))
    var bi_sun= bi_engine1.get_sun_vector(bi_engine1.get_julian_date(date))
    var max_sun_diff= Math.max(...bi_sun.map((v, i) => Math.abs(mpnrw_fp_to_number(v) - fp_sun[i])))
    assert.ok(max_sun_diff < 1e-3, `maxSunDiff=${max_sun_diff}`)

    var addr1= bi_engine1.get_universal_address(
        mpnrw_fp_from_decimal('51.477928'),
        0n,
        mpnrw_fp_from_decimal('46'),
        new Date('2000-01-01T12:00:00Z'),
        [0n, 0n, 0n]
    )
    var addr2= bi_engine2.get_universal_address(
        mpnrw_fp_from_decimal('51.477928'),
        0n,
        mpnrw_fp_from_decimal('46'),
        new Date('2000-01-01T12:00:00Z'),
        [0n, 0n, 0n]
    )
    assert.deepEqual(addr1, addr2)
})

// SOFA validation suite t_c2t06a, MJD 53736 TT = UT1.
// The equatorial Greenwich site selects the first row of the published C2T matrix.
for (const fixed of [false, true]) {
    test('celestial reference, quarter turn, milliseconds and solar offset, fixed=' + fixed, () => {
        const e= fixed ? new M_PNRW_BigInt() : new M_PNRW()
        const encode= fixed ? mpnrw_fp_from_number : Number
        const decode= fixed ? mpnrw_fp_to_number : Number
        const zero= encode(0)
        const origin= [zero, zero, zero]
        const orientation= { xp: encode(2.55060238e-7), yp: encode(1.860359247e-6) }
        const at= ms => e.get_universal_address(zero, zero, zero, new Date(ms), origin, orientation).map(decode)
        const site= at(1136073600000)
        const expected= [-0.181033213, 0.983476981, 0.000065556].map(v => v * 6378137)
        assert.ok(Math.hypot(...site.map((v, i) => v - expected[i])) < 10)
        const later= at(1136073600001)
        const step= Math.hypot(...later.map((v, i) => v - site[i]))
        assert.ok(step > 0.44 && step < 0.49)
        const quarter= encode(2451545 + (6 - 18.697374558) / 24.06570982441908)
        const rotated= e.multiply_matrix_vector(e.get_r_matrix(quarter), [encode(1), zero, zero]).map(decode)
        assert.ok(Math.hypot(rotated[0], rotated[1] - 1, rotated[2]) < 1e-8)
        const jd= encode(2451545)
        const offset= e.get_heliocentric_earth_m(jd)
        assert.deepEqual(offset, e.get_sun_vector(jd).map(v => v * (fixed ? -1000n : -1000)))
    })
}
