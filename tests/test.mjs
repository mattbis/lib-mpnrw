// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Matthew & Contributors — Free Software (GNU GPLv3)

// 1p
import test from 'node:test'
import assert from 'node:assert/strict'

// 2p
import { M_PNRW } from '../src-mjs.dev/M_PNRW.mjs'

// 3p
// (none)

test('M_PNRW full pipeline produces finite 3D vector', () => {
    var engine= new M_PNRW()
    var result= engine.get_universal_address(51.477928, 0, 0.046, new Date('2000-01-01T12:00:00Z'), [0, 0, 0])
    assert.equal(result.length, 3)
    assert.ok(result.every((v) => typeof v === 'number' && Number.isFinite(v)))
})

test('M_PNRW GMST sanity at J2000.0', () => {
    var engine= new M_PNRW()
    var r_mat= engine.get_r_matrix(2451545.0)
    var theta_rad= Math.atan2(r_mat[0][1], r_mat[0][0])
    if (theta_rad < 0) theta_rad= theta_rad + (2 * Math.PI)
    var gmst= theta_rad / (15 * (Math.PI / 180))
    assert.ok(Math.abs(gmst - 18.697) < 0.01)
})

test('M_PNRW precession matrix orthogonality', () => {
    var engine= new M_PNRW()
    var p= engine.get_p_matrix(0)
    for (let i= 0; i < 3; i= i + 1) {
        for (let j= 0; j < 3; j= j + 1) {
            var dot= 0
            for (let k= 0; k < 3; k= k + 1) {
                dot= dot + (p[i][k] * p[j][k])
            }
            var expected_val= (i === j)
                ? 1
                : 0
            assert.ok(Math.abs(dot - expected_val) <= 1e-6)
        }
    }
})

test('M_PNRW polar motion identity when xp=0, yp=0', () => {
    var engine= new M_PNRW()
    var w= engine.get_w_matrix(0, 0)
    var identity= [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for (let i= 0; i < 3; i= i + 1) {
        for (let j= 0; j < 3; j= j + 1) {
            assert.ok(Math.abs(w[i][j] - identity[i][j]) <= 1e-10)
        }
    }
})

test('M_PNRW mas_to_rad constant present', () => {
    var engine= new M_PNRW()
    assert.equal(typeof engine.mas_to_rad, 'number')
    assert.equal(typeof engine.ms_to_rad, 'undefined')
})
