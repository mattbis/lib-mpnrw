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

import { M_PNRW } from '../src.mjs/M_PNRW.mjs'


async function runTest() {
    console.log('--- M_PNRW Engine Test ---')

    var engine= new M_PNRW()
    var passed= 0
    var failed= 0

    // --- SHARED TEST INPUTS ---
    var lat=        51.477928   // Greenwich Observatory
    var lon=        0
    var height=     0.046       // km (46m)
    var date=       new Date('2000-01-01T12:00:00Z')
    var sun_vector= [0, 0, 0]   // Centered at barycenter for simplicity

    console.log(`\nCoordinates : Lat ${lat}, Lon ${lon}, Height ${height}km`)
    console.log(`Date        : ${date.toISOString()}`)


    // --- TEST 1: Full pipeline output ---
    console.log('\n[Test 1] Full pipeline (get_universal_address)')
    try {
        var result= engine.get_universal_address(lat, lon, height, date, sun_vector)

        console.log('  X:', result[0].toFixed(6))
        console.log('  Y:', result[1].toFixed(6))
        console.log('  Z:', result[2].toFixed(6))

        if (result.length === 3 && result.every(v => typeof v === 'number' && isFinite(v))) {
            console.log('  Result: PASS')
            passed++
        } else {
            console.log('  Result: FAIL (unexpected format or non-finite values)')
            failed++
        }
    } catch (error) {
        console.error('  Result: FAIL (threw)', error)
        failed++
    }


    // --- TEST 2: GMST sanity check at J2000.0 ---
    // At J2000.0 (JD 2451545.0), GMST should be ~18.697h (within 0.01h tolerance)
    console.log('\n[Test 2] GMST sanity at J2000.0 (get_r_matrix)')
    try {
        var jd_j2000= 2451545.0
        var r_mat = engine.get_r_matrix(jd_j2000)
        // Extract theta from R matriz: R[0][0] = cos(theta), R[0][1] = sin(theta)
        var theta_rad = Math.atan2(r_mat[0][1], r_mat[0][0])
        if (theta_rad < 0) theta_rad += 2 * Math.PI
        var gmst = theta_rad / (15 * (Math.PI / 180)) // convert back to hours
        
        var expected= 18.697
        var tolerance= 0.01

        console.log('  GMST (hours):', gmst.toFixed(6))
        console.log('  Expected ~  :', expected)

        if (Math.abs(gmst - expected) < tolerance) {
            console.log('  Result: PASS')
            passed++
        } else {
            console.log('  Result: FAIL (outside tolerance)')
            failed++
        }
    } catch (error) {
        console.error('  Result: FAIL (threw)', error)
        failed++
    }


    // --- TEST 3: Precession matrix orthogonality (get_p_matrix) ---
    // A valid rotation matrix M should satisfy M * M^T = I
    console.log('\n[Test 3] Precession matrix orthogonality (get_p_matrix)')
    try {
        var t= 0
        var p= engine.get_p_matrix(t)
        var ok= true

        for (var i= 0; i < 3; i++) {
            for (var j= 0; j < 3; j++) {
                var dot= 0
                for (var k= 0; k < 3; k++) dot += p[i][k] * p[j][k]
                var expected_val= (i === j) ? 1 : 0
                if (Math.abs(dot - expected_val) > 1e-6) ok= false
            }
        }

        console.log('  M * M^T = I:', ok)
        if (ok) {
            console.log('  Result: PASS')
            passed++
        } else {
            console.log('  Result: FAIL (not orthogonal)')
            failed++
        }
    } catch (error) {
        console.error('  Result: FAIL (threw)', error)
        failed++
    }


    // --- TEST 4: Polar motion passthrough (get_w_matrix) ---
    // With xp=0, yp=0 it should return identity
    console.log('\n[Test 4] Polar motion identity when xp=yp=0 (get_w_matrix)')
    try {
        var w=       engine.get_w_matrix(0, 0)
        var identity= [[1,0,0],[0,1,0],[0,0,1]]
        var match=   true

        for (var i= 0; i < 3; i++) {
            for (var j= 0; j < 3; j++) {
                if (Math.abs(w[i][j] - identity[i][j]) > 1e-10) match= false
            }
        }

        if (match) {
            console.log('  Result: PASS')
            passed++
        } else {
            console.log('  Result: FAIL (unexpected values)', w)
            failed++
        }
    } catch (error) {
        console.error('  Result: FAIL (threw)', error)
        failed++
    }


    // --- TEST 5: mas_to_rad rename check ---
    // Verify the old ms_to_rad is gone and mas_to_rad is present
    console.log('\n[Test 5] mas_to_rad constant present (rename check)')
    if (typeof engine.mas_to_rad === 'number' && typeof engine.ms_to_rad === 'undefined') {
        console.log('  mas_to_rad:', engine.mas_to_rad)
        console.log('  Result: PASS')
        passed++
    } else {
        console.log('  Result: FAIL (mas_to_rad missing or ms_to_rad still present)')
        failed++
    }


    // --- SUMMARY ---
    console.log('\n--- Summary ---')
    console.log(`Passed: ${passed} / ${passed + failed}`)
    if (failed > 0) process.exit(1)
}

runTest()