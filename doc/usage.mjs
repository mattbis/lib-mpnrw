// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Matthew & Contributors — Free Software (GNU GPLv3)

// 1p
// (none)

// 2p
import { M_PNRW } from '../src-mjs.dev/M_PNRW.mjs'
import { M_PNRW_BigInt } from '../src-mjs.dev/M_PNRW_BigInt.mjs'

// 3p
// (none)

var fp_engine= new M_PNRW()
var bi_engine= new M_PNRW_BigInt()

function display_universal_address(name, lat_str, lon_str, height_str) {
    var date= new Date('2026-03-20T12:00:00Z')
    var lat_num= Number(lat_str)
    var lon_num= Number(lon_str)
    var height_num= Number(height_str)

    var fp_addr= fp_engine.get_universal_address(lat_num, lon_num, height_num, date, [0, 0, 0])
    var bi_addr= bi_engine.get_universal_address(
        M_PNRW_BigInt.fromDecimal(lat_str),
        M_PNRW_BigInt.fromDecimal(lon_str),
        M_PNRW_BigInt.fromDecimal(height_str),
        date,
        [0n, 0n, 0n]
    )

    console.log(`[ ${name} ]`)
    console.log(`Geodetic Coordinates  : Lat ${lat_num.toFixed(4)}, Lon ${lon_num.toFixed(4)}, Height ${height_num}m`)
    console.log(`Universal Space Vector (f64 vs BigInt Q.80):`)
    console.log(`  X: ${fp_addr[0].toFixed(3)} m | ${M_PNRW_BigInt.toDecimal(bi_addr[0], 6)} m`)
    console.log(`  Y: ${fp_addr[1].toFixed(3)} m | ${M_PNRW_BigInt.toDecimal(bi_addr[1], 6)} m`)
    console.log(`  Z: ${fp_addr[2].toFixed(3)} m | ${M_PNRW_BigInt.toDecimal(bi_addr[2], 6)} m\n`)
}

console.log('=========================================')
console.log(' M_PNRW Celestial Transformation Usages  ')
console.log('=========================================\n')

display_universal_address('Charing Cross (London Center)', '51.5073', '-0.1277', '15')
display_universal_address('The London Stone (Roman Measure)', '51.5115', '-0.0898', '15')
display_universal_address('Stonehenge', '51.1789', '-1.8262', '103')
display_universal_address('Maeshowe Runic Inscriptions', '58.9965', '-3.1884', '20')
display_universal_address('Jelling Runestones (Denmark)', '55.7566', '9.4194', '50')
