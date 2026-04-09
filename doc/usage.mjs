import { M_PNRW } from '../src.mjs/M_PNRW.mjs'

var engine = new M_PNRW()

/**
 * A helper function to print the Universal Address for a specific location.
 * Uses a fixed date (Vernal Equinox 2026) and a simplified Sun vector for deterministic results.
 * In practice, the Sun vector would represent the ICRS coordinates of the Sun or Barycenter.
 * 
 * @param {string} name - The name of the location.
 * @param {number} lat - The latitude of the location in degrees.
 * @param {number} lon - The longitude of the location in degrees.
 * @param {number} height - The height of the location in meters.
 * @returns {void}
 */
function displayUniversalAddress(name, lat, lon, height) {
    // A fixed date (e.g. Vernal Equinox 2026)
    var date = new Date('2026-03-20T12:00:00Z')
    
    // For demonstration, assuming a 0-offset from the Solar System Barycenter.
    // In a real scenario, this would be the actual position of the Sun or Barycenter in ICRS.
    var sun_vector = [0, 0, 0] 

    try {
        var address = engine.get_universal_address(lat, lon, height, date, sun_vector)
        
        console.log(`[ ${name} ]`)
        console.log(`Geodetic Coordinates  : Lat ${lat.toFixed(4)}, Lon ${lon.toFixed(4)}, Height ${height}m`)
        console.log(`Universal Space Vector:`)
        console.log(`  X: ${address[0].toFixed(3)} m`)
        console.log(`  Y: ${address[1].toFixed(3)} m`)
        console.log(`  Z: ${address[2].toFixed(3)} m\n`)
    } catch (e) {
        console.error(`Failed to calculate for ${name}: ${e.message}\n`)
    }
}


console.log("=========================================")
console.log(" M_PNRW Celestial Transformation Usages  ")
console.log("=========================================\n")


// 1. Charing Cross, London 
// (Considered the exact center of London from which distances are measured)
displayUniversalAddress("Charing Cross (London Center)", 51.5073, -0.1277, 15)

// 2. The London Stone
// (An ancient Roman stone from which distances in Roman Britain were supposedly measured)
displayUniversalAddress("The London Stone (Roman Measure)", 51.5115, -0.0898, 15)

// 3. Stonehenge, Wiltshire
// (Mysterious prehistoric standing stones)
displayUniversalAddress("Stonehenge", 51.1789, -1.8262, 103)

// 4. Maeshowe, Orkney, Scotland
// (Contains the largest collection of runic inscriptions in the world)
displayUniversalAddress("Maeshowe Runic Inscriptions", 58.9965, -3.1884, 20)

// 5. The Jelling Runestones (Denmark)
// (Famous runestones from the Viking Age. Going a bit outside the UK for runestones!)
displayUniversalAddress("Jelling Runestones (Denmark)", 55.7566, 9.4194, 50)
