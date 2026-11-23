//! Integration tests mirroring `tests/test.mjs` and verifying parity with `src-mjs.dev/M_PNRW.mjs`.

use lib_mpnrw::{M_PNRW, MpnrwError, UtcDateTime};
use std::f64::consts::PI;

#[test]
fn test_full_pipeline_greenwich_j2000() {
    let engine = M_PNRW::new();
    let date = UtcDateTime::new(2000, 1, 1, 12, 0, 0).unwrap();
    let result = engine
        .get_universal_address(51.477928, 0.0, 0.046, date, [0.0, 0.0, 0.0])
        .unwrap();

    assert!(result.iter().all(|v| v.is_finite()));
    // Exact parity check against `node tests/test.mjs`
    assert!(
        (result[0] - 722_578.329_441).abs() < 1e-5,
        "X={}",
        result[0]
    );
    assert!(
        (result[1] - -3_914_552.214_023).abs() < 1e-5,
        "Y={}",
        result[1]
    );
    assert!(
        (result[2] - 4_966_742.813_298).abs() < 1e-5,
        "Z={}",
        result[2]
    );
}

#[test]
fn test_usage_parity_charing_cross() {
    let engine = M_PNRW::new();
    let date = UtcDateTime::new(2026, 3, 20, 12, 0, 0).unwrap();
    let result = engine
        .get_universal_address(51.5073, -0.1277, 15.0, date, [0.0, 0.0, 0.0])
        .unwrap();

    // Exact parity check against `node doc/usage.mjs`
    assert!((result[0] - 3_987_154.848).abs() < 1e-2, "X={}", result[0]);
    assert!((result[1] - -168_439.660).abs() < 1e-2, "Y={}", result[1]);
    assert!((result[2] - 4_958_694.494).abs() < 1e-2, "Z={}", result[2]);
}

#[test]
fn test_gmst_sanity_at_j2000() {
    let engine = M_PNRW::new();
    let r_mat = engine.get_r_matrix(2_451_545.0);

    let mut theta_rad = r_mat[1][0].atan2(r_mat[0][0]);
    if theta_rad < 0.0 {
        theta_rad += 2.0 * PI;
    }
    let gmst = theta_rad / (15.0 * (PI / 180.0));

    assert!((gmst - 18.697).abs() < 0.01, "gmst={gmst}");
}

#[test]
fn test_precession_matrix_orthogonality() {
    let engine = M_PNRW::new();
    for t in [0.0, 0.25, -0.5] {
        let p = engine.get_p_matrix(t);
        for (i, row_i) in p.iter().enumerate() {
            for (j, row_j) in p.iter().enumerate() {
                let dot: f64 = row_i.iter().zip(row_j.iter()).map(|(a, b)| a * b).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-6,
                    "t={t}, i={i}, j={j}, dot={dot}"
                );
            }
        }
    }
}

#[test]
fn test_polar_motion_identity() {
    let engine = M_PNRW::new();
    let w = engine.get_w_matrix(0.0, 0.0);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for i in 0..3 {
        for j in 0..3 {
            assert!((w[i][j] - identity[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn test_mas_to_rad_constant() {
    let engine = M_PNRW::new();
    assert!((engine.mas_to_rad - 4.848_136_811_095_36e-9).abs() < 1e-20);
}

#[test]
fn test_depack_vector() {
    let engine = M_PNRW::new();
    let mut buf = Vec::with_capacity(40);
    buf.extend_from_slice(&10i64.to_be_bytes());
    buf.extend_from_slice(&(-20i64).to_be_bytes());
    buf.extend_from_slice(&30i64.to_be_bytes());
    buf.extend_from_slice(&1_700_000_000i64.to_be_bytes());
    buf.extend_from_slice(&42i64.to_be_bytes());

    let depacked = engine.depack_vector(&buf).unwrap();
    assert_eq!(depacked.lat, 10);
    assert_eq!(depacked.lon, -20);
    assert_eq!(depacked.height, 30);
    assert_eq!(depacked.date, 1_700_000_000);
    assert_eq!(depacked.sun, 42);

    assert_eq!(
        engine.depack_vector(&buf[..32]),
        Err(MpnrwError::BufferTooShort {
            expected: 40,
            actual: 32
        })
    );
}

#[test]
fn test_unix_timestamp_conversion() {
    assert_eq!(
        UtcDateTime::from_unix_timestamp(0),
        UtcDateTime::new(1970, 1, 1, 0, 0, 0).unwrap()
    );
    assert_eq!(
        UtcDateTime::from_unix_timestamp(946_728_000),
        UtcDateTime::new(2000, 1, 1, 12, 0, 0).unwrap()
    );
}

#[test]
fn test_input_validation() {
    let engine = M_PNRW::new();
    let date = UtcDateTime::new(2000, 1, 1, 12, 0, 0).unwrap();

    assert_eq!(
        engine.get_universal_address(91.0, 0.0, 0.0, date, [0.0; 3]),
        Err(MpnrwError::InvalidLatitude)
    );
    assert_eq!(
        engine.get_universal_address(0.0, -181.0, 0.0, date, [0.0; 3]),
        Err(MpnrwError::InvalidLongitude)
    );
    assert_eq!(
        engine.get_universal_address(0.0, 0.0, -1.0, date, [0.0; 3]),
        Err(MpnrwError::InvalidHeight)
    );
}

#[test]
fn terrestrial_to_celestial_reference() {
    use lib_mpnrw::EarthOrientation;
    let e = M_PNRW::new();
    let orientation = EarthOrientation {
        xp: 2.55060238e-7,
        yp: 1.860359247e-6,
        ..Default::default()
    };
    // SOFA t_c2t06a, MJD 53736 TT = UT1. First row of C2T is column one of T2C.
    let expected = [-0.181033213, 0.983476981, 0.000065556].map(|v| v * 6378137.0);
    let site = e
        .get_universal_address_ms(0., 0., 0., 1136073600000, [0.; 3], orientation)
        .unwrap();
    assert!(e.get_distance(site, expected) < 10.);
    let later = e
        .get_universal_address_ms(0., 0., 0., 1136073600001, [0.; 3], orientation)
        .unwrap();
    assert!((0.44..0.49).contains(&e.get_distance(site, later)));
    let jd = 2451545.0;
    let offset = e.get_heliocentric_earth_m(jd);
    assert_eq!(offset, e.get_sun_vector(jd).map(|v| -1000. * v));
    let quarter = 2451545.0 + (6.0 - 18.697374558) / 24.06570982441908;
    let rotated = e.multiply_matrix_vector(e.get_r_matrix(quarter), [1., 0., 0.]);
    assert!(e.get_distance(rotated, [0., 1., 0.]) < 1e-8);
}
