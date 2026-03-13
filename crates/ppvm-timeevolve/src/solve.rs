pub struct SolverConfig {
    pub rtol: f64,
    pub atol: f64,
    pub h0: Option<f64>,
    pub hmin: f64,
    pub hmax: f64,
}

impl Default for SolverConfig {
    fn default() -> Self {
        SolverConfig {
            rtol: 1e-6,
            atol: 1e-9,
            h0: None,
            hmin: 1e-12,
            hmax: f64::INFINITY,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solver_config_defaults() {
        let c = SolverConfig::default();
        assert_eq!(c.rtol, 1e-6);
        assert_eq!(c.atol, 1e-9);
        assert_eq!(c.h0, None);
        assert_eq!(c.hmin, 1e-12);
        assert_eq!(c.hmax, f64::INFINITY);
    }
}
