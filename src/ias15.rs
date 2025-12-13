pub struct Ias15 {
    pub epsilon: f64,
    pub min_dt: f64,
    pub adaptive_mode: AdaptiveMode,
    pub iterations_max_exceeded: usize,
    pub n_allocated: usize,
}

pub enum AdaptiveMode {
    Individual,
    Global,
    Prs23,
    Aarseth85,
}
