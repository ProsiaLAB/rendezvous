use crate::{
    eos::Eos, gbs::Gbs, ias15::Ias15, janus::Janus, leapfrog::LeapFrog, mercurius::Mercurius,
    saba::Saba, sei::Sei, trace::Trace, whfast::WHFast,
};

pub enum Integrator {
    Ias15(Ias15),
    WHFast(WHFast),
    LeapFrog(LeapFrog),
    Gbs(Gbs),
    Sei(Sei),
    Janus(Janus),
    Mercurius(Mercurius),
    Saba(Saba),
    Eos(Eos),
    Trace(Trace),
    None,
}
