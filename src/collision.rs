pub enum Collision {
    None,
    Direct,
    Tree,
    Line,
    LineTree,
}

pub enum CollisionResolver {
    HardSphere,
    Halt,
    Merge,
}

pub enum ResolutionOutcome {
    RemoveNone,
    RemoveFirst,
    RemoveSecond,
    RemoveBoth,
}

impl CollisionResolver {
    pub fn resolve(&self) -> ResolutionOutcome {
        match self {
            CollisionResolver::HardSphere => {
                // Implement hard sphere collision resolution
                todo!()
            }
            CollisionResolver::Halt => {
                // Implement halt on collision
                todo!()
            }
            CollisionResolver::Merge => {
                // Implement merge on collision
                todo!()
            }
        }
    }
}
