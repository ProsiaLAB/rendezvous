#[cfg(feature = "quadrupole")]
pub type TreeParam = Quadrupole;

#[cfg(not(feature = "quadrupole"))]
pub type TreeParam = ();

pub type TreeType = Tree<TreeParam>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeId(pub usize);

pub struct Node<M = ()> {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
    pub m: f64,
    pub mx: f64,
    pub my: f64,
    pub mz: f64,
    pub moment: M,
    pub children: [Option<NodeId>; 8],
    pub kind: NodeKind,
    pub remote: bool,
}

pub enum NodeKind {
    Leaf {
        particle: usize,
    },
    Internal {
        count: usize, // total particles in subtree
    },
}

pub struct Tree<M> {
    nodes: Vec<Node<M>>,
}

impl<M> Tree<M> {
    pub fn alloc() -> Self {
        Self { nodes: Vec::new() }
    }

    pub fn new(size: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(size),
        }
    }

    pub fn get_node(&self, id: NodeId) -> &Node<M> {
        &self.nodes[id.0]
    }

    pub fn get_node_mut(&mut self, id: NodeId) -> &mut Node<M> {
        &mut self.nodes[id.0]
    }

    pub fn size(&self) -> usize {
        self.nodes.len()
    }
}

#[cfg(feature = "quadrupole")]
pub struct Quadrupole {
    pub mxx: f64,
    pub mxy: f64,
    pub mxz: f64,
    pub myy: f64,
    pub myz: f64,
    pub mzz: f64,
}
