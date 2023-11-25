/// Jeder DOF ist definiert über eine Liste von Vertizes
pub enum DOFObject {
    VERTEX(usize),
    EDGE(Vec<usize>),
    FACE(Vec<usize>),
}

pub trait FEMesh {
    fn get_elements();
}
