use nalgebra::{Dyn, OMatrix, SVector};

type GenMatrixf64 = OMatrix<f64, Dyn, Dyn>;
type ConnectivityMatrix = OMatrix<usize, Dyn, Dyn>;
pub type Vertex<const DIM: usize> = SVector<f64, DIM>;

/// A Curve in [`DIM`] Dimensions, that starts somewhere and ends somewhere.
/// Is considered a cycle if start and end coincide.
/// should be parametrized by a parameter t
pub trait Curve<const DIM: usize> {
    /// The explizit parametrization of the Curve
    /// t should only be defined in the range of [0,1]
    fn get_vertex_at(&self, t: f64) -> Vertex<DIM>;

    /// The curve has to start somewhere, even in the parametrization
    fn get_starting_boundary_vertex(&self) -> Vertex<DIM>;

    /// the curve has to end somewhere, even in the parametrization
    fn get_ending_boundary_vertex(&self) -> Vertex<DIM>;

    /// Whether the Curve is a loop or not.
    fn is_unbounded(&self) -> bool;
}

/// A face in [`DIM`] Dimensions, that is bounded by curves
/// Is considered a cycle if it has no boundary
/// should be parametrized by two parameters t and s
pub trait Face<const DIM: usize> {
    /// The explizit parametrization of the Face
    /// t and s should only be defined in the range of [0,1]
    fn get_vertex_at(&self, t: f64, s: f64) -> Vertex<DIM>;

    /// The curve has to start somewhere, even in the parametrization
    fn get_boundary_curves(&self) -> Vec<Box<dyn Curve<DIM>>>;
}

/// A face in [`DIM`] Dimensions, that is bounded by curves
/// Is considered a cycle if it has no boundary
/// should be parametrized by two parameters t,s and r
pub trait Volume<const DIM: usize> {
    /// The explizit parametrization of the Face
    /// t,s and r should only be defined in the range of [0,1]
    fn get_vertex_at(&self, t: f64, s: f64, r: f64) -> Vertex<DIM>;

    /// The curve has to start somewhere, even in the parametrization
    fn get_boundary_faces(&self) -> Vec<Box<dyn Face<DIM>>>;
}

pub struct StraightEdge<const DIM: usize> {
    start_vertex: Vertex<DIM>,
    end_vertex: Vertex<DIM>,
}

impl<const DIM: usize> StraightEdge<DIM> {
    fn new(first: Vertex<DIM>, second: Vertex<DIM>) -> Self {
        Self {
            start_vertex: first,
            end_vertex: second,
        }
    }
}

impl<const DIM: usize> Curve<DIM> for StraightEdge<DIM> {
    fn get_vertex_at(&self, t: f64) -> Vertex<DIM> {
        return self.start_vertex * (1.0 - t) + self.end_vertex * t;
    }
    fn is_unbounded(&self) -> bool {
        false
    }
    fn get_starting_boundary_vertex(&self) -> Vertex<DIM> {
        self.start_vertex
    }
    fn get_ending_boundary_vertex(&self) -> Vertex<DIM> {
        self.end_vertex
    }
}

pub struct SecondOrderEdge<const DIM: usize> {
    start_vertex: Vertex<DIM>,
    middle_vertex: Vertex<DIM>,
    end_vertex: Vertex<DIM>,
}

impl<const DIM: usize> SecondOrderEdge<DIM> {
    /// erster knoten
    fn interpol1(t: f64) -> f64 {
        let t = t * 2. - 1.;
        0.5 * t * (t - 1.)
    }
    /// zweiter Knoten
    fn interpol2(t: f64) -> f64 {
        let t = t * 2. - 1.;
        (t + 1.) * (t - 1.)
    }

    /// dritter Knoten
    fn interpol3(t: f64) -> f64 {
        let t = t * 2. - 1.;
        0.5 * t * (t + 1.)
    }
}

impl<const DIM: usize> Curve<DIM> for SecondOrderEdge<DIM> {
    fn get_vertex_at(&self, t: f64) -> Vertex<DIM> {
        return self.start_vertex * SecondOrderEdge::<DIM>::interpol1(t)
            + self.middle_vertex * SecondOrderEdge::<DIM>::interpol2(t)
            + self.end_vertex * SecondOrderEdge::<DIM>::interpol3(t);
    }
    fn is_unbounded(&self) -> bool {
        false
    }
    fn get_starting_boundary_vertex(&self) -> Vertex<DIM> {
        self.start_vertex
    }
    fn get_ending_boundary_vertex(&self) -> Vertex<DIM> {
        self.end_vertex
    }
}

/// A Quadrilaterial is made of 4 points in space, which are arranged in counter clockwise order
pub struct StraightQuadrilaterial<const DIM: usize> {
    vertices: [Vertex<DIM>; 4],
}

impl<const DIM: usize> StraightQuadrilaterial<DIM> {
    fn get_lin_interpol0(t: f64, s: f64) -> f64 {
        (1. - t) * (1. - s)
    }

    fn get_lin_interpol1(t: f64, s: f64) -> f64 {
        t * (1. - s)
    }

    fn get_lin_interpol2(t: f64, s: f64) -> f64 {
        t * s
    }

    fn get_lin_interpol3(t: f64, s: f64) -> f64 {
        (1. - t) * s
    }
}

impl<const DIM: usize> Face<DIM> for StraightQuadrilaterial<DIM> {
    fn get_vertex_at(&self, t: f64, s: f64) -> Vertex<DIM> {
        self.vertices[0] * StraightQuadrilaterial::<DIM>::get_lin_interpol0(t, s)
            + self.vertices[1] * StraightQuadrilaterial::<DIM>::get_lin_interpol1(t, s)
            + self.vertices[2] * StraightQuadrilaterial::<DIM>::get_lin_interpol2(t, s)
            + self.vertices[3] * StraightQuadrilaterial::<DIM>::get_lin_interpol3(t, s)
    }

    fn get_boundary_curves(&self) -> Vec<Box<dyn Curve<DIM>>> {
        let mut result: Vec<Box<dyn Curve<DIM>>> = Vec::new();

        result.push(Box::new(StraightEdge::new(
            self.vertices[0],
            self.vertices[1],
        )));
        result.push(Box::new(StraightEdge::new(
            self.vertices[1],
            self.vertices[2],
        )));
        result.push(Box::new(StraightEdge::new(
            self.vertices[2],
            self.vertices[3],
        )));
        result.push(Box::new(StraightEdge::new(
            self.vertices[3],
            self.vertices[0],
        )));

        return result;
    }
}

struct Domain2D<const DIM: usize> {
    pub constituends: Vec<Box<dyn Face<DIM>>>,
    pub connectivity: ConnectivityMatrix,
}
