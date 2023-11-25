use nalgebra::{DVector, Dyn, OMatrix, SVector};

type GenMatrixf64 = OMatrix<f64, Dyn, Dyn>;
type ConnectivityMatrix = OMatrix<usize, Dyn, Dyn>;
pub type Vertex<const DIM: usize> = SVector<f64, DIM>;

pub type Vertex2D = Vertex<2>;
pub type Vertex3D = Vertex<3>;

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
    fn new(start_vertex: Vertex<DIM>, middle_vertex: Vertex<DIM>, end_vertex: Vertex<DIM>) -> Self {
        Self {
            start_vertex,
            middle_vertex,
            end_vertex,
        }
    }
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

pub struct SecondOrderQuadrilaterial<const DIM: usize> {
    vertices: [Vertex<DIM>; 9],
}

impl<const DIM: usize> SecondOrderQuadrilaterial<DIM> {
    fn get_interpolant(index: usize, t: f64, s: f64) -> f64 {
        match index {
            0 => 0.5 * t * (t - 1.) * 0.5 * s * (s - 1.),
            1 => 0.5 * t * (t + 1.) * 0.5 * s * (s - 1.),
            2 => 0.5 * t * (t + 1.) * 0.5 * s * (s + 1.),
            3 => 0.5 * t * (t - 1.) * 0.5 * s * (s + 1.),
            4 => -(t + 1.) * (t - 1.) * 0.5 * s * (s - 1.),
            5 => -0.5 * t * (t + 1.) * (s + 1.) * (s - 1.),
            6 => -(t + 1.) * (t - 1.) * 0.5 * s * (s + 1.),
            7 => -0.5 * t * (t - 1.) * (s + 1.) * (s - 1.),
            8 => (t + 1.) * (t - 1.) * (s + 1.) * (s - 1.),
            _ => {
                panic!("Illegal Index given")
            }
        }
    }

    pub fn new_from_vector(vector: Vec<Vertex<DIM>>) -> Self {
        match vector.len() {
            0..=3 => {
                panic!("Insufficient Points given")
            }
            4 => {
                let g = StraightQuadrilaterial {
                    vertices: [vector[0], vector[1], vector[2], vector[3]],
                };
                let mut vertices = [Vertex::<DIM>::zeros(); 9];
                for (index, (xi, eta)) in [
                    (0.0, 0.0),
                    (1.0, 0.0),
                    (1.0, 1.0),
                    (0.0, 1.0),
                    (0.5, 0.0),
                    (1.0, 0.5),
                    (0.5, 1.0),
                    (0.0, 0.5),
                    (0.5, 0.5),
                ]
                .iter()
                .enumerate()
                {
                    vertices[index] = g.get_vertex_at(*xi, *eta);
                }
                Self { vertices }
            }
            5..=9 => {
                todo!();
            }
            _ => {
                panic!("Illegal Length given!")
            }
        }
    }
}

impl<const DIM: usize> Face<DIM> for SecondOrderQuadrilaterial<DIM> {
    fn get_vertex_at(&self, t: f64, s: f64) -> Vertex<DIM> {
        let mut res = Vertex::<DIM>::zeros();
        for i in 0..9 {
            res += self.vertices[i] * SecondOrderQuadrilaterial::<DIM>::get_interpolant(i, t, s);
        }
        res
    }

    fn get_boundary_curves(&self) -> Vec<Box<dyn Curve<DIM>>> {
        let mut result: Vec<Box<dyn Curve<DIM>>> = Vec::new();

        result.push(Box::new(SecondOrderEdge::new(
            self.vertices[0],
            self.vertices[4],
            self.vertices[1],
        )));
        result.push(Box::new(SecondOrderEdge::new(
            self.vertices[1],
            self.vertices[5],
            self.vertices[2],
        )));
        result.push(Box::new(SecondOrderEdge::new(
            self.vertices[2],
            self.vertices[6],
            self.vertices[3],
        )));
        result.push(Box::new(SecondOrderEdge::new(
            self.vertices[3],
            self.vertices[7],
            self.vertices[0],
        )));

        return result;
    }
}

struct Domain2D<const DIM: usize> {
    pub constituends: Vec<Box<dyn Face<DIM>>>,
    pub connectivity: ConnectivityMatrix,
}

impl<const DIM: usize> Domain2D<DIM> {
    pub fn get_FE_mesh(&self, n: usize, m: usize) {
        let m_f = 1. / (m as f64);
        let n_f = 1. / (n as f64);

        let mut xi_vec = DVector::<f64>::zeros(n);
        let mut eta_vec = DVector::<f64>::zeros(m);

        let mut cumulative = 0.0;
        for i in 0..n {
            xi_vec[i] = cumulative;
            cumulative += n_f;
        }
        cumulative = 0.0;
        for i in 0..m {
            eta_vec[i] = cumulative;
            cumulative += m_f;
        }
    }
}
