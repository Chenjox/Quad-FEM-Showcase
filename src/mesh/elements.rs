use nalgebra::{Const, Dyn, OMatrix, OVector, SMatrix, SVector};

use crate::DIM;

pub trait ReferenceElement<const DIM: usize> {
    fn get_shape_functions(&self, ef_coordinates: SVector<f64, DIM>) -> OVector<f64, Dyn>;

    fn get_shape_function_derivatives(
        &self,
        ref_coordinates: SVector<f64, DIM>,
    ) -> OMatrix<f64, Dyn, Const<DIM>>;

    fn get_jacobian_determinant(
        &self,
        nodal_coordinates: &OMatrix<f64, Const<DIM>, Dyn>,
        ref_coordinates: SVector<f64, DIM>,
    ) -> SMatrix<f64, DIM, DIM> {
        let derivatives = self.get_shape_function_derivatives(ref_coordinates);

        let jacobian = nodal_coordinates * derivatives;
        jacobian
    }

    fn get_edge_coordinate(
        &self,
        edge_index: usize,
        orientation: usize,
        ref_coordinates: &OVector<f64, Dyn>,
    ) -> SVector<f64, DIM>;

    fn get_edge_nodes(&self, edge_index: usize, orientation: usize) -> OVector<usize, Dyn>;
}

#[derive(Debug)]
pub struct Quad4Element {}

fn g1(coord: f64) -> f64 {
    0.5 - coord / 2.0
}
fn g1_d(_coord: f64) -> f64 {
    -0.5
}
fn g2(coord: f64) -> f64 {
    0.5 + coord / 2.0
}
fn g2_d(_coord: f64) -> f64 {
    0.5
}

impl ReferenceElement<2> for Quad4Element {
    fn get_shape_functions(&self, ref_coordinates: SVector<f64, 2>) -> OVector<f64, Dyn> {
        let mut result = OVector::<f64, Dyn>::zeros(4);

        result[0] = g1(ref_coordinates[0]) * g1(ref_coordinates[1]);
        result[1] = g2(ref_coordinates[0]) * g1(ref_coordinates[1]);
        result[2] = g2(ref_coordinates[0]) * g2(ref_coordinates[1]);
        result[3] = g1(ref_coordinates[0]) * g2(ref_coordinates[1]);

        result
    }

    fn get_shape_function_derivatives(
        &self,
        ref_coordinates: SVector<f64, 2>,
    ) -> OMatrix<f64, Dyn, Const<2>> {
        let mut result = OMatrix::<f64, Dyn, Const<2>>::zeros(4);

        result[(0, 0)] = g1_d(ref_coordinates[0]) * g1(ref_coordinates[1]);
        result[(1, 0)] = g2_d(ref_coordinates[0]) * g1(ref_coordinates[1]);
        result[(2, 0)] = g2_d(ref_coordinates[0]) * g2(ref_coordinates[1]);
        result[(3, 0)] = g1_d(ref_coordinates[0]) * g2(ref_coordinates[1]);

        result[(0, 1)] = g1(ref_coordinates[0]) * g1_d(ref_coordinates[1]);
        result[(1, 1)] = g2(ref_coordinates[0]) * g1_d(ref_coordinates[1]);
        result[(2, 1)] = g2(ref_coordinates[0]) * g2_d(ref_coordinates[1]);
        result[(3, 1)] = g1(ref_coordinates[0]) * g2_d(ref_coordinates[1]);

        result
    }

    fn get_edge_nodes(&self, edge_index: usize, orientation: usize) -> OVector<usize, Dyn> {
        if orientation != 0 {
            println!("Orientation is not implemented yet");
        }
        match edge_index {
            0 => {
                // xi is coord, eta is -1
                return OVector::<usize, Dyn>::from_row_slice(&[0, 1]);
            }
            1 => {
                // eta is coord, xi is 1
                return OVector::<usize, Dyn>::from_row_slice(&[1, 2]);
            }
            2 => {
                // xi is coord, eta is 1
                return OVector::<usize, Dyn>::from_row_slice(&[2, 3]);
            }
            3 => {
                // eta is coord, xi is -1
                return OVector::<usize, Dyn>::from_row_slice(&[3, 0]);
            }
            _ => {
                panic!("Illegal Edge Index given!")
            }
        }
    }

    fn get_edge_coordinate(
        &self,
        edge_index: usize,
        orientation: usize,
        ref_coordinates: &OVector<f64, Dyn>,
    ) -> SVector<f64, 2> {
        if orientation != 0 {
            println!("Orientation is not implemented yet");
        }
        // Check whether given ref_coordinates is okay
        if ref_coordinates.len() >= DIM {
            panic!("This should not have happened");
        }
        match edge_index {
            0 => {
                // xi is coord, eta is -1
                return SVector::<f64, 2>::new(ref_coordinates[0], -1.0);
            }
            1 => {
                // eta is coord, xi is 1
                return SVector::<f64, 2>::new(1.0, ref_coordinates[0]);
            }
            2 => {
                // xi is coord, eta is 1
                return SVector::<f64, 2>::new(ref_coordinates[0], 1.0);
            }
            3 => {
                // eta is coord, xi is -1
                return SVector::<f64, 2>::new(-1.0, ref_coordinates[0]);
            }
            _ => {
                panic!("Illegal Edge Index given!")
            }
        }
    }
}
