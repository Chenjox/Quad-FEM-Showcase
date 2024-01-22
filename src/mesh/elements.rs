use nalgebra::{Const, Dyn, OMatrix, OVector, SMatrix, SVector};

pub trait ReferenceElement<const DIM: usize> {
    fn get_shape_functions(&self, ef_coordinates: SVector<f64, DIM>) -> OVector<f64, Dyn>;

    fn get_shape_function_derivatives(
        &self,
        ref_coordinates: SVector<f64, DIM>,
    ) -> OMatrix<f64, Dyn, Const<DIM>>;
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
}
