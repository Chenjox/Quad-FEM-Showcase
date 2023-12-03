use nalgebra::Dyn;
use nalgebra::SVector;
use nalgebra::OMatrix;


pub type Point<const DIM: usize> = SVector<f64, DIM>;

pub struct LinearLagrangeGeom<const DIM: usize> {
  
}


impl<const DIM: usize> LinearLagrangeGeom<DIM> {
  pub fn get_basis_functions(point: Point<DIM>) -> OMatrix<f64,Dyn,Dyn> {
    todo!()
  }
}