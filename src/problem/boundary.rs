use nalgebra::{Const, Dyn, OMatrix, OVector, SVector};



pub trait WeakForm {
    fn num_dof_per_node(&self) -> usize;

    fn stiffness_term(
        &self,
        virt_node: usize,
        real_node: usize,
        element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OMatrix<f64, Dyn, Dyn>;

    fn right_hand_side_body(
        &self,
        virt_node: usize,
        element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OVector<f64, Dyn>;
}

pub trait NeumannBoundaryTerm {
    fn num_dof_per_node(&self) -> usize;

    fn eval_function(
        &self,
        boundary_marker: usize,
        current_global_coordinates: &SVector<f64, 2>,
        outward_normal_vector: &SVector<f64, 2>,
    ) -> OVector<f64, Dyn>;
}

pub trait DirichletBoundary {
    fn num_dof_per_node(&self) -> usize;

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool;

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64;
}