use nalgebra::{Const, Dyn, OMatrix, OVector, SVector};

use crate::mesh::elements::ReferenceElement;

use super::integration::get_gauss_rule;

pub trait LocalStiffnessAssembler {
    fn num_dof_per_node(&self) -> usize;

    fn assemble_local_stiffness_term(&self, 
        ref_element: &Box<dyn ReferenceElement<2>>, 
        element_nodes: &[usize],
        nodal_coordinates: &OMatrix<f64, Const<2>, Dyn>,
        local_stiffness_matrix: &mut OMatrix<f64,Dyn,Dyn>,
        local_rhs_vector: &mut OVector<f64,Dyn>,
    );
}

pub struct WeakFormAssembler<W: WeakForm> {
    pub weak_form: W
}

impl<W: WeakForm> LocalStiffnessAssembler for WeakFormAssembler<W> {
    fn num_dof_per_node(&self) -> usize {
        self.weak_form.num_dof_per_node()
    }

    fn assemble_local_stiffness_term(&self, 
        ref_element: &Box<dyn ReferenceElement<2>>, 
        element_nodes: &[usize],
        nodal_coordinates: &OMatrix<f64, Const<2>, Dyn>,
        local_stiffness_matrix: &mut OMatrix<f64,Dyn,Dyn>,
        local_rhs_vector: &mut OVector<f64,Dyn>,
    ){
        let num_dof_per_node = self.num_dof_per_node();
        let gauss = get_gauss_rule(2);

        for gauss_point in gauss.column_iter() {
            let xi_1 = gauss_point[0];
            let xi_2 = gauss_point[1];
            let weight = gauss_point[2];

            let ref_coordinates = SVector::<f64, 2>::new(xi_1, xi_2);

            //let normal_func = ref_element.get_shape_functions(ref_coordinates);

            //println!("{}",normal_func);

            let shape_functions = ref_element.get_shape_functions(ref_coordinates);
            let derivatives = ref_element.get_shape_function_derivatives(ref_coordinates);
            let jacobian =
                ref_element.get_jacobian_determinant(&nodal_coordinates, ref_coordinates);
            let inv_jacobian = jacobian.try_inverse().unwrap();
            let determinant =
                jacobian[(0, 0)] * jacobian[(1, 1)] - jacobian[(0, 1)] * jacobian[(1, 0)];

            //println!("{}",jacobian);

            // Transformation auf tatsächliche Elemente
            let derivatives = derivatives * inv_jacobian;

            for j in element_nodes.iter().enumerate() {
                let virt_node_number = j.0;
                for i in element_nodes.iter().enumerate() {
                    let real_node_number = i.0;

                    let result = weight
                        * self.weak_form.stiffness_term(
                            virt_node_number,
                            real_node_number,
                            element_nodes,
                            &shape_functions,
                            &derivatives,
                        )
                        * determinant;
                    // Offset
                    for i in 0..(num_dof_per_node * num_dof_per_node) {
                        let i_k = i % num_dof_per_node;
                        let j_k = i / num_dof_per_node;

                        local_stiffness_matrix[(
                            real_node_number * num_dof_per_node + i_k,
                            virt_node_number * num_dof_per_node + j_k,
                        )] += result[(i_k, j_k)];
                    }
                }

                let rhs_term = self.weak_form.right_hand_side_body(
                    virt_node_number,
                    &element_nodes,
                    &shape_functions,
                    &derivatives,
                );
                for i in 0..num_dof_per_node {
                    local_rhs_vector[virt_node_number * num_dof_per_node + i] +=
                        weight * rhs_term[i] * determinant;
                }
            }
        }
        println!("{:3.3}",local_stiffness_matrix);
    }
}

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