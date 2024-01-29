use nalgebra::{Const, Dyn, OMatrix, OVector, Rotation2, SMatrix, SVector};

use crate::{mesh::elements::ReferenceElement, problem::{boundary::{DirichletBoundary, LocalStiffnessAssembler, WeakForm}, integration::get_gauss_rule}};


pub struct ElasticityMixed {
    pub youngs_modulus: f64,
    pub poissons_ratio: f64
}

impl ElasticityMixed {

    fn get_shear_and_compressive(&self) -> (f64,f64) {
        let compression_modulus = self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poissons_ratio));
        let second_lame = self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio));

        return (second_lame, compression_modulus);
    }

    fn get_compressive_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poissons_ratio))
    }

    fn get_shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio))
    }

    fn first_weak_form(&self,
        virt_node: usize,
        real_node: usize,
        _element_nodes: &[usize],
        _shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OMatrix<f64, Dyn, Dyn>{
        let b_mat_i = {
            let mut result = SMatrix::<f64, 3, 2>::zeros();
            result[(0, 0)] = shape_derivatives[(virt_node, 0)]; // N1x
            result[(1, 1)] = shape_derivatives[(virt_node, 1)]; // N1y
            result[(2, 0)] = shape_derivatives[(virt_node, 1)]; // N1y
            result[(2, 1)] = shape_derivatives[(virt_node, 0)]; // N1x
            result
        };

        let b_mat_i = b_mat_i - SMatrix::<f64,3,3>::new(1.,1.,0.,1.,1.,0.,0.,0.,0.) * b_mat_i;

        let b_mat_j = {
            let mut result = SMatrix::<f64, 3, 2>::zeros();
            result[(0, 0)] = shape_derivatives[(real_node, 0)]; // N1x
            result[(1, 1)] = shape_derivatives[(real_node, 1)]; // N1y
            result[(2, 0)] = shape_derivatives[(real_node, 1)]; // N1y
            result[(2, 1)] = shape_derivatives[(real_node, 0)]; // N1x
            result
        };
        let b_mat_j = b_mat_j - SMatrix::<f64,3,3>::new(1.,1.,0.,1.,1.,0.,0.,0.,0.) * b_mat_j;


        let result = b_mat_i.transpose() * 2.0 * self.get_shear_modulus() * b_mat_j;
        //println!("{}",result);

        let mut resulting = OMatrix::<f64, Dyn, Dyn>::zeros(2, 2);
        for i in 0..(2 * 2) {
            let i_k = i % 2;
            let j_k = i / 2;

            resulting[(i_k, j_k)] = result[(i_k, j_k)];
        }
        resulting
    }

    fn second_weak_form(&self,
    virt_node: usize,
    real_node: usize,
    _element_nodes: &[usize],
    _shape_functions: &OMatrix<f64, Dyn, Const<1>>,
    shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
) -> (OMatrix<f64, Const<2>, Const<1>>,OMatrix<f64, Const<1>, Const<2>>) {
    let b_mat_i = {
        let mut result = SMatrix::<f64, 3, 2>::zeros();
        result[(0, 0)] = shape_derivatives[(virt_node, 0)]; // N1x
        result[(1, 1)] = shape_derivatives[(virt_node, 1)]; // N1y
        result[(2, 0)] = shape_derivatives[(virt_node, 1)]; // N1y
        result[(2, 1)] = shape_derivatives[(virt_node, 0)]; // N1x
        result
    };
    let discrete_div_i = SMatrix::<f64,1,3>::new(1.,1.,0.) * b_mat_i;

    let b_mat_j = {
        let mut result = SMatrix::<f64, 3, 2>::zeros();
        result[(0, 0)] = shape_derivatives[(real_node, 0)]; // N1x
        result[(1, 1)] = shape_derivatives[(real_node, 1)]; // N1y
        result[(2, 0)] = shape_derivatives[(real_node, 1)]; // N1y
        result[(2, 1)] = shape_derivatives[(real_node, 0)]; // N1x
        result
    };
    let discrete_div_j = SMatrix::<f64,1,3>::new(1.,1.,0.) * b_mat_j;

    return (discrete_div_j.transpose(), discrete_div_i);
}
}

impl LocalStiffnessAssembler for ElasticityMixed {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn assemble_local_stiffness_term(&self, 
        ref_element: &Box<dyn ReferenceElement<2>>, 
        element_nodes: &[usize],
        nodal_coordinates: &OMatrix<f64, Const<2>, Dyn>,
        local_stiffness_matrix: &mut OMatrix<f64,Dyn,Dyn>,
        local_rhs_vector: &mut OVector<f64,Dyn>,
    ) {
        let num_dof_per_node = self.num_dof_per_node();
        let gauss = get_gauss_rule(2);

        let mut pressure_correction_left = OMatrix::<f64,Dyn,Dyn>::zeros(local_stiffness_matrix.nrows(),1);
        let mut pressure_correction_right = OMatrix::<f64,Dyn,Dyn>::zeros(1,local_stiffness_matrix.ncols());

        let mut area = 0.;
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
            area += weight * determinant;

            // Transformation auf tatsächliche Elemente
            let derivatives = derivatives * inv_jacobian;

            for j in element_nodes.iter().enumerate() {
                let virt_node_number = j.0;
                for i in element_nodes.iter().enumerate() {
                    let real_node_number = i.0;

                    let result = weight
                        * self.first_weak_form(
                            virt_node_number,
                            real_node_number,
                            element_nodes,
                            &shape_functions,
                            &derivatives,
                        )
                        * determinant;

                    let (pressure_left,pressure_right) = self.second_weak_form(
                        virt_node_number,
                        real_node_number,
                        element_nodes,
                        &shape_functions,
                        &derivatives,
                    );
                    // Offset
                    for i in 0..(num_dof_per_node * num_dof_per_node) {
                        let i_k = i % num_dof_per_node;
                        let j_k = i / num_dof_per_node;

                        local_stiffness_matrix[(
                            real_node_number * num_dof_per_node + i_k,
                            virt_node_number * num_dof_per_node + j_k,
                        )] += result[(i_k, j_k)];

                        pressure_correction_left[real_node_number * num_dof_per_node + i_k] += weight*pressure_left[i_k]*determinant;
                        pressure_correction_right[real_node_number * num_dof_per_node + j_k] += weight*pressure_right[j_k]*determinant;
                    }
                }

                let rhs_term = SVector::<f64,2>::zeros();
                for i in 0..num_dof_per_node {
                    local_rhs_vector[virt_node_number * num_dof_per_node + i] +=
                        weight * rhs_term[i] * determinant;
                }
            }
        }// Gausspunkt ende
        let pressure_correction = self.get_compressive_modulus()/area * pressure_correction_left * pressure_correction_right;


        for i in 0..local_stiffness_matrix.ncols() {
            for j in 0..local_stiffness_matrix.nrows() {
                local_stiffness_matrix[(i,j)] += pressure_correction[(i,j)];
            }
        }
    }
}

pub struct Elasticity {
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
}

impl Elasticity {
    fn material_matrix(&self) -> SMatrix<f64, 3, 3> {
        let first_lame = self.youngs_modulus * self.poissons_ratio
            / ((1.0 + self.poissons_ratio) * (1.0 - 2.0 * self.poissons_ratio));
        let second_lame = self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio));
        return SMatrix::<f64, 3, 3>::new(
            2.0 * second_lame + first_lame,
            first_lame,
            0.0,
            first_lame,
            2.0 * second_lame + first_lame,
            0.0,
            0.0,
            0.0,
            second_lame,
        );
    }
}

impl WeakForm for Elasticity {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn stiffness_term(
        &self,
        virt_node: usize,
        real_node: usize,
        _element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OMatrix<f64, Dyn, Dyn> {
        let b_mat_i = {
            let mut result = SMatrix::<f64, 3, 2>::zeros();
            result[(0, 0)] = shape_derivatives[(virt_node, 0)]; // N1x
            result[(1, 1)] = shape_derivatives[(virt_node, 1)]; // N1y
            result[(2, 0)] = shape_derivatives[(virt_node, 1)]; // N1y
            result[(2, 1)] = shape_derivatives[(virt_node, 0)]; // N1x
            result
        };
        let b_mat_j = {
            let mut result = SMatrix::<f64, 3, 2>::zeros();
            result[(0, 0)] = shape_derivatives[(real_node, 0)]; // N1x
            result[(1, 1)] = shape_derivatives[(real_node, 1)]; // N1y
            result[(2, 0)] = shape_derivatives[(real_node, 1)]; // N1y
            result[(2, 1)] = shape_derivatives[(real_node, 0)]; // N1x
            result
        };

        let result = b_mat_j.transpose() * self.material_matrix() * b_mat_i;
        //println!("{},{},{}",virt_node,real_node,result);

        let mut resulting = OMatrix::<f64, Dyn, Dyn>::zeros(2, 2);
        for i in 0..(2 * 2) {
            let i_k = i % 2;
            let j_k = i / 2;

            resulting[(i_k, j_k)] = result[(i_k, j_k)];
        }
        resulting
    }

    fn right_hand_side_body(
        &self,
        _virt_node: usize,
        _element_nodes: &[usize],
        _shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        _shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OVector<f64, Dyn> {
        return OVector::<f64, Dyn>::zeros(2);
    }
}

pub struct StressSmootherRHS {
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub solution_vec: OVector<f64, Dyn>,
    pub num_nodes_solution: usize,
    pub num_dofs_per_node: usize,
}

pub struct StressSmootherStiffness {
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub solution_vec: OVector<f64, Dyn>,
    pub num_nodes_solution: usize,
    pub num_dofs_per_node: usize,
}

impl StressSmootherRHS {
    fn material_matrix(&self) -> SMatrix<f64, 3, 3> {
        let first_lame = self.youngs_modulus * self.poissons_ratio
            / ((1.0 + self.poissons_ratio) * (1.0 - 2.0 * self.poissons_ratio));
        let second_lame = self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio));
        return SMatrix::<f64, 3, 3>::new(
            2.0 * second_lame + first_lame,
            first_lame,
            0.0,
            first_lame,
            2.0 * second_lame + first_lame,
            0.0,
            0.0,
            0.0,
            second_lame,
        );
    }
}

impl WeakForm for StressSmootherStiffness {
    fn num_dof_per_node(&self) -> usize {
        1 // sigma_x, sigma_y, sigma_xy
    }

    fn stiffness_term(
        &self,
        virt_node: usize,
        real_node: usize,
        _element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OMatrix<f64, Dyn, Dyn> {
        let shape_1 = OVector::<f64, Dyn>::from_element(1, shape_functions[virt_node]);
        let shape_2 = OVector::<f64, Dyn>::from_element(1, shape_functions[real_node]);

        let result = shape_1 * shape_2.transpose();

        return result;
    }

    fn right_hand_side_body(
        &self,
        _virt_node: usize,
        _element_nodes: &[usize],
        _shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        _shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OVector<f64, Dyn> {
        return OVector::<f64, Dyn>::zeros(1);
    }
}

impl WeakForm for StressSmootherRHS {
    fn num_dof_per_node(&self) -> usize {
        3 // sigma_x, sigma_y, sigma_xy
    }

    fn stiffness_term(
        &self,
        virt_node: usize,
        real_node: usize,
        _element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        _shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OMatrix<f64, Dyn, Dyn> {
        let shape_1 = OVector::<f64, Dyn>::from_element(3, shape_functions[virt_node]);
        let shape_2 = OVector::<f64, Dyn>::from_element(3, shape_functions[real_node]);

        let result = shape_1 * shape_2.transpose();

        return result;
    }

    fn right_hand_side_body(
        &self,
        virt_node: usize,
        element_nodes: &[usize],
        shape_functions: &OMatrix<f64, Dyn, Const<1>>,
        shape_derivatives: &OMatrix<f64, Dyn, Const<2>>,
    ) -> OVector<f64, Dyn> {
        let mut residual = OVector::<f64, Dyn>::zeros(3);

        for i in 0..shape_derivatives.nrows() {
            let b_mat_j = {
                // alle knoten durchlaufen!
                let mut result = SMatrix::<f64, 3, 2>::zeros();
                result[(0, 0)] = shape_derivatives[(i, 0)]; // N1x
                result[(1, 1)] = shape_derivatives[(i, 1)]; // N1y
                result[(2, 0)] = shape_derivatives[(i, 1)]; // N1y
                result[(2, 1)] = shape_derivatives[(i, 0)]; // N1x
                result
            };

            let disp = self
                .solution_vec
                .rows(self.num_dofs_per_node * element_nodes[i], 2);
            let stress = shape_functions[virt_node] * self.material_matrix() * b_mat_j * disp;

            for j in 0..3 {
                residual[j] += stress[j]
            }
        }
        //println!("Residual for Node {}, {:?}", virt_node, residual.as_slice());
        // get correct displacements
        //println!("{}",virt_node);

        //println!("{},{:?}",virt_node,disp.as_slice());

        //println!("{}",stress);

        return residual;
    }
}

pub struct YValueDirichlet {
    pub y_value: f64,
}

pub struct XValueDirichlet {
    pub x_value: f64,
}

pub struct XValueCoordinateDirichlet {
    pub x_stretch: f64,
}

pub struct YValueRotationDirichlet {
    pub rot: f64,
}

pub struct XValueClamped {
    pub x_coord: f64
}

impl DirichletBoundary for XValueClamped {
    fn num_dof_per_node(&self) -> usize {
        2
    }


    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        (coords[0] - self.x_coord).abs() < 1e-10 
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        0.0
    }
}

impl DirichletBoundary for YValueRotationDirichlet {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        true
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        //if coords[1] < 1e-10 && dof_num == 1 {
        //    // bei x==0 soll y festgehalten werden
        //    return coords[0] * self.rot;
        //}
        //if coords[0] < 1e-10 && coords[1] < 1e-10 && dof_num == 0 {
        //    return 0.0;
        //}
        //0.0
        let normal_vec = SVector::<f64, 2>::new(coords[1], -coords[0]);
        let rotation = Rotation2::new(self.rot);
        (coords - rotation * coords)[dof_num]
    }
}

impl DirichletBoundary for YValueDirichlet {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        // ist y 0, dann wird y-richtung festgehalten
        if coords[1] < 1e-10 && dof_num == 1 {
            // bei x==0 soll y festgehalten werden
            //println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }
        if coords[0] < 1e-10 && coords[1] < 1e-10 && dof_num == 0 {
            // bei x== 0 und y== 0 soll x auch festgehalten werden
            //println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }

        false
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        if dof_num == 1 {
            // bei x==0 soll y festgehalten werden
            return self.y_value;
        }
        if dof_num == 0 {
            return 0.0;
        }
        0.0
    }
}

impl DirichletBoundary for XValueDirichlet {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        //println!("{},{}",node,coords);
        if coords[0] < 1e-10 && dof_num == 0 {
            // bei x==0 soll x festgehalten werden
            //println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }
        if coords[0] < 1e-10 && coords[1] < 1e-10 && dof_num == 1 {
            // bei x== 0 und y== 0 soll y auch festgehalten werden
            //println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }

        false
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        if dof_num == 0 {
            // bei x==0 soll x festgehalten werden
            return self.x_value;
        }
        if dof_num == 1 {
            return 0.0;
        }
        0.0
    }
}

impl DirichletBoundary for XValueCoordinateDirichlet {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        return dof_num == 0;
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        
        return self.x_stretch * coords[0];
        
    }
}