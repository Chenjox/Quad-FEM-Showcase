// Einlesen eines Meshes
// connectivität des Meshes

pub mod element;
pub mod mesh;

use core::num;
use std::{collections::HashMap, fs::File, io::Write, path::PathBuf};

use crate::mesh::{
    elements::{Quad4Element, ReferenceElement},
    femesh::FEMesh,
};
use mshio::ElementType;
use nalgebra::{
    coordinates, Const, Dim, Dyn, Matrix, MatrixView2x1, OMatrix, OVector, SMatrix, SVector,
    Storage,
};
use nalgebra_sparse::{coo, SparseEntryMut::NonZero};
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use vtkio::{
    model::{
        Attribute, Attributes, CellType, Cells, DataArray, DataSet, MetaData,
        UnstructuredGridPiece, Version, VertexNumbers,
    },
    Vtk,
};

// Connectivität zu Quad Elementen

// Assemblierung Elemente

// Assemblierung Randbedingungen

// Penalty Verfahren für festgesetzte Randbedingungen.

// Lösen

// Rückwärtseinsetzen

// Visualisierung

fn get_gauss_rule(order: usize) -> OMatrix<f64, Const<3>, Dyn> {
    match order {
        1 => {
            let mut m = OMatrix::<f64, Const<3>, Dyn>::zeros(1);
            m[(0, 0)] = 0.;
            m[(1, 0)] = 0.;
            m[(2, 0)] = 2.;
            m
        }
        2 => {
            let mut m = OMatrix::<f64, Const<3>, Dyn>::zeros(4);
            m[(0, 0)] = (1. / 3.0_f64).sqrt();
            m[(1, 0)] = (1. / 3.0_f64).sqrt();
            m[(2, 0)] = 1.;

            m[(0, 1)] = -(1. / 3.0_f64).sqrt();
            m[(1, 1)] = (1. / 3.0_f64).sqrt();
            m[(2, 1)] = 1.;

            m[(0, 2)] = (1. / 3.0_f64).sqrt();
            m[(1, 2)] = -(1. / 3.0_f64).sqrt();
            m[(2, 2)] = 1.;

            m[(0, 3)] = -(1. / 3.0_f64).sqrt();
            m[(1, 3)] = -(1. / 3.0_f64).sqrt();
            m[(2, 3)] = 1.;

            m
        }
        _ => {
            panic!()
        }
    }
}

fn get_1d_gauss_rule(order: usize) -> OMatrix<f64, Const<2>, Dyn> {
    match order {
        1 => {
            let mut m = OMatrix::<f64, Const<2>, Dyn>::zeros(1);
            m[(0, 0)] = 0.;
            m[(1, 0)] = 2.;
            m
        }
        2 => {
            let mut m = OMatrix::<f64, Const<2>, Dyn>::zeros(2);
            m[(0, 0)] = (1. / 3.0_f64).sqrt();
            m[(1, 0)] = 1.;

            m[(0, 1)] = -(1. / 3.0_f64).sqrt();
            m[(1, 1)] = 1.;

            m
        }
        _ => {
            panic!()
        }
    }
}

const DIM: usize = 2;

fn compute_sparsity_pattern<const DIM: usize>(
    mesh: &FEMesh<DIM>,
    num_dof_per_node: usize,
) -> CsrMatrix<f64> {
    let num_dofs = num_dof_per_node * mesh.num_nodes();

    let mut coo: CooMatrix<f64> = CooMatrix::new(num_dofs, num_dofs);

    for element in mesh.elements.column_iter() {
        for node_i in 0..element.nrows() {
            for node_j in 0..element.nrows() {
                let pos_i = element[node_i] * num_dof_per_node;
                let pos_j = element[node_j] * num_dof_per_node;

                for i in 0..(num_dof_per_node * num_dof_per_node) {
                    let i_k = i % num_dof_per_node;
                    let j_k = i / num_dof_per_node;

                    //println!("{},{},{}",i,i_k,j_k)
                    coo.push(pos_i + i_k, pos_j + j_k, 0.0);
                }
            }
        }
    }

    let csr = CsrMatrix::from(&coo);
    return csr;
}

struct Elasticity {
    youngs_modulus: f64,
    poissons_ratio: f64,
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

        let result = b_mat_i.transpose() * self.material_matrix() * b_mat_j;
        //println!("{}",result);

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

struct StressSmootherRHS {
    youngs_modulus: f64,
    poissons_ratio: f64,
    solution_vec: OVector<f64, Dyn>,
    num_nodes_solution: usize,
    num_dofs_per_node: usize,
}

struct StressSmootherStiffness {
    youngs_modulus: f64,
    poissons_ratio: f64,
    solution_vec: OVector<f64, Dyn>,
    num_nodes_solution: usize,
    num_dofs_per_node: usize,
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
        println!("Residual for Node {}, {:?}",virt_node,residual.as_slice());
        // get correct displacements
        //println!("{}",virt_node);

        //println!("{},{:?}",virt_node,disp.as_slice());

        //println!("{}",stress);

        return residual;
    }
}

// Diese gibt eine Dimension des Meshs vor, noch nicht implementiert.
trait WeakForm {
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

fn assemble_stiffness_matrix<W>(
    mesh: &FEMesh<2>,
    weakform: &W,
    rhs: &mut OVector<f64, Dyn>,
) -> CsrMatrix<f64>
where
    W: WeakForm,
{
    let num_dof_per_node = weakform.num_dof_per_node();

    let mut stiffness = compute_sparsity_pattern(&mesh, num_dof_per_node);

    for element in mesh
        .elements
        .column_iter()
        .enumerate()
        .map(|f| (mesh.ref_element_index[f.0], f.1, f.0))
    {
        //println!("{},{}", element.0, element.1);

        let current_element_index = element.2;
        let ref_element = &mesh.ref_elements[element.0];
        let num_element_nodes = element.1.nrows();
        let gauss = get_gauss_rule(2);
        let nodal_coordinates = mesh.get_nodal_coordinates(element.1.as_slice());

        // Lokale Steifigkeitsmatrix
        let mut k = OMatrix::<f64, Dyn, Dyn>::zeros(
            num_dof_per_node * num_element_nodes,
            num_dof_per_node * num_element_nodes,
        );

        let mut rhs_local = OVector::<f64,Dyn>::zeros(
            num_dof_per_node * num_element_nodes
        );

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

            for j in element.1.row_iter().enumerate() {
                let virt_node_number = j.0;
                for i in element.1.row_iter().enumerate() {
                    let real_node_number = i.0;

                    let result = weight
                        * weakform.stiffness_term(
                            virt_node_number,
                            real_node_number,
                            element.1.as_slice(),
                            &shape_functions,
                            &derivatives,
                        )
                        * determinant;
                    // Offset
                    for i in 0..(num_dof_per_node * num_dof_per_node) {
                        let i_k = i % num_dof_per_node;
                        let j_k = i / num_dof_per_node;

                        k[(
                            real_node_number * num_dof_per_node + i_k,
                            virt_node_number * num_dof_per_node + j_k,
                        )] += result[(i_k, j_k)];
                    }
                }

                let rhs_term = weakform.right_hand_side_body(
                    virt_node_number,
                    &element.1.as_slice(),
                    &shape_functions,
                    &derivatives,
                );
                for i in 0..num_dof_per_node {
                    rhs_local[virt_node_number * num_dof_per_node + i] += weight * rhs_term[i] * determinant;
                }
            }
        }
        // Rechte Seite Vektor:

        // Wenn das Element tatsächlich einen Rand besitzt

        // Einbauen in die Stiffness Matrix
        for node_i in 0..element.1.nrows() {
            let pos_i = element.1[node_i] * num_dof_per_node;
            for node_j in 0..element.1.nrows() {
                let pos_j = element.1[node_j] * num_dof_per_node;

                for i in 0..(num_dof_per_node * num_dof_per_node) {
                    let i_k = i % num_dof_per_node;
                    let j_k = i / num_dof_per_node;

                    //println!("{},{},{}",i,i_k,j_k)
                    if let NonZero(entry) = stiffness.index_entry_mut(pos_i + i_k, pos_j + j_k) {
                        *entry += k[(
                            node_i * num_dof_per_node + i_k,
                            node_j * num_dof_per_node + j_k,
                        )];
                    };
                }
            }
            for i in 0..num_dof_per_node {
                rhs[pos_i + i] += rhs_local[node_i * num_dof_per_node + i];
            }
        }
    }
    return stiffness;
}

trait NeumannBoundaryTerm {
    fn num_dof_per_node(&self) -> usize;

    fn eval_function(
        &self,
        boundary_marker: usize,
        current_global_coordinates: &SVector<f64, 2>,
        outward_normal_vector: &SVector<f64, 2>,
    ) -> OVector<f64, Dyn>;
}

struct ConstantTractionForces {
    map: HashMap<usize, [f64; 2]>,
}

impl NeumannBoundaryTerm for ConstantTractionForces {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn eval_function(
        &self,
        boundary_marker: usize,
        current_global_coordinates: &SVector<f64, 2>,
        outward_normal_vector: &SVector<f64, 2>,
    ) -> OVector<f64, Dyn> {
        if let Some(component_slice) = self.map.get(&boundary_marker) {
            return OVector::<f64, Dyn>::from_column_slice(component_slice);
        } else {
            return OVector::<f64, Dyn>::zeros(2);
        }
    }
}

fn assemble_rhs_vector<W>(mesh: &FEMesh<2>, boundary_term: &W, rhs: &mut OVector<f64, Dyn>)
where
    W: NeumannBoundaryTerm,
{
    let num_dof_per_node = boundary_term.num_dof_per_node();
    // Wenn das Element tatsächlich einen Rand besitzt
    for (current_element_index, boundary_list) in mesh.element_boundary_map.iter() {
        //
        //
        let ref_element = &mesh.ref_elements[mesh.ref_element_index[*current_element_index]];
        let gauss_points = get_1d_gauss_rule(2);
        let element_nodes = mesh.elements.column(*current_element_index);
        let nodal_coordinates = mesh.get_nodal_coordinates(element_nodes.as_slice());

        for boundary_element_index in boundary_list {
            // Lokale Kantennummer und Orientierung
            let edge_index_orientation =
                mesh.boundary_elements_local.column(*boundary_element_index);
            let edge_index = edge_index_orientation[0];
            let edge_orientation = edge_index_orientation[1];

            let boundary_marker = mesh.boundary_type[(0, *boundary_element_index)];
            //println!("{}", boundary_marker);
            // Globale Knotennummern (auch orientiert)
            let global_edge_nodes = mesh.boundary_elements.column(*boundary_element_index);
            // Lokale Knotennummern
            let local_edge_nodes = ref_element.get_edge_nodes(edge_index, edge_orientation);

            //println!("{}",local_edge_nodes);
            // Randintegralsache
            for gauss in gauss_points.column_iter() {
                let xi_1 = gauss[0];
                let co = OVector::<f64, Dyn>::from(vec![xi_1]);
                let weight = gauss[1];

                // Was wäre denn das in Referenzcoordinaten?
                let ref_coords = ref_element.get_edge_coordinate(edge_index, edge_orientation, &co);

                let jacobian = ref_element.get_jacobian_determinant(&nodal_coordinates, ref_coords);
                let differential = if edge_index % 2 == 0 {
                    jacobian.column(0)
                } else {
                    jacobian.column(1)
                };
                let transformation_factor = differential.norm();

                // Welche Shape Functions sind aktiv?

                let shape_function = ref_element.get_shape_functions(ref_coords);
                let coordinates = &nodal_coordinates * &shape_function;

                let neumann_value =
                    boundary_term.eval_function(boundary_marker, &coordinates, &coordinates);
                // Hier bräucht ich dann die tatsächliche Randfunktion
                for i in 0..global_edge_nodes.len() {
                    for j in 0..num_dof_per_node {
                        // Neumannfunktion in Richtung f(0) und f(1)

                        rhs[global_edge_nodes[i] * num_dof_per_node + j] += //func[j] * 
                          neumann_value[j] * weight * &shape_function[local_edge_nodes[i]] * transformation_factor;
                    }
                }
            }
        }
        //println!("{:?}",boundary_list)
    }
}

trait DirichletBoundary {
    fn num_dof_per_node(&self) -> usize;

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool;

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64;
}

struct PointWiseDirichlet {}

impl DirichletBoundary for PointWiseDirichlet {
    fn num_dof_per_node(&self) -> usize {
        2
    }

    fn is_constrained(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> bool {
        //println!("{},{}",node,coords);
        if coords[1] < 1e-10 && dof_num == 1 {
            // bei x==0 soll y festgehalten werden
            println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }
        if coords[0] < 1e-10 && coords[1] < 1e-10 && dof_num == 0 {
            // bei x== 0 und y== 0 soll x auch festgehalten werden
            println!("Constraining Node {} in Direction {}", node, dof_num);
            return true;
        }

        false
    }

    fn get_constrained_value(&self, node: usize, dof_num: usize, coords: SVector<f64, 2>) -> f64 {
        if coords[0] < 1e-10 && dof_num == 1 {
            // bei x==0 soll y festgehalten werden
            return 0.0;
        }
        if coords[0] < 1e-10 && coords[1] < 1e-10 && dof_num == 0 {
            return 0.0;
        }
        0.0
    }
}

fn get_dirichlet_vector_and_map<D: DirichletBoundary>(
    mesh: &FEMesh<2>,
    dirichlet: &D,
) -> (usize, Vec<usize>, OVector<f64, Dyn>) {
    let num_dof_per_node = dirichlet.num_dof_per_node();
    let num_nodes = mesh.num_nodes();
    let num_dofs = num_dof_per_node * mesh.num_nodes();

    let mut marked_dofs = HashMap::new();
    for (node_index, coords) in mesh.coordinates.column_iter().enumerate() {
        // Decision Rule
        let coords = SVector::<f64, 2>::from(coords);
        let mut is_constrained_slice = [false, false];
        let mut is_value_slice = [0.0, 0.0];
        for j in 0..num_dof_per_node {
            if dirichlet.is_constrained(node_index, j, coords) {
                is_constrained_slice[j] = true;
                is_value_slice[j] =
                    dirichlet.get_constrained_value(node_index, num_dof_per_node, coords);
                //if coords[0] < 1e-10 {
                //    marked_dofs.insert(node_index, ([true, true], [1., 0.]));
                //} else {
                //    marked_dofs.insert(node_index, ([true, false], [1., 0.]));
                //}
            }
        }
        marked_dofs.insert(node_index, (is_constrained_slice, is_value_slice));
    }

    let mut dirichlet_vector = OVector::<f64, Dyn>::zeros(num_dofs);

    let mut constraint_marker = Vec::new();

    let mut num_constrained = 0;
    for (marked_node, (is_constrained, value)) in marked_dofs {
        for i in 0..num_dof_per_node {
            if is_constrained[i] {
                dirichlet_vector[marked_node * num_dof_per_node + i] += value[i];
                constraint_marker.push(marked_node * num_dof_per_node + i);
                num_constrained += 1;
            }
        }
    }
    let dirichlet_vector = dirichlet_vector;

    return (num_constrained, constraint_marker, dirichlet_vector);
}

fn incorporate_dirichlet_boundary<D: DirichletBoundary>(
    mesh: &FEMesh<2>,
    stiffness: &CsrMatrix<f64>,
    rhs: &OVector<f64, Dyn>,
    dirichlet: &D,
) -> (CsrMatrix<f64>, OVector<f64, Dyn>) {
    let num_dof_per_node = dirichlet.num_dof_per_node();
    let num_nodes = mesh.num_nodes();
    let num_dofs = num_dof_per_node * mesh.num_nodes();

    let (num_constrained, constraint_marker, dirichlet_vector) =
        get_dirichlet_vector_and_map(mesh, dirichlet);

    let splitting_correction = stiffness * &dirichlet_vector;

    let rhs = rhs - splitting_correction;

    let mut reduced_system_matrix =
        CooMatrix::new(num_dofs - num_constrained, num_dofs - num_constrained);
    let mut reduced_rhs = OVector::<f64, Dyn>::zeros(num_dofs - num_constrained);

    //println!("{:?}", constraint_marker);

    //let mut dense_stiffness =
    //OMatrix::<f64, Dyn, Dyn>::zeros(num_dofs, num_dofs);

    for (row, column, value) in stiffness.triplet_iter() {
        if !constraint_marker.contains(&row) && !constraint_marker.contains(&column) {
            let offset_column: usize = constraint_marker
                .iter()
                .map(|f| if f < &column { 1 } else { 0 })
                .sum();
            let offset_row: usize = constraint_marker
                .iter()
                .map(|f| if f < &row { 1 } else { 0 })
                .sum();
            //println!("{}", column);
            //println!("{},{}",column,offset_column);
            //dense_stiffness[(row,column)] = *value;
            //println!("{},{},{}",row-offset_row,column-offset_column,*value);
            reduced_system_matrix.push(row - offset_row, column - offset_column, *value)
        }
    }

    //println!("{:3.3}", dense_stiffness);
    for i in 0..num_dofs {
        if !constraint_marker.contains(&i) {
            let offset: usize = constraint_marker
                .iter()
                .map(|f| if f < &i { 1 } else { 0 })
                .sum();
            reduced_rhs[i - offset] = rhs[i];
        }
    }

    let reduced_stiffness = CsrMatrix::from(&reduced_system_matrix);

    return (reduced_stiffness, reduced_rhs);
}

fn dirichlet_backsubstitution<D: DirichletBoundary>(
    mesh: &FEMesh<2>,
    solution: &OVector<f64, Dyn>,
    dirichlet: &D,
) -> OVector<f64, Dyn> {
    let num_dof_per_node = dirichlet.num_dof_per_node();
    let num_nodes = mesh.num_nodes();
    let num_dofs = num_dof_per_node * mesh.num_nodes();

    let (num_constrained, constraint_marker, dirichlet_vector) =
        get_dirichlet_vector_and_map(mesh, dirichlet);

    let mut correct_rhs = OVector::<f64, Dyn>::zeros(num_dofs);

    let mut offset_counter = 0;
    for i in 0..num_nodes {
        for j in 0..num_dof_per_node {
            if constraint_marker.contains(&(i * num_dof_per_node + j)) {
                correct_rhs[i * num_dof_per_node + j] = dirichlet_vector[i * num_dof_per_node + j];
                offset_counter += 1;
            } else {
                correct_rhs[i * num_dof_per_node + j] =
                    solution[i * num_dof_per_node + j - offset_counter];
            }
        }
    }

    return correct_rhs;
}

fn main() {
    let m = FEMesh::<DIM>::read_from_gmsh(
        "meshStraigthQuad.msh4",
        HashMap::from([(ElementType::Qua4, 0), (ElementType::Lin2, 1)]),
        vec![Box::new(Quad4Element {})],
    )
    .unwrap();

    // Wenn ich 2 DOF pro Knoten habe, besitzt jeder Knoten eine untersteifigkeit von num_dof_per_node x num_dof_per_node größe

    let youngs = 1.0;
    let poisson = 0.1;

    let elast = Elasticity {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
    };

    let traction = ConstantTractionForces {
        map: HashMap::from([(3, [0., -1.])]),
    };

    let dirichlet = PointWiseDirichlet {};

    let num_dof_per_node = elast.num_dof_per_node();
    let num_nodes = m.num_nodes();
    let num_dofs = num_dof_per_node * m.num_nodes();

    let mut rhs = OVector::<f64, Dyn>::zeros(num_dofs);

    //let mut stiffness = compute_sparsity_pattern(&m, num_dof_per_node);

    println!("Setting up Stiffness Matrix");
    let stiffness = assemble_stiffness_matrix(&m, &elast, &mut rhs);

    println!("Assembling Neumann Terms for RHS");
    assemble_rhs_vector(&m, &traction, &mut rhs);

    println!("{}", rhs);
    // Dirichlet

    println!("Incorporating Dirichlet Terms");
    let (reduced_stiffness, reduced_rhs) =
        incorporate_dirichlet_boundary(&m, &stiffness, &rhs, &dirichlet);

    let num_free_dofs = reduced_stiffness.ncols();
    let mut dense_stiffness = OMatrix::<f64, Dyn, Dyn>::zeros(num_free_dofs, num_free_dofs);
    for values in reduced_stiffness.triplet_iter() {
        dense_stiffness[(values.0, values.1)] = *values.2;
    }

    println!("{:3.3}", dense_stiffness);
    println!("{}", reduced_rhs);

    println!("Solving System");

    let b = dense_stiffness.full_piv_lu();
    let k = b.solve(&reduced_rhs).unwrap();

    println!("{}", k);

    // jetzt die Rückwärtsrolle
    let correct_rhs = dirichlet_backsubstitution(&m, &k, &dirichlet);

    // Postprocessing
    let stress_smoother = StressSmootherRHS {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
        solution_vec: correct_rhs.clone(),
        num_nodes_solution: num_nodes,
        num_dofs_per_node: 2,
    };

    println!("Smothing Stresses");

    let mut stress_rhs =
        OVector::<f64, Dyn>::zeros(stress_smoother.num_dof_per_node() * m.num_nodes());

    let _stressness = assemble_stiffness_matrix(&m, &stress_smoother, &mut stress_rhs);

    let stress_smoother_stiffness = StressSmootherStiffness {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
        solution_vec: correct_rhs.clone(),
        num_nodes_solution: num_nodes,
        num_dofs_per_node: 2,
    };

    let mut dummy_stress_rhs = OVector::<f64, Dyn>::zeros(m.num_nodes());

    let stressness =
        assemble_stiffness_matrix(&m, &stress_smoother_stiffness, &mut dummy_stress_rhs);

    let mut dense_stressness =
        OMatrix::<f64, Dyn, Dyn>::zeros(stressness.nrows(), stressness.ncols());
    for values in stressness.triplet_iter() {
        dense_stressness[(values.0, values.1)] = *values.2;
    }

    println!("{:3.3}", dense_stressness);
    println!("{:3.3}", stress_rhs);

    let chol_stress = dense_stressness.cholesky().unwrap();

    let mut stress_solution = OVector::<f64, Dyn>::zeros(stress_rhs.nrows());

    for i in 0..3 {
        let stress_component_rhs = stress_rhs.rows_with_step(i, stressness.nrows(), 2);
        let solution = chol_stress.solve(&stress_component_rhs);
        for j in 0..solution.nrows() {
            stress_solution[j * 3 + i] = solution[j]
        }
        println!("{}", solution);
        println!("{}", stress_solution);
    }

    write_vkt_from_mesh(m, correct_rhs, stress_solution);

    // CG Verfahren

    /*

    let mut precon = reduced_stiffness.diagonal_as_csr();
    precon.triplet_iter_mut().for_each(|f| *f.2 = 1.0 / (*f.2));

    let mut start_vec = OVector::<f64, Dyn>::zeros(num_dofs - num_constrained);

    //println!("{},{}", reduced_stiffness.nrows(), reduced_rhs.nrows());
    let mut residual = &reduced_rhs - (&reduced_stiffness * &start_vec);
    let mut h0 = &precon * &residual;
    let mut d0 = h0.clone();


    loop {
        let z0 = &reduced_stiffness * &d0;

        let alpha = &residual.dot(&h0) / (&d0.dot(&z0));

        start_vec = &start_vec + alpha * &d0;
        let new_residual = &residual - alpha * z0;
        let h1 = &precon * &residual;

        let beta = (&new_residual.dot(&h1)) / &residual.dot(&h0);

        d0 = &h1 + beta * d0;
        residual = new_residual;
        h0 = h1;

        println!("{}",residual.norm());
        if residual.norm() < 1e-6 {
            break;
        }
    }
    */

    // Dirichlet Rand
    /*
    let mut dense_stiffness = OMatrix::<f64,Dyn,Dyn>::zeros(num_dofs,num_dofs);
    for values in stiffness.triplet_iter() {
        dense_stiffness[(values.0,values.1)] = *values.2;
    }
    println!("{:3.3}",dense_stiffness)
    */
}

fn write_vkt_from_mesh(
    mesh: FEMesh<2>,
    solution: OMatrix<f64, Dyn, Const<1>>,
    stress_solution: OMatrix<f64, Dyn, Const<1>>,
) {
    let path = PathBuf::from(r"\test\test.vtu");

    let mut points_vec = Vec::new();
    for point in mesh.coordinates.column_iter() {
        points_vec.push(point[0]);
        points_vec.push(point[1]);
        points_vec.push(0.0);
    }

    let num_elems = mesh.elements.ncols();
    let mut connectivity_vec = Vec::new();
    let mut offset_vec = Vec::new();

    for elems in mesh.elements.column_iter() {
        for vertex in elems.row_iter() {
            connectivity_vec.push(vertex[0] as u64);
        }
        offset_vec.push(connectivity_vec.len() as u64);
    }

    let mut result_vector = Vec::new();
    for (index, result) in solution.row_iter().enumerate() {
        result_vector.push(result[0]);
        if index % 2 == 1 {
            result_vector.push(0.0)
        }
    }
    let mut stress_vector = Vec::new();
    for (index, result) in stress_solution.row_iter().enumerate() {
        stress_vector.push(result[0]);
    }

    println!("{:?}", result_vector);

    let sols = vec![
        Attribute::DataArray(DataArray {
            name: String::from("Displacement"),
            elem: vtkio::model::ElementType::Vectors,
            data: vtkio::IOBuffer::F64(result_vector),
        }),
        Attribute::DataArray(DataArray {
            name: String::from("Stresses"),
            elem: vtkio::model::ElementType::Generic(3),
            data: vtkio::IOBuffer::F64(stress_vector),
        }),
    ];

    let vt = Vtk {
        version: Version::new((0, 1)),
        title: String::from("Displacements"),
        file_path: Some(path),
        byte_order: vtkio::model::ByteOrder::LittleEndian,
        data: DataSet::inline(UnstructuredGridPiece {
            points: points_vec.into(),
            cells: Cells {
                cell_verts: VertexNumbers::XML {
                    offsets: offset_vec,
                    connectivity: connectivity_vec,
                },
                types: vec![vec![CellType::Quad; num_elems]]
                    .into_iter()
                    .flatten()
                    .collect::<Vec<CellType>>(),
            },
            data: Attributes {
                point: sols,
                cell: vec![],
            },
        })
        .into(),
    };

    let mut f = File::create("test.vtu").unwrap();

    vt.write_xml(&mut f).unwrap();
}
