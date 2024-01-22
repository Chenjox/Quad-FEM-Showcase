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
use nalgebra::{coordinates, Const, Dim, Dyn, MatrixView2x1, OMatrix, OVector, SMatrix, SVector};
use nalgebra_sparse::SparseEntryMut::NonZero;
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

    let mut csr = CsrMatrix::from(&coo);
    return csr;
}

fn main() {
    let m = FEMesh::<DIM>::read_from_gmsh(
        "test.msh4",
        HashMap::from([(ElementType::Qua4, 0), (ElementType::Lin2, 1)]),
        vec![Box::new(Quad4Element {})],
    )
    .unwrap();

    // Wenn ich 2 DOF pro Knoten habe, besitzt jeder Knoten eine untersteifigkeit von num_dof_per_node x num_dof_per_node größe
    let num_dof_per_node = 2;

    let num_dofs = num_dof_per_node * m.num_nodes();

    let mut rhs = OVector::<f64, Dyn>::zeros(num_dofs);

    let mut stiffness = compute_sparsity_pattern(&m, num_dof_per_node);

    //let mut dense_stiffness = OMatrix::<f64,Dyn,Dyn>::zeros(num_dofs,num_dofs);
    //for values in csr.triplet_iter() {
    //    dense_stiffness[(values.0,values.1)] = 1.0;
    //}
    //println!("{}",dense_stiffness)

    // Assemblierung der Steifigkeitsmatrix
    // Jedes Element im Mesh
    for element in m
        .elements
        .column_iter()
        .enumerate()
        .map(|f| (m.ref_element_index[f.0], f.1))
    {
        //println!("{},{}", element.0, element.1);

        let ref_element = &m.ref_elements[element.0];
        let num_element_nodes = element.1.nrows();
        let gauss = get_gauss_rule(2);

        // Lokale Steifigkeitsmatrix
        let mut k = OMatrix::<f64, Dyn, Dyn>::zeros(
            num_dof_per_node * num_element_nodes,
            num_dof_per_node * num_element_nodes,
        );

        for gauss_point in gauss.column_iter() {
            let xi_1 = gauss_point[0];
            let xi_2 = gauss_point[1];
            let weight = gauss_point[2];

            let ref_coordinates = SVector::<f64, 2>::new(xi_1, xi_2);

            //let normal_func = ref_element.get_shape_functions(ref_coordinates);

            //println!("{}",normal_func);

            let derivatives = ref_element.get_shape_function_derivatives(ref_coordinates);

            //println!("{}",derivatives);

            let mut jacobian = SMatrix::<f64, DIM, DIM>::zeros();
            for (local_index, point_index) in element.1.iter().enumerate() {
                let point_index = *point_index;
                let coords = m.get_node_coordinates(point_index);

                // Grandienten einer shape function in beide richtungen
                let gradient = derivatives.row(local_index);
                let loc_jacobian = coords * gradient;

                jacobian = jacobian + loc_jacobian;
            }
            let jacobian = jacobian;
            let inv_jacobian = jacobian.try_inverse().unwrap();

            // Transformation auf tatsächliche Elemente
            // println!("{}", jacobian);
            let derivatives = (inv_jacobian * derivatives.transpose()).transpose();

            for j in element.1.row_iter().enumerate() {
                let virt_node_number = j.0;
                for i in element.1.row_iter().enumerate() {
                    let real_node_number = i.0;

                    // mit den virtuellen Funktionen wird getestet!
                    let B_mat_i = {
                        let mut result = SMatrix::<f64, 3, 2>::zeros();
                        result[(0, 0)] = derivatives[(virt_node_number, 0)]; // N1x
                        result[(1, 1)] = derivatives[(virt_node_number, 1)]; // N1y
                        result[(2, 0)] = derivatives[(virt_node_number, 1)]; // N1y
                        result[(2, 1)] = derivatives[(virt_node_number, 0)]; // N1x
                        result
                    };
                    let B_mat_j = {
                        let mut result = SMatrix::<f64, 3, 2>::zeros();
                        result[(0, 0)] = derivatives[(real_node_number, 0)]; // N1x
                        result[(1, 1)] = derivatives[(real_node_number, 1)]; // N1y
                        result[(2, 0)] = derivatives[(real_node_number, 1)]; // N1y
                        result[(2, 1)] = derivatives[(real_node_number, 0)]; // N1x
                        result
                    };

                    let result = weight * B_mat_i.transpose() * B_mat_j;
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
            }
        }
        //println!("{:3.3}", k);

        for node_i in 0..element.1.nrows() {
            for node_j in 0..element.1.nrows() {
                let pos_i = element.1[node_i] * num_dof_per_node;
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
        }
    }

    stiffness

    // Jetzt die Ränder
    // Neumann Rand

    // Dirichlet Rand
    /*
    let mut dense_stiffness = OMatrix::<f64,Dyn,Dyn>::zeros(num_dofs,num_dofs);
    for values in stiffness.triplet_iter() {
        dense_stiffness[(values.0,values.1)] = *values.2;
    }
    println!("{:3.3}",dense_stiffness)
    */
}

/*
fn main2() {
    let m = get_FE_mesh(30.0, 3.0, 30, 10);
    write_elements_to_file(&m);

    let solution_dim = 1;
    let mut global_matrix = OMatrix::<f64, Dyn, Dyn>::zeros(
        solution_dim * m.points.ncols(),
        solution_dim * m.points.ncols(),
    );
    //println!("{:?}", global_matrix.shape());

    // Hier beginnt die Assemblierung der Steifigkeitsmatrix
    // Für jedes Element im Mesh
    for (elem_index, elem) in m.elements.column_iter().enumerate() {
        //println!("{}",elem);
        //println!("Elem: {}", elem_index);
        // Gauss-Regel,
        let gauss = get_gauss_rule(2);

        // Gauss-Points
        //Elementsteifigkeit
        let mut k = OMatrix::<f64, Const<4>, Const<4>>::zeros();
        for gauss_point in gauss.column_iter() {
            let xi_1 = gauss_point[0];
            let xi_2 = gauss_point[1];
            let weight = gauss_point[2];
            //println!("Gauss: {:2.2},{:2.2}", xi_1, xi_2);

            // Alle Elementkoordinaten des Geometrischen Elements
            let shape_deriv_1 = g2_dxi1_vec(xi_1, xi_2);
            let shape_deriv_2 = g2_dxi2_vec(xi_1, xi_2);

            let mut jacobian: Matrix2x2 = Matrix2x2::zeros();
            for (local_index, point_index) in elem.iter().enumerate() {
            let point_index = *point_index;
            let coords = m.get_point(point_index);

                let gradient = OMatrix::<f64, Const<2>, Const<1>>::new(
                    shape_deriv_1[local_index],
                    shape_deriv_2[local_index],
                );
                let loc_jacobian = coords * gradient.transpose();

                jacobian = jacobian + loc_jacobian;
            }
            let jacobian = jacobian;
            let inv_jacobian = jacobian.try_inverse().unwrap();
            //println!("{:2.2},{:2.2},{}", jacobian, inv_jacobian,jacobian.determinant());

            // Jetzt alles des Tatsächlichen Elements

            // Hier brauchen wir die Elementsteifigkeitsmatrix!

            for (local_index_virt, point_index_virt) in elem.iter().enumerate() {
                for (local_index, point_index) in elem.iter().enumerate() {
                    //let point_index = *point_index;
                    //let coords = m.get_point(point_index);

                    let shape_deriv =
                        inv_jacobian * shape_functions_2d_grad(local_index, xi_1, xi_2);

                    let test_deriv =
                        inv_jacobian * shape_functions_2d_grad(local_index_virt, xi_1, xi_2);

                    let value = weight * shape_deriv.dot(&test_deriv);

                    k[(local_index, local_index_virt)] += value;

                    //println!("{:2.2},{:2.2}", point_index,point_index_virt);

                    // So könnte man globale Koordinaten errechnen, aber das brauchen wir garnicht, da alles in lokalen ausgedrückt ist.
                    //let global_coordinate = coords * shape_geom[local_index];
                    //println!("{}",global_coordinate);
                }
            }
        }
        // Den Teil hier könnnte man sich sparen, wenn man das einpackt in ne größere Routine
        for (local_p, global_p) in elem.iter().enumerate() {
            for (local_q, global_q) in elem.iter().enumerate() {
                global_matrix[(*global_p, *global_q)] += k[(local_p, local_q)]
            }
        }
    }
    //println!("{:2.4}", global_matrix);

    let A_cop = global_matrix.clone();
    let A_cop_inv = A_cop.try_inverse().unwrap();
    // Inhomogen Dirichlet
    let mut projection = OMatrix::<f64, Dyn, Const<1>>::zeros(solution_dim * m.points.ncols());
    let mut constrained =
        OMatrix::<bool, Dyn, Const<1>>::from_element(solution_dim * m.points.ncols(), false);

    // Hier wird jetzt der Rand hardgecoded
    for i in m.bounds.bottom.iter() {
        projection[*i] = 100.;
        //constrained[*i] = true;
    }
    for i in m.bounds.top.iter() {
        projection[*i] = 0.;
        constrained[*i] = true;
    }
    let mut rhs: OMatrix<f64, Dyn, Const<1>> =
        OMatrix::<f64, Dyn, Const<1>>::zeros(m.points.ncols());

    for i in 0..m.points.ncols() {
        if constrained[i] {
            for j in 0..m.points.ncols() {
                global_matrix[(i, j)] = 0.;
                global_matrix[(j, i)] = 0.;
            }
            global_matrix[(i, i)] = 1.;
        }
        rhs[i] = projection[i];
    }

    //println!("{:2.4}", global_matrix);
    // Solving the system

    // Hier landen die Unbekannten
    //let x_vec = OMatrix::<f64,Dyn,Const<1>>::zeros(solution_dim*m.points.ncols());

    let A = global_matrix;
    let b = rhs;

    //println!("{:2.4}",b);

    let chol = A.cholesky().unwrap();
    let inv_a = chol.inverse();

    //let x = inv_a * b;
    let x = chol.solve(&b);

    //let sol = A_cop * x;

    //println!("{}",x);

    write_vkt_from_mesh(m, x);
}

*/
