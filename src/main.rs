// Einlesen eines Meshes
// connectivität des Meshes

pub mod element;
pub mod mesh;

use std::{collections::HashMap, fs::File, io::Write, path::PathBuf};

use mshio::ElementType;
use nalgebra::{coordinates, Const, Dyn, MatrixView2x1, OMatrix, OVector, SVector};
use vtkio::{
    model::{
        Attribute, Attributes, CellType, Cells, DataArray, DataSet, MetaData,
        UnstructuredGridPiece, Version, VertexNumbers,
    },
    Vtk,
};

use crate::mesh::{
    elements::{Quad4Element, ReferenceElement},
    femesh::FEMesh,
};

// Connectivität zu Quad Elementen

// Lineare Shape Functions (Locking)
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
//

fn shape_functions_2d(index: usize, xi_1: f64, xi_2: f64) -> f64 {
    return g2_vec(xi_1, xi_2)[index];
}

fn shape_functions_2d_grad(index: usize, xi_1: f64, xi_2: f64) -> OMatrix<f64, Const<2>, Const<1>> {
    return OMatrix::<f64, Const<2>, Const<1>>::new(
        g2_dxi1_vec(xi_1, xi_2)[index],
        g2_dxi2_vec(xi_1, xi_2)[index],
    );
}

fn g2_vec(xi_1: f64, xi_2: f64) -> OVector<f64, Const<4>> {
    let mut m = OVector::<f64, Const<4>>::zeros();
    m[0] = g1(xi_1) * g1(xi_2);
    m[1] = g2(xi_1) * g1(xi_2);
    m[2] = g2(xi_1) * g2(xi_2);
    m[3] = g1(xi_1) * g2(xi_2);
    return m;
}

fn g2_dxi1_vec(xi_1: f64, xi_2: f64) -> OVector<f64, Const<4>> {
    let mut m = OVector::<f64, Const<4>>::zeros();
    m[0] = g1_d(xi_1) * g1(xi_2);
    m[1] = g2_d(xi_1) * g1(xi_2);
    m[2] = g2_d(xi_1) * g2(xi_2);
    m[3] = g1_d(xi_1) * g2(xi_2);
    return m;
}

fn g2_dxi2_vec(xi_1: f64, xi_2: f64) -> OVector<f64, Const<4>> {
    let mut m = OVector::<f64, Const<4>>::zeros();
    m[0] = g1(xi_1) * g1_d(xi_2);
    m[1] = g2(xi_1) * g1_d(xi_2);
    m[2] = g2(xi_1) * g2_d(xi_2);
    m[3] = g1(xi_1) * g2_d(xi_2);
    return m;
}

// Vector der Formfunktionen
fn g_vec(coord: f64) -> OVector<f64, Const<2>> {
    let mut m = OVector::<f64, Const<2>>::zeros();
    m[0] = g1(coord);
    m[1] = g2(coord);
    return m;
}

fn gd_vec(coord: f64) -> OVector<f64, Const<2>> {
    let mut m = OVector::<f64, Const<2>>::zeros();
    m[0] = g1_d(coord);
    m[1] = g2_d(coord);
    return m;
}

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

fn main() {
    let m = FEMesh::<DIM>::read_from_gmsh(
        "test.msh4",
        HashMap::from([(ElementType::Qua4, 0), (ElementType::Lin2, 1)]),
        vec![Box::new(Quad4Element {})],
    )
    .unwrap();

    //if let Some(entities) = &msh.data.entities {
    //    println!("{:?}",entities.surfaces);
    //}
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
