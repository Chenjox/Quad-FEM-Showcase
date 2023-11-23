// Einlesen eines Meshes
// connectivität des Meshes

use std::{fs::File, io::Write, path::PathBuf};

use nalgebra::{Const, Dyn, MatrixView2x1, OMatrix, OVector};
use vtkio::{Vtk, model::{Version, DataSet, MetaData, UnstructuredGridPiece, Cells, VertexNumbers, CellType, Attributes, DataArray, Attribute}};

type Points2D = OMatrix<f64, Const<2>, Dyn>;
type Elements = OMatrix<usize, Const<4>, Dyn>;
type BoundList = OVector<usize, Dyn>;
type Matrix3x3 = OMatrix<f64, Const<3>, Const<3>>;
type Matrix2x2 = OMatrix<f64, Const<2>, Const<2>>;
// TODO: Boundaries abspeichern!
struct Boundaries {
    left: BoundList,
    right: BoundList,
    top: BoundList,
    bottom: BoundList,
}

struct FEMesh {
    points: Points2D,
    elements: Elements,
    bounds: Boundaries,
}

impl FEMesh {
    pub fn get_point(&self, index: usize) -> MatrixView2x1<f64> {
        return self.points.column(index);
    }
}

// FIXME: Boundaries abspeichern!
fn get_FE_mesh(length: f64, height: f64, length_n: usize, height_n: usize) -> FEMesh {
    // Get increase
    let l_increase = length / (length_n as f64);
    let h_increase = height / (height_n as f64);

    let mut p = Points2D::zeros((length_n + 1) * (height_n + 1));
    let mut count = 0;
    for i in 0..=length_n {
        for j in 0..=height_n {
            p[(0, count)] = (i as f64) * l_increase;
            p[(1, count)] = (j as f64) * h_increase;
            count += 1;
        }
    }

    //println!("{}", p);
    count = 0;
    let mut elem = Elements::zeros(length_n * height_n);

    for i in 0..length_n {
        for j in 0..height_n {
            let row = i * (height_n + 1);
            elem[(0, count)] = row + j;
            elem[(1, count)] = row + j + (height_n + 1);
            elem[(2, count)] = row + j + 1 + (height_n + 1);
            elem[(3, count)] = row + j + 1;
            count += 1;
        }
    }
    // Boundaries
    let mut left = BoundList::zeros(height_n + 1);
    let mut right = BoundList::zeros(height_n + 1);
    let mut top = BoundList::zeros(length_n + 1);
    let mut bottom = BoundList::zeros(length_n + 1);

    // Linker und rechter Rand
    for j in 0..=height_n {
        left[(height_n - j)] = j;
        right[(j)] = j + (length_n + 1) * (height_n + 1) - (height_n + 1);
    }

    // Oberer und unterer Rand
    for i in 0..=length_n {
        bottom[(i)] = i * (height_n + 1);
        top[(length_n - i)] = i * (height_n + 1) + height_n;
    }

    // FIXME: Orientierung ist wichtig! Immer gegen den Uhrzeigersinn

    // println!("{}", left);
    // println!("{}", bottom);
    // println!("{}", right);
    // println!("{}", top);

    // println!("{}", elem);

    return FEMesh {
        points: p,
        elements: elem,
        bounds: Boundaries {
            left,
            right,
            top,
            bottom,
        },
    };
}

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

// Stoffmatrix
fn C_tensor(emodul: f64, thickness: f64, nu: f64) -> Matrix3x3 {
    let d = emodul * thickness / (1.0 - nu.powi(2));
    let mut m = Matrix3x3::zeros();
    m[(0, 0)] = d;
    m[(1, 1)] = d;
    m[(2, 2)] = 0.5 * (1.0 - nu) * d;
    m[(0, 1)] = nu * d;
    m[(1, 0)] = nu * d;
    return m;
}

// Assemblierung Elemente

// Assemblierung Randbedingungen

// Penalty Verfahren für festgesetzte Randbedingungen.

// Lösen

// Rückwärtseinsetzen

// Visualisierung

fn write_elements_to_file(elem: &FEMesh) {
    let p = elem.points.shape().1;
    let mut f = File::create("output.csv").expect("Unable to create file");
    for i in 0..p {
        write!(
            &mut f,
            "{} {} {}\n",
            elem.points[(0, i)],
            elem.points[(1, i)],
            i
        )
        .unwrap();
    }
    let q = elem.elements.shape().1;
    for i in 0..q {
        //println!("{},{}/{}", p, q, i);
        let cen_x = elem.points[(0, elem.elements[(0, i)])] * 0.25
            + elem.points[(0, elem.elements[(1, i)])] * 0.25
            + elem.points[(0, elem.elements[(2, i)])] * 0.25
            + elem.points[(0, elem.elements[(3, i)])] * 0.25;
        let cen_y = elem.points[(1, elem.elements[(0, i)])] * 0.25
            + elem.points[(1, elem.elements[(1, i)])] * 0.25
            + elem.points[(1, elem.elements[(2, i)])] * 0.25
            + elem.points[(1, elem.elements[(3, i)])] * 0.25;
        write!(&mut f, "{} {} {}\n", cen_x, cen_y, i).unwrap();
    }
}

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

fn main() {
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
    let mut projection = OMatrix::<f64,Dyn,Const<1>>::zeros(solution_dim*m.points.ncols());
    let mut constrained = OMatrix::<bool,Dyn,Const<1>>::from_element(solution_dim*m.points.ncols(), false);
    
    // Hier wird jetzt der Rand hardgecoded
    for i in m.bounds.bottom.iter() {
        projection[*i] = 100.;
        //constrained[*i] = true;
    }
    for i in m.bounds.top.iter() {
        projection[*i] = 0.;
        constrained[*i] = true;
    }
    let mut rhs : OMatrix<f64,Dyn,Const<1>> = OMatrix::<f64,Dyn,Const<1>>::zeros(m.points.ncols());

    for i in 0..m.points.ncols() {
        if constrained[i] {
            for j in 0..m.points.ncols() {
                global_matrix[(i,j)] = 0.;
                global_matrix[(j,i)] = 0.;
            }
            global_matrix[(i,i)] = 1.;
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

fn write_vkt_from_mesh(mesh: FEMesh, solution: OMatrix<f64,Dyn,Const<1>>) {

    let path = PathBuf::from(r"\test\test.vtu");

    let mut points_vec = Vec::new();
    for point in mesh.points.column_iter() {
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
    for result in solution.row_iter() {
        result_vector.push(result[0])
    }

    let sols = Attribute::DataArray(
        DataArray {
            name: String::from("Temperature"),
            elem: vtkio::model::ElementType::Scalars { num_comp: 1, lookup_table: None },
            data: vtkio::IOBuffer::F64(result_vector)
        }
    );

    let vt = Vtk {
        version: Version::new((0, 1)),
        title: String::from("Heatmap"),
        file_path: Some(path),
        byte_order: vtkio::model::ByteOrder::LittleEndian,
        data: DataSet::inline(
            UnstructuredGridPiece {
                points: points_vec.into(),
                cells: Cells {
                    cell_verts: VertexNumbers::XML { 
                        offsets: offset_vec, 
                        connectivity: connectivity_vec, 
                    },
                    types: vec![
                        vec![CellType::Quad; num_elems]
                    ].into_iter().flatten().collect::<Vec<CellType>>()
                },
                data: Attributes {
                    point: vec![sols],
                    cell: vec![]
                }
            }
        ).into()
    };

    let mut f = File::create("test.vtu").unwrap();

    vt.write_xml(&mut f).unwrap();
}
