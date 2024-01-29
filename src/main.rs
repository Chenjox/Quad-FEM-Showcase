// Einlesen eines Meshes
// connectivität des Meshes

pub mod element;
pub mod mesh;
pub mod problem;
pub mod elasticity;

use core::num;
use std::{collections::HashMap, fs::File, io::Write, path::PathBuf};

use crate::{elasticity::{ElasticityMixed, XValueClamped, XValueCoordinateDirichlet, XValueDirichlet, YValueDirichlet, YValueRotationDirichlet}, mesh::{
    elements::{Quad4Element, ReferenceElement},
    femesh::FEMesh,
}};
use elasticity::{Elasticity, StressSmootherRHS, StressSmootherStiffness};
use mshio::ElementType;
use nalgebra::{
    coordinates, Const, Dim, Dyn, Matrix, MatrixView2x1, OMatrix, OVector, Rotation2, SMatrix,
    SVector, Storage,
};
use nalgebra_sparse::{coo, factorization::CscCholesky, SparseEntryMut::NonZero};
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use problem::{boundary::{DirichletBoundary, LocalStiffnessAssembler, WeakForm, WeakFormAssembler}, integration::get_1d_gauss_rule};
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



// Diese gibt eine Dimension des Meshs vor, noch nicht implementiert.


fn assemble_stiffness_matrix<A>(
    mesh: &FEMesh<2>,
    assembler: &A,
    rhs: &mut OVector<f64, Dyn>,
) -> CsrMatrix<f64>
where
    A: LocalStiffnessAssembler,
{
    let num_dof_per_node = assembler.num_dof_per_node();

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
        let element_nodes = element.1.as_slice();
        let nodal_coordinates = mesh.get_nodal_coordinates(element.1.as_slice());
        
        // Lokale Steifigkeitsmatrix
        let mut k = OMatrix::<f64, Dyn, Dyn>::zeros(
            num_dof_per_node * num_element_nodes,
            num_dof_per_node * num_element_nodes,
        );
        
        let mut rhs_local = OVector::<f64, Dyn>::zeros(num_dof_per_node * num_element_nodes);
        
        // Ab hier wird der LocalStiffnessAssembler eingeführt
        assembler.assemble_local_stiffness_term(
            &ref_element, 
            element_nodes, 
            &nodal_coordinates, 
            &mut k, 
            &mut rhs_local);

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



fn get_dirichlet_vector_and_map(
    mesh: &FEMesh<2>,
    dirichlet: &Vec<Box<dyn DirichletBoundary>>,
) -> (usize, Vec<usize>, OVector<f64, Dyn>) {
    let mut marked_dofs = HashMap::new();
    let num_dof_per_node = dirichlet[0].num_dof_per_node();
    let num_nodes = mesh.num_nodes();
    let num_dofs = num_dof_per_node * mesh.num_nodes();

    for dirichlet in dirichlet {
        for (node_index, coords) in mesh.coordinates.column_iter().enumerate() {
            // Decision Rule
            let coords = SVector::<f64, 2>::from(coords);
            let mut is_constrained_slice = [false, false];
            let mut is_value_slice = [0.0, 0.0];
            let mut has_constrained = false;
            for j in 0..num_dof_per_node {
                if dirichlet.is_constrained(node_index, j, coords) {
                    has_constrained = true;
                    is_constrained_slice[j] = true;
                    is_value_slice[j] = dirichlet.get_constrained_value(node_index, j, coords);
                    //if coords[0] < 1e-10 {
                    //    marked_dofs.insert(node_index, ([true, true], [1., 0.]));
                    //} else {
                    //    marked_dofs.insert(node_index, ([true, false], [1., 0.]));
                    //}
                }
            };
            if has_constrained { marked_dofs.insert(node_index, (is_constrained_slice, is_value_slice)); };
        }
    }
    //println!("{:?}",marked_dofs);
    let mut dirichlet_vector = OVector::<f64, Dyn>::zeros(num_dofs);
    let mut dirichlet_marker = OVector::<bool, Dyn>::from_element(num_dofs,false);
    let mut constraint_marker = Vec::new();

    let mut num_constrained = 0;
    for (marked_node, (is_constrained, value)) in marked_dofs {
        for i in 0..num_dof_per_node {
            if is_constrained[i] {
                dirichlet_vector[marked_node * num_dof_per_node + i] += value[i];
                dirichlet_marker[marked_node * num_dof_per_node + i] = true;
                constraint_marker.push(marked_node * num_dof_per_node + i);
                num_constrained += 1;
            }
        }
    }
    //println!("{}",dirichlet_marker.transpose());
    let dirichlet_vector = dirichlet_vector;

    //println!("{:3.3}", dirichlet_vector.transpose());

    return (num_constrained, constraint_marker, dirichlet_vector);
}

fn incorporate_dirichlet_boundary(
    mesh: &FEMesh<2>,
    stiffness: &CsrMatrix<f64>,
    rhs: &OVector<f64, Dyn>,
    dirichlet: &Vec<Box<dyn DirichletBoundary>>,
) -> (CsrMatrix<f64>, OVector<f64, Dyn>) {
    let num_dof_per_node = dirichlet[0].num_dof_per_node();
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

fn dirichlet_backsubstitution(
    mesh: &FEMesh<2>,
    solution: &OVector<f64, Dyn>,
    dirichlet: &Vec<Box<dyn DirichletBoundary>>,
) -> OVector<f64, Dyn> {
    let num_dof_per_node = dirichlet[0].num_dof_per_node();
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
    
    run_mixed_form(&"CooksMembrane-006", 2.1, 0.49);
    run_mixed_form(&"CooksMembrane-012", 2.1, 0.49);
    run_mixed_form(&"CooksMembrane-024", 2.1, 0.49);
    run_mixed_form(&"CooksMembrane-048", 2.1, 0.49);
    //run_mixed_form(&"CooksMembrane-096", 2.1, 0.49);
}

fn patch_tests() {
    let emodulus = 2.1; //e5;
    let poisson = 0.2;

    run_patch_test_for_file(&"patchA", emodulus, poisson);
    run_patch_test_for_file(&"patchB", emodulus, poisson);
    run_patch_test_for_file(&"patchC", emodulus, poisson);
    run_patch_test_for_file(&"patchD", emodulus, poisson);
    run_patch_test_for_file(&"patchE", emodulus, poisson);
    run_patch_test_for_file(&"patchF", emodulus, poisson);
}

fn run_mixed_form(file_name: &str, youngs: f64, poisson: f64) {
    println!("Running {}", file_name);
    let m = FEMesh::<DIM>::read_from_gmsh(
        &format!("{}.msh4", file_name),
        HashMap::from([(ElementType::Qua4, 0), (ElementType::Lin2, 1)]),
        vec![Box::new(Quad4Element {})],
    )
    .unwrap();

    let elast_form = ElasticityMixed {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
    };

    let traction =  ConstantTractionForces {
        map: HashMap::from([(2, [0., -1./16.])]),
    };

    let dirichlet: Vec<Box<dyn DirichletBoundary>> = vec![
        Box::new(XValueClamped { x_coord: 0.0 })
    ];

    let num_dof_per_node = elast_form.num_dof_per_node();
    let num_nodes = m.num_nodes();
    let num_dofs = num_dof_per_node * m.num_nodes();

    let mut rhs = OVector::<f64, Dyn>::zeros(num_dofs);

    //let mut stiffness = compute_sparsity_pattern(&m, num_dof_per_node);

    //println!("Setting up Stiffness Matrix");
    let stiffness = assemble_stiffness_matrix(&m, &elast_form, &mut rhs);

    //println!("Assembling Neumann Terms for RHS");
    assemble_rhs_vector(&m, &traction, &mut rhs);

    //println!("{}", rhs);
    // Dirichlet

    //println!("Incorporating Dirichlet Terms");
    let (reduced_stiffness, reduced_rhs) =
        incorporate_dirichlet_boundary(&m, &stiffness, &rhs, &dirichlet);

    let num_free_dofs = reduced_stiffness.ncols();
    let mut dense_stiffness = OMatrix::<f64, Dyn, Dyn>::zeros(num_free_dofs, num_free_dofs);
    for values in reduced_stiffness.triplet_iter() {
        dense_stiffness[(values.0, values.1)] = *values.2;
    }

    //println!("Solving System");

    let reduced_stiffness = reduced_stiffness.transpose_as_csc();

    let choles = CscCholesky::factor(&reduced_stiffness).unwrap();

    let sol_k = choles.solve(&reduced_rhs);

    let mut k = OMatrix::<f64, Dyn, Const<1>>::zeros(sol_k.nrows());
    for i in 0..sol_k.nrows() {
        k[i] = sol_k[i]
    }

    //println!("{}",k);

    // jetzt die Rückwärtsrolle
    let correct_rhs = dirichlet_backsubstitution(&m, &k, &dirichlet);

    println!("{}",&m.get_nodal_coordinates(&[1]));

    write_vkt_from_mesh_disp(file_name, &m, correct_rhs);
}

fn run_patch_test_for_file(file_name: &str, youngs: f64, poisson: f64) {
    println!("Running {}", file_name);
    let m = FEMesh::<DIM>::read_from_gmsh(
        &format!("{}.msh4", file_name),
        HashMap::from([(ElementType::Qua4, 0), (ElementType::Lin2, 1)]),
        vec![Box::new(Quad4Element {})],
    )
    .unwrap();

    println!("Homogen Y");
    let traction = ConstantTractionForces {
        map: HashMap::from([(3, [0., -1.])]),
    };
    let dirichlet: Vec<Box<dyn DirichletBoundary>> =
        vec![Box::new(YValueDirichlet { y_value: 0.0 })];

    patch_test_run(
        &format!("{}-Y-Loading", file_name),
        &m,
        youngs,
        poisson,
        &dirichlet,
        traction,
    );

    println!("Homogen X");
    let traction = ConstantTractionForces {
        map: HashMap::from([(2, [0., 0.])]),
    };
    let dirichlet: Vec<Box<dyn DirichletBoundary>> = vec![
        Box::new(XValueCoordinateDirichlet {
            x_stretch: 1.0
        }),
        Box::new(XValueDirichlet { x_value: 0.0 }),
    ];
    patch_test_run(
        &format!("{}-X-Loading", file_name),
        &m,
        youngs,
        poisson,
        &dirichlet,
        traction,
    );

    println!("X Verschiebung");
    let traction = ConstantTractionForces {
        map: HashMap::new(),
    };
    let dirichlet: Vec<Box<dyn DirichletBoundary>> =
        vec![Box::new(XValueDirichlet { x_value: 1.0 })];
    patch_test_run(
        &format!("{}-X-Disp", file_name),
        &m,
        youngs,
        poisson,
        &dirichlet,
        traction,
    );

    println!("Y Verschiebung");
    let traction = ConstantTractionForces {
        map: HashMap::new(),
    };
    let dirichlet: Vec<Box<dyn DirichletBoundary>> =
        vec![Box::new(YValueDirichlet { y_value: 1.0 })];
    patch_test_run(
        &format!("{}-Y-Disp", file_name),
        &m,
        youngs,
        poisson,
        &dirichlet,
        traction,
    );

    // Rotation test
    println!("Rotationspatch");
    let traction = ConstantTractionForces {
        map: HashMap::new(),
    };
    let dirichlet: Vec<Box<dyn DirichletBoundary>> =
        vec![Box::new(YValueRotationDirichlet { rot: 0.1 })]; // pi/16 ~ noch in der kleinwinkelnäherung
    patch_test_run(
        &format!("{}-XY-Rot", file_name),
        &m,
        youngs,
        poisson,
        &dirichlet,
        traction,
    );
}

fn patch_test_run<N: NeumannBoundaryTerm>(
    name: &str,
    mesh: &FEMesh<2>,
    youngs: f64,
    poisson: f64,
    dirichlet: &Vec<Box<dyn DirichletBoundary>>,
    neumann: N,
) {
    let m = mesh;

    // Wenn ich 2 DOF pro Knoten habe, besitzt jeder Knoten eine untersteifigkeit von num_dof_per_node x num_dof_per_node größe

    let elast = Elasticity {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
    };

    let elast_form = WeakFormAssembler {
        weak_form: elast
    };

    // Am Rand 3
    let traction = neumann;

    let dirichlet = dirichlet;

    let num_dof_per_node = elast_form.num_dof_per_node();
    let num_nodes = m.num_nodes();
    let num_dofs = num_dof_per_node * m.num_nodes();

    let mut rhs = OVector::<f64, Dyn>::zeros(num_dofs);

    //let mut stiffness = compute_sparsity_pattern(&m, num_dof_per_node);

    //println!("Setting up Stiffness Matrix");
    let stiffness = assemble_stiffness_matrix(&m, &elast_form, &mut rhs);

    //println!("Assembling Neumann Terms for RHS");
    assemble_rhs_vector(&m, &traction, &mut rhs);

    //println!("{}", rhs);
    // Dirichlet

    //println!("Incorporating Dirichlet Terms");
    let (reduced_stiffness, reduced_rhs) =
        incorporate_dirichlet_boundary(&m, &stiffness, &rhs, &dirichlet);

    let num_free_dofs = reduced_stiffness.ncols();
    let mut dense_stiffness = OMatrix::<f64, Dyn, Dyn>::zeros(num_free_dofs, num_free_dofs);
    for values in reduced_stiffness.triplet_iter() {
        dense_stiffness[(values.0, values.1)] = *values.2;
    }

    //println!("Solving System");

    let reduced_stiffness = reduced_stiffness.transpose_as_csc();

    let choles = CscCholesky::factor(&reduced_stiffness).unwrap();

    let sol_k = choles.solve(&reduced_rhs);

    let mut k = OMatrix::<f64, Dyn, Const<1>>::zeros(sol_k.nrows());
    for i in 0..sol_k.nrows() {
        k[i] = sol_k[i]
    }

    //println!("{}",k);

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

    let stress_smoother_form = WeakFormAssembler {
        weak_form: stress_smoother
    };

    //println!("Smothing Stresses");

    let mut stress_rhs =
        OVector::<f64, Dyn>::zeros(stress_smoother_form.num_dof_per_node() * m.num_nodes());

    let _stressness = assemble_stiffness_matrix(&m, &stress_smoother_form, &mut stress_rhs);

    let stress_smoother_stiffness = StressSmootherStiffness {
        youngs_modulus: youngs,
        poissons_ratio: poisson,
        solution_vec: correct_rhs.clone(),
        num_nodes_solution: num_nodes,
        num_dofs_per_node: 2,
    };

    let stress_stiff_smoother_form = WeakFormAssembler {
        weak_form: stress_smoother_stiffness
    };

    let mut dummy_stress_rhs = OVector::<f64, Dyn>::zeros(m.num_nodes());

    let stressness =
        assemble_stiffness_matrix(&m, &stress_stiff_smoother_form, &mut dummy_stress_rhs);

    let mut dense_stressness =
        OMatrix::<f64, Dyn, Dyn>::zeros(stressness.nrows(), stressness.ncols());
    for values in stressness.triplet_iter() {
        dense_stressness[(values.0, values.1)] = *values.2;
    }

    //println!("{:3.3}", dense_stressness);
    //println!("{:3.3}", stress_rhs);

    let chol_stress = dense_stressness.cholesky().unwrap();

    let mut stress_solution = OVector::<f64, Dyn>::zeros(stress_rhs.nrows());

    for i in 0..3 {
        let stress_component_rhs = stress_rhs.rows_with_step(i, stressness.nrows(), 2);
        let solution = chol_stress.solve(&stress_component_rhs);
        for j in 0..solution.nrows() {
            stress_solution[j * 3 + i] = solution[j]
        }
        //println!("{}", solution);
        //println!("{}", stress_solution);
    }

    write_vkt_from_mesh(name, &m, correct_rhs, stress_solution);
}

fn write_vkt_from_mesh(
    file: &str,
    mesh: &FEMesh<2>,
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

    //println!("{:?}", result_vector);

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

    let mut f = File::create(format!("{}.vtu", file)).unwrap();

    vt.write_xml(&mut f).unwrap();
}

fn write_vkt_from_mesh_disp(
    file: &str,
    mesh: &FEMesh<2>,
    solution: OMatrix<f64, Dyn, Const<1>>
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

    //println!("{:?}", result_vector);

    let sols = vec![
        Attribute::DataArray(DataArray {
            name: String::from("Displacement"),
            elem: vtkio::model::ElementType::Vectors,
            data: vtkio::IOBuffer::F64(result_vector),
        })
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

    let mut f = File::create(format!("{}.vtu", file)).unwrap();

    vt.write_xml(&mut f).unwrap();
}
