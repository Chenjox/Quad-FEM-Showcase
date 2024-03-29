use core::slice;
use std::{
    collections::HashMap,
    error::{self, Error},
    slice::Iter,
};

use mshio::ElementType;
use nalgebra::{Const, Dyn, Matrix, OMatrix, OVector, SMatrix, SVector, ViewStorage};

use crate::mesh::elements::Quad4Element;

use super::elements::ReferenceElement;

/// Ein FE Mesh ist eine Ansammlung aus:
///   - Geometrischer Information (Knotenkoordinaten)
///   - Topologischer Information (Wo sind welche Elemente, was ist Rand)
///   - Weiteren problemspezifischer Information (Bspw.: was ist der E Modul im Element)
///   
/// Dieses FE Mesh geht von einer Isoparametrischen Elementierung aus.
/// Es kümmert sich Ebenfalls um die Zuweisung der Referenzelemente zu den Knotensammlungen
/// Mittels der Referenzelemente wird die Dimension des Problems bestimmt.
/// Ränder sind (n-1) Elemente (oder vllt wird das ander numerisch verarbeitet)
pub struct FEMesh<const DIM: usize> {
    pub coordinates: OMatrix<f64, Const<DIM>, Dyn>, // Der Index der Coordinaten ist der Index des Knotens
    pub ref_elements: Vec<Box<dyn ReferenceElement<DIM>>>, // Die referenzelemente mit der Bilinearform
    pub ref_element_index: OVector<usize, Dyn>, // Every Element gets a reference element via this index
    pub elements: OMatrix<usize, Dyn, Dyn>,     // All the elements
    pub boundary_type: OMatrix<usize, Const<2>, Dyn>, // the type of the reference element of this specific boundary, and its physical tag (multiple elements may or may not belong to the same boundary)
    pub boundary_elements: OMatrix<usize, Dyn, Dyn>,  // All the Boundaries (as Elements...)
    pub boundary_elements_local: OMatrix<usize, Const<2>, Dyn>, // All the Boundaries in local edge coordinates, the second row denotes the orientation
    pub element_boundary_map: HashMap<usize, Vec<usize>>,       // Maps Element in Domain to Rand
}

fn find2(haystack: &[usize], needle: &[usize]) -> Option<usize> {
    //https://stackoverflow.com/questions/57118537/rust-how-can-i-search-a-vec-for-a-subset-and-find-the-start-index-of-the-sub
    for i in 0..haystack.len() - needle.len() + 1 {
        if haystack[i..i + needle.len()] == needle[..] {
            return Some(i);
        }
    }
    None
}

impl<const DIM: usize> FEMesh<DIM> {
    pub fn num_nodes(&self) -> usize {
        return self.coordinates.ncols();
    }

    pub fn get_node_coordinates(
        &self,
        node_index: usize,
    ) -> Matrix<
        f64,
        Const<DIM>,
        Const<1>,
        ViewStorage<'_, f64, Const<DIM>, Const<1>, Const<1>, Const<DIM>>,
    > {
        return self.coordinates.column(node_index);
    }

    pub fn get_nodal_coordinates(&self, node_index_vec: &[usize]) -> OMatrix<f64, Const<DIM>, Dyn> {
        let num = node_index_vec.len();
        let mut b = OMatrix::<f64, Const<DIM>, Dyn>::zeros(num);
        for i in 0..num {
            b.set_column(i, &self.get_node_coordinates(node_index_vec[i]))
        }
        b
    }

    pub fn get_nodal_coordinates_const_size<const NUM: usize>(
        &self,
        node_index_vec: OMatrix<usize, Dyn, Const<1>>,
    ) -> SMatrix<f64, DIM, NUM> {
        let mut b = SMatrix::<f64, DIM, NUM>::zeros();
        //let num = node_index_vec.ncols();
        for i in 0..NUM {
            b.set_column(i, &self.get_node_coordinates(node_index_vec[i]))
        }
        b
    }

    pub fn is_element_on_boundary(&self, element_index: usize) -> bool {
        return self.element_boundary_map.contains_key(&element_index);
    }
}

impl<const DIM: usize> FEMesh<DIM> {
    pub fn read_from_gmsh(
        file_path: &str,
        mapper: HashMap<ElementType, usize>,
        ref_elements: Vec<Box<dyn ReferenceElement<DIM>>>,
    ) -> Result<FEMesh<DIM>, Box<dyn Error>> {
        let msh_bytes = std::fs::read(file_path)?;
        let parser_result = mshio::parse_msh_bytes(msh_bytes.as_slice());

        // Note that the a parser error cannot be propagated directly using the ?-operator, as it
        // contains a reference into the u8 slice where the error occurred.
        let msh = parser_result.map_err(|e| format!("Error while parsing:\n{}", e))?;

        // println!("{:?}",msh.count_element_types());

        let mut coordinates;

        let test = &msh.data.nodes.ok_or("No coordinates found")?;
        {
            coordinates = OMatrix::<f64, Const<DIM>, Dyn>::zeros(test.num_nodes as usize);
            let mut counter = 0;
            for block in &test.node_blocks {
                for node in &block.nodes {
                    let mut coords: SVector<f64, DIM> = SVector::<f64, DIM>::zeros();
                    if DIM == 2 {
                        coords[0] = node.x;
                        coords[1] = node.y;
                    }
                    if DIM == 3 {
                        coords[0] = node.x;
                        coords[1] = node.y;
                        coords[2] = node.z;
                    }
                    coordinates.set_column(counter, &coords);
                    counter += 1;
                }
            }
        }

        //println!("{:2.2}", &coordinates);

        // Wie komme ich an die Randinformation?
        // Mit den physical_tags!
        let mut entity_physical_map: HashMap<i32, i32> = HashMap::new();
        if let Some(metadata) = msh.data.entities {
            for curve in &metadata.curves {
                //println!("{:?}",curve);
                // der Erste index des Physical Tags ist der Randindex
                if let Some(boundary_index) = curve.physical_tags.get(0) {
                    entity_physical_map.insert(curve.tag, *boundary_index);
                };
            }
        }

        // Hier sind alle Randentitäten auf die physischen Tags gemappt
        // entity-Tag -> physical Tag
        //println!("{:?}",entity_physical_map);

        let mapper = mapper;

        /*
        {
            let mut map =  HashMap::new();
            map.insert(ElementType::Qua4, 0);
            map.insert(ElementType::Lin2, 1);
            map
        };
        */

        let mut elem_type_store = Vec::<usize>::new();
        let mut node_store = Vec::<Vec<usize>>::new();

        let mut boundary_elem_type_store = Vec::<(usize, i32)>::new();
        let mut boundary_node_store = Vec::<Vec<usize>>::new();
        // Alle Elemente
        let elements = &msh.data.elements.ok_or("No elements found")?;
        {
            // Eingeteilt in Blöcke
            for block in &elements.element_blocks {
                //println!("Block {:?}", block.element_type);
                // Ist es vllt ein Rand?
                if let Some(physical_tag) = entity_physical_map.get(&block.entity_tag) {
                    // wenn ja, packe es in die Randliste
                    // println!("Boundary detected {:?}", block.element_type);
                    if let Some(elem_type) = mapper.get(&block.element_type) {
                        for element in &block.elements {
                            // Index shift
                            let nodes: Vec<usize> = element
                                .nodes
                                .iter()
                                .map(|index| (index - 1) as usize)
                                .collect();
                            //println!("{:?}",nodes);
                            boundary_elem_type_store.push((*elem_type, *physical_tag));
                            boundary_node_store.push(nodes);
                        }
                    }
                } else if let Some(elem_type) = mapper.get(&block.element_type) {
                    for element in &block.elements {
                        // Index shift
                        let nodes: Vec<usize> = element
                            .nodes
                            .iter()
                            .map(|index| (index - 1) as usize)
                            .collect();
                        //println!("{:?}",nodes);
                        elem_type_store.push(*elem_type);
                        node_store.push(nodes);
                    }
                }
            }
        }

        let num_elements = elem_type_store.len();
        let max_num_nodes = node_store.iter().fold(0, |agg, elem| agg.max(elem.len()));
        // Die Typen zum rekonstruieren von FEMesh
        let mut elements: OMatrix<usize, Dyn, Dyn> =
            OMatrix::<usize, Dyn, Dyn>::zeros(max_num_nodes, num_elements);
        let mut element_type: OVector<usize, Dyn> = OVector::<usize, Dyn>::zeros(num_elements);

        for i in 0..num_elements {
            let elem = &node_store[i];
            for j in 0..elem.len() {
                elements[(j, i)] = elem[j];
            }
            element_type[i] = elem_type_store[i];
        }

        let num_elements = boundary_elem_type_store.len();
        let max_num_nodes = boundary_node_store
            .iter()
            .fold(0, |agg, elem| agg.max(elem.len()));
        // Die Typen zum rekonstruieren von FEMesh
        let mut boundary_elements: OMatrix<usize, Dyn, Dyn> =
            OMatrix::<usize, Dyn, Dyn>::zeros(max_num_nodes, num_elements);
        let mut boundary_element_type: OMatrix<usize, Const<2>, Dyn> =
            OMatrix::<usize, Const<2>, Dyn>::zeros(num_elements);

        for i in 0..num_elements {
            let elem = &boundary_node_store[i];
            for j in 0..elem.len() {
                boundary_elements[(j, i)] = elem[j];
            }
            boundary_element_type[(0, i)] = boundary_elem_type_store[i].1 as usize;
            boundary_element_type[(1, i)] = boundary_elem_type_store[i].0;
        }

        //println!("{},{}", elements, element_type);
        //println!("{},{}", boundary_elements, boundary_element_type);

        // Jetzt das Boundary Mapping
        let mut boundary_containing_elements = HashMap::<usize, Vec<usize>>::new();
        let shape = boundary_elements.shape();
        let mut boundary_local_mapping = OMatrix::<usize, Const<2>, Dyn>::zeros(shape.1);

        // Abbildung:
        for boundary_curve in boundary_elements.column_iter().enumerate() {
            let (boundary_curve_index, boundary_curve) = boundary_curve;
            for element in elements.column_iter().enumerate() {
                let (element_index, element) = element;
                let slice = element.as_slice();
                let mut slice = [slice, &[slice[0]]].concat();
                //println!("{} in {}?",boundary_curve.transpose(), element.transpose());
                if let Some(index) = find2(&slice, boundary_curve.as_slice()) {
                    //println!("{},{}",element_index,index);

                    // Alle Elemente die einen Rand besitzen
                    if let Some(saved_boundaries_indezes) =
                        boundary_containing_elements.get_mut(&element_index)
                    {
                        saved_boundaries_indezes.push(boundary_curve_index);
                    } else {
                        boundary_containing_elements
                            .insert(element_index, vec![boundary_curve_index]);
                    };
                    //
                    boundary_local_mapping[(0, boundary_curve_index)] = index % 4;

                    //println!("{}",boundary_local_mapping);
                    //boundary_local_mapping.push((index,(index+1) % 4));
                    //
                    break;
                };
                slice.reverse(); // maybe the element is backwards?
                if let Some(index) = find2(&slice, boundary_curve.as_slice()) {
                    if let Some(saved_boundaries_indezes) =
                        boundary_containing_elements.get_mut(&element_index)
                    {
                        saved_boundaries_indezes.push(boundary_curve_index);
                    } else {
                        boundary_containing_elements
                            .insert(element_index, vec![boundary_curve_index]);
                    };
                    //
                    boundary_local_mapping[(0, boundary_curve_index)] = index % 4;
                    boundary_local_mapping[(1, boundary_curve_index)] = 1;
                    break;
                };
            }
        }

        //println!("{}", boundary_element_type);

        let final_mesh = FEMesh {
            coordinates: coordinates,
            ref_elements: ref_elements,
            ref_element_index: element_type,
            elements: elements,
            boundary_type: boundary_element_type,
            boundary_elements: boundary_elements,
            boundary_elements_local: boundary_local_mapping,
            element_boundary_map: boundary_containing_elements,
        };

        return Ok(final_mesh);
    }
}
