
// Einlesen eines Meshes
// connectivität des Meshes

use std::{fs::File, io::Write};

use nalgebra::{OMatrix, Dyn, Const, OVector};

type Points2D = OMatrix<f64,Const<2>,Dyn>;
type Elements = OMatrix<usize,Const<4>,Dyn>;
type BoundList = OVector<usize,Dyn>;
// TODO: Boundaries abspeichern!
struct Boundaries {
    left: BoundList,
    right: BoundList,
    top: BoundList,
    bottom: BoundList
}

struct FEMesh {
    points: Points2D,
    elements: Elements
}

// FIXME: Boundaries abspeichern!
fn get_FE_mesh(length: f64, height: f64, length_n: usize, height_n: usize) -> FEMesh {
    
    // Get increase
    let l_increase = length/(length_n as f64);
    let h_increase = height/(height_n as f64);

    let mut p = Points2D::zeros((length_n+1)*(height_n+1));
    let mut count = 0;
    for i in 0..=length_n {
        for j in 0..=height_n {
            p[(0,count)] = (i as f64)*l_increase;
            p[(1,count)] = (j as f64)*h_increase;
            count+=1;
        }
    }

    println!("{}",p);
    count = 0;
    let mut elem = Elements::zeros(length_n*height_n);

    for i in 0..length_n {
        for j in 0..height_n {
            let row = i*(length_n-1);
            elem[(0,count)] = row+j;
            elem[(1,count)] = row+j+(length_n-1);
            elem[(2,count)] = row+j+1+(length_n-1);
            elem[(3,count)] = row+j+1;
            count+=1;
        }
    }
    // Boundaries
    let mut left = BoundList::zeros(height_n+1);
    let mut right = BoundList::zeros(height_n+1);
    let mut top = BoundList::zeros(length_n+1);
    let mut bottom = BoundList::zeros(length_n+1);

    // Linker und rechter Rand 
    for j in 0..=height_n {
        left[(j)] = j;
        right[(j)] = j + (length_n+1)*(height_n+1) - (height_n+1);
    }

    // Oberer und unterer Rand
    for i in 0..=length_n {
        bottom[(i)] = i*(height_n+1);
        top[(i)] = i*(height_n+1) + height_n;
    }

    // FIXME: Orientierung ist wichtig! Immer gegen den Uhrzeigersinn

    println!("{}",bottom);
    println!("{}",top);

    println!("{}",elem);

    return FEMesh {
        points: p,
        elements: elem
    };
}

// Connectivität zu Quad Elementen

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
        write!(&mut f,"{} {} {}\n",elem.points[(0,i)],elem.points[(1,i)],i).unwrap();
    }
    let q = elem.elements.shape().1;
    for i in 0..q {
        let cen_x = 
            elem.points[(0,elem.elements[(0,i)])]*0.25 + 
            elem.points[(0,elem.elements[(1,i)])]*0.25 +
            elem.points[(0,elem.elements[(2,i)])]*0.25 +
            elem.points[(0,elem.elements[(3,i)])]*0.25;
        let cen_y = 
            elem.points[(1,elem.elements[(0,i)])]*0.25 + 
            elem.points[(1,elem.elements[(1,i)])]*0.25 +
            elem.points[(1,elem.elements[(2,i)])]*0.25 +
            elem.points[(1,elem.elements[(3,i)])]*0.25;
        write!(&mut f,"{} {} {}\n",cen_x,cen_y,i).unwrap();
    }
}

fn main() {
    let m = get_FE_mesh(30.0, 3.0, 5, 3);
    write_elements_to_file(&m);
}
