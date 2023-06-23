
// Einlesen eines Meshes
// connectivität des Meshes

use std::{fs::File, io::Write};

use nalgebra::{OMatrix, Dyn, Const, OVector};

type Points2D = OMatrix<f64,Const<2>,Dyn>;
type Elements = OMatrix<usize,Const<4>,Dyn>;
type BoundList = OVector<usize,Dyn>;
type Matrix3x3 = OMatrix<f64,Const<3>,Const<3>>;
// TODO: Boundaries abspeichern!
struct Boundaries {
    left: BoundList,
    right: BoundList,
    top: BoundList,
    bottom: BoundList
}

struct FEMesh {
    points: Points2D,
    elements: Elements,
    bounds: Boundaries
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
        left[(height_n-j)] = j;
        right[(j)] = j + (length_n+1)*(height_n+1) - (height_n+1);
    }

    // Oberer und unterer Rand
    for i in 0..=length_n {
        bottom[(i)] = i*(height_n+1);
        top[(length_n-i)] = i*(height_n+1) + height_n;
    }

    // FIXME: Orientierung ist wichtig! Immer gegen den Uhrzeigersinn

    println!("{}",left);
    println!("{}",bottom);
    println!("{}",right);
    println!("{}",top);

    println!("{}",elem);

    return FEMesh {
        points: p,
        elements: elem,
        bounds: Boundaries { 
            left, 
            right, 
            top, 
            bottom 
        }
    };
}

// Connectivität zu Quad Elementen

// Lineare Shape Functions (Locking)
fn g1(coord: f64) -> f64 {
    0.5 - coord/2.0
}
fn g1_d(_coord: f64) -> f64 {
    - 0.5
}
fn g2(coord: f64) -> f64 {
    0.5 + coord/2.0
}
fn g2_d(_coord: f64) -> f64 {
    0.5
}

// Vector der Formfunktionen
fn g_vec(coord: f64) -> OVector<f64,Const<2>> {
    let mut m = OVector::<f64,Const<2>>::zeros();
    m[0] = g1(coord);
    m[1] = g2(coord);
    return m;
}

// Stoffmatrix 
fn C_tensor(emodul: f64, thickness: f64, nu: f64) -> Matrix3x3 {
    let d = emodul*thickness/(1.0-nu.powi(2));
    let mut m = Matrix3x3::zeros();
    m[(0,0)] = d;
    m[(1,1)] = d;
    m[(2,2)] = 0.5*(1.0 - nu)*d;
    m[(0,1)] = nu*d;
    m[(1,0)] = nu*d;
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
