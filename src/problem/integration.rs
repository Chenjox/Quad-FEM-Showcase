use nalgebra::{Const, Dyn, OMatrix};



pub fn get_gauss_rule(order: usize) -> OMatrix<f64, Const<3>, Dyn> {
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

pub fn get_1d_gauss_rule(order: usize) -> OMatrix<f64, Const<2>, Dyn> {
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