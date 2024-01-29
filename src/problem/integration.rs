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
        3 => {
            let mut m = OMatrix::<f64, Const<3>, Dyn>::zeros(9);
            let weights_1D = [5.0/9.0,8.0/9.0,5.0/9.0];
            let pos_1D = [(3. / 5.0_f64).sqrt(), 0.0, -(3. / 5.0_f64).sqrt()];
            for k in 0..9 {
                let i_k = k % 3;
                let j_k = k / 3;

                m[(0,k)] = pos_1D[i_k];
                m[(1,k)] = pos_1D[j_k];
                m[(2,k)] = weights_1D[i_k] * weights_1D[j_k];
            }
            
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