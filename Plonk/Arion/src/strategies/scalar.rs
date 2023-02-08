use crate::{constants::MATRIX, WIDTH};

use super::Strategy;
use dusk_bls12_381::BlsScalar;

/// Implements a Arion strategy for `BlsScalar` as input input.
#[derive(Default)]
pub struct ScalarStrategy {}

impl ScalarStrategy {
    /// Constructs a new `ScalarStrategy`.
    pub fn new() -> Self {
        Default::default()
    }
}

impl Strategy<BlsScalar> for ScalarStrategy {
    fn quinti_s_box(&mut self, input: &mut BlsScalar) {
        *input = input.square().square() * *input;
    }

    fn gtds(&mut self, input: &mut [BlsScalar], constants_g: &[[BlsScalar; 2]; WIDTH - 1], constants_h: &[BlsScalar; WIDTH - 1]) {
        let mut output = [BlsScalar::zero(); WIDTH];
        output.copy_from_slice(&input);
         
        // High degree inverse S-Box is inverse of x^{257}
        output[WIDTH - 1].pow_vartime(&[
                                            8469711284772863745,
                                            1214928404647555091,
                                            15849274830579433833,
                                            3867970841541904563,
                                            ]);
        
        let mut s = input[WIDTH - 1].clone();
        s += output[WIDTH - 1];
        
        for i in (0..WIDTH - 1).rev() {   
            // Quintic S-Box
            self.quinti_s_box(&mut output[i]);
            
            // Evaluate g and h
            let tmp = s.square();
            // add linear term
            let mut g = s.clone();
            g *= constants_g[i][0];
            let mut h = s.clone();
            h *= constants_h[i];
            // add quadratic term
            g += tmp;
            h += tmp;
            // add constant_term
            g += constants_g[i][1];
            
            // Multply g and add h
            output[i] *= g;
            output[i] += h;

            s += output[i];
            s += input[i];
        }

        input.copy_from_slice(&output);
    }

    fn mul_matrix(&mut self, input: &mut [BlsScalar]) {
        let mut result = [BlsScalar::zero(); WIDTH];
        for row in 0..WIDTH {
            for col in 0..WIDTH {
                result[row] += MATRIX[row][col] * input[col];
            }
        }
        input.copy_from_slice(&result);
    }

    fn affine_layer(&mut self, input: &mut [BlsScalar], constants_aff: &[BlsScalar; WIDTH]) {
        self.mul_matrix(input);
        let mut result = [BlsScalar::zero(); WIDTH];
        for i in 0..WIDTH {
            result[i] += input[i] + constants_aff[i];
        }
        input.copy_from_slice(&result);
    }
}

#[cfg(test)]
mod tests {
    use dusk_bls12_381::BlsScalar;

    // High degree inverse S-Box is inverse of `x^{257}`
    fn high_degree_permutation(samples: usize) {
        for _i in 1..samples {
            let val = BlsScalar::random(&mut rand::thread_rng());
            let val_2 = val.pow_vartime(&[
                                                        8469711284772863745,
                                                        1214928404647555091,
                                                        15849274830579433833,
                                                        3867970841541904563,
                                                        ]);
            
            assert_eq!(val_2.pow_vartime(&[257, 0, 0, 0]), val);
        }  
    }

    #[test]
    fn high_degree_permutation_1000() {
        high_degree_permutation(1000);    
    }
}

