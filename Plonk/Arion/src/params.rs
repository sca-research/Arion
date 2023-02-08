use crate::{WIDTH, ROUNDS, constants::{CONSTANTS_G, CONSTANTS_H, CONSTANTS_AFF}};
use dusk_bls12_381::BlsScalar;

#[derive(Clone, Debug)]
pub struct ArionParams {
    pub(crate) n: usize,
    pub(crate) rounds: usize,
    pub(crate) constants_g: Vec<Vec<[BlsScalar; 2]>>,
    pub(crate) constants_h: Vec<Vec<BlsScalar>>,
    pub(crate) constants_aff: Vec<Vec<BlsScalar>>,
    pub(crate) mat: Vec<Vec<BlsScalar>>,
}

impl ArionParams {
    pub fn new(n: usize, rounds: usize) -> Self {
        assert!(n <= WIDTH);
        assert!(rounds <= ROUNDS);

        let constants_g = Self::collect_constants_g(n, rounds);
        let constants_h = Self::collect_constants_h(n, rounds);
        let constants_aff = Self::collect_constants_aff(n, rounds);
        let mat = Self::instantiate_matrix(n);

        ArionParams {
            n,
            rounds,
            constants_g,
            constants_h,
            constants_aff,
            mat,
        }
    }

    fn collect_constants_g(n: usize, rounds: usize) -> Vec<Vec<[BlsScalar; 2]>> {
        let mut constants_g: Vec<Vec<[BlsScalar; 2]>> = Vec::new();
        for i in 0..rounds {
            let mut current_round: Vec<[BlsScalar; 2]> = Vec::new();
            for j in 0..(n - 1) {
                current_round.push(CONSTANTS_G[i][j]);
            }
            constants_g.push(current_round);
        }
        constants_g
    }

    fn collect_constants_h(n: usize, rounds: usize) -> Vec<Vec<BlsScalar>> {
        let mut constants_h: Vec<Vec<BlsScalar>> = Vec::new();
        for i in 0..rounds {
            let mut current_round: Vec<BlsScalar> = Vec::new();
            for j in 0..(n - 1) {
                current_round.push(CONSTANTS_H[i][j]);
            }
            constants_h.push(current_round);
        }
        constants_h
    }

    fn collect_constants_aff(n: usize, rounds: usize) -> Vec<Vec<BlsScalar>> {
        let mut constants_aff: Vec<Vec<BlsScalar>> = Vec::new();
        for i in 0..rounds {
            let mut current_round: Vec<BlsScalar> = Vec::new();
            for j in 0..n {
                current_round.push(CONSTANTS_AFF[i][j]);
            }
            constants_aff.push(current_round);
        }
        constants_aff
    }

    fn circ_mat(row: &[BlsScalar]) -> Vec<Vec<BlsScalar>> {
        let n = row.len();
        let mut mat: Vec<Vec<BlsScalar>> = Vec::with_capacity(n);
        let mut rot = row.to_owned();
        mat.push(rot.clone());
        for _ in 1..n {
            rot.rotate_right(1);
            mat.push(rot.clone());
        }
        mat
    }

    fn instantiate_matrix(n: usize) -> Vec<Vec<BlsScalar>> {
        let mut row: Vec<BlsScalar> = vec![BlsScalar::zero(); n];
        for i in 0..n {
            let j = i as u64;
            row[i] = BlsScalar::from(j);
        }
        Self::circ_mat(&row)
    }
}

mod tests {
    use super::*;
    use crate::{WIDTH, ROUNDS};

    pub const MATRIX: [[BlsScalar; WIDTH]; WIDTH] = {
        let mut matrix = [[BlsScalar::zero(); WIDTH]; WIDTH];
    
        let one = BlsScalar::one();
    
        let mut i = 1;
        while i < WIDTH {
            matrix[0][i] = matrix[0][i  - 1].add(&one);
            i += 1;
        }
    
        i = 1;
        let mut j = 1;
        while i < WIDTH {
            matrix[i][0] = matrix[i - 1][WIDTH - 1];
            while j < WIDTH {
                matrix[i][j] = matrix[i - 1][j - 1];
                j += 1;
            }
            j = 1;
            i += 1;
        }
    
        matrix
    };

    #[test]
    fn check_constants_g_generation() {
        let params = ArionParams::new(WIDTH, ROUNDS);
        for i in 0..ROUNDS {
            for j in 0..(WIDTH - 1) {
                assert_eq!(params.constants_g[i][j], CONSTANTS_G[i][j]);
            }
        }
    }

    #[test]
    fn check_constants_h_generation() {
        let params = ArionParams::new(WIDTH, ROUNDS);
        for i in 0..ROUNDS {
            for j in 0..(WIDTH - 1) {
                assert_eq!(params.constants_h[i][j], CONSTANTS_H[i][j]);
            }
        }
    }

    #[test]
    fn check_constants_aff_generation() {
        let params = ArionParams::new(WIDTH, ROUNDS);
        for i in 0..ROUNDS {
            for j in 0..WIDTH {
                assert_eq!(params.constants_aff[i][j], CONSTANTS_AFF[i][j]);
            }
        }
    }

    #[test]
    fn check_matrix_generation() {
        let params = ArionParams::new(WIDTH, ROUNDS);
        for i in 0..ROUNDS {
            for j in 0..WIDTH {
                assert_eq!(params.mat[i][j], MATRIX[i][j]);
            }
        }
    }
}
