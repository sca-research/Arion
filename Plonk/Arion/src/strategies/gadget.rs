use crate::{constants::MATRIX, WIDTH};

use super::Strategy;
use dusk_plonk::prelude::*;

/// Implements a Arion strategy for `Witness` as input values.
/// Requires a reference to a `ConstraintSystem`.
pub struct GadgetStrategy<'a, C> {
    /// A reference to the constraint system used by the gadgets
    cs: &'a mut C,
}

impl<'a, C> GadgetStrategy<'a, C>
where
    C: Composer,
{
    /// Constructs a new `GadgetStrategy` with the constraint system.
    pub fn new(cs: &'a mut C) -> Self {
        GadgetStrategy { cs }
    }

    /// Perform the arion permutation on a plonk circuit
    pub fn gadget(composer: &'a mut C, x: &mut [Witness]) {
        let mut strategy = GadgetStrategy::new(composer);

        strategy.perm(x);
    }
}

impl<C> AsMut<C> for GadgetStrategy<'_, C> {
    fn as_mut(&mut self) -> &mut C {
        self.cs
    }
}

impl<'a, C> Strategy<Witness> for GadgetStrategy<'a, C>
where
    C: Composer,
{
    fn quinti_s_box(&mut self, input: &mut Witness) {
        let constraint = Constraint::new().mult(1).a(*input).b(*input);
        let v2 = self.cs.gate_mul(constraint);

        let constraint = Constraint::new().mult(1).a(v2).b(v2);
        let v4 = self.cs.gate_mul(constraint);

        let constraint = Constraint::new().mult(1).a(v4).b(*input);
        *input = self.cs.gate_mul(constraint);
    }

    fn gtds(
        &mut self,
        input: &mut [Witness],
        constants_g: &[[BlsScalar; 2]; WIDTH - 1],
        constants_h: &[BlsScalar; WIDTH - 1],
    ) {
        let mut output = [C::ZERO; WIDTH];
        output.copy_from_slice(&input);

        // High degree inverse S-Box is inverse of x^{257}
        let mut tmp = BlsScalar::zero();
        tmp.clone_from(&self.cs.index(output[WIDTH - 1]));
        tmp = tmp.pow_vartime(&[
            8469711284772863745,
            1214928404647555091,
            15849274830579433833,
            3867970841541904563,
        ]);
        let wit = self.cs.append_witness(tmp);
        // y^2
        let constraint = Constraint::new().mult(1).a(wit).b(wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^4
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^8
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^16
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^32
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^64
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^128
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^256
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(tmp_wit);
        let tmp_wit = self.cs.gate_mul(constraint);
        // y^257
        let constraint = Constraint::new().mult(1).a(tmp_wit).b(wit);
        output[WIDTH - 1] = self.cs.gate_mul(constraint);

        let constraint = Constraint::new()
            .left(1)
            .a(output[WIDTH - 1])
            .right(1)
            .b(input[WIDTH - 1]);
        let mut s = self.cs.gate_add(constraint);

        // (WIDTH - 1)-th round to 2-nd component
        for i in (1..WIDTH - 1).rev() {
            // Quintic S-Box
            self.quinti_s_box(&mut output[i]);

            // Evaluate g
            let constraint = Constraint::new()
                .mult(1)
                .a(s)
                .b(s)
                .right(constants_g[i][0])
                .constant(constants_g[i][1]);
            let g = self.cs.gate_add(constraint);

            // Evaluate h
            let constraint = Constraint::new().mult(1).a(s).b(s).right(constants_h[i]);
            let h = self.cs.gate_add(constraint);

            // Evaluate x^d*g + h
            let constraint = Constraint::new().mult(1).a(output[i]).b(g).d(h).fourth(1);
            output[i] = self.cs.gate_add(constraint);

            let constraint = Constraint::new()
                .left(1)
                .a(s)
                .right(1)
                .b(input[i])
                .fourth(1)
                .d(output[i]);
            s = self.cs.gate_add(constraint);
        }

        // 1-st component
        // Quintic S-Box
        // Quintic S-Box
        self.quinti_s_box(&mut output[0]);

        // Evaluate g
        let constraint = Constraint::new()
            .mult(1)
            .a(s)
            .b(s)
            .right(constants_g[0][0])
            .constant(constants_g[0][1]);
        let g = self.cs.gate_add(constraint);

        // Evaluate h
        let constraint = Constraint::new().mult(1).a(s).b(s).right(constants_h[0]);
        let h = self.cs.gate_add(constraint);

        // Evaluate x^d*g + h
        let constraint = Constraint::new().mult(1).a(output[0]).b(g).d(h).fourth(1);
        output[0] = self.cs.gate_add(constraint);

        input.copy_from_slice(&output);
    }

    fn mul_matrix(&mut self, input: &mut [Witness]) {
        let mut output = [C::ZERO; WIDTH];

        for row in 0..WIDTH {
            for col in 0..WIDTH {
                let constraint = Constraint::new().left(1).a(output[row]).right(MATRIX[row][col]).b(input[col]);
                output[row] = self.cs.gate_add(constraint);
            }
        }

        input.copy_from_slice(&output);
    }

    fn affine_layer(&mut self, input: &mut [Witness], constants_aff: &[BlsScalar; WIDTH]) {
        let mut output = [C::ZERO; WIDTH];
        output.copy_from_slice(&input);

        self.mul_matrix(&mut output);

        for col in 0..WIDTH {
            let constraint = Constraint::new().left(1).a(output[col]).constant(constants_aff[col]);
            output[col] = self.cs.gate_add(constraint);
        }
        
        input.copy_from_slice(&output);
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        strategies::{GadgetStrategy, ScalarStrategy, Strategy},
        WIDTH,
    };
    use core::result::Result;
    use dusk_plonk::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[derive(Default)]
    struct TestCircuit {
        i: [BlsScalar; WIDTH],
        o: [BlsScalar; WIDTH],
    }

    impl Circuit for TestCircuit {
        fn circuit<C>(&self, composer: &mut C) -> Result<(), Error>
        where
            C: Composer,
        {
            let zero = C::ZERO;

            let mut perm: [Witness; WIDTH] = [zero; WIDTH];

            let mut i_var: [Witness; WIDTH] = [zero; WIDTH];
            self.i.iter().zip(i_var.iter_mut()).for_each(|(i, v)| {
                *v = composer.append_witness(*i);
            });

            let mut o_var: [Witness; WIDTH] = [zero; WIDTH];
            self.o.iter().zip(o_var.iter_mut()).for_each(|(o, v)| {
                *v = composer.append_witness(*o);
            });

            // Apply Arion gadget strategy.
            GadgetStrategy::gadget(composer, &mut i_var);

            // Copy the result of the permutation into the perm.
            perm.copy_from_slice(&i_var);

            // Check that the Gadget perm results = BlsScalar perm results
            i_var.iter().zip(o_var.iter()).for_each(|(p, o)| {
                composer.assert_equal(*p, *o);
            });

            Ok(())
        }
    }

    /// Generate a random input and perform a permutation
    fn arion() -> ([BlsScalar; WIDTH], [BlsScalar; WIDTH]) {
        let mut input = [BlsScalar::zero(); WIDTH];

        input
            .iter_mut()
            .for_each(|s| *s = BlsScalar::random(&mut rand::thread_rng()));

        let mut output = [BlsScalar::zero(); WIDTH];

        output.copy_from_slice(&input);
        ScalarStrategy::new().perm(&mut output);

        (input, output)
    }

    /// Setup the test circuit prover and verifier
    fn setup() -> Result<(Prover<TestCircuit>, Verifier<TestCircuit>), Error> {
        const CAPACITY: usize = 1 << 10;

        let pp = PublicParameters::setup(CAPACITY, &mut rand::thread_rng())?;
        let label = b"arion_gadget_tester";

        Compiler::compile::<TestCircuit>(&pp, label)
    }

    #[test]
    fn preimage() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        let (i, o) = arion();

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving
        let (proof, public_inputs) = prover.prove(&mut rng, &circuit)?;

        // Verifying
        verifier.verify(&proof, &public_inputs)?;

        Ok(())
    }

    #[test]
    fn preimage_constant() -> Result<(), Error> {
        let (prover, verifier) = setup()?;

        // Prepare input & output
        let i = [BlsScalar::from(5000u64); WIDTH];
        let mut o = [BlsScalar::from(5000u64); WIDTH];
        ScalarStrategy::new().perm(&mut o);

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving
        let (proof, public_inputs) = prover.prove(&mut rng, &circuit)?;

        // Verifying
        verifier.verify(&proof, &public_inputs)?;

        Ok(())
    }

    #[test]
    fn preimage_fails() -> Result<(), Error> {
        let (prover, _) = setup()?;

        // Generate [31, 0, 0, 0, 0] as real input to the perm but build the
        // proof with [31, 31, 31, 31, 31]. This should fail on verification
        // since the Proof contains incorrect statements.
        let x_scalar = BlsScalar::from(31u64);

        let mut i = [BlsScalar::zero(); WIDTH];
        i[1] = x_scalar;

        let mut o = [BlsScalar::from(31u64); WIDTH];
        ScalarStrategy::new().perm(&mut o);

        let circuit = TestCircuit { i, o };
        let mut rng = StdRng::seed_from_u64(0xbeef);

        // Proving should fail
        assert!(
            prover.prove(&mut rng, &circuit).is_err(),
            "proving should fail since the circuit is invalid"
        );

        Ok(())
    }
}
