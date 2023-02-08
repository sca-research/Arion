use crate::tree::ArionLeaf;

use core::borrow::Borrow;

use dusk_bls12_381::BlsScalar;
use crate::{strategies::{ScalarStrategy, Strategy}, WIDTH};

use nstack::annotation::{Cardinality, Keyed, MaxKey};
use nstack::NStack;
use ranno::Annotation;

/// Annotation holding the root, cardinality, and the maximum value of a generic
/// key.
#[derive(Debug, Clone, Default)]
pub struct ArionAnnotation<K> {
    arion_root: BlsScalar,
    cardinality: Cardinality,
    max_key: MaxKey<K>,
}

impl<K> Borrow<BlsScalar> for ArionAnnotation<K> {
    fn borrow(&self) -> &BlsScalar {
        &self.arion_root
    }
}

impl<K> Borrow<Cardinality> for ArionAnnotation<K> {
    fn borrow(&self) -> &Cardinality {
        &self.cardinality
    }
}

impl<K> Borrow<MaxKey<K>> for ArionAnnotation<K> {
    fn borrow(&self) -> &MaxKey<K> {
        &self.max_key
    }
}

impl<L, K> Annotation<NStack<L, ArionAnnotation<K>>>
    for ArionAnnotation<K>
where
    L: ArionLeaf + Keyed<K>,
    K: Clone + PartialOrd,
{
    fn from_child(stack: &NStack<L, ArionAnnotation<K>>) -> Self {
        let mut perm = [BlsScalar::zero(); WIDTH];
        let mut flag = 1;
        let mut mask = 0;
        let mut cardinality = 0;
        let mut max_key = MaxKey::<K>::NegativeInfinity;

        match stack {
            NStack::Leaf(leaf) => {
                for (i, l) in leaf.iter().enumerate() {
                    if let Some(l) = l {
                        mask |= flag;
                        perm[i + 1] = l.arion_hash();
                        cardinality += 1;

                        let key = l.key();
                        if &max_key < key {
                            max_key = MaxKey::Maximum(key.clone());
                        }
                    }
                    flag <<= 1;
                }
            }
            NStack::Node(node) => {
                for (i, n) in node.iter().enumerate() {
                    if let Some(annotated) = n {
                        mask |= flag;
                        let anno = annotated.anno();
                        let anno = &*anno;

                        perm[i + 1] = anno.arion_root;
                        cardinality += *anno.cardinality;

                        if max_key < anno.max_key {
                            max_key = anno.max_key.clone();
                        }
                    }
                    flag <<= 1;
                }
            }
        }

        perm[0] = BlsScalar::from(mask);
        ScalarStrategy::new().perm(&mut perm);
        let arion_root = perm[1];

        Self {
            cardinality: cardinality.into(),
            arion_root,
            max_key,
        }
    }
}
