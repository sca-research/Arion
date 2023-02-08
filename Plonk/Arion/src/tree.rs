mod annotation;
mod branch;
mod leaf;
mod zk;

pub use annotation::ArionAnnotation;

#[cfg(feature = "rkyv-impl")]
pub use branch::{
    ArchivedArionBranch, ArchivedArionLevel, ArionBranchResolver,
    ArionLevelResolver,
};
pub use branch::{ArionBranch, ArionLevel};

pub use leaf::ArionLeaf;
pub use zk::merkle_opening;

use core::borrow::Borrow;

use dusk_bls12_381::BlsScalar;
use microkelvin::{Branch, Walker};
use nstack::annotation::{Cardinality, Keyed};
use nstack::NStack;
use ranno::Annotation;

/// Represents a Merkle Tree with a given depth that will be calculated using
/// the Arion Hash technique.
#[derive(Debug, Default)]
pub struct ArionTree<L, K, const DEPTH: usize> {
    inner: NStack<L, ArionAnnotation<K>>,
}

impl<L, K, const DEPTH: usize> Clone for ArionTree<L, K, DEPTH>
where
    L: Clone + ArionLeaf + Keyed<K>,
    K: Clone + PartialOrd,
{
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }
}

impl<L, K, const DEPTH: usize> AsRef<NStack<L, ArionAnnotation<K>>>
    for ArionTree<L, K, DEPTH>
{
    fn as_ref(&self) -> &NStack<L, ArionAnnotation<K>> {
        &self.inner
    }
}

impl<L, K, const DEPTH: usize> AsMut<NStack<L, ArionAnnotation<K>>>
    for ArionTree<L, K, DEPTH>
{
    fn as_mut(&mut self) -> &mut NStack<L, ArionAnnotation<K>> {
        &mut self.inner
    }
}

impl<L, K, const DEPTH: usize> ArionTree<L, K, DEPTH> {
    /// Creates a new arion tree
    pub const fn new() -> Self {
        Self {
            inner: NStack::new(),
        }
    }
}

impl<L, K, const DEPTH: usize> ArionTree<L, K, DEPTH>
where
    L: ArionLeaf + Keyed<K>,
    K: Clone + PartialOrd,
{
    /// Append a leaf to the tree. Return the index of the appended leaf.
    pub fn push(&mut self, mut leaf: L) -> u64 {
        let anno = ArionAnnotation::from_child(&self.inner);
        let cardinality: &Cardinality = anno.borrow();

        let pos = **cardinality;

        leaf.set_pos(pos);
        self.inner.push(leaf);

        pos
    }

    /// Fetch, remove and return the last inserted leaf, if present.
    pub fn pop(&mut self) -> Option<L> {
        self.inner.pop()
    }

    /// Fetch a leaf on a provided index.
    pub fn get(&self, n: u64) -> Option<L>
    where
        L: Clone,
    {
        self.inner.nth(n).map(|b| (*b).clone())
    }

    /// Return a full merkle opening for this arion tree for a given index.
    pub fn branch(&self, n: u64) -> Option<ArionBranch<DEPTH>> {
        self.inner.nth(n).as_ref().map(ArionBranch::from)
    }

    /// Return the current root/state of the tree.
    pub fn root(&self) -> BlsScalar {
        self.branch(0).map(|b| *b.root()).unwrap_or_default()
    }

    /// Provides an iterator over the leaves of the tree from a provided
    /// starting point. To iterate the entire tree, simply provide `0` as
    /// `start`.
    pub fn iter_walk(
        &self,
        start: u64,
    ) -> Option<impl IntoIterator<Item = &L>> {
        self.inner.nth(start)
    }

    /// Provides an iterator over the leaves of the tree which have been
    /// previously annotated via a custom `Walker` passed as argument.
    ///
    /// # Note
    /// This is only useful if annotate the tree is going to make the iteration
    /// perform sub-linearly.
    pub fn annotated_iter_walk<W>(
        &self,
        walker: W,
    ) -> Option<impl IntoIterator<Item = &L>>
    where
        W: Walker<NStack<L, ArionAnnotation<K>>, ArionAnnotation<K>>,
    {
        Branch::walk(&self.inner, walker)
    }
}
