use crate::backend_extfield::DietMacAndCheeseExtField;
use crate::circuit_ir::ConvGate;
use crate::circuit_ir::{
    CircInputs, CompiledInfo, FieldInputs, FunId, FunStore, FuncDecl, GateM, TypeSpecification,
    TypeStore, WireCount, WireId, WireRange,
};
use crate::dora::DoraState;
use crate::edabits::{Conv, Edabits};
use crate::fields::SieveIrDeserialize;
use crate::homcom::FCom;
use crate::mac::{Mac, MacT};
use crate::memory::Memory;
use crate::plaintext::DietMacAndCheesePlaintext;
use crate::plugins::{
    DisjunctionBody, PluginExecution, PluginType, Ram, RamArithV0, RamArithV1, RamBoolV0,
    RamBoolV1, RamOp, RamVersion,
};
use crate::ram::ArithmeticRam;
use crate::sieveir_reader_fbs::BufRelation;
use crate::sieveir_reader_text::TextRelation;
use crate::svole_thread::SvoleAtomicRoundRobin;
use crate::svole_trait::{Svole, SvoleStopSignal, SvoleT};
use crate::{backend_trait::BackendT, circuit_ir::FunctionBody};
use crate::{
    dora::{Disjunction, Dora},
    gadgets::less_than_eq_with_public,
};
use crate::{mapping_lpn_size, mapping_lpn_size_large_field, LpnSize};
use crate::{number_to_u64, DietMacAndCheese};
use eyre::{bail, ensure, OptionExt, Result};
use generic_array::typenum::Unsigned;
use log::{debug, info, warn};
use mac_n_cheese_sieve_parser::text_parser::RelationReader;
use mac_n_cheese_sieve_parser::{Number, PluginTypeArg};
use ocelot::svole::LpnParams;
use scuttlebutt::AbstractChannel;
use scuttlebutt::AesRng;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Debug;
use std::io::{Read, Seek};
use std::iter;
use std::marker::PhantomData;
use std::path::PathBuf;
use swanky_field::{FiniteField, FiniteRing, PrimeFiniteField, StatisticallySecureField};
use swanky_field_binary::{F40b, F2};
use swanky_field_f61p::F61p;
use swanky_field_ff_primes::{F127p, F128p, F384p, F384q, Secp256k1, Secp256k1order};
use swanky_party::private::{ProverPrivate, ProverPrivateCopy};
use swanky_party::{IsParty, Party, Prover, WhichParty};

// This file implements IR0+ support for diet-mac-n-cheese and is broken up into the following components:
//
// 0)   Assuming `DietMacAndCheeseProver/Verifier` and `BackendT` which provides the interface and implementation of
//         primitive arithmetic gates for a single field
// I)   Field Conversion. Extend BackendT with interface for field-switching conversions
// II)  Circuit. Instance/Witness/Relation/FunStore
// III) Memory. Structure for Wires and stack to support function calls
// IV)  EvaluatorSingle. An evaluator for a single field
// V)   EvaluatorMulti. An evaluator holding multiple single-field evaluators.

// NOTES: optimizations to consider
//
// * use `*mut X` instead of AbsoluteAddr
// * retry the StreamGate but without the Arc, that might be the reason of the slow down.
//   Profiler was showing a lot of time on `drop Gate`, which implies not building the intermediate gate
// * Reduce to one round the `mult_check` and `check_zero`.

// ## Design Choices Integrating Edabits Conversion Check
//
// The protocol for conversions is designed to check a batch of conversions at once, as opposed to individual conversions.
// It is secure when the batch contains more than N conversions to check, depending on a second parameter B.
// According to the protocol, the number of voles underlying a conversion check is proportional to N*B,
// and the memory required is also proportional to N*B.
// According to Theorem 1 of the Appenzeller to Brie paper, if N is the batch size, B is a parameter and
// s be the security parameter, then the protocol is secure when
// N >= 2^{s/(B-1)} for B=3 or 4 or 5.
// So the lower B, the higher becomes N for the protocol to be secure.
// Working out the equation we get the following constants for N and B with s=40.
// For B = 5, N >= 1024
// For B = 4, N >= 10_321;
// For B = 3, N >= 1_048_576;
// There is a memory/time tradeoff for the different set of parameters. We wont use B=3 because even though
// it would be the most time efficient, it could require several GB of memory. Therefore we are going to
// limit ourselves to using B=4 and B=5, where we call "safe" the later.
// Since we do not know ahead of time how many conversions are in a circuit, we decide to
// presume that the number of conversions will be larger than N=10_321 with B=4.
// During the execution of a circuit we maintain a bucket of conversions to check.
// When the bucket becomes larger than 2*N, we perform a conversion check for N of them,
// and leave the other N unchecked in the bucket. Once the entire circuit is executed, and we reach finalize(),
// we analyze the number of conversions that remain to be checked and may use the B=5 parameter when appropriate.
// And if there are fewer than 1024 conversions in the bucket then a warning is emitted to indicate that the
// conversion check might be insecure.
const CONVERSION_PARAM_B: usize = 4;
const CONVERSION_BATCH_SIZE: usize = 10_321;
const CONVERSION_PARAM_B_SAFE: usize = 5;
const CONVERSION_BATCH_SIZE_SAFE: usize = 1_024;

/// This trait extends the [`BackendT`] trait with `assert_conv_*`
/// functions to go to bits.
pub trait BackendConvT<P: Party>: BackendT {
    // Convert a wire to bits in lower-endian
    fn assert_conv_to_bits(&mut self, w: &Self::Wire) -> Result<Vec<Mac<P, F2, F40b>>>;
    // convert bits in lower-endian to a wire
    fn assert_conv_from_bits(&mut self, x: &[Mac<P, F2, F40b>]) -> Result<Self::Wire>;

    // Finalize the field switching conversions, by running edabits conversion checks
    fn finalize_conv(&mut self) -> Result<()>;
}

pub trait BackendLiftT: BackendT {
    type LiftedBackend: BackendT<Wire = <Self::Wire as MacT>::LiftedMac>;
    fn lift(&mut self) -> &mut Self::LiftedBackend;
}

impl<
        P: Party,
        T: PrimeFiniteField,
        C: AbstractChannel + Clone,
        SVOLE1: SvoleT<P, F2, F40b>,
        SVOLE2: SvoleT<P, T, T>,
    > BackendLiftT for DietMacAndCheeseConv<P, T, C, SVOLE1, SVOLE2>
{
    type LiftedBackend = Self;

    fn lift(&mut self) -> &mut Self::LiftedBackend {
        self
    }
}

pub trait BackendDisjunctionT: BackendT {
    // finalize the disjunctions, by running the final Dora checks
    fn finalize_disj(&mut self) -> Result<()>;

    // execute a disjunction on the given inputs
    fn disjunction(
        &mut self,
        inswit: &mut FieldInputs,
        fun_store: &FunStore,
        inputs: &[Self::Wire],
        disj: &DisjunctionBody,
    ) -> Result<Vec<Self::Wire>>;
}

/// Identifiers for RAMs.
///
/// Per underlying field type, RAMs are globally scoped, and referred to by a
/// unique ID. Wires of RAM type carry these unique IDs so RAMs may be correctly
/// accessed from within any function scope.
pub type RamId = usize;

/// A backend type that supports operations on RAMs.
///
/// Backends implementing this trait may create, read from, and write to RAMs.
pub trait BackendRamT: BackendT {
    /// Create a new RAW with `size` cells, `addr_count`-wide addresses, and
    /// `value_count`-wide values. Cells will be initialized to `init_value`.
    ///
    /// For arithmetic field backends, `addr_count` and `value_count` must be 1.
    fn init_ram(
        &mut self,
        size: usize,
        addr_count: usize,
        value_count: usize,
        init_value: &[Self::Wire],
    ) -> Result<RamId>;

    /// Read and return the value at `addr` in `ram`.
    fn ram_read(&mut self, ram: RamId, addr: &[Self::Wire]) -> Result<Vec<Self::Wire>>;

    /// Write `new` to `addr` in `ram`.
    fn ram_write(&mut self, ram: RamId, addr: &[Self::Wire], new: &[Self::Wire]) -> Result<()>;

    /// Finalize all operations on all RAMs.
    fn finalize_rams(&mut self) -> Result<()>;
}

impl<P: Party, C: AbstractChannel + Clone, SVOLE: SvoleT<P, F2, F40b>> BackendConvT<P>
    for DietMacAndCheese<P, F2, F40b, C, SVOLE>
{
    fn assert_conv_to_bits(&mut self, w: &Self::Wire) -> Result<Vec<Mac<P, F2, F40b>>> {
        debug!("CONV_TO_BITS {:?}", w);
        Ok(vec![*w])
    }

    fn assert_conv_from_bits(&mut self, x: &[Mac<P, F2, F40b>]) -> Result<Self::Wire> {
        Ok(x[0])
    }

    fn finalize_conv(&mut self) -> Result<()> {
        // We dont need to finalize the conversion
        // for the binary functionality because they are for free.
        Ok(())
    }
}

// this structure is for grouping the edabits with the same number of bits.
// This is necessary for example when F1 -> F2 and F3 -> F2, with F1 and F3
// requiring a different number of bits.
// NOTE: We use a BTreeMap instead of a HashMap so that the iterator is sorted over the keys, which
// is essential during `finalize`. The keys must be sorted deterministically so that the prover and verifier are in sync while finalizing
// the conversion for every bit width.
struct EdabitsMap<E>(BTreeMap<usize, Vec<E>>);

impl<E> EdabitsMap<E> {
    fn new() -> Self {
        EdabitsMap(BTreeMap::new())
    }

    fn push_elem(&mut self, bit_width: usize, e: E) {
        self.0
            .entry(bit_width)
            .and_modify(|v| v.push(e))
            .or_default();
    }

    fn get_edabits(&mut self, bit_width: usize) -> Option<&mut Vec<E>> {
        self.0.get_mut(&bit_width)
    }

    fn set_edabits(&mut self, bit_width: usize, edabits: Vec<E>) {
        self.0.insert(bit_width, edabits);
    }
}

pub(crate) struct DietMacAndCheeseConv<
    P: Party,
    FE: FiniteField,
    C: AbstractChannel + Clone,
    SvoleF2: SvoleT<P, F2, F40b>,
    SvoleFE: SvoleT<P, FE, FE>,
> {
    dmc: DietMacAndCheese<P, FE, FE, C, SvoleFE>,
    conv: Conv<P, FE, SvoleF2, SvoleFE>,
    dora_states: HashMap<usize, DoraState<P, FE, FE, C, SvoleFE>>,
    ram_states: Vec<ArithmeticRam<P, FE, FE, C, SvoleFE>>,
    edabits_map: EdabitsMap<Edabits<P, FE>>,
    dmc_f2: DietMacAndCheese<P, F2, F40b, C, SvoleF2>,
    no_batching: bool,
}

impl<
        P: Party,
        FE: PrimeFiniteField,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleFE: SvoleT<P, FE, FE>,
    > DietMacAndCheeseConv<P, FE, C, SvoleF2, SvoleFE>
{
    pub fn init(
        channel: &mut C,
        mut rng: AesRng,
        fcom_f2: &FCom<P, F2, F40b, SvoleF2>,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
        no_batching: bool,
    ) -> Result<Self> {
        let rng2 = rng.fork();
        let dmc = DietMacAndCheese::<P, FE, FE, C, SvoleFE>::init(
            channel,
            rng,
            lpn_setup,
            lpn_extend,
            no_batching,
        )?;
        let conv = Conv::init_with_fcoms(fcom_f2, &dmc.fcom)?;
        Ok(DietMacAndCheeseConv {
            dmc,
            conv,
            dora_states: Default::default(),
            ram_states: Default::default(),
            edabits_map: EdabitsMap::new(),
            dmc_f2: DietMacAndCheese::<P, F2, F40b, C, SvoleF2>::init_with_fcom(
                channel,
                rng2,
                fcom_f2,
                no_batching,
            )?,
            no_batching,
        })
    }

    pub fn init_with_svole(
        channel: &mut C,
        mut rng: AesRng,
        fcom_f2: &FCom<P, F2, F40b, SvoleF2>,
        svole2: SvoleFE,
        no_batching: bool,
    ) -> Result<Self> {
        let rng2 = rng.fork();

        let fcom = FCom::init_with_vole(svole2)?;
        let dmc = DietMacAndCheese::<P, FE, FE, C, SvoleFE>::init_with_fcom(
            channel,
            rng,
            &fcom,
            no_batching,
        )?;
        let conv = Conv::init_with_fcoms(fcom_f2, &dmc.fcom)?;
        Ok(DietMacAndCheeseConv {
            dmc,
            conv,
            dora_states: Default::default(),
            ram_states: Default::default(),
            edabits_map: EdabitsMap::new(),
            dmc_f2: DietMacAndCheese::<P, F2, F40b, C, SvoleF2>::init_with_fcom(
                channel,
                rng2,
                fcom_f2,
                no_batching,
            )?,
            no_batching,
        })
    }

    fn maybe_do_conversion_check(&mut self, bit_width: usize) -> Result<()> {
        let edabits = self.edabits_map.get_edabits(bit_width).unwrap();
        let num = edabits.len();
        if self.no_batching {
            self.conv.conv(
                &mut self.dmc.channel,
                &mut self.dmc.rng,
                CONVERSION_PARAM_B_SAFE,
                CONVERSION_PARAM_B_SAFE,
                edabits,
                None,
            )?;
            self.edabits_map.set_edabits(bit_width, vec![]);
        } else if num >= 2 * CONVERSION_BATCH_SIZE {
            // If there is more than twice the conversion_batch then lets do half of them, and keep the other
            // half for finalize so that it is still safe
            let index_to_split = edabits.len() - CONVERSION_BATCH_SIZE;
            let (_, edabits_to_process) = edabits.split_at(index_to_split);
            self.conv.conv(
                &mut self.dmc.channel,
                &mut self.dmc.rng,
                CONVERSION_PARAM_B,
                CONVERSION_PARAM_B,
                edabits_to_process,
                None,
            )?;
            edabits.truncate(index_to_split);
            assert_eq!(
                edabits.len(),
                self.edabits_map.get_edabits(bit_width).unwrap().len()
            );
        }

        Ok(())
    }
}

impl<
        P: Party,
        FE: PrimeFiniteField,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleFE: SvoleT<P, FE, FE>,
    > BackendT for DietMacAndCheeseConv<P, FE, C, SvoleF2, SvoleFE>
{
    type Wire = <DietMacAndCheese<P, FE, FE, C, SvoleFE> as BackendT>::Wire;
    type FieldElement = <DietMacAndCheese<P, FE, FE, C, SvoleFE> as BackendT>::FieldElement;

    fn wire_value(&self, wire: &Self::Wire) -> Option<Self::FieldElement> {
        self.dmc.wire_value(wire)
    }
    fn random(&mut self) -> Result<Self::FieldElement> {
        self.dmc.random()
    }
    fn copy(&mut self, wire: &Self::Wire) -> Result<Self::Wire> {
        self.dmc.copy(wire)
    }
    fn constant(&mut self, val: Self::FieldElement) -> Result<Self::Wire> {
        self.dmc.constant(val)
    }
    fn assert_zero(&mut self, wire: &Self::Wire) -> Result<()> {
        self.dmc.assert_zero(wire)
    }
    fn add(&mut self, a: &Self::Wire, b: &Self::Wire) -> Result<Self::Wire> {
        self.dmc.add(a, b)
    }
    fn sub(&mut self, a: &Self::Wire, b: &Self::Wire) -> Result<Self::Wire> {
        self.dmc.sub(a, b)
    }
    fn mul(&mut self, a: &Self::Wire, b: &Self::Wire) -> Result<Self::Wire> {
        self.dmc.mul(a, b)
    }
    fn add_constant(&mut self, a: &Self::Wire, b: Self::FieldElement) -> Result<Self::Wire> {
        self.dmc.add_constant(a, b)
    }
    fn mul_constant(&mut self, a: &Self::Wire, b: Self::FieldElement) -> Result<Self::Wire> {
        self.dmc.mul_constant(a, b)
    }

    fn input_public(&mut self, val: Self::FieldElement) -> Result<Self::Wire> {
        self.dmc.input_public(val)
    }
    fn input_private(&mut self, val: Option<Self::FieldElement>) -> Result<Self::Wire> {
        self.dmc.input_private(val)
    }
    fn finalize(&mut self) -> Result<()> {
        self.dmc.finalize()?;
        self.dmc_f2.finalize()?;
        Ok(())
    }
}

// Note: The restriction to a prime field is not caused by Dora
// This should be expanded in the future to allow disjunctions over extension fields.
impl<
        P: Party,
        FP: PrimeFiniteField + SieveIrDeserialize,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleFP: SvoleT<P, FP, FP>,
    > BackendDisjunctionT for DietMacAndCheeseConv<P, FP, C, SvoleF2, SvoleFP>
{
    fn finalize_disj(&mut self) -> Result<()> {
        for (_, disj) in std::mem::take(&mut self.dora_states) {
            disj.dora.finalize(&mut self.dmc)?;
        }
        Ok(())
    }

    fn disjunction(
        &mut self,
        inswit: &mut FieldInputs,
        fun_store: &FunStore,
        inputs: &[Self::Wire],
        disj: &DisjunctionBody,
    ) -> Result<Vec<Self::Wire>> {
        fn execute_branch<
            I: Iterator<Item = F>,
            P: Party,
            F: FiniteField,
            C: AbstractChannel + Clone,
            SvoleF: SvoleT<P, F, F>,
        >(
            ev: IsParty<P, Prover>,
            dmc: &mut DietMacAndCheese<P, F, F, C, SvoleF>,
            wit_tape: I,
            inputs: &[<DietMacAndCheese<P, F, F, C, SvoleF> as BackendT>::Wire],
            cond: usize,
            st: &mut DoraState<P, F, F, C, SvoleF>,
        ) -> Result<Vec<<DietMacAndCheese<P, F, F, C, SvoleF> as BackendT>::Wire>> {
            // currently only support 1 field element switch
            debug_assert_eq!(cond, 1);

            // so the guard is the last input
            let guard_val = inputs[inputs.len() - 1].value().into_inner(ev);

            // lookup the clause based on the guard
            let opt = *st
                .clause_resolver
                .as_ref()
                .into_inner(ev)
                .get(&guard_val)
                .expect("no clause guard is satisfied");

            st.dora
                .mux(dmc, wit_tape, inputs, ProverPrivateCopy::new(opt))
        }

        match self.dora_states.entry(disj.id()) {
            Entry::Occupied(mut entry) => match P::WHICH {
                WhichParty::Prover(ev) => execute_branch(
                    ev,
                    &mut self.dmc,
                    inswit.wit_iter::<FP::PrimeField>(),
                    inputs,
                    disj.cond() as usize,
                    entry.get_mut(),
                ),
                WhichParty::Verifier(ev) => entry.get_mut().dora.mux(
                    &mut self.dmc,
                    iter::empty(),
                    inputs,
                    ProverPrivateCopy::empty(ev),
                ),
            },
            Entry::Vacant(entry) => {
                // compile disjunction to the field
                let disjunction = Disjunction::compile(disj, disj.cond(), fun_store);

                let mut resolver: ProverPrivate<P, HashMap<FP, _>> =
                    ProverPrivate::new(Default::default());
                if let WhichParty::Prover(ev) = P::WHICH {
                    for (i, guard) in disj.guards().enumerate() {
                        let guard = FP::try_from_int(*guard).unwrap();
                        resolver.as_mut().into_inner(ev).insert(guard, i);
                    }
                }

                // create new Dora instance
                let dora = entry.insert(DoraState {
                    dora: Dora::new(disjunction),
                    clause_resolver: resolver,
                });

                // compute opt
                match P::WHICH {
                    WhichParty::Prover(ev) => execute_branch(
                        ev,
                        &mut self.dmc,
                        inswit.wit_iter::<FP::PrimeField>(),
                        inputs,
                        disj.cond() as usize,
                        dora,
                    ),
                    WhichParty::Verifier(ev) => dora.dora.mux(
                        &mut self.dmc,
                        iter::empty(),
                        inputs,
                        ProverPrivateCopy::empty(ev),
                    ),
                }
            }
        }
    }
}

impl<
        P: Party,
        FE: PrimeFiniteField,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleFE: SvoleT<P, FE, FE>,
    > BackendConvT<P> for DietMacAndCheeseConv<P, FE, C, SvoleF2, SvoleFE>
{
    fn assert_conv_to_bits(&mut self, a: &Self::Wire) -> Result<Vec<Mac<P, F2, F40b>>> {
        debug!("CONV_TO_BITS {:?}", a);
        let mut v;

        match P::WHICH {
            WhichParty::Prover(ev) => {
                let bits = a.value().into_inner(ev).bit_decomposition();

                v = Vec::with_capacity(bits.len());
                for b in bits {
                    let b2 = F2::from(b);
                    let mac = self.conv.fcom_f2.input1_prover(
                        ev,
                        &mut self.dmc.channel,
                        &mut self.dmc.rng,
                        b2,
                    )?;
                    v.push(Mac::new(ProverPrivateCopy::new(b2), mac));
                }
            }
            WhichParty::Verifier(ev) => {
                v = Vec::with_capacity(FE::NumberOfBitsInBitDecomposition::to_usize());
                for _ in 0..FE::NumberOfBitsInBitDecomposition::to_usize() {
                    let mac = self.conv.fcom_f2.input1_verifier(
                        ev,
                        &mut self.dmc.channel,
                        &mut self.dmc.rng,
                    )?;
                    v.push(mac);
                }
            }
        }

        less_than_eq_with_public(
            &mut self.dmc_f2,
            &v,
            (-FE::ONE)
                .bit_decomposition()
                .into_iter()
                .map(F2::from)
                .collect::<Vec<_>>()
                .as_slice(),
        )?;

        let r = v.clone();

        let bit_width = v.len();
        self.edabits_map
            .push_elem(bit_width, Edabits::<P, FE> { bits: v, value: *a });
        self.maybe_do_conversion_check(bit_width)?;

        Ok(r)
    }

    fn assert_conv_from_bits(&mut self, x: &[Mac<P, F2, F40b>]) -> Result<Self::Wire> {
        let mut power_twos = ProverPrivateCopy::new(FE::ONE);
        let mut recomposed_value = ProverPrivateCopy::new(FE::ZERO);

        let mut bits = Vec::with_capacity(x.len());

        for m in x {
            if let WhichParty::Prover(ev) = P::WHICH {
                *recomposed_value.as_mut().into_inner(ev) += (if m.value().into_inner(ev) == F2::ONE
                {
                    FE::ONE
                } else {
                    FE::ZERO
                }) * power_twos.into_inner(ev);
                power_twos
                    .as_mut()
                    .map(|power_twos| *power_twos += *power_twos);
            }

            bits.push(*m);
        }

        if let WhichParty::Prover(ev) = P::WHICH {
            debug!("CONV_FROM_BITS {:?}", recomposed_value.into_inner(ev));
        }

        let mac = <DietMacAndCheese<P, FE, FE, C, SvoleFE>>::input_private(
            &mut self.dmc,
            match P::WHICH {
                WhichParty::Prover(ev) => Some(recomposed_value.into_inner(ev)),
                WhichParty::Verifier(_) => None,
            },
        )?;

        let bit_width = bits.len();
        self.edabits_map
            .push_elem(bit_width, Edabits::<P, FE> { bits, value: mac });
        self.maybe_do_conversion_check(bit_width)?;

        Ok(mac)
    }

    fn finalize_conv(&mut self) -> Result<()> {
        // The keys must be sorted deterministically so that the prover and verifier are in sync.
        // This is the reason why the EdabitsMap is using a BTreeMap instead of a HashMap.
        for (_key, edabits) in self.edabits_map.0.iter() {
            // because they are periodically executed by maybe_do_conversion
            assert!(edabits.len() < 2 * CONVERSION_BATCH_SIZE);
            if edabits.len() < CONVERSION_BATCH_SIZE_SAFE {
                warn!(
                    "Insecure conversion check in finalize() because there are only {}, less than {}",
                    edabits.len(),
                    CONVERSION_BATCH_SIZE_SAFE,
                );
                self.conv.conv(
                    &mut self.dmc.channel,
                    &mut self.dmc.rng,
                    CONVERSION_PARAM_B_SAFE,
                    CONVERSION_PARAM_B_SAFE,
                    edabits,
                    None,
                )?;
            } else if edabits.len() >= CONVERSION_BATCH_SIZE_SAFE
                && edabits.len() < CONVERSION_BATCH_SIZE
            {
                self.conv.conv(
                    &mut self.dmc.channel,
                    &mut self.dmc.rng,
                    CONVERSION_PARAM_B_SAFE,
                    CONVERSION_PARAM_B_SAFE,
                    edabits,
                    None,
                )?;
            } else if edabits.len() >= CONVERSION_BATCH_SIZE
                && edabits.len() < (CONVERSION_BATCH_SIZE + CONVERSION_BATCH_SIZE_SAFE)
            {
                self.conv.conv(
                    &mut self.dmc.channel,
                    &mut self.dmc.rng,
                    CONVERSION_PARAM_B,
                    CONVERSION_PARAM_B,
                    edabits,
                    None,
                )?;
            } else {
                let index_to_split = CONVERSION_BATCH_SIZE;
                let (edabits1, edabits2) = edabits.split_at(index_to_split);
                self.conv.conv(
                    &mut self.dmc.channel,
                    &mut self.dmc.rng,
                    CONVERSION_PARAM_B,
                    CONVERSION_PARAM_B,
                    edabits1,
                    None,
                )?;

                // TODO: maybe split this in small chunks
                self.conv.conv(
                    &mut self.dmc.channel,
                    &mut self.dmc.rng,
                    CONVERSION_PARAM_B_SAFE,
                    CONVERSION_PARAM_B_SAFE,
                    edabits2,
                    None,
                )?;
            }
        }
        self.dmc.channel.flush()?;
        Ok(())
    }
}

impl<
        P: Party,
        FE: PrimeFiniteField,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleFE: SvoleT<P, FE, FE>,
    > BackendRamT for DietMacAndCheeseConv<P, FE, C, SvoleF2, SvoleFE>
{
    fn init_ram(
        &mut self,
        size: usize,
        addr_count: usize,
        value_count: usize,
        init_value: &[Self::Wire],
    ) -> Result<RamId> {
        debug_assert_eq!(addr_count, 1);
        debug_assert_eq!(value_count, 1);
        debug_assert_eq!(init_value.len(), 1);

        let ram_id = self.ram_states.len();
        self.ram_states
            .push(ArithmeticRam::new(size, init_value[0]));

        Ok(ram_id)
    }

    fn ram_read(&mut self, ram: RamId, addr: &[Self::Wire]) -> Result<Vec<Self::Wire>> {
        debug_assert!(ram < self.ram_states.len());
        debug_assert_eq!(addr.len(), 1);

        Ok(vec![self.ram_states[ram].read(&mut self.dmc, &addr[0])?])
    }

    fn ram_write(&mut self, ram: RamId, addr: &[Self::Wire], new: &[Self::Wire]) -> Result<()> {
        debug_assert!(ram < self.ram_states.len());
        debug_assert_eq!(addr.len(), 1);
        debug_assert_eq!(new.len(), 1);

        self.ram_states[ram].write(&mut self.dmc, &addr[0], &new[0])
    }

    fn finalize_rams(&mut self) -> Result<()> {
        self.ram_states
            .iter_mut()
            .try_for_each(|ram| ram.finalize(&mut self.dmc))
    }
}

// II) Instance/Witness/Relation/Gates/FunStore
// See circuit_ir.rs

// III Memory layout
// See memory.rs

// IV Evaluator for single field

/// A trait for evaluating circuits on a single field.
trait EvaluatorT<P: Party> {
    /// Evaluate a [`GateM`] alongside an optional instance and witness value.
    fn evaluate_gate(
        &mut self,
        gate: &GateM,
        instances: Option<Vec<Number>>,
        witnesses: Option<Vec<Number>>,
    ) -> Result<()>;

    /// Start the conversion for a [`ConvGate`].
    fn conv_gate_get(&mut self, gate: &ConvGate) -> Result<Vec<Mac<P, F2, F40b>>>;
    /// Finish the conversion for a [`ConvGate`].
    fn conv_gate_set(&mut self, gate: &ConvGate, bits: &[Mac<P, F2, F40b>]) -> Result<()>;

    /// Get the [`RamId`] at `wire_id`.
    fn get_ram_id(&mut self, wire_id: WireId) -> Result<RamId>;

    /// Set the [`RamId`] at `wire_id`.
    fn set_ram_id(&mut self, wire_id: WireId, ram_id: RamId) -> Result<()>;

    /// Evaluate a plugin call.
    ///
    /// `ram_id` is `Some(...)` iff `plugin` is a RAM read/write execution,
    /// and the return value is `Ok(Some(...))` iff `plugin` is a RAM init
    /// execution. Otherwise-successful evaluations should return `Ok(None)`.
    ///
    /// NOTE: This interface is a bit clunky, but then again, plugins are a bit
    /// clunky. A better design in anticipation of additional need for auxiliary
    /// inputs/outputs of plugins would be a type PluginExtras that is
    /// Option-like, and starting with a variant for RamIds. Since we don't
    /// anticipate additional new plugins being proposed, we go with the design
    /// as documented above.
    fn plugin_call_gate(
        &mut self,
        inswit: &mut FieldInputs,
        ram_id: Option<RamId>,
        fun_store: &FunStore,
        outputs: &[WireRange],
        inputs: &[WireRange],
        plugin: &PluginExecution,
    ) -> Result<Option<RamId>>;

    fn push_frame(&mut self, compiled_info: &CompiledInfo);
    fn pop_frame(&mut self);
    // TODO: Make allocate_slice return a result in case the operation violate some memory management
    fn allocate_slice(
        &mut self,
        src_first: WireId,
        src_last: WireId,
        start: WireId,
        count: WireId,
        allow_allocation: bool,
    );

    fn finalize(&mut self) -> Result<()>;
}

/// An [`EvaluatorT`] for RAM types.
///
/// This evaluator does nothing but manage the memory for RAM-typed wires. The
/// SIEVE IR RAM specification does not allow evaluation of any gates (besides
/// function calls) on RAM-typed wires, and RAM type conversions are undefined.
///
/// The evaluation of RAM operations and finalization is handled by the
/// underlying address/value field's [`EvaluatorSingle`], so `plugin_call_gate`
/// and `finalize` are also (counterintuitively) considered undefined for this
/// evaluator.
struct EvaluatorRam(Memory<RamId>);

impl<P: Party> EvaluatorT<P> for EvaluatorRam {
    fn evaluate_gate(
        &mut self,
        gate: &GateM,
        _instances: Option<Vec<Number>>,
        _witnesses: Option<Vec<Number>>,
    ) -> Result<()> {
        use GateM::*;

        match gate {
            New(_, first, last) => {
                self.0.allocation_new(*first, *last);
            }
            Delete(_, first, last) => {
                self.0.allocation_delete(*first, *last);
            }
            _ => unimplemented!("EvaluatorRam cannot evaluate gates other than new and delete."),
        }

        Ok(())
    }

    fn conv_gate_get(&mut self, _gate: &ConvGate) -> Result<Vec<Mac<P, F2, F40b>>> {
        unimplemented!("EvaluatorRam cannot convert between field types.")
    }

    fn conv_gate_set(&mut self, _gate: &ConvGate, _bits: &[Mac<P, F2, F40b>]) -> Result<()> {
        unimplemented!("EvaluatorRam cannot convert between field types.")
    }

    fn get_ram_id(&mut self, wire_id: WireId) -> Result<RamId> {
        Ok(*self.0.get(wire_id))
    }

    fn set_ram_id(&mut self, wire_id: WireId, ram_id: RamId) -> Result<()> {
        self.0.set(wire_id, &ram_id);
        Ok(())
    }

    fn plugin_call_gate(
        &mut self,
        _inswit: &mut FieldInputs,
        _ram_id: Option<RamId>,
        _fun_store: &FunStore,
        _outputs: &[WireRange],
        _inputs: &[WireRange],
        _plugin: &PluginExecution,
    ) -> Result<Option<RamId>> {
        unimplemented!("RAM operations are handled by the underlying address/value field.")
    }

    fn push_frame(&mut self, compiled_info: &CompiledInfo) {
        self.0.push_frame(compiled_info);
    }

    fn pop_frame(&mut self) {
        self.0.pop_frame();
    }

    fn allocate_slice(
        &mut self,
        src_first: WireId,
        src_last: WireId,
        start: WireId,
        count: WireId,
        allow_allocation: bool,
    ) {
        // RAMs may only be allocated one wire at a time.
        debug_assert_eq!(src_first, src_last);
        self.0
            .allocate_slice(src_first, src_last, start, count, allow_allocation);
    }

    fn finalize(&mut self) -> Result<()> {
        // RAM finalization is handled by the underlying address/value field.
        Ok(())
    }
}

/// A circuit evaluator for a single [`BackendT`].
///
/// The evaluator uses [`BackendT`] to evaluate the circuit, and uses `Memory`
/// to manage memory for the evaluation.
pub struct EvaluatorSingle<B: BackendT> {
    memory: Memory<<B as BackendT>::Wire>,
    backend: B,
    is_boolean: bool,
}

impl<B: BackendT> EvaluatorSingle<B> {
    fn new(backend: B, is_boolean: bool) -> Self {
        let memory = Memory::new();
        EvaluatorSingle {
            memory,
            backend,
            is_boolean,
        }
    }
}

impl<P: Party, B: BackendConvT<P> + BackendDisjunctionT + BackendLiftT + BackendRamT> EvaluatorT<P>
    for EvaluatorSingle<B>
where
    B::FieldElement: SieveIrDeserialize,
{
    // TODO: Revisit the type of instances / witnesses when we implement
    // streaming for these values. They should probably be
    // Option<&[Number]>!
    #[inline]
    fn evaluate_gate(
        &mut self,
        gate: &GateM,
        instances: Option<Vec<Number>>,
        witnesses: Option<Vec<Number>>,
    ) -> Result<()> {
        use GateM::*;

        match gate {
            Constant(_, out, value) => {
                let v = self
                    .backend
                    .constant(B::FieldElement::from_number(value)?)?;
                self.memory.set(*out, &v);
            }

            AssertZero(_, inp) => {
                let wire = self.memory.get(*inp);
                debug!("AssertZero wire: {wire:?}");
                if self.backend.assert_zero(wire).is_err() {
                    bail!("Assert zero fails on wire {}", *inp);
                }
            }

            Copy(_, out, inp) => {
                let mut curr_out = out.0;
                for curr_inp_range in inp.iter() {
                    for curr_inp in curr_inp_range.0..=curr_inp_range.1 {
                        let in_wire = self.memory.get(curr_inp);
                        let out_wire = self.backend.copy(in_wire)?;
                        self.memory.set(curr_out, &out_wire);
                        curr_out += 1;
                    }
                }
            }

            Add(_, out, left, right) => {
                let l = self.memory.get(*left);
                let r = self.memory.get(*right);
                let v = self.backend.add(l, r)?;
                self.memory.set(*out, &v);
            }

            Sub(_, out, left, right) => {
                let l = self.memory.get(*left);
                let r = self.memory.get(*right);
                let v = self.backend.sub(l, r)?;
                self.memory.set(*out, &v);
            }

            Mul(_, out, left, right) => {
                let l = self.memory.get(*left);
                let r = self.memory.get(*right);
                let v = self.backend.mul(l, r)?;
                self.memory.set(*out, &v);
            }

            AddConstant(_, out, inp, constant) => {
                let l = self.memory.get(*inp);
                let r = constant;
                let v = self
                    .backend
                    .add_constant(l, B::FieldElement::from_number(r)?)?;
                self.memory.set(*out, &v);
            }

            MulConstant(_, out, inp, constant) => {
                let l = self.memory.get(*inp);
                let r = constant;
                let v = self
                    .backend
                    .mul_constant(l, B::FieldElement::from_number(r)?)?;
                self.memory.set(*out, &v);
            }

            Instance(_, out) => {
                self.memory.allocate_possibly(out.0, out.1);
                let mut curr_out = out.0;
                for instance in instances.unwrap() {
                    let v = self
                        .backend
                        .input_public(B::FieldElement::from_number(&instance)?)?;
                    self.memory.set(curr_out, &v);
                    curr_out += 1;
                }
            }

            Witness(_, out) => {
                self.memory.allocate_possibly(out.0, out.1);
                for (i, curr_out) in (out.0..=out.1).enumerate() {
                    let w = witnesses
                        .as_ref()
                        .and_then(|wv| B::FieldElement::from_number(&wv[i]).ok());
                    let v = self.backend.input_private(w)?;
                    self.memory.set(curr_out, &v);
                }
            }
            New(_, first, last) => {
                self.memory.allocation_new(*first, *last);
            }
            Delete(_, first, last) => {
                self.memory.allocation_delete(*first, *last);
            }
            Call(_) => {
                panic!("Call should be intercepted earlier")
            }
            Conv(_) => {
                panic!("Conv should be intercepted earlier")
            }
            Challenge(_, out) => {
                let v = self.backend.random()?;
                let v = self.backend.input_public(v)?;
                self.memory.set(*out, &v);
            }
            Comment(_) => {
                panic!("Comment should be intercepted earlier")
            }
        }
        Ok(())
    }

    fn plugin_call_gate(
        &mut self,
        inswit: &mut FieldInputs,
        ram_id: Option<RamId>,
        fun_store: &FunStore,
        outputs: &[WireRange],
        inputs: &[WireRange],
        plugin: &PluginExecution,
    ) -> Result<Option<RamId>> {
        fn copy_mem<W>(mem: &Memory<W>, range: WireRange) -> impl Iterator<Item = &W>
        where
            W: Copy + Clone + Debug + Default,
        {
            let (start, end) = range;
            (start..=end).map(|i| mem.get(i))
        }

        match plugin {
            PluginExecution::PermutationCheck(plugin) => {
                assert_eq!(outputs.len(), 0);
                assert_eq!(inputs.len(), 2);
                let xs = copy_mem(&self.memory, inputs[0]).copied();
                let ys = copy_mem(&self.memory, inputs[1]).copied();
                plugin.execute::<_>(xs, ys, &mut self.backend)?
            }
            PluginExecution::Disjunction(disj) => {
                assert!(!inputs.is_empty(), "must provide condition");

                // retrieve input wires
                let mut wires = Vec::with_capacity(disj.inputs() as usize + disj.cond() as usize);

                // copy enviroment / inputs
                for range in inputs[1..].iter() {
                    wires.extend(copy_mem(&self.memory, *range));
                }

                // copy condition
                wires.extend(copy_mem(&self.memory, inputs[0]));

                // sanity check
                debug_assert_eq!(wires.len() as WireCount, disj.inputs() + disj.cond());

                // invoke disjunction implement on the backend
                let wires = self
                    .backend
                    .disjunction(inswit, fun_store, &wires[..], disj)?;
                debug_assert_eq!(wires.len() as u64, disj.outputs());

                // write back output wires
                let mut wires = wires.into_iter();
                for range in outputs {
                    for w in (range.0)..=(range.1) {
                        self.memory.set(w, &wires.next().unwrap())
                    }
                }
                debug_assert!(wires.next().is_none());
            }
            PluginExecution::Mux(plugin) => {
                plugin.execute::<P, B>(&mut self.backend, &mut self.memory)?
            }
            PluginExecution::Ram(plugin) => match plugin {
                RamVersion::RamBoolV0(RamBoolV0(Ram {
                    addr_count,
                    value_count,
                    op,
                    ..
                }))
                | RamVersion::RamBoolV1(RamBoolV1(Ram {
                    addr_count,
                    value_count,
                    op,
                    ..
                }))
                | RamVersion::RamArithV0(RamArithV0(Ram {
                    addr_count,
                    value_count,
                    op,
                    ..
                }))
                | RamVersion::RamArithV1(RamArithV1(Ram {
                    addr_count,
                    value_count,
                    op,
                    ..
                })) => match op {
                    RamOp::Init(size) => {
                        assert_eq!(inputs.len(), 1);
                        assert_eq!((inputs[0].1 - inputs[0].0 + 1) as usize, *value_count);
                        assert_eq!(outputs.len(), 1);
                        assert_eq!(outputs[0].0, outputs[0].1);

                        let init_value: Vec<_> =
                            copy_mem(&self.memory, inputs[0]).copied().collect();

                        let ram_id =
                            self.backend
                                .init_ram(*size, *addr_count, *value_count, &init_value)?;

                        return Ok(Some(ram_id));
                    }
                    RamOp::Read => {
                        assert_eq!(inputs.len(), 2);
                        assert_eq!(inputs[0].0, inputs[0].1);
                        assert_eq!((inputs[1].1 - inputs[1].0 + 1) as usize, *addr_count);
                        assert_eq!(outputs.len(), 1);
                        assert_eq!((outputs[0].1 - outputs[0].0 + 1) as usize, *value_count);

                        let ram_id = ram_id.ok_or_eyre("ram_id should be Some for Read/Write")?;
                        let addr: Vec<_> = copy_mem(&self.memory, inputs[1]).copied().collect();

                        let values = self.backend.ram_read(ram_id, &addr)?;
                        debug_assert_eq!(values.len(), *value_count);

                        for (w, value) in (outputs[0].0..=outputs[0].1).zip(values.into_iter()) {
                            self.memory.set(w, &value);
                        }
                    }
                    RamOp::Write => {
                        assert_eq!(inputs.len(), 3);
                        assert_eq!(inputs[0].0, inputs[0].1);
                        assert_eq!((inputs[1].1 - inputs[1].0 + 1) as usize, *addr_count);
                        assert_eq!((inputs[2].1 - inputs[2].0 + 1) as usize, *value_count);
                        assert_eq!(outputs.len(), 0);

                        let ram_id = ram_id.ok_or_eyre("ram_id should be Some for Read/Write")?;
                        let addr: Vec<_> = copy_mem(&self.memory, inputs[1]).copied().collect();
                        let new: Vec<_> = copy_mem(&self.memory, inputs[2]).copied().collect();

                        self.backend.ram_write(ram_id, &addr, &new)?
                    }
                },
            },
            _ => bail!("Plugin {plugin:?} is unsupported"),
        };
        Ok(None)
    }

    // The cases covered for field switching are:
    // 1) b <- x
    // 2) x <- b
    // 3) b0..b_n <- x   with n = log2(X)
    // 4) x <- b0..b_n   with n = log2(X)
    // 5) b0..b_n <- x   with n < log2(X)
    // 6) x <- b0..b_n   with n < log2(X)
    // 7) y <- x         with Y > X
    // 8) x <- y         with Y > X
    fn conv_gate_get(
        &mut self,
        (_, _, _, (start, end)): &ConvGate,
    ) -> Result<Vec<Mac<P, F2, F40b>>> {
        if *start != *end {
            if self.is_boolean {
                let mut v = Vec::with_capacity((end + 1 - start).try_into()?);
                for inp in *start..(*end + 1) {
                    let in_wire = self.memory.get(inp);
                    debug!("CONV GET {:?}", in_wire);
                    let bits = self.backend.assert_conv_to_bits(in_wire)?;
                    assert_eq!(bits.len(), 1);
                    v.push(bits[0]);
                }
                Ok(v.into_iter().rev().collect())
                // NOTE: Without reverse in case conversation gates are little-endian instead of big-endian
                //return Ok(v);
            } else {
                bail!("field switching from multiple wires on non-boolean field is not supported");
            }
        } else {
            let in_wire = self.memory.get(*start);
            debug!("CONV GET {:?}", in_wire);
            let bits = self.backend.assert_conv_to_bits(in_wire)?;
            debug!("CONV GET bits {:?}", bits);
            Ok(bits)
        }
    }

    fn conv_gate_set(
        &mut self,
        (_, (start, end), _, _): &ConvGate,
        bits: &[Mac<P, F2, F40b>],
    ) -> Result<()> {
        if *start != *end {
            if self.is_boolean {
                assert!((*end - *start + 1) as usize <= bits.len());

                for (i, _) in (*start..(*end + 1)).enumerate() {
                    let v = self.backend.assert_conv_from_bits(&[bits[i]])?;
                    debug!("CONV SET {:?}", v);
                    let out_wire = end - (i as WireId);
                    // NOTE: Without reverse in case conversation gates are little-endian instead of big-endian
                    // let out_wire = out1 + i as WireId;
                    self.memory.set(out_wire, &v);
                }
                Ok(())
            } else {
                bail!("field switching to multiple wires on non-boolean field is not supported");
            }
        } else {
            let v = self.backend.assert_conv_from_bits(bits)?;
            debug!("CONV SET {:?}", v);
            self.memory.set(*start, &v);
            Ok(())
        }
    }

    fn get_ram_id(&mut self, _wire_id: WireId) -> Result<RamId> {
        unimplemented!("This method should only be called on EvaluatorRams")
    }

    fn set_ram_id(&mut self, _wire_id: WireId, _ram_id: RamId) -> Result<()> {
        unimplemented!("This method should only be called on EvaluatorRams")
    }

    fn push_frame(&mut self, compiled_info: &CompiledInfo) {
        self.memory.push_frame(compiled_info);
    }

    fn pop_frame(&mut self) {
        self.memory.pop_frame();
    }

    fn allocate_slice(
        &mut self,
        src_first: WireId,
        src_last: WireId,
        start: WireId,
        count: WireId,
        allow_allocation: bool,
    ) {
        self.memory
            .allocate_slice(src_first, src_last, start, count, allow_allocation);
    }

    fn finalize(&mut self) -> Result<()> {
        debug!("Finalize in EvaluatorSingle");
        self.backend.finalize_conv()?;
        self.backend.finalize_disj()?;
        self.backend.finalize_rams()?;
        self.backend.finalize()?;
        Ok(())
    }
}

// V) Evaluator for multiple fields

/// A SIEVE IR circuit evaluator.
///
/// `EvaluatorCirc` provides everything needed to execute a zero-knowledge
/// proof described by SIEVE IR, provided either on disk via text or flatbuffers
/// or as in-memory resources as defined by [`crate::circuit_ir`]. A plaintext
/// evaluation mode (with limited functionality) is available.
///
/// A full example of using `EvaluatorCirc` as a library can be found in the
/// crate's `examples` directory. Briefly, the process is:
///
/// 1. Create an `EvaluatorCirc` (using the appropriate `new_*` method for your
///    situation)
/// 2. Call the corresponding `load_backends_*` method
/// 3. Evaluate a relation
///    - The [`crate::sieveir_reader_fbs`] and [`crate::sieveir_reader_text`]
///      modules can read SIEVE IR resources from disk (to provide the types,
///      function definitions, and gates)
///    - The module [`crate::circuit_ir`] can otherwise be used to construct
///      SIEVE IR resources. **Caution!** `EvaluatorCirc` expects all resources
///      to be syntactically valid and resource valid as defined by the SIEVE IR
///      specification. It is an unchecked runtime error to invoke
///      `EvaluatorCirc` methods on invalid SIEVE IR resources
/// 4. Call [`Self::terminate`] to stop and reset  any sVOLE functionality
pub struct EvaluatorCirc<
    P: Party,
    // See use of phantom for details on `C`
    C: AbstractChannel + Clone + 'static,
    SvoleF2: SvoleT<P, F2, F40b>,
    SvoleF2Ext: SvoleT<P, F40b, F40b>,
> {
    inputs: CircInputs,
    // Fcom <F2,F40b> functionality shared between the F2 backend and for field switching.
    // This is an optional to allow using none in the plaintext mode.
    // NOTE: If not in plaintext mode, in the current implementation it is always used,
    // which could be a TODO to only use it if an F2 backend is necessary or if any field switching is used.
    fcom_f2: Option<FCom<P, F2, F40b, SvoleF2>>,
    // Fcom functionality for full field voles on F40b. This is used for the permutation check on F2,
    // and for the disjunction plugin on F2.
    fcom_f2_ext: Option<FCom<P, F40b, F40b, SvoleF2Ext>>,
    type_store: TypeStore,
    eval: Vec<Box<dyn EvaluatorT<P>>>,
    rng: AesRng,
    multithreaded_voles: Vec<Box<dyn SvoleStopSignal>>,
    // Helper array used in `callframe_start` to avoid allocating a new array everytime the function is called.
    // This is an optimization.
    helper_callframe_arr: [WireId; 256],
    no_batching: bool,
    phantom: PhantomData<C>, // Storing the channel type here, is almost a convenience to prevent adding the type paramter to
                             // a handful of functions in the `impl`, maybe a TODO to eliminate the phatom and bring the type down to the
                             // functions. One benefit would be to not have to pass a dummy channel for the plaintext mode
}

impl<
        P: Party,
        C: AbstractChannel + Clone + 'static,
        SvoleF2: SvoleT<P, F2, F40b> + 'static,
        SvoleF2Ext: SvoleT<P, F40b, F40b> + 'static,
    > EvaluatorCirc<P, C, SvoleF2, SvoleF2Ext>
{
    /// Initialize a new (single-threaded) `EvaluatorCirc`.
    ///
    /// This requires an [`AbstractChannel`] and an [`AesRng`] to initialize the
    /// homomorphic commitment scheme for [`F2`] and its extension field
    /// [`F40b`]. This is suboptimal; ideally, these protocols would only start
    /// in circuits where they are strictly required (i.e., Boolean circuits or
    /// circuits where field conversion is used).
    ///
    /// The [`CircInputs`] paramater are instance and witness values as defined
    /// by the SIEVE IR; the [`TypeStore`] is the type environment as defined by
    /// the relation's header. The [`LpnSize`] configures the sVOLE protocol,
    /// and `no_batching`, if `true`, disables conversion check batching.
    ///
    /// For plaintext evaluation, see [`Self::new_plaintext`]. For multithreaded
    /// evaluation, see [`Self::new_multithreaded`].
    pub fn new(
        channel: &mut C,
        mut rng: AesRng,
        inputs: CircInputs,
        type_store: TypeStore,
        lpn_size: LpnSize,
        no_batching: bool,
    ) -> Result<Self> {
        let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
        let fcom_f2 = FCom::<_, F2, F40b, SvoleF2>::init(channel, &mut rng, lpn_setup, lpn_extend)?;

        let fcom_f2_ext = match P::WHICH {
            // For the verifier Fcom is initialized with delta
            WhichParty::Verifier(ev) => FCom::<_, F40b, F40b, SvoleF2Ext>::init_with_delta(
                channel,
                &mut rng,
                lpn_setup,
                lpn_extend,
                fcom_f2.get_delta().into_inner(ev),
            )?,
            WhichParty::Prover(_) => {
                FCom::<_, F40b, F40b, SvoleF2Ext>::init(channel, &mut rng, lpn_setup, lpn_extend)?
            }
        };

        Ok(EvaluatorCirc {
            inputs,
            fcom_f2: Some(fcom_f2),
            fcom_f2_ext: Some(fcom_f2_ext),
            type_store,
            eval: Vec::new(),
            rng,
            multithreaded_voles: vec![],
            helper_callframe_arr: [0; 256],
            no_batching,
            phantom: PhantomData,
        })
    }

    /// Initialize a new (plaintext) `EvaluatorCirc`.
    ///
    /// An alternative constructor, specifically curated for plaintext SIEVE IR
    /// evaluation (in particular: sVOLE initialization and instance/witness
    /// values are not required).
    ///
    /// For single-threaded zero-knowledge evaluation, use [`Self::new`]. For
    /// multithreaded zero-knowledge evaluation, use
    /// [`Self::new_multithreaded`].
    pub fn new_plaintext(inputs: CircInputs, type_store: TypeStore) -> Result<Self> {
        Ok(EvaluatorCirc {
            inputs,
            fcom_f2: None,
            fcom_f2_ext: None,
            type_store,
            eval: Vec::new(),
            rng: AesRng::new(),
            multithreaded_voles: vec![],
            helper_callframe_arr: [0; 256],
            no_batching: false, // unused
            phantom: PhantomData,
        })
    }

    /// Initialize a new (multithreaded) `EvaluatorCirc`.
    ///
    /// sVOLE tends to be where most time is spent during evaluation. To improve
    /// performance, `new_multithreaded` is provided as an alternative to
    /// [`Self::new`], the primary difference being that the F2 and F40b VOLEs
    /// are handled on (potentially many) separate threads.
    ///
    /// The returned `JoinHandle`s should, after [`Self::terminate`] is
    /// called, be `join`ed.
    pub fn new_multithreaded<C2: AbstractChannel + Clone + 'static + Send, I>(
        channels_vole: &mut I,
        threads_per_field: usize,
        mut rng: AesRng,
        inputs: CircInputs,
        type_store: TypeStore,
        no_batching: bool,
        lpn_size: LpnSize,
    ) -> Result<(
        EvaluatorCirc<
            P,
            C,
            SvoleAtomicRoundRobin<P, F2, F40b>,
            SvoleAtomicRoundRobin<P, F40b, F40b>,
        >,
        Vec<std::thread::JoinHandle<Result<()>>>,
    )>
    where
        I: Iterator<Item = C2>,
    {
        let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);

        let delta =
            // The verified generates a random delta for voles.
            match P::WHICH {
                WhichParty::Verifier(_ev) => Some(F40b::random(&mut rng)),
                WhichParty::Prover(_) => None,
            }
        ;

        let mut multithreaded_voles: Vec<Box<dyn SvoleStopSignal>> = vec![];

        // Setting up all the threads for F2 voles
        let (svole_round_robin, svoles_atomics, threads_f2) =
            SvoleAtomicRoundRobin::<P, F2, F40b>::create_and_spawn_svole_threads::<C2, I>(
                channels_vole,
                threads_per_field,
                lpn_setup,
                lpn_extend,
                delta,
            )?;
        let multithreaded_voles_f2: Vec<Box<dyn SvoleStopSignal>> = svoles_atomics
            .into_iter()
            .map::<Box<dyn SvoleStopSignal>, _>(|s| Box::new(s))
            .collect();
        let fcom_f2 = FCom::init_with_vole(svole_round_robin)?;

        // Setting up all the threads for full field F40b voles
        let (svole_round_robin, svoles_atomics, threads_f2_ext) =
            SvoleAtomicRoundRobin::<P, F40b, F40b>::create_and_spawn_svole_threads::<C2, I>(
                channels_vole,
                threads_per_field,
                lpn_setup,
                lpn_extend,
                delta,
            )?;
        let multithreaded_voles_f2_ext: Vec<Box<dyn SvoleStopSignal>> = svoles_atomics
            .into_iter()
            .map::<Box<dyn SvoleStopSignal>, _>(|s| Box::new(s))
            .collect();
        let fcom_f2_ext = FCom::init_with_vole(svole_round_robin)?;

        // Combining components from F2 voles and their full field F40b version
        multithreaded_voles.extend(multithreaded_voles_f2);
        multithreaded_voles.extend(multithreaded_voles_f2_ext);
        let mut threads = threads_f2;
        threads.extend(threads_f2_ext);

        Ok((
            EvaluatorCirc {
                inputs,
                fcom_f2: Some(fcom_f2),
                fcom_f2_ext: Some(fcom_f2_ext),
                type_store,
                eval: Vec::new(),
                rng,
                multithreaded_voles,
                helper_callframe_arr: [0; 256],
                no_batching,
                phantom: PhantomData,
            },
            threads,
        ))
    }

    /// Load the necessary backends for a (single-threaded) `EvaluatorCirc`.
    ///
    /// Using the `EvaluatorCirc`'s [`TypeStore`], prepare an appropriate set of
    /// backend evaluators to handle incoming relations. In most cases,
    /// `lpn_size` should be the same as what was passed to [`Self::new`].
    ///
    /// This method should be preferred over [`Self::load_backend`], which
    /// requires explicit use of [`std::any::TypeId`].
    ///
    /// For plaintext evaluation, use [`Self::load_backends_plaintext`]. For
    /// multithreaded evaluation, use [`Self::load_backends_multithreaded`].
    pub fn load_backends(&mut self, channel: &mut C, lpn_size: LpnSize) -> Result<()> {
        let type_store = self.type_store.clone();
        for (idx, spec) in type_store.iter() {
            let rng = self.rng.fork();
            match spec {
                TypeSpecification::Field(field) => {
                    self.load_backend(channel, rng, *field, *idx as usize, lpn_size)?;
                }
                TypeSpecification::Plugin(PluginType {
                    name,
                    operation,
                    params,
                }) => self.load_backend_plugin(&type_store, name, operation, params)?,
            }
        }
        Ok(())
    }

    fn load_backend_plugin(
        &mut self,
        type_store: &TypeStore,
        name: &str,
        operation: &str,
        params: &[PluginTypeArg],
    ) -> Result<()> {
        match name {
            "ram_bool_v0" => match operation {
                "ram" => {
                    ensure!(
                        params.len() == 6,
                        "v0 Boolean RAM types must specify six parameters, but {} were given.",
                        params.len()
                    );

                    // Validate type index refers to F2
                    let PluginTypeArg::Number(field_id) = params[0] else {
                        bail!("The v0 Boolean RAM type expects a number as its first parameter, but a string was found")
                    };
                    let field_id = u8::try_from(number_to_u64(&field_id)?)?;
                    let &TypeSpecification::Field(field_rust_id) = type_store.get(&field_id)?
                    else {
                        bail!("The first v0 Boolean RAM type parameter must refer to a field")
                    };
                    ensure!(
                        field_rust_id == std::any::TypeId::of::<F2>(),
                        "v0 Boolean RAMs only support Boolean addresses/values"
                    );

                    // Just check that the other parameters are numeric
                    for (i, p) in params.iter().enumerate().skip(1) {
                        if let PluginTypeArg::String(_) = p {
                            bail!("The v0 Boolean RAM type expects a number for parameter {}, but a string was found", i)
                        }
                    }

                    self.eval.push(Box::new(EvaluatorRam(Memory::new())));

                    Ok(())
                }
                _ => bail!("{} is not a type defined by {}", operation, name),
            },
            "ram_bool_v1" => match operation {
                "ram" => {
                    ensure!(
                        params.len() == 3,
                        "Boolean RAM types must specify three parameters, but {} were given",
                        params.len()
                    );

                    // Validate type index refers to F2
                    let PluginTypeArg::Number(field_id) = params[0] else {
                        bail!("The Boolean RAM type expects a number as its first parameter, but a string was found")
                    };
                    let field_id = u8::try_from(number_to_u64(&field_id)?)?;
                    let &TypeSpecification::Field(field_rust_id) = type_store.get(&field_id)?
                    else {
                        bail!("The first Boolean RAM type parameter must refer to a field")
                    };
                    ensure!(
                        field_rust_id == std::any::TypeId::of::<F2>(),
                        "Boolean RAMs only support Boolean addresses/values"
                    );

                    // Just check that the other parameters are numeric
                    if let PluginTypeArg::String(_) = params[1] {
                        bail!("The Boolean RAM type expects a number as its second parameter, but a string was found")
                    }

                    if let PluginTypeArg::String(_) = params[2] {
                        bail!("The Boolean RAM type expects a number as its third parameter, but a string was found")
                    }

                    self.eval.push(Box::new(EvaluatorRam(Memory::new())));

                    Ok(())
                }
                _ => bail!("{} is not a type defined by {}", operation, name),
            },
            "ram_arith_v0" => match operation {
                "ram" => {
                    ensure!(
                        params.len() == 4,
                        "v0 Arithmetic RAM types must specify four parameters, but {} were given",
                        params.len()
                    );

                    // Validate type index refers to a field
                    let PluginTypeArg::Number(field_id) = params[0] else {
                        bail!("The v0 arithmetic RAM type expects a number as its first parameter, but a string was found")
                    };
                    let field_id = u8::try_from(number_to_u64(&field_id)?)?;
                    if let TypeSpecification::Plugin(_) = type_store.get(&field_id)? {
                        bail!("The first v0 arithmetic RAM type parameter must refer to a field")
                    }

                    // Just check that the other parameters are numeric
                    for (i, p) in params.iter().enumerate().skip(1) {
                        if let PluginTypeArg::String(_) = p {
                            bail!("The v0 arithmetic RAM type expects a number for parameter {}, but a string was found", i)
                        }
                    }

                    self.eval.push(Box::new(EvaluatorRam(Memory::new())));

                    Ok(())
                }
                _ => bail!("{} is not a type defined by {}", operation, name),
            },
            "ram_arith_v1" => match operation {
                "ram" => {
                    ensure!(
                        params.len() == 1,
                        "Arithmetic RAM types must specify one parameter, but {} were given",
                        params.len()
                    );

                    // Validate type index refers to a field
                    let PluginTypeArg::Number(field_id) = params[0] else {
                        bail!("The arithmetic RAM type expects a number as its first parameter, but a string was found")
                    };
                    let field_id = u8::try_from(number_to_u64(&field_id)?)?;
                    if let TypeSpecification::Plugin(_) = type_store.get(&field_id)? {
                        bail!("The first arithmetic RAM type parameter must refer to a field")
                    }

                    self.eval.push(Box::new(EvaluatorRam(Memory::new())));

                    Ok(())
                }
                _ => bail!("{} is not a type defined by {}", operation, name),
            },
            _ => bail!("Plugin type not supported yet: {}:{}", name, operation),
        }
    }

    // All of these parameters are required to load a backend
    #[allow(clippy::too_many_arguments)]
    fn load_backend_fe<FE: PrimeFiniteField + StatisticallySecureField + SieveIrDeserialize>(
        &mut self,
        channel: &mut C,
        rng: AesRng,
        idx: usize,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
    ) -> Result<()> {
        assert!(idx == self.eval.len());
        let back: Box<dyn EvaluatorT<P>>;
        let dmc = DietMacAndCheeseConv::<P, FE, _, _, Svole<_, FE, FE>>::init(
            channel,
            rng,
            self.fcom_f2.as_ref().unwrap(),
            lpn_setup,
            lpn_extend,
            self.no_batching,
        )?;
        back = Box::new(EvaluatorSingle::new(dmc, false));
        self.eval.push(back);
        Ok(())
    }

    /// Load a backend for a (single-threaded) `EvaluatorCirc`.
    ///
    /// **CAUTION!** This method should be used with care, as it requires the
    /// use of [`std::any::TypeId`] to associate SIEVE IR type IDs with Rust
    /// types. If possible, [`Self::load_backends`] should be preferred over
    /// this method. If you must use this method, consider using
    /// [`crate::modulus_to_type_id`].
    ///
    /// For `EvaluatorCirc` created with [`Self::new_plaintext`], see
    /// [`Self::load_backend_plaintext`].
    ///
    /// DMC does not currently support client loading of individual backends
    /// when multithreading - use [`Self::load_backends_multithreaded`] instead.
    pub fn load_backend(
        &mut self,
        channel: &mut C,
        rng: AesRng,
        field: std::any::TypeId,
        idx: usize,
        lpn_size: LpnSize,
    ) -> Result<()> {
        // Loading the backends in order
        let back: Box<dyn EvaluatorT<P>>;
        if field == std::any::TypeId::of::<F2>() {
            info!("loading field F2");
            assert_eq!(idx, self.eval.len());

            // Note for F2 we do not use the backend with Conv but the one that can lift values to its extension field
            let dmc = DietMacAndCheeseExtField::<P, F40b, _, SvoleF2, SvoleF2Ext>::init_with_fcom(
                channel,
                rng,
                self.fcom_f2.as_ref().unwrap(),
                self.fcom_f2_ext.as_ref().unwrap(),
                self.no_batching,
            )?;
            back = Box::new(EvaluatorSingle::new(dmc, true));
            self.eval.push(back);
            Ok(())
        } else if field == std::any::TypeId::of::<F61p>() {
            info!("loading field F61p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_fe::<F61p>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F127p>() {
            info!("loading field F127p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_fe::<F127p>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F128p>() {
            info!("loading field F128p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_fe::<F128p>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<Secp256k1>() {
            info!("loading field Secp256k1");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_fe::<Secp256k1>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<Secp256k1order>() {
            info!("loading field Secp256k1order");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_fe::<Secp256k1order>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F384p>() {
            info!("loading field F384p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_fe::<F384p>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F384q>() {
            info!("loading field F384q");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_fe::<F384q>(channel, rng, idx, lpn_setup, lpn_extend)?;
            Ok(())
        } else {
            bail!("Unknown or unsupported field {:?}", field);
        }
    }

    /// Load the necessary backends for a (plaintext) `EvaluatorCirc`.
    ///
    /// Analogous to [`Self::load_backends`] for `EvaluatorCirc`s that have been
    /// created with [`Self::new_plaintext`].
    ///
    /// For single-threaded zero-knowledge loading, use [`Self::load_backends`].
    /// For multithreaded zero-knowedge loading, use
    /// [`Self::load_backends_multithreaded`].
    pub fn load_backends_plaintext(&mut self) -> Result<()> {
        let type_store = self.type_store.clone();
        for (idx, spec) in type_store.iter() {
            match spec {
                TypeSpecification::Field(field) => {
                    self.load_backend_plaintext(*field, *idx as usize)?;
                }
                _ => {
                    bail!("Type not supported yet: {:?}", spec);
                }
            }
        }
        Ok(())
    }

    fn load_backend_fe_plaintext<
        FE: PrimeFiniteField + StatisticallySecureField + SieveIrDeserialize,
    >(
        &mut self,
        idx: usize,
    ) -> Result<()> {
        assert!(idx == self.eval.len());
        let dmc = DietMacAndCheesePlaintext::<FE, FE>::new()?;
        let back = Box::new(EvaluatorSingle::new(dmc, false));
        self.eval.push(back);
        Ok(())
    }

    /// Load a backend for a (plaintext) `EvaluatorCirc`.
    ///
    /// Analogous to [`Self::load_backend`] for `EvaluatorCirc`s that have been
    /// created with [`Self::new_plaintext`].
    ///
    /// For `EvaluatorCirc` created with [`Self::new`], see [`Self::load_backend`].
    pub fn load_backend_plaintext(&mut self, field: std::any::TypeId, idx: usize) -> Result<()> {
        // Loading the backends in order
        let back: Box<dyn EvaluatorT<P>>;
        if field == std::any::TypeId::of::<F2>() {
            info!("loading field F2");
            assert_eq!(idx, self.eval.len());

            // For F2, it requires the extension field backend for the permutation_check
            let mut dmc = DietMacAndCheesePlaintext::<F2, F40b>::new()?;
            let ext_backend = DietMacAndCheesePlaintext::<F40b, F40b>::new()?;
            dmc.set_extfield_backend(ext_backend);
            back = Box::new(EvaluatorSingle::new(dmc, true));
            self.eval.push(back);
            Ok(())
        } else if field == std::any::TypeId::of::<F61p>() {
            info!("loading field F61p");
            self.load_backend_fe_plaintext::<F61p>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F127p>() {
            info!("loading field F127p");
            self.load_backend_fe_plaintext::<F127p>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F128p>() {
            info!("loading field F128p");
            self.load_backend_fe_plaintext::<F128p>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<Secp256k1>() {
            info!("loading field Secp256k1");
            self.load_backend_fe_plaintext::<Secp256k1>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<Secp256k1order>() {
            info!("loading field Secp256k1order");
            self.load_backend_fe_plaintext::<Secp256k1order>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F384p>() {
            info!("loading field F384p");
            self.load_backend_fe_plaintext::<F384p>(idx)?;
            Ok(())
        } else if field == std::any::TypeId::of::<F384q>() {
            info!("loading field F384q");
            self.load_backend_fe_plaintext::<F384q>(idx)?;
            Ok(())
        } else {
            bail!("Unknown or unsupported field {:?}", field);
        }
    }

    fn load_backend_multithreaded_f2(
        &mut self,
        channel: &mut C,
        rng: AesRng,
        _idx: usize,
    ) -> Result<()> {
        info!("loading field F2");

        let back: Box<dyn EvaluatorT<P>> = {
            let dmc = DietMacAndCheeseExtField::<P, F40b, _, SvoleF2, SvoleF2Ext>::init_with_fcom(
                channel,
                rng,
                self.fcom_f2.as_ref().unwrap(),
                self.fcom_f2_ext.as_ref().unwrap(),
                self.no_batching,
            )?;
            Box::new(EvaluatorSingle::new(dmc, true))
        };
        self.eval.push(back);
        Ok(())
    }

    // All of these parameters are required to load a backend
    #[allow(clippy::too_many_arguments)]
    fn load_backend_multithreaded_fe<
        FE: PrimeFiniteField + StatisticallySecureField + SieveIrDeserialize,
        C2: AbstractChannel + Clone + 'static + Send,
        I,
    >(
        &mut self,
        channel: &mut C,
        channels_vole: &mut I,
        threads_per_field: usize,
        mut rng: AesRng,
        idx: usize,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
    ) -> Result<Vec<std::thread::JoinHandle<Result<()>>>>
    where
        I: Iterator<Item = C2>,
    {
        assert!(idx == self.eval.len());
        let back: Box<dyn EvaluatorT<P>>;

        let delta =
            // The verifier generates a random delta for the voles.
            match P::WHICH {
                WhichParty::Verifier(_ev) => Some(FE::random(&mut rng)),
                WhichParty::Prover(_) => None,
            }        ;
        let (svole_round_robin, svoles_atomics, threads) =
            SvoleAtomicRoundRobin::<P, FE, FE>::create_and_spawn_svole_threads(
                channels_vole,
                threads_per_field,
                lpn_setup,
                lpn_extend,
                delta,
            )?;
        self.multithreaded_voles.extend(
            svoles_atomics
                .into_iter()
                .map::<Box<dyn SvoleStopSignal>, _>(|s| Box::new(s))
                .collect::<Vec<_>>(),
        );

        debug!("Starting DietMacAndCheese...");
        let dmc =
            DietMacAndCheeseConv::<P, FE, _, SvoleF2, SvoleAtomicRoundRobin<P, FE, FE>>::init_with_svole(
                channel,
                rng,
                self.fcom_f2.as_ref().unwrap(),
                svole_round_robin,
                self.no_batching,
            )?;
        back = Box::new(EvaluatorSingle::new(dmc, false));
        self.eval.push(back);
        Ok(threads)
    }

    #[allow(clippy::too_many_arguments)]
    fn load_backend_multi_any<C2: AbstractChannel + Clone + 'static + Send, I>(
        &mut self,
        channel: &mut C,
        channels_vole: &mut I,
        threads_per_field: usize,
        rng: AesRng,
        field: std::any::TypeId,
        idx: usize,
        lpn_size: LpnSize,
    ) -> Result<Vec<std::thread::JoinHandle<Result<()>>>>
    where
        I: Iterator<Item = C2>,
    {
        if field == std::any::TypeId::of::<F61p>() {
            info!("loading field F161p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_multithreaded_fe::<F61p, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<F127p>() {
            info!("loading field F127p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_multithreaded_fe::<F127p, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<F128p>() {
            info!("loading field F128p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size(lpn_size);
            self.load_backend_multithreaded_fe::<F128p, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<Secp256k1>() {
            info!("loading field Secp256k1");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_multithreaded_fe::<Secp256k1, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<Secp256k1order>() {
            info!("loading field Secp256k1order");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_multithreaded_fe::<Secp256k1order, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<F384p>() {
            info!("loading field F384p");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_multithreaded_fe::<F384p, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else if field == std::any::TypeId::of::<F384q>() {
            info!("loading field F384q");
            let (lpn_setup, lpn_extend) = mapping_lpn_size_large_field(lpn_size);
            self.load_backend_multithreaded_fe::<F384q, C2, I>(
                channel,
                channels_vole,
                threads_per_field,
                rng,
                idx,
                lpn_setup,
                lpn_extend,
            )
        } else {
            bail!("Unknown or unsupported field {:?}", field);
        }
    }

    /// Load the necessary backends for a (multithreaded) `EvaluatorCirc`.
    ///
    /// Analogous to [`Self::load_backends`] for `EvaluatorCirc`s that have been
    /// created with [`Self::new_multithreaded`].
    ///
    /// For single-threaded zero-knowledge loading, use [`Self::load_backends`].
    /// For single-threaded plaintext loading, use
    ///[`Self::load_backends_plaintext`].
    pub fn load_backends_multithreaded<C2: AbstractChannel + Clone + 'static + Send, I>(
        &mut self,
        channel: &mut C,
        channels_vole: &mut I,
        lpn_size: LpnSize,
        threads_per_field: usize,
    ) -> Result<Vec<std::thread::JoinHandle<Result<()>>>>
    where
        I: Iterator<Item = C2>,
    {
        let type_store = self.type_store.clone();
        let mut handles = vec![];
        for (idx, spec) in type_store.iter() {
            let rng = self.rng.fork();
            match spec {
                TypeSpecification::Field(field) => {
                    if *field == std::any::TypeId::of::<F2>() {
                        self.load_backend_multithreaded_f2(channel, rng, *idx as usize)?;
                    } else {
                        let more_handles = self.load_backend_multi_any(
                            channel,
                            channels_vole,
                            threads_per_field,
                            rng,
                            *field,
                            *idx as usize,
                            lpn_size,
                        )?;
                        handles.extend(more_handles);
                    }
                }
                TypeSpecification::Plugin(PluginType {
                    name,
                    operation,
                    params,
                }) => self.load_backend_plugin(&type_store, name, operation, params)?,
            }
        }
        Ok(handles)
    }

    /// Finalize all evaluation backends loaded by this `EvaluatorCirc`.
    ///
    /// This method only needs to be called if using
    /// [`Self::evaluate_gates_with_inputs`]; other evaluation methods take care
    /// of finalization.
    ///
    /// **CAUTION!** `finish` is different than [`Self::terminate`], which must
    /// be called to cleanly reset the `EvaluatorCirc`'s sVOLE functionalities!
    pub fn finish(&mut self) -> Result<()> {
        for i in 0..self.eval.len() {
            self.eval[i].finalize()?;
        }
        Ok(())
    }

    /// Evaluate a relation, given as a slice of [`GateM`].
    ///
    /// This method assumes that the `EvaluatorCirc` has been created with
    /// a `TypeStore`, `FunStore`, and `CircInputs`  consistent with the gates.
    ///
    /// Calls [`Self::finish`] after evaluation is complete.
    pub fn evaluate_gates(&mut self, gates: &[GateM], fun_store: &FunStore) -> Result<()> {
        self.evaluate_gates_passed(gates, fun_store)?;
        self.finish()
    }

    fn evaluate_gates_passed(&mut self, gates: &[GateM], fun_store: &FunStore) -> Result<()> {
        for gate in gates.iter() {
            self.eval_gate(gate, fun_store)?;
        }
        Ok(())
    }

    /// Evaluate a relation, given as a slice of [`GateM`].
    ///
    /// This is almost the same as `evaluate_gates`, except it uses `CircInputs`
    /// provided by the caller rather than that owned by the `EvaluatorCirc` (to
    /// allow for dynamic updates).
    ///
    /// **CAUTION!* Does **not** call [`Self::finish`]! Make sure to do so, else
    /// the protocols will continue waiting to execute additional steps.
    pub fn evaluate_gates_with_inputs(
        &mut self,
        gates: &[GateM],
        fun_store: &FunStore,
        inputs: &mut CircInputs,
    ) -> Result<()> {
        for gate in gates.iter() {
            self.eval_gate_with_inputs(gate, fun_store, inputs)?;
        }
        Ok(())
    }

    /// Evaluate a relation,  given as a path to a SIEVE IR flatbuffer file.
    ///
    /// Calls [`Self::finish`] after evaluation is complete.
    pub fn evaluate_relation(&mut self, path: &PathBuf) -> Result<()> {
        let mut buf_rel = BufRelation::new(path, &self.type_store)?;

        loop {
            let r = buf_rel.read_next();
            match r {
                None => {
                    break;
                }
                Some(()) => {
                    self.evaluate_gates_passed(&buf_rel.gates, &buf_rel.fun_store)?;
                }
            }
        }
        self.finish()
    }

    /// Evaluate a relation, provided as a text relation via `T: Read + Seek`.
    ///
    /// Calls [`Self::finish`] after evaluation is complete.
    pub fn evaluate_relation_text<T: Read + Seek>(&mut self, rel: T) -> Result<()> {
        let rel = RelationReader::new(rel)?;

        let mut buf_rel = TextRelation::new(self.type_store.clone(), FunStore::default());

        rel.read(&mut buf_rel)?;
        self.evaluate_gates_passed(&buf_rel.gates, &buf_rel.fun_store)?;

        self.finish()
    }

    fn callframe_start(
        &mut self,
        func: &FuncDecl,
        out_ranges: &[WireRange],
        in_ranges: &[WireRange],
    ) -> Result<()> {
        // 2)
        // We use the analysis on function body to find the types used in the body and only push a frame to those field backends.
        // TODO: currently push the size of args or vec without differentiating based on type.
        for ty in func.compiled_info.type_ids.iter() {
            self.eval[*ty as usize].push_frame(&func.compiled_info);
        }

        // For every type appearing in the function, initialize the counter to 0
        for i in func.compiled_info.type_ids.iter() {
            self.helper_callframe_arr[*i as usize] = 0;
        }

        let output_counts = func.output_counts();
        ensure!(
            out_ranges.len() == output_counts.len(),
            "Output range does not match output counts: {} != {}",
            out_ranges.len(),
            output_counts.len()
        );
        for i in 0..output_counts.len() {
            let (field_idx, count) = output_counts[i];
            let (src_first, src_last) = out_ranges[i];
            self.eval[field_idx as usize].allocate_slice(
                src_first,
                src_last,
                self.helper_callframe_arr[field_idx as usize],
                count,
                true,
            );
            self.helper_callframe_arr[field_idx as usize] += count;
        }

        let input_counts = func.input_counts();
        ensure!(
            in_ranges.len() == input_counts.len(),
            "Input range does not match input counts: {} != {}",
            in_ranges.len(),
            input_counts.len()
        );
        for i in 0..input_counts.len() {
            let (field_idx, count) = input_counts[i];
            let (src_first, src_last) = in_ranges[i];

            self.eval[field_idx as usize].allocate_slice(
                src_first,
                src_last,
                self.helper_callframe_arr[field_idx as usize],
                count,
                false,
            );
            self.helper_callframe_arr[field_idx as usize] += count;
        }
        Ok(())
    }

    fn callframe_end(&mut self, func: &FuncDecl) {
        // 4)
        // TODO: dont do the push blindly on all backends
        for ty in func.compiled_info.type_ids.iter() {
            self.eval[*ty as usize].pop_frame();
        }
    }

    #[inline]
    fn evaluate_call_gate(
        &mut self,
        fun_id: FunId,
        out_ranges: &[WireRange],
        in_ranges: &[WireRange],
        fun_store: &FunStore,
    ) -> Result<()> {
        let func = fun_store.get_func(fun_id)?;
        match &func.body() {
            FunctionBody::Gates(body) => {
                self.callframe_start(func, out_ranges, in_ranges)?;
                self.evaluate_gates_passed(body.gates(), fun_store)?;
                self.callframe_end(func);
            }
            FunctionBody::Plugin(body) => match body {
                PluginExecution::Gates(body) => {
                    self.callframe_start(func, out_ranges, in_ranges)?;
                    self.evaluate_gates_passed(body.gates(), fun_store)?;
                    self.callframe_end(func);
                }
                PluginExecution::PermutationCheck(plugin) => {
                    let type_id = plugin.type_id() as usize;
                    // The permutation plugin does not need to execute `callframe_start` or `callframe_end`
                    self.eval[type_id].plugin_call_gate(
                        self.inputs.get(type_id),
                        None,
                        fun_store,
                        out_ranges,
                        in_ranges,
                        body,
                    )?;
                }
                PluginExecution::Disjunction(plugin) => {
                    let type_id = plugin.field() as usize;
                    // disjunction does not use a callframe:
                    // since the inputs/outputs must be flattened to an R1CS witness.
                    self.eval[type_id].plugin_call_gate(
                        self.inputs.get(type_id),
                        None,
                        fun_store,
                        out_ranges,
                        in_ranges,
                        body,
                    )?;
                }
                PluginExecution::Mux(plugin) => {
                    let type_id = plugin.type_id() as usize;
                    self.callframe_start(func, out_ranges, in_ranges)?;
                    self.eval[type_id].plugin_call_gate(
                        self.inputs.get(type_id),
                        None,
                        fun_store,
                        out_ranges,
                        in_ranges,
                        body,
                    )?;
                    self.callframe_end(func);
                }
                PluginExecution::Ram(plugin) => {
                    // No callframes: RAMs are all 'global', and RAM-typed
                    // wires are to be thought of as (mutable) references.
                    match plugin {
                        RamVersion::RamArithV0(RamArithV0(Ram {
                            ram_type_id,
                            field_id,
                            op,
                            ..
                        }))
                        | RamVersion::RamArithV1(RamArithV1(Ram {
                            ram_type_id,
                            field_id,
                            op,
                            ..
                        }))
                        | RamVersion::RamBoolV0(RamBoolV0(Ram {
                            ram_type_id,
                            field_id,
                            op,
                            ..
                        }))
                        | RamVersion::RamBoolV1(RamBoolV1(Ram {
                            ram_type_id,
                            field_id,
                            op,
                            ..
                        })) => match op {
                            RamOp::Init(_) => {
                                let ram_id = self.eval[*field_id as usize]
                                    .plugin_call_gate(
                                        self.inputs.get(*field_id as usize),
                                        None,
                                        fun_store,
                                        out_ranges,
                                        in_ranges,
                                        body,
                                    )?
                                    .ok_or_eyre("RamInit should return a RamId")?;

                                // Safety: The above will fail if out_ranges is
                                // not exactly one wire range of length 1
                                self.eval[*ram_type_id as usize]
                                    .set_ram_id(out_ranges[0].0, ram_id)?;
                            }
                            RamOp::Read | RamOp::Write => {
                                assert!(!in_ranges.is_empty());

                                // Safety: in_ranges is not empty, and if any
                                // other invariants are violated, the subsequent
                                // call to plugin_call_gate will fail
                                let ram_id =
                                    self.eval[*ram_type_id as usize].get_ram_id(in_ranges[0].0)?;

                                self.eval[*field_id as usize].plugin_call_gate(
                                    self.inputs.get(*field_id as usize),
                                    Some(ram_id),
                                    fun_store,
                                    out_ranges,
                                    in_ranges,
                                    body,
                                )?;
                            }
                        },
                    }
                }
            },
        };

        Ok(())
    }

    fn eval_gate(&mut self, gate: &GateM, fun_store: &FunStore) -> Result<()> {
        debug!("GATE: {:?}", gate);
        match gate {
            GateM::Conv(gate) => {
                debug!("CONV IN");
                let (ty1, _, ty2, _) = gate.as_ref();
                // First we get the bits from the input and then we convert to the output.
                let bits = self.eval[*ty2 as usize].conv_gate_get(gate.as_ref())?;
                // then we convert the bits to the out field.
                self.eval[*ty1 as usize].conv_gate_set(gate.as_ref(), &bits)?;
                debug!("CONV OUT");
            }
            GateM::Instance(ty, out) => {
                let i = *ty as usize;
                self.eval[i].evaluate_gate(
                    gate,
                    self.inputs.pop_instances(i, out.1 - out.0 + 1),
                    None,
                )?;
            }
            GateM::Witness(ty, out) => {
                let i = *ty as usize;
                self.eval[i].evaluate_gate(
                    gate,
                    None,
                    self.inputs.pop_witnesses(i, out.1 - out.0 + 1),
                )?;
            }
            GateM::Call(arg) => {
                let (fun_id, out_ranges, in_ranges) = arg.as_ref();
                self.evaluate_call_gate(*fun_id, out_ranges, in_ranges, fun_store)?;
            }
            GateM::Comment(str) => {
                debug!("Comment: {:?}", str);
            }
            _ => {
                let ty = gate.type_id();
                self.eval[ty as usize].evaluate_gate(gate, None, None)?;
            }
        }
        Ok(())
    }

    // This function is a copy of `eval_gate` (added for Cybernetica/ZKSC) where the inputs
    // are passed because they could be dynamically updated.
    fn eval_gate_with_inputs(
        &mut self,
        gate: &GateM,
        fun_store: &FunStore,
        inputs: &mut CircInputs,
    ) -> Result<()> {
        debug!("GATE: {:?}", gate);
        match gate {
            GateM::Conv(gate) => {
                debug!("CONV IN");
                let (ty1, _, ty2, _) = gate.as_ref();
                // First we get the bits from the input and then we convert to the output.
                let bits = self.eval[*ty2 as usize].conv_gate_get(gate.as_ref())?;
                // then we convert the bits to the out field.
                self.eval[*ty1 as usize].conv_gate_set(gate.as_ref(), &bits)?;
                debug!("CONV OUT");
            }
            GateM::Instance(ty, out) => {
                let i = *ty as usize;
                self.eval[i].evaluate_gate(
                    gate,
                    inputs.pop_instances(i, out.1 - out.0 + 1),
                    None,
                )?;
            }
            GateM::Witness(ty, out) => {
                let i = *ty as usize;
                self.eval[i].evaluate_gate(
                    gate,
                    None,
                    inputs.pop_witnesses(i, out.1 - out.0 + 1),
                )?;
            }
            GateM::Call(arg) => {
                let (fun_id, out_ranges, in_ranges) = arg.as_ref();
                self.evaluate_call_gate(*fun_id, out_ranges, in_ranges, fun_store)?;
            }
            GateM::Comment(str) => {
                debug!("Comment: {:?}", str);
            }
            _ => {
                let ty = gate.type_id();
                self.eval[ty as usize].evaluate_gate(gate, None, None)?;
            }
        }
        Ok(())
    }

    /// Terminate the evaluator.
    ///
    /// This functions sends stop signals to all the registered sVOLE functionalities.
    pub fn terminate(&mut self) -> Result<()> {
        for e in self.multithreaded_voles.iter_mut() {
            info!("Sending stop signal");
            e.send_stop_signal()?;
        }
        self.multithreaded_voles = vec![];
        Ok(())
    }
}

impl<
        P: Party,
        C: AbstractChannel + Clone,
        SvoleF2: SvoleT<P, F2, F40b>,
        SvoleF2Ext: SvoleT<P, F40b, F40b>,
    > Drop for EvaluatorCirc<P, C, SvoleF2, SvoleF2Ext>
{
    fn drop(&mut self) {
        if !self.multithreaded_voles.is_empty() {
            warn!(
                "Need to call terminate()? There are {} multithreaded voles not terminated.",
                self.multithreaded_voles.len()
            );
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::TypeStore;
    use crate::svole_trait::Svole;
    use crate::LpnSize;
    use crate::{
        backend_multifield::EvaluatorCirc,
        fields::{F2_MODULUS, F61P_MODULUS, SECP256K1ORDER_MODULUS, SECP256K1_MODULUS},
    };
    use crate::{
        circuit_ir::{CircInputs, FunStore, FuncDecl, GateM, WireId, WireRange},
        fields::{F384P_MODULUS, F384Q_MODULUS},
    };
    use mac_n_cheese_sieve_parser::Number;
    use rand::SeedableRng;
    use scuttlebutt::field::F2;
    use scuttlebutt::field::{F384p, F384q, PrimeFiniteField};
    use scuttlebutt::field::{Secp256k1, Secp256k1order};
    use scuttlebutt::ring::FiniteRing;
    use scuttlebutt::SyncChannel;
    use scuttlebutt::{field::F61p, AesRng, Channel};
    use std::env;
    use std::net::TcpStream;
    use std::{collections::VecDeque, thread::JoinHandle};
    use std::{
        io::{BufReader, BufWriter},
        os::unix::net::UnixStream,
    };
    use swanky_party::{Prover, Verifier};

    pub(crate) const FF0: u8 = 0;
    const FF1: u8 = 1;
    const FF2: u8 = 2;
    const FF3: u8 = 3;

    pub(crate) fn zero<FE: PrimeFiniteField>() -> Number {
        FE::ZERO.as_int()
    }
    pub(crate) fn one<FE: PrimeFiniteField>() -> Number {
        FE::ONE.as_int()
    }
    pub(crate) fn two<FE: PrimeFiniteField>() -> Number {
        (FE::ONE + FE::ONE).as_int()
    }
    pub(crate) fn minus_one<FE: PrimeFiniteField>() -> Number {
        (-FE::ONE).as_int()
    }
    pub(crate) fn minus_two<FE: PrimeFiniteField>() -> Number {
        (-(FE::ONE + FE::ONE)).as_int()
    }
    pub(crate) fn three<FE: PrimeFiniteField>() -> Number {
        (FE::ONE + FE::ONE + FE::ONE).as_int()
    }
    pub(crate) fn minus_three<FE: PrimeFiniteField>() -> Number {
        (-(FE::ONE + FE::ONE + FE::ONE)).as_int()
    }
    pub(crate) fn four<FE: PrimeFiniteField>() -> Number {
        (FE::ONE + FE::ONE + FE::ONE + FE::ONE).as_int()
    }
    pub(crate) fn minus_four<FE: PrimeFiniteField>() -> Number {
        (-(FE::ONE + FE::ONE + FE::ONE + FE::ONE)).as_int()
    }
    pub(crate) fn minus_five<FE: PrimeFiniteField>() -> Number {
        (-(FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE)).as_int()
    }
    pub(crate) fn minus_nine<FE: PrimeFiniteField>() -> Number {
        (-(FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE + FE::ONE))
            .as_int()
    }

    fn wr(w: WireId) -> WireRange {
        (w, w)
    }

    // This function is useful when debugging tests.
    #[cfg(test)]
    #[allow(dead_code)]
    fn setup_logger() {
        // if log-level `RUST_LOG` not already set, then set to info
        match env::var("RUST_LOG") {
            Ok(val) => println!("loglvl: {}", val),
            Err(_) => env::set_var("RUST_LOG", "info"),
        };

        pretty_env_logger::init_timed();
    }

    pub(crate) fn test_circuit(
        fields: Vec<Number>,
        func_store: FunStore,
        gates: Vec<GateM>,
        instances: Vec<Vec<Number>>,
        witnesses: Vec<Vec<Number>>,
    ) -> eyre::Result<()> {
        let func_store_prover = func_store.clone();
        let gates_prover = gates.clone();
        let ins_prover = instances.clone();
        let wit_prover = witnesses;
        let type_store = TypeStore::try_from(fields.clone())?;
        let type_store_prover = type_store.clone();
        let (sender, receiver) = UnixStream::pair()?;
        let handle: JoinHandle<eyre::Result<()>> = std::thread::spawn(move || {
            let rng = AesRng::from_seed(Default::default());
            let reader = BufReader::new(sender.try_clone().unwrap());
            let writer = BufWriter::new(sender);
            let mut channel = Channel::new(reader, writer);

            let mut inputs = CircInputs::default();

            for (id, instances) in ins_prover.into_iter().enumerate() {
                inputs.ingest_instances(id, VecDeque::from(instances));
            }

            for (id, witnesses) in wit_prover.into_iter().enumerate() {
                inputs.ingest_witnesses(id, VecDeque::from(witnesses));
            }

            let mut eval = EvaluatorCirc::<Prover, _, Svole<_, _, _>, Svole<_, _, _>>::new(
                &mut channel,
                rng,
                inputs,
                type_store_prover,
                LpnSize::Small,
                false,
            )?;
            eval.load_backends(&mut channel, LpnSize::Small)?;
            eval.evaluate_gates(&gates_prover, &func_store_prover)?;
            eyre::Result::Ok(())
        });

        let rng = AesRng::from_seed(Default::default());
        let reader = BufReader::new(receiver.try_clone().unwrap());
        let writer = BufWriter::new(receiver);
        let mut channel = Channel::new(reader, writer);

        let mut inputs = CircInputs::default();

        for (id, inst) in instances.into_iter().enumerate() {
            inputs.ingest_instances(id, VecDeque::from(inst));
        }

        let mut eval = EvaluatorCirc::<Verifier, _, Svole<_, _, _>, Svole<_, _, _>>::new(
            &mut channel,
            rng,
            inputs,
            type_store,
            LpnSize::Small,
            false,
        )
        .unwrap();
        eval.load_backends(&mut channel, LpnSize::Small)?;
        eval.evaluate_gates(&gates, &func_store)?;

        handle.join().unwrap()
    }

    pub(crate) fn test_circuit_plaintext(
        fields: Vec<Number>,
        func_store: FunStore,
        gates: Vec<GateM>,
        instances: Vec<Vec<Number>>,
        witnesses: Vec<Vec<Number>>,
    ) -> eyre::Result<()> {
        let type_store = TypeStore::try_from(fields.clone())?;

        let mut inputs = CircInputs::default();

        for (id, instances) in instances.into_iter().enumerate() {
            inputs.ingest_instances(id, VecDeque::from(instances));
        }

        for (id, witnesses) in witnesses.into_iter().enumerate() {
            inputs.ingest_witnesses(id, VecDeque::from(witnesses));
        }

        let mut eval = EvaluatorCirc::<
            Prover,
            SyncChannel<BufReader<TcpStream>, BufWriter<TcpStream>>, // unnecessary type
            Svole<_, _, _>,
            Svole<_, _, _>,
        >::new_plaintext(inputs, type_store)?;
        eval.load_backends_plaintext()?;
        eval.evaluate_gates(&gates, &func_store)?;
        eyre::Result::Ok(())
    }

    fn test_conv_00() {
        // Test simple conversion from F61p to F2
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, wr(0)),
            GateM::Witness(FF0, wr(1)),
            GateM::Mul(FF0, 2, 0, 1),
            GateM::Conv(Box::new((FF1, wr(3), FF0, wr(2)))),
            GateM::AssertZero(FF1, 3),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![zero::<F61p>(), one::<F61p>()]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_01() {
        // Test simple conversion from F2 to F61p
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF1, wr(0)),
            GateM::Witness(FF1, wr(1)),
            GateM::Add(FF1, 2, 0, 1),
            GateM::Conv(Box::new((FF0, wr(3), FF1, wr(2)))),
            GateM::AddConstant(FF0, 4, 3, Box::from(minus_one::<F61p>())),
            GateM::AssertZero(FF0, 4),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![], vec![zero::<F2>(), one::<F2>()]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_02_twoway() {
        // Test that convert from F61p to F2 and from F2 to F61p works
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, (0, 1)),
            GateM::Mul(FF0, 2, 0, 1),
            GateM::Conv(Box::new((FF1, wr(3), FF0, wr(2)))),
            GateM::AssertZero(FF1, 3),
            GateM::Witness(FF1, (4, 5)),
            GateM::Add(FF1, 6, 5, 4),
            GateM::Conv(Box::new((FF0, wr(7), FF1, wr(6)))),
            GateM::AddConstant(FF0, 8, 7, Box::from(minus_one::<F61p>())),
            GateM::AssertZero(FF0, 8),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![
            vec![zero::<F61p>(), one::<F61p>()],
            vec![zero::<F2>(), one::<F2>()],
        ];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_binary_to_field() {
        // Test conversion from 2 bits to F61p
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF1, (0, 1)),
            GateM::Conv(Box::new((FF0, wr(3), FF1, (0, 1)))),
            GateM::AddConstant(FF0, 4, 3, Box::from(minus_three::<F61p>())),
            GateM::AssertZero(FF0, 4),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![], vec![one::<F2>(), one::<F2>()]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    // Same as `test_conv_binary_to_field()` but with plaintext backend
    fn test_conv_binary_to_field_plaintext() {
        // Test conversion from 2 bits to F61p
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF1, (0, 1)),
            GateM::Conv(Box::new((FF0, wr(3), FF1, (0, 1)))),
            GateM::AddConstant(FF0, 4, 3, Box::from(minus_three::<F61p>())),
            GateM::AssertZero(FF0, 4),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![], vec![one::<F2>(), one::<F2>()]];

        test_circuit_plaintext(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_field_to_binary() {
        // Test conversion from F61p to a vec of F2
        // 3 bit decomposition is 11000 on 5 bits, 00011
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, wr(0)),
            GateM::Conv(Box::new((FF1, (1, 5), FF0, wr(0)))),
            GateM::AssertZero(FF1, 1),
            GateM::AssertZero(FF1, 2),
            GateM::AssertZero(FF1, 3),
            GateM::AddConstant(FF1, 6, 4, Box::from(one::<F2>())),
            GateM::AddConstant(FF1, 7, 5, Box::from(one::<F2>())),
            GateM::AssertZero(FF1, 6),
            GateM::AssertZero(FF1, 7),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![three::<F61p>()], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    // Same as `test_conv_field_to_binary` but with plaintext evaluator
    fn test_conv_field_to_binary_plaintext() {
        // Test conversion from F61p to a vec of F2
        // 3 bit decomposition is 11000 on 5 bits, 00011
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, wr(0)),
            GateM::Conv(Box::new((FF1, (1, 5), FF0, wr(0)))),
            GateM::AssertZero(FF1, 1),
            GateM::AssertZero(FF1, 2),
            GateM::AssertZero(FF1, 3),
            GateM::AddConstant(FF1, 6, 4, Box::from(one::<F2>())),
            GateM::AddConstant(FF1, 7, 5, Box::from(one::<F2>())),
            GateM::AssertZero(FF1, 6),
            GateM::AssertZero(FF1, 7),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![three::<F61p>()], vec![]];

        test_circuit_plaintext(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_publics() {
        // Test conversion from F61p to a vec of F2 on public values
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Instance(FF1, (0, 3)),
            GateM::Conv(Box::new((FF0, wr(4), FF1, (0, 3)))),
            GateM::AddConstant(FF0, 5, 4, Box::from(minus_five::<F61p>())),
            GateM::AssertZero(FF0, 5),
        ];

        let instances = vec![
            vec![],
            vec![
                F2::ZERO.as_int(),
                F2::ONE.as_int(),
                F2::ZERO.as_int(),
                F2::ONE.as_int(),
            ],
        ];
        let witnesses = vec![vec![], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_shift() {
        // Test conversion and shift
        // 2 = 010000..., shifted as 10+010000...]= 10010000...] = 9, with truncation
        let fields = vec![F61P_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let mut gates = vec![
            GateM::New(FF0, 0, 0),
            GateM::Witness(FF0, wr(0)),
            GateM::New(FF1, 1, 61),
            GateM::Conv(Box::new((FF1, (1, 61), FF0, wr(0)))),
            GateM::New(FF1, 62, 122),
        ];
        gates.push(GateM::Copy(FF1, (62, 120), Box::new(vec![(3, 61)])));
        gates.push(GateM::Constant(FF1, 121, Box::new(zero::<F2>())));
        gates.push(GateM::Constant(FF1, 122, Box::new(one::<F2>())));
        gates.push(GateM::New(FF0, 123, 124));
        gates.push(GateM::Conv(Box::new((FF0, wr(123), FF1, (100, 122))))); // Beware!! truncate here, but that's only the zero upper bits
        gates.push(GateM::AddConstant(
            FF0,
            124,
            123,
            Box::from(minus_nine::<F61p>()),
        ));
        gates.push(GateM::AssertZero(FF0, 124));

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![two::<F61p>()], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    #[test]
    fn test_conv_ff_1() {
        let fields = vec![F61P_MODULUS, F384P_MODULUS, F384Q_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, (0, 1)),
            GateM::Mul(FF0, 2, 0, 1),
            GateM::Conv(Box::new((FF1, wr(3), FF0, wr(2)))),
            GateM::AssertZero(FF1, 3),
        ];

        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![vec![zero::<F61p>(), one::<F61p>()], vec![], vec![], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    #[test]
    fn test_conv_ff_2() {
        let fields = vec![F61P_MODULUS, F384P_MODULUS, F384Q_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            //GateM::New(FF3, 4, 4),
            //GateM::New(FF2, 5, 7),
            GateM::Witness(FF3, wr(4)),
            GateM::Conv(Box::new((FF2, wr(5), FF3, wr(4)))),
            GateM::Constant(FF2, 6, Box::from(minus_one::<F384q>())),
            GateM::Add(FF2, 7, 5, 6),
            GateM::AssertZero(FF2, 7),
        ];

        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![vec![], vec![], vec![], vec![one::<F2>()]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    #[test]
    fn test_conv_ff_3() {
        // tests that conversions from big fields to bools
        let fields = vec![F61P_MODULUS, F384P_MODULUS, F384Q_MODULUS, F2_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF2, wr(4)),
            GateM::Conv(Box::new((FF3, wr(5), FF2, wr(4)))),
            GateM::Witness(FF1, wr(1)),
            GateM::Conv(Box::new((FF3, wr(2), FF1, wr(1)))),
            GateM::Constant(FF3, 6, Box::from(minus_one::<F2>())),
            GateM::Add(FF3, 7, 5, 6),
            GateM::AssertZero(FF3, 7),
            GateM::AssertZero(FF3, 2),
        ];

        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![
            vec![],
            vec![F384p::ZERO.as_int()],
            vec![F384q::ONE.as_int()],
            vec![],
        ];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    #[test]
    fn test_conv_ff_4() {
        // test conversion from large field to smaller field
        let fields = vec![F61P_MODULUS, F384P_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF1, (0, 1)),
            GateM::Mul(FF1, 2, 0, 1),
            GateM::Conv(Box::new((FF0, wr(3), FF1, wr(2)))),
            GateM::AssertZero(FF0, 3),
            GateM::Add(FF1, 3, 1, 1),
            GateM::Add(FF1, 4, 3, 1),
            GateM::Conv(Box::new((FF0, wr(5), FF1, wr(4)))),
            GateM::AddConstant(FF0, 6, 5, Box::from(minus_three::<F61p>())),
            GateM::AssertZero(FF0, 6),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![vec![], vec![zero::<F61p>(), one::<F61p>()]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test_conv_ff_5() {
        // tests that conversions from big fields secp
        let fields = vec![SECP256K1_MODULUS, SECP256K1ORDER_MODULUS];
        let func_store = FunStore::default();

        let gates = vec![
            GateM::Witness(FF0, wr(0)),
            GateM::Conv(Box::new((FF1, wr(1), FF0, wr(0)))),
            GateM::Witness(FF1, wr(2)),
            GateM::Conv(Box::new((FF0, wr(3), FF1, wr(2)))),
            GateM::Constant(FF1, 4, Box::from(zero::<Secp256k1order>())),
            GateM::Add(FF1, 5, 1, 4),
            GateM::AssertZero(FF1, 5),
            GateM::Constant(FF0, 6, Box::from(minus_one::<Secp256k1>())),
            GateM::Add(FF0, 7, 3, 6),
            GateM::AssertZero(FF0, 7),
        ];

        let instances = vec![vec![], vec![]];
        let witnesses = vec![
            vec![Secp256k1::ZERO.as_int()],
            vec![Secp256k1order::ONE.as_int()],
            vec![],
        ];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test4_simple_fun() {
        // tests the simplest function

        let fields = vec![F61P_MODULUS];
        let mut func_store = FunStore::default();

        let gates_func = vec![GateM::Add(FF0, 0, 2, 4), GateM::Add(FF0, 1, 3, 5)];

        let mut func = FuncDecl::new_function(
            gates_func,
            vec![(FF0, 1), (FF0, 1)],
            vec![(FF0, 2), (FF0, 2)],
        );

        // The following instruction disable the vector optimization
        func.compiled_info.body_max = None;

        let fun_id = func_store.insert("myadd".into(), func).unwrap();

        let gates = vec![
            GateM::New(FF0, 0, 7), // TODO: Test when not all the New is done
            GateM::Witness(FF0, (0, 3)),
            GateM::Call(Box::new((
                fun_id,
                vec![(4, 4), (5, 5)],
                vec![(0, 1), (2, 3)],
            ))),
            GateM::Add(FF0, 6, 4, 5),
            GateM::AddConstant(
                FF0,
                7,
                6,
                Box::from((-(F61p::ONE + F61p::ONE + F61p::ONE + F61p::ONE)).as_int()),
            ),
            GateM::AssertZero(FF0, 7),
        ];

        let one = one::<F61p>();
        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![vec![one, one, one, one], vec![], vec![], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test4_simple_fun_plaintext() {
        // tests the simplest function

        let fields = vec![F61P_MODULUS];
        let mut func_store = FunStore::default();

        let gates_func = vec![GateM::Add(FF0, 0, 2, 4), GateM::Add(FF0, 1, 3, 5)];

        let mut func = FuncDecl::new_function(
            gates_func,
            vec![(FF0, 1), (FF0, 1)],
            vec![(FF0, 2), (FF0, 2)],
        );

        // The following instruction disable the vector optimization
        func.compiled_info.body_max = None;

        let fun_id = func_store.insert("myadd".into(), func).unwrap();

        let gates = vec![
            GateM::New(FF0, 0, 7), // TODO: Test when not all the New is done
            GateM::Witness(FF0, (0, 3)),
            GateM::Call(Box::new((
                fun_id,
                vec![(4, 4), (5, 5)],
                vec![(0, 1), (2, 3)],
            ))),
            GateM::Add(FF0, 6, 4, 5),
            GateM::AddConstant(
                FF0,
                7,
                6,
                Box::from((-(F61p::ONE + F61p::ONE + F61p::ONE + F61p::ONE)).as_int()),
            ),
            GateM::AssertZero(FF0, 7),
        ];

        let one = one::<F61p>();
        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![vec![one, one, one, one], vec![], vec![], vec![]];

        test_circuit_plaintext(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test5_simple_fun_with_vec() {
        // tests the simplest function with vec

        let fields = vec![F61P_MODULUS];
        let mut func_store = FunStore::default();

        let gates_fun = vec![
            GateM::Add(FF0, 6, 2, 4),
            GateM::AddConstant(FF0, 0, 6, Box::from(zero::<F61p>())),
            GateM::Add(FF0, 1, 3, 5),
        ];

        let func = FuncDecl::new_function(
            gates_fun,
            vec![(FF0, 1), (FF0, 1)],
            vec![(FF0, 2), (FF0, 2)],
        );
        let fun_id = func_store.insert("myadd".into(), func).unwrap();

        let gates = vec![
            GateM::New(FF0, 0, 7), // TODO: Test when not all the New is done
            GateM::Witness(FF0, (0, 3)),
            GateM::Call(Box::new((
                fun_id,
                vec![(4, 4), (5, 5)],
                vec![(0, 1), (2, 3)],
            ))),
            GateM::Add(FF0, 6, 4, 5),
            GateM::AddConstant(
                FF0,
                7,
                6,
                Box::from((-(F61p::ONE + F61p::ONE + F61p::ONE + F61p::ONE)).as_int()),
            ),
            GateM::AssertZero(FF0, 7),
        ];

        let one = one::<F61p>();
        let instances = vec![vec![], vec![], vec![], vec![]];
        let witnesses = vec![vec![one, one, one, one], vec![], vec![], vec![]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    fn test6_fun_slice_and_unallocated() {
        // tests a simple function passing instances in allocated slice and unallocated wire

        let fields = vec![F61P_MODULUS];
        let mut func_store = FunStore::default();

        let gates_func = vec![
            GateM::Copy(FF0, wr(0), Box::new(vec![wr(3)])),
            GateM::Add(FF0, 1, 4, 6),
            GateM::Copy(FF0, wr(2), Box::new(vec![wr(5)])),
        ];

        let mut func = FuncDecl::new_function(
            gates_func,
            vec![(FF0, 2), (FF0, 1)],
            vec![(FF0, 3), (FF0, 1)],
        );

        // The following instruction disable the vector optimization
        func.compiled_info.body_max = None;
        let fun_id = func_store.insert("myfun".into(), func).unwrap();

        let two = (F61p::ONE + F61p::ONE).as_int();
        let minus_four = (-(F61p::ONE + F61p::ONE + F61p::ONE + F61p::ONE)).as_int();
        let gates = vec![
            // New(0,2)
            // New(3,3)
            // Witness(0)  2
            // Witness(1)  2
            // Instance(2) 2
            // Witness(3)  2
            // 4..5, 6 <- Call(f, 0..2, 3)
            // AddConstant(7, 4, -2)
            // AddConstant(8, 5, -4)
            // AddConstant(9, 6, -2)
            // AssertZero(7)
            // AssertZero(8)
            // AssertZero(9)
            GateM::New(FF0, 0, 2),
            GateM::New(FF0, 3, 3),
            GateM::Witness(FF0, (0, 1)),
            GateM::Instance(FF0, wr(2)),
            GateM::Witness(FF0, wr(3)),
            GateM::Call(Box::new((
                fun_id,
                vec![(4, 5), (6, 6)],
                vec![(0, 2), (3, 3)],
            ))),
            GateM::AddConstant(FF0, 7, 4, Box::from(minus_two::<F61p>())),
            GateM::AddConstant(FF0, 8, 5, Box::from(minus_four)),
            GateM::AddConstant(FF0, 9, 6, Box::from(minus_two::<F61p>())),
            GateM::AssertZero(FF0, 7),
            GateM::AssertZero(FF0, 8),
            GateM::AssertZero(FF0, 9),
        ];

        let instances = vec![vec![two]];
        let witnesses = vec![vec![two, two, two]];

        test_circuit(fields, func_store, gates, instances, witnesses).unwrap();
    }

    #[test]
    fn test_multifield_conv() {
        test_conv_00();
        test_conv_01();
        test_conv_02_twoway();
        test_conv_binary_to_field();
        test_conv_field_to_binary();
        test_conv_publics();
        test_conv_shift();
    }

    #[test]
    fn test_multifield_conv_plaintext() {
        test_conv_binary_to_field_plaintext();
        test_conv_field_to_binary_plaintext();
    }

    #[test]
    fn test_multifield_ff_secp256() {
        test_conv_ff_5();
    }

    #[test]
    fn test_func() {
        test4_simple_fun();
        test4_simple_fun_plaintext();
        test5_simple_fun_with_vec();
        test6_fun_slice_and_unallocated()
    }
}
