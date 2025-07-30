//! Create split k-mers from sequences by rolling through input
//!
//! The [`SplitKmer`] struct stores a reference to then input sequence, the
//! necessary values to extract k-mers, and the current split k-mer:middle base
//! (as well as its reverse complement)
//!
//! Use [`SplitKmer::new()`] to initialise with sequence input from [`needletail`],
//! and [`SplitKmer::get_next_kmer()`] to advance through the sequence.
//!
//! (This could be made into a rust iterator, but I didn't do this when I wrote
//! it as I didn't know how yet.)

use std::borrow::Cow;
use std::cmp::Ordering;

use super::super::QualFilter;

use super::bit_encoding::*;
use super::nthash::NtHashIterator;
use crate::ska_dict::AscMinima;

use std::collections::VecDeque;
use std::process::exit;

/*
#[derive(Debug)]
pub struct AscMinima {
    /// Ascending minima FIFO queue
    mfifo : VecDeque<(u64, usize)>,
    hvec  : Vec<u64>,
}

impl AscMinima {
    /// Quality score is at least minimum.
    #[inline(always)]
    pub fn new(k: usize, g: usize) -> Self {
        let mut tmphvec = Vec::with_capacity(k - g - 1);
        for _ in 0..(k-g-1) {tmphvec.push(0)};
        Self {
            mfifo : VecDeque::<(u64, usize)>::with_capacity(k - g - 1),
//             hvec  : Vec::with_capacity(k - g - 1),
            hvec  : tmphvec,
        }
    }

//     pub fn clear(mut self) {
//         // We only clear the asc_minima FIFO queue, we don't care about the other.
//         self.asc_minima.clear();
//     }
}*/


/// Struct to generate all split k-mers from an input sequence
///
/// Holds reference to input sequence, current encoded k-mer, and other
/// information (k-mer size, masks etc.)
#[derive(Debug)]
pub struct SplitKmer<'a, IntT> {
    /// K-mer size
    k: usize,
    /// Mask to extract upper k-mer
    upper_mask: IntT,
    /// Mask to extract lower k-mer
    lower_mask: IntT,
    /// Reference to input sequence
    seq: Cow<'a, [u8]>,
    /// Size of seq
    seq_len: usize,
    /// Reference to sequence quality scores
    qual: Option<&'a [u8]>,
    /// How to filter on quality scores
    qual_filter: QualFilter,
    /// Minimum quality score to allow in a split k-mers
    min_qual: u8,
    /// Current index in input sequence
    index: usize,
    /// Current upper part of split k-mer
    upper: IntT,
    /// Current lower part of split k-mer
    lower: IntT,
    /// Current middle base
    middle_base: u8,
    /// Whether reverse complements are being used
    rc: bool,
    /// Current upper part of split k-mer, reverse complemented
    rc_upper: IntT,
    /// Current lower part of split k-mer, reverse complemented
    rc_lower: IntT,
    /// Current middle base, reverse complemented
    rc_middle_base: u8,
    /// Hash generator for reads
    hash_gen: Option<NtHashIterator>,
    /// Hash generator for extracting minimisers
    hash_min_gen: Option<NtHashIterator>,
    /// Minimiser size
    g : usize,
//     /// Ascending minima FIFO queue
//     asc_minima: VecDeque<(u64, usize)>,
    /// Ascending minima FIFO queue and utility vector
    asc_minima: &'a mut AscMinima,
}

impl<'a, IntT: for<'b> UInt<'b>> SplitKmer<'a, IntT> {
    /// Quality score is at least minimum.
    #[inline(always)]
    fn valid_qual(idx: usize, qual: Option<&'a [u8]>, min_qual: u8) -> bool {
        match qual {
            Some(qual_seq) => (qual_seq[idx] - 33) > min_qual, // ASCII encoding starts from b'!' = 33
            None => true,
        }
    }

    /// Build a new split k-mer at the given index.
    ///
    /// Called when initialised, or after skipping unknown bases. Returns
    /// [`None`] if the end of the input has been reached.
    #[allow(clippy::too_many_arguments)]
    fn build(
        seq: &[u8],
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        idx: &mut usize,
        qual_filter: &QualFilter,
        min_qual: u8,
        is_reads: bool,
        rc: bool,
        g: usize, // TEMP
//     ) -> Option<(IntT, IntT, u8, Option<NtHashIterator>, Option<NtHashIterator>, VecDeque<(u64, usize)>)> {
    ) -> Option<(IntT, IntT, u8, Option<NtHashIterator>, Option<NtHashIterator>)> {
        if *idx + k >= seq_len {
            return None;
        }
        let mut upper = IntT::zero_init();
        let mut lower = IntT::zero_init();
        let mut middle_base: u8 = 0;
        let middle_idx = (k + 1) / 2 - 1;
        let mut i = 0;
        while i < k {
            if valid_base(seq[i + *idx])
                && (*qual_filter != QualFilter::Strict
                    || Self::valid_qual(i + *idx, qual, min_qual))
            {
                // Checks for N or n
                let next_base = encode_base(seq[i + *idx]);
                match i.cmp(&middle_idx) {
                    Ordering::Greater => {
                        lower <<= 2;
                        lower |= IntT::from_encoded_base(next_base);
                    }
                    Ordering::Less => {
                        upper <<= 2;
                        upper |= IntT::from_encoded_base(next_base) << (middle_idx * 2);
                    }
                    Ordering::Equal => {
                        middle_base = next_base;
                    }
                }
                i += 1;
            } else {
                // Start again, skipping over N
                *idx += i + 1;
                if *idx + k >= seq_len {
                    return None;
                }
                upper = IntT::zero_init();
                lower = IntT::zero_init();
                middle_base = 0;
                i = 0;
            }
        }

        // For reads, start a rolling hash
        let hash_gen = if is_reads {
            Some(NtHashIterator::new(&seq[*idx..(*idx + k)], k, rc))
        } else {
            None
        };

//         println!("{:?}", &seq[*idx..(*idx + k)]);
//         println!("{:?}", &seq[(*idx + 1)..(*idx + 1 + g)]);
//         exit(1);

        // For reads, start a rolling hash for the minimisers TEMP
        let hash_min_gen = if is_reads {
            Some(NtHashIterator::new(&seq[(*idx + 1)..(*idx + 1 + g)], g, rc))
        } else {
            None
        };

        // ...and initialise an empty ring buffer for the ascending minima
//         let asc_minima = VecDeque::<(u64, usize)>::with_capacity(k - g - 1);

        *idx += k - 1;
//         Some((upper, lower, middle_base, hash_gen, hash_min_gen, asc_minima))
        Some((upper, lower, middle_base, hash_gen, hash_min_gen))
    }

    /// Checks if the split k-mer arms are palindromes, i.e. k-mer is its own reverse complement
    /// In this case the middle base needs ambiguity with its rc added.
    pub fn self_palindrome(&mut self) -> bool {
        self.rc && self.upper == self.rc_upper && self.lower == self.rc_lower
    }

    /// Update the stored reverse complement using the stored split-k and middle base
    fn update_rc(&mut self) {
        self.rc_upper = self.lower.rev_comp(self.k - 1) & self.upper_mask;
        self.rc_middle_base = rc_base(self.middle_base);
        self.rc_lower = self.upper.rev_comp(self.k - 1) & self.lower_mask;
    }

    /// Move forward to the next valid split k-mer
    ///
    /// Usually the next base, but if an N skips over it.
    /// If end of sequence encountered then returns `false`.
    fn roll_fwd(&mut self) -> bool {
        let mut success = false;
        self.index += 1;
        if self.index >= self.seq_len {
            return success;
        }
        let base = self.seq[self.index];
        if !valid_base(base)
            || (self.qual_filter == QualFilter::Strict
                && !Self::valid_qual(self.index, self.qual, self.min_qual))
        {
            let new_kmer = Self::build(
                &self.seq,
                self.seq_len,
                self.qual,
                self.k,
                &mut self.index,
                &self.qual_filter,
                self.min_qual,
                self.hash_gen.is_some(),
                self.rc,
                self.g,
            );
            if let Some(kmer_tuple) = new_kmer {
//                 (self.upper, self.lower, self.middle_base, self.hash_gen, self.hash_min_gen, self.asc_minima) = kmer_tuple;
                (self.upper, self.lower, self.middle_base, self.hash_gen, self.hash_min_gen) = kmer_tuple;
                if self.rc {
                    self.update_rc();
                }
                success = true;
            }
            self.asc_minima.mfifo.clear();
        } else {
            let half_k: usize = (self.k - 1) / 2;
            let new_base = encode_base(base);

            // Update the hash, if needed (i.e. if reads)
            if let Some(ref mut roll_hash) = self.hash_gen {
                let old_base = self.upper >> ((self.k - 2) * 2);
                roll_hash.roll_fwd(old_base.as_u8(), new_base);
            }

            // Update the hash for calculating minimisers if needed (i.e. if reads) // NOT OK NOW, TEMP
//             if self.hash_gen.is_some() { ///// OLD
//                 self.hash_min_gen = Some(NtHashIterator::new(&self.seq[(self.index - self.k + 1)..(self.index - self.k + 1 + self.g)], self.g, self.rc));
//             }
            if let Some(ref mut roll_hash) = self.hash_min_gen {
//                 println!("Before: {:?}", self.asc_minima.mfifo);
//                 println!("{:#066b}", self.lower);
//                 println!("{:#066b}", (self.lower >> (self.g * 2)).lsb_u8());
//                 println!("{:#066b}", self.lower.lsb_u8());
//                 exit(1);

                roll_hash.roll_fwd((self.lower >> (self.g * 2)).lsb_u8(), self.lower.lsb_u8());

                // After updating the rolling hash, now update the minima FIFO queue
                // Details in https://richardhartersworld.com/slidingmin/

//                 if self.asc_minima.front().expect("Empty ascending minima FIFO queue").1 == 0 {
//                     self.asc_minima.pop_front();
//                 }
//                 let newhash = roll_hash.curr_hash();
//                 let mut iq = self.asc_minima.len();
//                 let mut found : bool = false;
//                 while iq != 0 {
//                     iq -= 1;
//                     if !found && (newhash < self.asc_minima.back().expect("Empty ascending minima FIFO queue").0) {
//                         self.asc_minima.pop_back();
//                     } else {
//                         self.asc_minima.get_mut(iq).expect("Empty ascending minima FIFO queue").1 -= 1;
//                         found = true;
//                     }
//                 }
//                 self.asc_minima.push_back( (newhash, self.k - self.g - 1 - 1) );
// //                 println!("After:  {:?}", self.asc_minima);
// //                 exit(1);



//                 // NEW IMPL (with struct) BUT OLD (without ng)
//                 if self.asc_minima.mfifo.front().expect("Empty ascending minima FIFO queue").1 == 0 {
//                     self.asc_minima.mfifo.pop_front();
//                 }
//                 let newhash = roll_hash.curr_hash();
//                 let mut iq = self.asc_minima.mfifo.len();
//                 let mut found : bool = false;
//                 while iq != 0 {
//                     iq -= 1;
//                     if !found && (newhash < self.asc_minima.mfifo.back().expect("Empty ascending minima FIFO queue").0) {
//                         self.asc_minima.mfifo.pop_back();
//                     } else {
//                         self.asc_minima.mfifo.get_mut(iq).expect("Empty ascending minima FIFO queue").1 -= 1;
//                         found = true;
//                     }
//                 }
//                 self.asc_minima.mfifo.push_back( (newhash, self.k - self.g - 1 - 1) );
//
//             }

                // NEW IMPL (with struct) BUT NEW (with ng)
                if self.asc_minima.mfifo.front().expect("Empty ascending minima FIFO queue").1 == self.asc_minima.ng {
                    self.asc_minima.mfifo.pop_front();
                }
                let newhash = roll_hash.curr_hash();
//                 let mut iq = self.asc_minima.mfifo.len();
//                 while iq != 0 {
//                     iq -= 1;
//                     if newhash < self.asc_minima.mfifo.back().expect("Empty ascending minima FIFO queue").0 {
//                         self.asc_minima.mfifo.pop_back();
//                     } else {
//                         break;
//                     }
//                 }
                while (!self.asc_minima.mfifo.is_empty()) && (newhash < self.asc_minima.mfifo.back().unwrap().0) {
                    self.asc_minima.mfifo.pop_back();
                }
                self.asc_minima.mfifo.push_back( (newhash, self.asc_minima.ng + (self.k - self.g - 1)) );
                self.asc_minima.ng += 1;
//                 println!("After:  {:?}", self.asc_minima.mfifo);
            }
//             exit(1);

            // Update the k-mer
            self.upper = (self.upper << 2
                | (IntT::from_encoded_base(self.middle_base) << (half_k * 2)))
                & self.upper_mask;
            self.middle_base = (self.lower >> (2 * (half_k - 1))).as_u8();
            self.lower =
                ((self.lower << 2) | (IntT::from_encoded_base(new_base))) & self.lower_mask;
            if self.rc {
                self.rc_lower = (self.rc_lower >> 2
                    | ((IntT::from_encoded_base(self.rc_middle_base)) << (2 * (half_k - 1))))
                    & self.lower_mask;
                self.rc_middle_base = rc_base(self.middle_base);
                self.rc_upper = (self.rc_upper >> 2
                    | (IntT::from_encoded_base(rc_base(new_base))) << (2 * ((half_k * 2) - 1)))
                    & self.upper_mask;
            }
            success = true;
        }
        success
    }

    /// Create a [`SplitKmer`] iterator given reference to sequence input.
    ///
    /// Sequence, length and quality come from [`needletail`].
    ///
    /// Returns [`None`] if no valid split k-mers found in input (e.g. too short,
    /// no sequence, too many Ns).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        seq: Cow<'a, [u8]>,
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        rc: bool,
        min_qual: u8,
        qual_filter: QualFilter,
        is_reads: bool,
        g: usize, // TEMP
        asc_minima: &'a mut AscMinima, // TEMP
    ) -> Option<Self> {
        let mut index = 0;
        let first_kmer = Self::build(
            &seq,
            seq_len,
            qual,
            k,
            &mut index,
            &qual_filter,
            min_qual,
            is_reads,
            rc,
            g, //TEMP
        );

//         let asc_minima = AscMinima::new(k, g);

//         if let Some((upper, lower, middle_base, hash_gen, hash_min_gen, asc_minima)) = first_kmer {
        if let Some((upper, lower, middle_base, hash_gen, hash_min_gen)) = first_kmer {
            let (lower_mask, upper_mask) = IntT::generate_masks(k);
            let mut split_kmer = Self {
                k,
                upper_mask,
                lower_mask,
                seq_len,
                seq,
                qual,
                qual_filter,
                min_qual,
                upper,
                lower,
                middle_base,
                rc,
                rc_upper: IntT::zero_init(),
                rc_lower: IntT::zero_init(),
                rc_middle_base: 0,
                index,
                hash_gen,
                hash_min_gen, // TEMP
                g, // TEMP
                asc_minima, // TEMP
            };
            if rc {
                split_kmer.update_rc();
            }
            Some(split_kmer)
        } else {
            None
        }
    }

    /// Get the current split k-mer
    ///
    /// Returns split k-mer, middle base, and whether the reverse complement
    pub fn get_curr_kmer(&self) -> (IntT, u8, bool) {
        let split_kmer = self.upper | self.lower;
        // Some of the most useful prints for debugging left as comments here
        //let (upper, lower) = decode_kmer(self.k, split_kmer, self.upper_mask, self.lower_mask);
        //println!("{} {} {}", upper, lower, self.middle_base);
        if self.rc {
            let rc_split_kmer = self.rc_upper | self.rc_lower;
            //let (upper_rc, lower_rc) = decode_kmer(self.k, rc_split_kmer, self.upper_mask, self.lower_mask);
            //println!("{} {} {}", upper_rc, lower_rc, self.rc_middle_base);
            if split_kmer > rc_split_kmer {
                return (rc_split_kmer, self.rc_middle_base, true);
            }
        }
        (split_kmer, self.middle_base, false)
    }

    /// Get a `u64` hash of the current k-mer using [`NtHashIterator`]
    ///
    /// # Panics
    ///
    /// If called after creating with `is_reads = false`
    pub fn get_hash(&self) -> u64 {
        self.hash_gen
            .as_ref()
            .expect("Trying to get unitialised hash")
            .curr_hash()
    }

    /// Get the `u64` minimiser of the current k-mer
    ///
    ///
    ///
    ///
    pub fn get_minimiser(&mut self) -> u64 {
        if self.asc_minima.mfifo.is_empty() {
//             let thebits = ((self.upper >> (self.k - 3) | (IntT::from_encoded_base(self.middle_base))) << (self.k - 1)) | self.lower;
//             println!("Jijijiji");

            // Original debugging
    //         println!("{}, {}", self.get_curr_kmer().0, self.get_curr_kmer().1);

//             let (upper, lower) = decode_kmer(self.k, self.upper | self.lower, self.upper_mask, self.lower_mask);
//             println!("{} {} {}", upper, decode_base(self.middle_base) as char, lower);

    //         println!("{:#034b} {:#010b} {:#034b}", self.upper, self.middle_base, self.lower);
    //         println!("{:#034b} {:#010b} {:#034b}", self.upper >> 28, self.middle_base, self.lower);

    //
    //         println!("\nOriginal (upper) one\t\t {:#066b}",         self.upper);
    //         println!("With space for the middle one\t {:#066b}",    self.upper >> 28);
    //         println!("...PLUS the middle one\t\t {:#066b}",         self.upper >> 28 | (IntT::from_encoded_base(self.middle_base)));
    //         println!("Back again\t\t\t {:#066b}",                  (self.upper >> 28 | (IntT::from_encoded_base(self.middle_base))) << 28);
    //         println!("And this is the lower\t\t {:#066b}",          self.lower);
    //         println!("two to the left in the upper!\t {:#066b}",   (self.upper >> 28 | (IntT::from_encoded_base(self.middle_base))) << 30);
    //         println!("Now, if we sum...\t\t {:#066b}",            ((self.upper >> 28 | (IntT::from_encoded_base(self.middle_base))) << 30) | self.lower);
    //
    //         println!("{:#066b}", self.upper | (IntT::from_encoded_base(self.middle_base) << 28));
    //         println!("{:#066b}", self.upper | self.lower,);


            //quick (kinda) print
            //                                k-3                                                 k-1
    //         let mut thebits = ((self.upper >> 28 | (IntT::from_encoded_base(self.middle_base))) << 30) | self.lower;
//             let thebits = ((self.upper >> (self.k - 3) | (IntT::from_encoded_base(self.middle_base))) << (self.k - 1)) | self.lower;
    //         let mut thekmer = String::with_capacity(self.k); // k
    //         for _idx in 0..self.k {
    //             let base = decode_base(thebits.lsb_u8());
    //             thekmer.push(base as char);
    //             thebits >>= 2;
    //         }
    //         let thekmerstr: String = thekmer.chars().rev().collect::<String>();
    //         println!("Mine:\t\t{}", thekmerstr);
    //         println!("Original:\t{}{}{}", upper, decode_base(self.middle_base) as char, lower);


            ////// WITH STRINGS!!! TODO: bitwise, [[[ascending minima]]]
    /*
            let size_minimiser = 4;
            let mut min = String::from("ZZZZZZZ");

            for i in 0..(self.k - size_minimiser) {
                if thekmerstr[i..(i + size_minimiser)] < *min {
                    min = thekmerstr[i..(i + size_minimiser)].to_string();
                }
            }
    //         println!("{}", min);
            let minvec = min.into_bytes();
            let mut minint = 0u64;
    //         println!("{:?}", minvec);

            for i in 0..size_minimiser {
    //             println!("{}", minvec[i]);
    //             println!("{:#066b}", minvec[i]);
    //             println!("{:#066b}", encode_base(minvec[i]));
                minint <<= 2;
                minint |= encode_base(minvec[i]) as u64;
            }

            println!("{}", minint);
    //         println!("{:#066b}", minint);*/

            /*
            // BITWISE IMPLEMENTATION WITH HASH TODO: ascending minima /////OLLDDDDDDD
            println!("{:#066b}", thebits);
            println!("{:#066b}", (thebits >> (self.k - 2 - 1)          * 2).lsb_u8());
            println!("{:#066b}", (thebits >> (self.k - 2 - 1 - self.g) * 2).lsb_u8());
            let roll_hash  = &mut self.hash_min_gen.as_mut().unwrap();
            let mut minint = roll_hash.curr_hash();
    //         let mut dnbase : u8 = (thebits >> (self.k - 2 - 1)          * 2).lsb_u8();
    //         let mut upbase : u8 = (thebits >> (self.k - 2 - 1 - self.g) * 2).lsb_u8();
    //         println!("dnbase: {} upbase: {}", dnbase, upbase);
    //         println!("1st:{}\tNew:{}", decode_base(dnbase) as char, decode_base(upbase) as char);
    //         roll_hash.roll_fwd(dnbase, upbase);
            println!("{}", minint);
            println!("start");
            for ig in 0..(self.k - self.g - 1 - 1) {
                let dnbase = (thebits >> (self.k - 2          - ig) * 2).lsb_u8();
                let upbase = (thebits >> (self.k - 2 - self.g - ig) * 2).lsb_u8();
                roll_hash.roll_fwd(dnbase, upbase);
                let pothash = roll_hash.curr_hash();
                if pothash < minint {
                    minint = pothash;
                }
                println!("i:{}\t{}\t{}\tOld:{}\tNew:{}", ig, minint, pothash, decode_base(dnbase) as char, decode_base(upbase) as char);
            }


            println!("{}", minint);
            println!("{:#066b}", minint);
    //         exit(1);
            */



/*
            // IMPLEMENTATION WITH ASCENDING MINIMA ---- OLD
            // Review use of VecDeque!!!!!
//             println!("{:#066b}", thebits);
//             println!("{:#066b}", (thebits >> (self.k - 2 - 1)          * 2).lsb_u8());
//             println!("{:#066b}", (thebits >> (self.k - 2 - 1 - self.g) * 2).lsb_u8());
            let roll_hash  = &mut self.hash_min_gen.as_mut().expect("Not initialised rolling hash for minimisers");
            let mut minint = roll_hash.curr_hash();
            let mut ind    : usize = 0;
            let mut subind : usize = 0;
//             println!("jojo");
//             println!("{}", minint);
//             println!("k={}\tg={}\tk-1-1-g={}", self.k, self.g, self.k - self.g - 1 - 1);
//             println!("start");
//             let mut old_base;
//             let mut new_base;

//             let mut hash_vec : Vec<u64> = Vec::with_capacity(self.k - self.g - 1);
//             hash_vec.push(minint);

            self.asc_minima.hvec[subind] = minint;

//             println!("JOJOJO");

            let mut pothash;
            subind += 1;
//             for ig in 0..(self.k - self.g - 1 - 1) {
            while subind < self.k - self.g - 1 {
//                 let old_base = (thebits >> (self.k - 2          - ig) * 2).lsb_u8();
//                 let new_base = (thebits >> (self.k - 2 - self.g - ig) * 2).lsb_u8();
//                 roll_hash.roll_fwd(old_base, new_base);

                roll_hash.roll_fwd((thebits >> (self.k          - 1 - subind) * 2).lsb_u8(),
                                   (thebits >> (self.k - self.g - 1 - subind) * 2).lsb_u8());
                pothash = roll_hash.curr_hash();

                if pothash < minint {
                    minint = pothash;
                    ind    = subind;
                }
//                 hash_vec.push(pothash);
                self.asc_minima.hvec[subind] = pothash;
                subind += 1;

//                 println!("ig:{}\t{}\t{}\tOld:{}\tNew:{}", ig, minint, pothash, decode_base(old_base) as char, decode_base(new_base) as char);
//                 println!("abs ind:{}\t{}\t{}\tOld:{}\tNew:{}", ig + 1, minint, pothash,
//                 decode_base(old_base) as char, decode_base(new_base) as char);
            }
//             exit(1);

            // Construct asc. minima vector
            self.asc_minima.mfifo.push_back( (minint, ind) );
//             println!("asc_minima:\t{:?}", self.asc_minima);
            if ind != (self.k - self.g - 1 - 1) {
                ind += 1;
                let mut tmpmin;
                let mut minind;
                while ind < (self.k - self.g - 1) { // probably there is a smarter way of doing ALL of this
                    if ind == (self.k - self.g - 1 - 1) {
//                         println!("END, ind ={}", ind);
                        self.asc_minima.mfifo.push_back( (self.asc_minima.hvec[ind], ind) );
                        break;
                    } else {
                        tmpmin = self.asc_minima.hvec[ind];
                        minind = ind;
                        subind = ind + 1;
                        while subind < (self.k - self.g - 1) {
                            if self.asc_minima.hvec[subind] < tmpmin {
                                tmpmin = self.asc_minima.hvec[subind];
                                minind = subind;
                            }
                            subind += 1;
                        }
                        self.asc_minima.mfifo.push_back( (tmpmin, minind) );

                        ind = minind + 1;
                    }
                }
            }
//             println!("asc_minima:\t{:?}", self.asc_minima);
//             println!("{}", minint);
//             println!("{}", self.asc_minima.front().expect("Empty ascending minima FIFO queue").0);
//             println!("{:#066b}", minint);*/


            // IMPLEMENTATION WITH ASCENDING MINIMA ---- NEW, USE THE SAME ALG AS THE ONE USED TO UPDATE THE WINDOW
//             let thebits = ((self.upper >> (self.k - 3) | (IntT::from_encoded_base(self.middle_base))) << (self.k - 1)) | self.lower; // There's no difference in using this, and later rolling over them, and using directly the sequence, so let's use the sequence

//             println!("{:#066b}", thebits);
//             println!("{:#066b}", (thebits >> (self.k - 2 - 1)          * 2).lsb_u8());
//             println!("{:#066b}", (thebits >> (self.k - 2 - 1 - self.g) * 2).lsb_u8());
            let roll_hash  = &mut self.hash_min_gen.as_mut().expect("Not initialised rolling hash for minimisers");
            let mut minint = roll_hash.curr_hash();
            let ngs        = self.k - self.g - 1;
            let indseq     = self.index - self.k + 1;
            let mut subind : usize = 0;
//             println!("jojo");
//             println!("{}", minint);
//             println!("k={}\tg={}\tk-1-1-g={}\tindex={}", self.k, self.g, self.k - self.g - 1 - 1, self.index);
//             println!("start");
            self.asc_minima.mfifo.push_back( (minint, ngs) );

//             println!("JOJOJO");
//             println!("{:?}", self.asc_minima.mfifo);

            subind += 1;
            while subind < ngs {
//                 let old_base = (thebits >> (self.k - 1 - subind) * 2).lsb_u8();
//                 let new_base = (thebits >> (ngs - subind) * 2).lsb_u8();
//                 roll_hash.roll_fwd((thebits >> (self.k - 1 - subind) * 2).lsb_u8(),
//                                    (thebits >> (ngs        - subind) * 2).lsb_u8());
//                 println!("subind={}", subind);
                roll_hash.roll_fwd(encode_base(self.seq[indseq + subind]),
                                   encode_base(self.seq[indseq + subind + self.g]));
                minint = roll_hash.curr_hash();

//                 let mut iq = self.asc_minima.mfifo.len();
//                 while iq != 0 && minint < self.asc_minima.mfifo.back().unwrap().0 {
                while (!self.asc_minima.mfifo.is_empty()) && (minint < self.asc_minima.mfifo.back().unwrap().0) {
//                     iq -= 1;
                    self.asc_minima.mfifo.pop_back();
                }

//                 self.asc_minima.mfifo.push_back( (minint, subind) );
                self.asc_minima.mfifo.push_back( (minint, ngs + subind) );

//                 println!("ind:{}\t{}\tOld:{}\tNew:{}", subind, minint, decode_base(old_base) as char, decode_base(new_base) as char);
//                 println!("ind:{}\t{}\tOld:{}\tNew:{}", subind, minint, self.seq[indseq + subind] as char, self.seq[indseq + subind + self.g] as char);
//                 println!("{:?}", self.asc_minima.mfifo);

                subind += 1;
            }
            self.asc_minima.ng = ngs;
//             exit(1);
        }

    // In any case, the minimum for the k-mer will always be the first element of the fifo queue
    self.asc_minima.mfifo.front().expect("Empty ascending minima FIFO queue").0
    }

    /// Get the next split k-mer in the sequence
    ///
    /// Returns split k-mer, middle base, and whether the reverse complement
    /// or [`None`] if no next k-mer
    pub fn get_next_kmer(&mut self) -> Option<(IntT, u8, bool)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmer()),
            false => None,
        }
    }

    /// Get the index in the sequence of the current middle base
    pub fn get_middle_pos(&self) -> usize {
        let middle_idx = (self.k + 1) / 2 - 1;
        self.index - middle_idx
    }

    /// Returns true if the middle base passes the provided quality filtering criteria
    pub fn middle_base_qual(&self) -> bool {
        if self.qual.is_some() {
            match self.qual_filter {
                QualFilter::Middle | QualFilter::Strict => {
                    Self::valid_qual(self.get_middle_pos(), self.qual, self.min_qual)
                }
                QualFilter::NoFilter => true,
            }
        } else {
            true
        }
    }
}
