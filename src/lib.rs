use ndarray::prelude::*;
use std::mem::MaybeUninit;
use std::cmp::max;

pub struct GotohParameters {
    match_: isize,
    mismatch: isize,
    gap_opening: isize,
    gap_enlargement: isize
}

pub enum GotohMatrix {
    P,
    D,
    Q
}

pub struct GotohMatrices<F: Fn(usize, usize) -> bool> {
    params: GotohParameters,
    are_equal: F,
    p: Array2<isize>,
    d: Array2<isize>,
    q: Array2<isize>
}

pub struct Alignment {
    steps: Vec<(Option<usize>, Option<usize>)>
}

impl GotohParameters {
    pub fn new(match_: isize, mismatch: isize, gap_opening: isize, gap_enlargement: isize)
            -> Self {
        Self { match_, mismatch, gap_opening, gap_enlargement }
    }
}

impl<F: Fn(usize, usize) -> bool> GotohMatrices<F> {

    pub fn new(params: GotohParameters, seq1_len: usize, seq2_len: usize, are_equal: F) -> GotohMatrices<F> {
        let ncolumns = seq1_len;
        let nrows = seq2_len;

        let mut p_matrix = Array2::<isize>::uninit((nrows, ncolumns));
        let mut q_matrix = Array2::<isize>::uninit((nrows, ncolumns));
        let mut d_matrix = Array2::<isize>::uninit((nrows + 1, ncolumns + 1));

        d_matrix[[0, 0]] = MaybeUninit::new(0);
        for i in 1..nrows+1 {
            d_matrix[[i, 0]] = MaybeUninit::new(params.gap_opening + i as isize * params.gap_enlargement)
        }

        for j in 1..ncolumns+1 {
            d_matrix[[0, j]] = MaybeUninit::new(params.gap_opening + j as isize * params.gap_enlargement)
        }

        for i in 0..nrows {
            for j in 0..ncolumns {
                let new_p_value = {
                    let d_value = unsafe { d_matrix[[i, j+1]].assume_init() } + params.gap_opening + params.gap_enlargement;
                    if i == 0 { d_value } else {
                        let p_value = unsafe { p_matrix[[i-1, j]].assume_init() } + params.gap_enlargement;
                        max(d_value, p_value)
                    }
                };
                p_matrix[[i, j]] = MaybeUninit::new(new_p_value);

                let new_q_value = {
                    let d_value = unsafe { d_matrix[[i+1, j]].assume_init() }  + params.gap_opening + params.gap_enlargement;
                    if j == 0 { d_value } else {
                        let q_value = unsafe { q_matrix[[i, j-1]].assume_init() } + params.gap_enlargement;
                        max(d_value, q_value)
                    }
                };
                q_matrix[[i, j]] = MaybeUninit::new(new_q_value);

                let new_d_value = {
                    let d_value = unsafe { d_matrix[[i, j]].assume_init() };
                    let score = if are_equal(j, i) { params.match_ } else { params.mismatch };
                    max(max(d_value + score, new_p_value), new_q_value)
                };
                d_matrix[[i+1, j+1]] = MaybeUninit::new(new_d_value);
            }
        }

        unsafe {
            GotohMatrices {
                params,
                are_equal,
                p: p_matrix.assume_init(),
                d: d_matrix.assume_init(),
                q: q_matrix.assume_init(),
            }
        }
    }

    pub fn backtrack(&self) -> Alignment {
        let (p, d, q) = (&self.p, &self.d, &self.q);
        let score = |x, y| if (self.are_equal)(x, y) { self.params.match_ } else { self.params.mismatch };

        let (mut i, mut j) = self.p.dim();
        let mut matrix = GotohMatrix::D;
        let mut steps = vec![];

        while i != 0 && j != 0 {
            match matrix {
                GotohMatrix::D => {
                    if d[[i, j]] == d[[i-1, j-1]] + score(j-1, i-1) {
                        steps.push((Some(j-1), Some(i-1)));
                        i -= 1;
                        j -= 1;
                    } else if d[[i, j]] == p[[i-1, j-1]] {
                        matrix = GotohMatrix::P
                    } else if d[[i, j]] == q[[i-1, j-1]] {
                        matrix = GotohMatrix::Q
                    } else {
                        unreachable!()
                    }
                },
                GotohMatrix::P => {
                    if p[[i-1, j-1]] == d[[i-1, j]] + self.params.gap_opening + self.params.gap_enlargement {
                        matrix = GotohMatrix::D
                    } else if p[[i-1, j-1]] == p[[i-2, j-1]] + self.params.gap_enlargement {
                        matrix = GotohMatrix::P
                    } else {
                        unreachable!()
                    }
                    steps.push((None, Some(i-1)));
                    i -= 1
                },
                GotohMatrix::Q => {
                    if q[[i-1, j-1]] == d[[i, j-1]] + self.params.gap_opening + self.params.gap_enlargement {
                        matrix = GotohMatrix::D
                    } else if q[[i-1, j-1]] == q[[i-1, j-2]] + self.params.gap_enlargement {
                        matrix = GotohMatrix::Q
                    } else {
                        unreachable!()
                    }
                    steps.push((Some(j-1), None));
                    j -= 1
                }
            }
        }

        while j != 0 {
            steps.push((Some(j-1), None));
            j -= 1
        }

        while i != 0 {
            steps.push((None, Some(i-1)));
            i -= 1
        }

        Alignment { steps }
    }

    /*
    Worse than the sequential equivalent :/

    pub fn par_matrices<E: PartialEq + Sync, V: AsRef<[E]>>(&self, seq1: V, seq2: V)
            -> (Array2<isize>, Array2<isize>, Array2<isize>) {
        let seq1 = seq1.as_ref();
        let seq2 = seq2.as_ref();

        let ncolumns = seq1.len();
        let nrows = seq2.len();

        let mut p_matrix = Array2::<isize>::uninit((nrows, ncolumns));
        let mut q_matrix = Array2::<isize>::uninit((nrows, ncolumns));
        let mut d_matrix = Array2::<isize>::uninit((nrows + 1, ncolumns + 1));

        d_matrix[[0, 0]] = MaybeUninit::new(0);
        for i in 1..nrows+1 {
            d_matrix[[i, 0]] = MaybeUninit::new(self.gap_opening + i as isize * self.gap_enlargement)
        }

        for j in 1..ncolumns+1 {
            d_matrix[[0, j]] = MaybeUninit::new(self.gap_opening + j as isize * self.gap_enlargement)
        }

        let mut diags_buff: Vec<(isize, isize, isize)> = Vec::with_capacity(max(nrows, ncolumns));

        for l in 1..nrows + ncolumns {

            let (first_column, last_row) = if l <= nrows {
                (0, l)
            } else {
                (l - nrows, nrows)
            };

            let (first_row, last_column) = if l <= ncolumns {
                (0, l)
            } else {
                (l - ncolumns, ncolumns)
            };

            let columns_slice = first_column..last_column;
            let rows_slice = first_row..last_row;

            columns_slice.clone().into_par_iter().map(|j| {
                let i = l - 1 - j;
                let d_value = unsafe { d_matrix[[i, j+1]].assume_init() };
                let new_p_value = {
                    if i == 0 {
                        d_value + self.gap_opening + self.gap_enlargement
                    } else {
                        let p_value = unsafe { p_matrix[[i-1, j]].assume_init() };
                        max(d_value + self.gap_opening + self.gap_enlargement, p_value + self.gap_enlargement)
                    }
                };
                let new_q_value = {
                    let d_value = unsafe { d_matrix[[i+1, j]].assume_init() };
                    if j == 0 {
                        d_value + self.gap_opening + self.gap_enlargement
                    } else {
                        let q_value = unsafe { q_matrix[[i, j-1]].assume_init() };
                        max(d_value + self.gap_opening + self.gap_enlargement, q_value + self.gap_enlargement)
                    }
                };
                let new_d_value = {
                    let d_value = unsafe { d_matrix[[i, j]].assume_init() };
                    let score = if seq1[j] == seq2[i] { self.match_ } else { self.mismatch };
                    max(max(d_value + score, new_p_value), new_q_value)
                };
                (new_p_value, new_q_value, new_d_value)
            }).collect_into_vec(&mut diags_buff);

            let mut m = p_matrix.slice_mut(s![rows_slice.clone(), columns_slice.clone()]);
            m.invert_axis(m.max_stride_axis());
            let p_diag: Vec<_> = m.diag_mut().into_par_iter().collect();
            p_diag.into_par_iter().enumerate().for_each(|(k, v)| {
                *v = MaybeUninit::new(diags_buff[k].0)
            });

            let mut m = q_matrix.slice_mut(s![rows_slice.clone(), columns_slice.clone()]);
            m.invert_axis(m.max_stride_axis());
            let q_diag: Vec<_> = m.diag_mut().into_par_iter().collect();
            q_diag.into_par_iter().enumerate().for_each(|(k, v)| {
                *v = MaybeUninit::new(diags_buff[k].1)
            });

            let mut d_matrix_view = d_matrix.slice_mut(s![1.., 1..]);
            let mut m = d_matrix_view.slice_mut(s![rows_slice, columns_slice]);
            m.invert_axis(m.max_stride_axis());
            let d_diag: Vec<_> = m.diag_mut().into_par_iter().collect();
            d_diag.into_par_iter().enumerate().for_each(|(k, v)| {
                *v = MaybeUninit::new(diags_buff[k].2)
            })
        }

        unsafe {
            (p_matrix.assume_init(), q_matrix.assume_init(), d_matrix.assume_init())
        }
    }*/
}

impl Alignment {
    pub fn new<P: PartialEq, V: AsRef<[P]>>(params: GotohParameters, seq1: V, seq2: V) -> Alignment {
        let seq1 = seq1.as_ref();
        let seq2 = seq2.as_ref();
        let are_equal = |x, y| seq1[x] == seq2[y];
        Self::new_with_equality_fn(params, seq1.len(), seq2.len(), are_equal)
    }

    pub fn new_with_equality_fn<F: Fn(usize, usize) -> bool>(params: GotohParameters, seq1_len: usize, seq2_len: usize, are_equal: F) -> Alignment {
        let matrices: GotohMatrices<_> = GotohMatrices::new(params, seq1_len, seq2_len, are_equal);
        matrices.backtrack()
    }

    pub fn steps(&self) -> &Vec<(Option<usize>, Option<usize>)> {
        &self.steps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrices() {
        let params = GotohParameters::new(1, -1, -3, -1);
        let seq1 = vec!['C', 'C', 'G', 'A'];
        let seq2 = vec!['C', 'G'];
        let are_equal = |x, y| seq1[x] == seq2[y];
        let matrices: GotohMatrices<_> = GotohMatrices::new(params, seq1.len(), seq2.len(), are_equal);
        println!("{}", matrices.p);
        println!("{}", matrices.d);
        println!("{}", matrices.q);

        let p = array![[-8, -9, -10, -11],
                       [-3, -7,  -8,  -9]];
        let d = array![[0, -4, -5, -6, -7],
                       [-4, 1, -3, -4, -5],
                       [-5, -3, 0, -2, -5]];
        let q = array![[-8, -3, -4, -5],
                       [-9, -7, -4, -5]];

        assert_eq!(p, matrices.p);
        assert_eq!(d, matrices.d);
        assert_eq!(q, matrices.q);
    }

    #[test]
    fn alignment() {
        let params = GotohParameters::new(1, -3, -3, -1);
        let seq1 = vec!['C', 'C', 'G', 'A'];
        let seq2 = vec!['C', 'G'];
        let alignment = Alignment::new(params, &seq1, &seq2);
        let (s1, s2): (Vec<_>, Vec<_>) = alignment.steps.into_iter().rev().unzip();

        let s1: String = s1.into_iter().map(|s| match s {
            Some(i) => seq1[i],
            None => '-'
        }).collect();
        println!("{}", s1);

        let s2: String = s2.into_iter().map(|s| match s {
            Some(i) => seq2[i],
            None => '-'
        }).collect();
        println!("{}", s2);

        assert_eq!(s1, "CCGA");
        assert_eq!(s2, "-CG-");
    }
}
