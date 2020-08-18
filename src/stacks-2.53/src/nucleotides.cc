#include "constants.h"
#include "nucleotides.h"

const Nt4 Nt4::$ (0);  // 0000 (0)
const Nt4 Nt4::a (1);  // 0001 (1)
const Nt4 Nt4::c (2);  // 0010 (2)
const Nt4 Nt4::g (4);  // 0100 (4)
const Nt4 Nt4::t (8);  // 1000 (8)
const Nt4 Nt4::n (15); // 1111 (15)
const array<Nt4,5> Nt4::all = {{a, c, g, t, n}};

const Nt4 Nt4::from_ch[256] = {
    0,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n, // 0x
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n, // 1x
    n,n,n,n, 0,n,n,n, n,n,n,n, n,n,n,n, // 2x
    n,n,n,n, n,n,n,n, n,n,n,n, n,0,n,n, // 3x rem. 0 for 0x3D ('=') because htslib does it...
    n,a,n,c, n,n,n,g, n,n,n,n, n,n,n,n, // 4x
    n,n,n,n, t,n,n,n, n,n,n,n, n,n,n,n, // 5x
    n,a,n,c, n,n,n,g, n,n,n,n, n,n,n,n, // 6x
    n,n,n,n, t,n,n,n, n,n,n,n, n,n,n,n, // 7x

    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n
};

const Nt4 Nt4::from_nt2[4] = {
    a,c,g,t
};

const char Nt4::to_ch[16] {
    '=','A','C','?','G','?','?','?','T','?','?','?','?','?','?','N'
};

const Nt4 Nt4::rev_compl_[16] {
    n,t,g,n,c,n,n,n,a,n,n,n,n,n,n,n
};

const size_t Nt4::to_index[16] {
    0,1,2,0,3,0,0,0,4,0,0,0,0,0,0,5
};

const Nt2 Nt2::a (0);
const Nt2 Nt2::c (1);
const Nt2 Nt2::g (2);
const Nt2 Nt2::t (3);
const array<Nt2,4> Nt2::all = {{a, c, g, t}};

const Nt2 Nt2::from_ch[256] = {
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 0x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 1x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 2x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 3x
    0,a,0,c, 0,0,0,g, 0,0,0,0, 0,0,0,0, // 4x
    0,0,0,0, t,0,0,0, 0,0,0,0, 0,0,0,0, // 5x
    0,a,0,c, 0,0,0,g, 0,0,0,0, 0,0,0,0, // 6x
    0,0,0,0, t,0,0,0, 0,0,0,0, 0,0,0,0, // 7x

    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

const Nt2 Nt2::from_nt4[16] = {
    0,a,c,0,g,0,0,0,t,0,0,0,0,0,0,0
};

const char Nt2::to_ch[4] {
    'A','C','G','T'
};

const Nt2 Nt2::rev_compl_[4] {
    t,g,c,a
};
