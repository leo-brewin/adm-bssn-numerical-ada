package body ADMBase.Time_Derivs is

   dxdx, dydy, dzdz : Real;
   two_dx, two_dy, two_dz : Real;
   four_dxdy, four_dxdz, four_dydz : Real;

   procedure set_time_derivatives (point : GridPoint) is

      -- these are local variables used by the basic codes such as dot-gab.ad, dot-kab.ad etc.
      -- renames are good as they don't involve a memory copy

      ----------------------------------------------------------------------
      -- basic local data

      N       : Real;                                     -- lapse
      gab     : MetricPointArray;                         -- 3-metric g_{ab}
      iab     : MetricPointArray;                         -- 3-metric g^{ab}
      Kab     : ExtcurvPointArray;                        -- extrinsic curvature K_{ab}

      ----------------------------------------------------------------------
      -- time derivatives of local data

      dot_N   : Real;
      dot_gab : MetricPointArray;
      dot_Kab : ExtcurvPointArray;

      ----------------------------------------------------------------------
      -- partial derivatives of local data

      d1N     : Array (1..3) of Real;                     -- 1st & 2nd derivs of N
      d2N     : Array (symmetric) of Real;

      d1iab   : Array (1..3) of MetricPointArray;         -- 1st derivs of iab

      d1gab   : Array (1..3) of MetricPointArray;         -- 1st & 2nd derivs of N
      d2gab   : Array (symmetric) of MetricPointArray;

      d1Kab   : Array (1..3) of ExtcurvPointArray;        -- 1st derivs of Kab

      ----------------------------------------------------------------------
      -- miscellaneous data

      trK     : Real;

      R       : Real;                                     -- Ricci scalar

      Rab     : Array (symmetric) of Real;                -- Ricci tensor
      Hess    : Array (symmetric) of Real;                -- Hessian, N_{;ab}

      ----------------------------------------------------------------------
      -- constraints

      Ham    : Real;                                      -- Hamiltonian constraint
      Mom    : MomConstraintPointArray;                   -- Momentum constraints

      ----------------------------------------------------------------------
      -- mappings between the Cadabra and Ada objects
      -- LCB: the following renames could be eliminated by using a set of sed-scripts
      --      to substitute the symbols in the foo.ad files

      -- \partial_{a}{N} --> dN_{a} --> Na --> d1N (a)

      Nx : Real renames d1N (1);
      Ny : Real renames d1N (2);
      Nz : Real renames d1N (3);

      -- \partial_{a b}{N} --> dN_{a b} --> Nab --> d2N (ab)

      Nxx : Real renames d2N (xx);
      Nxy : Real renames d2N (xy);
      Nxz : Real renames d2N (xz);
      Nyy : Real renames d2N (yy);
      Nyz : Real renames d2N (yz);
      Nzz : Real renames d2N (zz);

      -- Rab_{p q} --> Rpq --> Rab (pq)

      Rxx : Real renames Rab (xx);
      Rxy : Real renames Rab (xy);
      Rxz : Real renames Rab (xz);
      Ryy : Real renames Rab (yy);
      Ryz : Real renames Rab (yz);
      Rzz : Real renames Rab (zz);

      -- Kab_{p q} --> Kpq --> Kab (pq)

      Kxx : Real Renames Kab (xx);
      Kxy : Real Renames Kab (xy);
      Kxz : Real Renames Kab (xz);
      Kyy : Real Renames Kab (yy);
      Kyz : Real Renames Kab (yz);
      Kzz : Real Renames Kab (zz);

      -- g_{p q} --> gpq --> gab (pq)

      gxx : Real Renames gab (xx);
      gxy : Real Renames gab (xy);
      gxz : Real Renames gab (xz);
      gyy : Real Renames gab (yy);
      gyz : Real Renames gab (yz);
      gzz : Real Renames gab (zz);

      -- g^{p q} --> ipq --> iab (pq)

      ixx : Real renames iab (xx);
      ixy : Real renames iab (xy);
      ixz : Real renames iab (xz);
      iyy : Real renames iab (yy);
      iyz : Real renames iab (yz);
      izz : Real renames iab (zz);

      -- \partial_{p}{Kab_{q r}} --> dKab_{q r p} --> Kqrp --> d1Kab (p)(qr)

      Kxxx : Real renames d1Kab (1)(xx);
      Kxxy : Real renames d1Kab (2)(xx);
      Kxxz : Real renames d1Kab (3)(xx);

      Kxyx : Real renames d1Kab (1)(xy);
      Kxyy : Real renames d1Kab (2)(xy);
      Kxyz : Real renames d1Kab (3)(xy);

      Kxzx : Real renames d1Kab (1)(xz);
      Kxzy : Real renames d1Kab (2)(xz);
      Kxzz : Real renames d1Kab (3)(xz);

      Kyyx : Real renames d1Kab (1)(yy);
      Kyyy : Real renames d1Kab (2)(yy);
      Kyyz : Real renames d1Kab (3)(yy);

      Kyzx : Real renames d1Kab (1)(yz);
      Kyzy : Real renames d1Kab (2)(yz);
      Kyzz : Real renames d1Kab (3)(yz);

      Kzzx : Real renames d1Kab (1)(zz);
      Kzzy : Real renames d1Kab (2)(zz);
      Kzzz : Real renames d1Kab (3)(zz);

      -- \partial_{p}{gab_{q r}} --> d1gab_{q r p} --> gqrp --> d1gab (p)(qr)

      gxxx : Real renames d1gab (1)(xx);
      gxxy : Real renames d1gab (2)(xx);
      gxxz : Real renames d1gab (3)(xx);

      gxyx : Real renames d1gab (1)(xy);
      gxyy : Real renames d1gab (2)(xy);
      gxyz : Real renames d1gab (3)(xy);

      gxzx : Real renames d1gab (1)(xz);
      gxzy : Real renames d1gab (2)(xz);
      gxzz : Real renames d1gab (3)(xz);

      gyyx : Real renames d1gab (1)(yy);
      gyyy : Real renames d1gab (2)(yy);
      gyyz : Real renames d1gab (3)(yy);

      gyzx : Real renames d1gab (1)(yz);
      gyzy : Real renames d1gab (2)(yz);
      gyzz : Real renames d1gab (3)(yz);

      gzzx : Real renames d1gab (1)(zz);
      gzzy : Real renames d1gab (2)(zz);
      gzzz : Real renames d1gab (3)(zz);

      -- \partial_{p}{gab^{q r}} --> d1gab^{q r}_{p} --> iqrp --> d1iab (p)(qr)

      ixxx : Real renames d1iab (1)(xx);
      ixxy : Real renames d1iab (2)(xx);
      ixxz : Real renames d1iab (3)(xx);

      ixyx : Real renames d1iab (1)(xy);
      ixyy : Real renames d1iab (2)(xy);
      ixyz : Real renames d1iab (3)(xy);

      ixzx : Real renames d1iab (1)(xz);
      ixzy : Real renames d1iab (2)(xz);
      ixzz : Real renames d1iab (3)(xz);

      iyyx : Real renames d1iab (1)(yy);
      iyyy : Real renames d1iab (2)(yy);
      iyyz : Real renames d1iab (3)(yy);

      iyzx : Real renames d1iab (1)(yz);
      iyzy : Real renames d1iab (2)(yz);
      iyzz : Real renames d1iab (3)(yz);

      izzx : Real renames d1iab (1)(zz);
      izzy : Real renames d1iab (2)(zz);
      izzz : Real renames d1iab (3)(zz);

      -- \partial_{p q}{gab_{r s}} --> d1gab_{r s p q} --> grspq --> d2gab (pq)(rs)

      gxxxx : Real renames d2gab (xx)(xx);
      gxyxx : Real renames d2gab (xx)(xy);
      gxzxx : Real renames d2gab (xx)(xz);
      gyyxx : Real renames d2gab (xx)(yy);
      gyzxx : Real renames d2gab (xx)(yz);
      gzzxx : Real renames d2gab (xx)(zz);

      gxxxy : Real renames d2gab (xy)(xx);
      gxyxy : Real renames d2gab (xy)(xy);
      gxzxy : Real renames d2gab (xy)(xz);
      gyyxy : Real renames d2gab (xy)(yy);
      gyzxy : Real renames d2gab (xy)(yz);
      gzzxy : Real renames d2gab (xy)(zz);

      gxxxz : Real renames d2gab (xz)(xx);
      gxyxz : Real renames d2gab (xz)(xy);
      gxzxz : Real renames d2gab (xz)(xz);
      gyyxz : Real renames d2gab (xz)(yy);
      gyzxz : Real renames d2gab (xz)(yz);
      gzzxz : Real renames d2gab (xz)(zz);

      gxxyy : Real renames d2gab (yy)(xx);
      gxyyy : Real renames d2gab (yy)(xy);
      gxzyy : Real renames d2gab (yy)(xz);
      gyyyy : Real renames d2gab (yy)(yy);
      gyzyy : Real renames d2gab (yy)(yz);
      gzzyy : Real renames d2gab (yy)(zz);

      gxxyz : Real renames d2gab (yz)(xx);
      gxyyz : Real renames d2gab (yz)(xy);
      gxzyz : Real renames d2gab (yz)(xz);
      gyyyz : Real renames d2gab (yz)(yy);
      gyzyz : Real renames d2gab (yz)(yz);
      gzzyz : Real renames d2gab (yz)(zz);

      gxxzz : Real renames d2gab (zz)(xx);
      gxyzz : Real renames d2gab (zz)(xy);
      gxzzz : Real renames d2gab (zz)(xz);
      gyyzz : Real renames d2gab (zz)(yy);
      gyzzz : Real renames d2gab (zz)(yz);
      gzzzz : Real renames d2gab (zz)(zz);

      -- Hess_{a b} --> hessab --> hess (ab)

      Hessxx : Real renames hess (xx);
      Hessxy : Real renames hess (xy);
      Hessxz : Real renames hess (xz);
      Hessyy : Real renames hess (yy);
      Hessyz : Real renames hess (yz);
      Hesszz : Real renames hess (zz);

      Procedure set_3d_hessian is
         x0,  x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,  x10, x11, x12, x13, x14, x15, x16,
         x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32 : Real;
      begin
         x0  := Nx*ixx;
         x1  := (1.0/2.0)*gxxx;
         x2  := Ny*ixy;
         x3  := Nz*ixz;
         x4  := (1.0/2.0)*gxxy;
         x5  := gxyx - x4;
         x6  := Nx*ixy;
         x7  := (1.0/2.0)*gxxz;
         x8  := gxzx - x7;
         x9  := Nx*ixz;
         x10 := Ny*iyy;
         x11 := Ny*iyz;
         x12 := Nz*iyz;
         x13 := Nz*izz;
         x14 := (1.0/2.0)*gyyx;
         x15 := (1.0/2.0)*gxyz;
         x16 := (1.0/2.0)*gxzy;
         x17 := (1.0/2.0)*gyzx;
         x18 := -x15 + x16 + x17;
         x19 := Nxy - x0*x4 - x10*x14 - x11*x18 - x12*x14 - x13*x18 - x14*x6 - x18*x9 - x2*x4 - x3*x4;
         x20 := (1.0/2.0)*gzzx;
         x21 := x15 - x16 + x17;
         x22 := Nxz - x0*x7 - x10*x21 - x11*x20 - x12*x21 - x13*x20 - x2*x7 - x20*x9 - x21*x6 - x3*x7;
         x23 := (1.0/2.0)*gyyy;
         x24 := gxyy - x14;
         x25 := (1.0/2.0)*gyyz;
         x26 := gyzy - x25;
         x27 := (1.0/2.0)*gzzy;
         x28 := x15 + x16 - x17;
         x29 := Nyz - x0*x28 - x10*x25 - x11*x27 - x12*x25 - x13*x27 - x2*x28 - x25*x6 - x27*x9 - x28*x3;
         x30 := (1.0/2.0)*gzzz;
         x31 := gxzz - x20;
         x32 := gyzz - x27;
         Hess (xx) := Nxx - x0*x1 - x1*x2 - x1*x3 - x10*x5 - x11*x8 - x12*x5 - x13*x8 - x5*x6 - x8*x9;
         Hess (xy) := x19;
         Hess (xz) := x22;
         Hess (xy) := x19;
         Hess (yy) := Nyy - x0*x24 - x10*x23 - x11*x26 - x12*x23 - x13*x26 - x2*x24 - x23*x6 - x24*x3 - x26*x9;
         Hess (yz) := x29;
         Hess (xz) := x22;
         Hess (yz) := x29;
         Hess (zz) := Nzz - x0*x31 - x10*x32 - x11*x30 - x12*x32 - x13*x30 - x2*x31 - x3*x31 - x30*x9 - x32*x6;
      end set_3d_hessian;

      Procedure set_3d_ricci is
         x0,   x1,   x2,   x3,   x4,   x5,   x6,   x7,   x8,   x9,   x10,  x11,  x12,  x13,
         x14,  x15,  x16,  x17,  x18,  x19,  x20,  x21,  x22,  x23,  x24,  x25,  x26,  x27,
         x28,  x29,  x30,  x31,  x32,  x33,  x34,  x35,  x36,  x37,  x38,  x39,  x40,  x41,
         x42,  x43,  x44,  x45,  x46,  x47,  x48,  x49,  x50,  x51,  x52,  x53,  x54,  x55,
         x56,  x57,  x58,  x59,  x60,  x61,  x62,  x63,  x64,  x65,  x66,  x67,  x68,  x69,
         x70,  x71,  x72,  x73,  x74,  x75,  x76,  x77,  x78,  x79,  x80,  x81,  x82,  x83,
         x84,  x85,  x86,  x87,  x88,  x89,  x90,  x91,  x92,  x93,  x94,  x95,  x96,  x97,
         x98,  x99,  x100, x101, x102, x103, x104, x105, x106, x107, x108, x109, x110, x111,
         x112, x113, x114, x115, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125,
         x126, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137, x138, x139,
         x140, x141, x142, x143, x144, x145, x146, x147, x148, x149, x150, x151, x152, x153,
         x154, x155, x156, x157, x158, x159, x160, x161, x162, x163, x164, x165, x166, x167,
         x168, x169, x170, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181,
         x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192, x193, x194, x195,
         x196, x197, x198, x199, x200, x201, x202, x203, x204, x205, x206, x207, x208, x209,
         x210, x211, x212, x213, x214, x215, x216, x217, x218, x219, x220, x221, x222, x223,
         x224, x225, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236, x237,
         x238, x239, x240, x241, x242, x243, x244, x245, x246, x247, x248, x249, x250 : Real;
         x251, x252, x253, x254, x255, x256, x257, x258, x259, x260, x261, x262, x263, x264,
         x265, x266, x267, x268, x269, x270, x271, x272, x273, x274, x275, x276, x277, x278,
         x279, x280, x281, x282, x283, x284, x285, x286, x287 : Real;
      begin
         x0   := (1.0/4.0)*gxxx;
         x1   := (1.0/2.0)*iyy;
         x2   := (1.0/2.0)*izz;
         x3   := (1.0/2.0)*gxyx;
         x4   := (1.0/2.0)*gxzx;
         x5   := (1.0/4.0)*gyyx;
         x6   := (1.0/2.0)*gyzx;
         x7   := (1.0/4.0)*gzzx;
         x8   := gxxx*ixx;
         x9   := gxyx*ixy;
         x10  := gxzx*ixz;
         x11  := gxyy*x9;
         x12  := gxyx*ixz;
         x13  := gxyz*iyy;
         x14  := gxyx*iyz;
         x15  := iyz*x9;
         x16  := gxzz*iyz;
         x17  := gyzy*iyy;
         x18  := gxzx*ixy;
         x19  := gxyy*iyz;
         x20  := iyz*x10;
         x21  := gxzy*izz;
         x22  := gxzz*izz;
         x23  := iyz*izz;
         x24  := gyzz*x23;
         x25  := ixyx + iyyy + iyzz;
         x26  := ixzx + iyzy + izzz;
         x27  := (iyz ** 2);
         x28  := gxyx*gyzz;
         x29  := gxyz*gxzy;
         x30  := gxzx*gyzy;
         x31  := (1.0/2.0)*ixy;
         x32  := ixz*x31;
         x33  := gxxx*x32;
         x34  := ixx*iyy;
         x35  := gyyx*x0;
         x36  := gyyy*iyy;
         x37  := ixy*x0;
         x38  := ixz*x0;
         x39  := gyyz*x38;
         x40  := iyz*x6;
         x41  := iyz*x31;
         x42  := gyzy*x41;
         x43  := (1.0/2.0)*ixz;
         x44  := iyz*x43;
         x45  := gyzz*x44;
         x46  := ixx*izz;
         x47  := gzzx*x46;
         x48  := gzzy*izz;
         x49  := gzzz*ixz;
         x50  := izz*x49;
         x51  := gxxz*ixx;
         x52  := (1.0/2.0)*iyz;
         x53  := gxxy*x52;
         x54  := gxxy*ixz;
         x55  := gxxz*ixy;
         x56  := ixx*x1;
         x57  := gxxy*x56;
         x58  := gxyy*ixy;
         x59  := x1*x58;
         x60  := gxxy*iyz;
         x61  := x31*x60;
         x62  := ixx*x4;
         x63  := x4*x54;
         x64  := gxzy*ixy;
         x65  := gxzy*x54;
         x66  := x16*x43;
         x67  := gxxy*ixy;
         x68  := iyy*x5;
         x69  := iyy*iyz;
         x70  := (1.0/4.0)*gxxy;
         x71  := gyyz*x70;
         x72  := ixy*x6;
         x73  := gyzy*x1;
         x74  := izz*x7;
         x75  := iyy*izz;
         x76  := gzzy*x70;
         x77  := gzzz*x23;
         x78  := iyz*x3;
         x79  := ixz*x3;
         x80  := gxxz*x31;
         x81  := x2*x55;
         x82  := ixz*iyz;
         x83  := gxyz*x82;
         x84  := ixx*x2;
         x85  := gxzx*x84;
         x86  := gxxz*x44;
         x87  := gxzz*ixz;
         x88  := gxxz*x2;
         x89  := gxxz*ixz;
         x90  := (1.0/4.0)*gyyy;
         x91  := x69*x90;
         x92  := (1.0/4.0)*gxxz;
         x93  := gyyz*x92;
         x94  := x6*x82;
         x95  := gyzz*iyz;
         x96  := (1.0/4.0)*gzzy;
         x97  := x23*x96;
         x98  := gyyx*x1;
         x99  := gyyz*x1;
         x100 := gzzx*x2;
         x101 := izz*x1;
         x102 := gzzy*x101;
         x103 := gzzz*x2;
         x104 := iyz*x1;
         x105 := gyyz*x101;
         x106 := gzzy*iyz;
         x107 := x106*x2;
         x108 := ixxx + ixyy + ixzz;
         x109 := (1.0/2.0)*x108;
         x110 := (1.0/2.0)*x25;
         x111 := (1.0/2.0)*x26;
         x112 := (ixy ** 2);
         x113 := (1.0/2.0)*x112;
         x114 := gxyy*x113;
         x115 := (ixz ** 2);
         x116 := (1.0/2.0)*x115;
         x117 := gxzz*x116;
         x118 := gxxy*x113;
         x119 := (iyy ** 2);
         x120 := gyyy*x119;
         x121 := (1.0/2.0)*x27;
         x122 := gyzz*x121;
         x123 := (gxxy ** 2);
         x124 := (1.0/4.0)*ixx;
         x125 := x123*x124;
         x126 := gxxz*x116;
         x127 := gyzy*x121;
         x128 := (izz ** 2);
         x129 := gzzz*x128;
         x130 := (gxxz ** 2);
         x131 := x124*x130;
         x132 := (gxyx ** 2);
         x133 := (gxyz ** 2);
         x134 := (gxzx ** 2);
         x135 := (gxzy ** 2);
         x136 := (ixx ** 2);
         x137 := gxzx*x116;
         x138 := gyyz*x121;
         x139 := x27*x6;
         x140 := (gyzx ** 2);
         x141 := (1.0/4.0)*gzzz;
         x142 := x128*x141;
         x143 := x0*x136;
         x144 := (gyyx ** 2);
         x145 := iyy*x144;
         x146 := ixx*x3;
         x147 := gxyy*x51;
         x148 := gxxz*x32;
         x149 := gxyx*x56;
         x150 := gyzx*ixx;
         x151 := gyzy*ixy;
         x152 := gxyy*ixz;
         x153 := ixy*x4;
         x154 := x1*x152;
         x155 := gyyx*iyz;
         x156 := x4*x82;
         x157 := ixy*x2;
         x158 := gxzy*x2;
         x159 := gxzy*gyyx;
         x160 := gxzy*gyzx;
         x161 := gxzy*x82;
         x162 := gyzy*x101;
         x163 := gyzx*x2;
         x164 := ixy*x155;
         x165 := gxyz*ixx;
         x166 := gyyz*ixy;
         x167 := gxyz*x2;
         x168 := x18*x2;
         x169 := ixz*iyy;
         x170 := x5*x51;
         x171 := gyzx*x46;
         x172 := iyz*x13;
         x173 := ixz*x74;
         x174 := (1.0/4.0)*gyyz;
         x175 := iyy*x21;
         x176 := gxzy*ixx;
         x177 := ixz*x51;
         x178 := iyz*x70;
         x179 := iyz*x176;
         x180 := x36*x70;
         x181 := ixy*izz;
         x182 := x165*x92;
         x183 := ixx*x21;
         x184 := x13*x5;
         x185 := izz*x13;
         x186 := ixz*x68;
         x187 := x21*x7;
         x188 := iyz*x21;
         x189 := gyyz*x69;
         x190 := gyzx*iyy;
         x191 := ixz*x190;
         x192 := gzzx*x5;
         x193 := gzzy*x75;
         x194 := gyzx*x174;
         x195 := gxxy*x109 + gxxy*x117 + gxxy*x143 + gxxy*x46*x7 + gxxyy*x31 + gxxyz*x43 + gxyx*x114 - gxyxy*ixy - gxyxz*x43 - gxyy*x149 + gxyy*x57 - gxyyz*x52 - gxyz*x111 - gxyz*x126 - gxyz*x138 - gxyz*x142 - gxyz*x154 - gxyz*x162 - gxyz*x173 - gxyz*x85 - gxyz*x97
               - gxyzz*x2 - gxzxy*x43 + gxzy*x111 + gxzy*x137 + gxzy*x138 + gxzy*x142 + gxzy*x154 + gxzy*x162 + gxzy*x186 + gxzy*x91 + gxzyy*x52 + gxzyz*x2 + gyyx*x110 + gyyx*x118 + gyyx*x122 + gyyx*x148 + gyyx*x149 + gyyx*x66 + gyyxx*x31 + gyyxz*x52 + gyyz*x156
               - gyyz*x168 + gyyz*x61 + gyyz*x81 - gyyz*x86 + gyzx*x111 + gyzx*x126 + gyzx*x142 + gyzx*x173 + gyzx*x85 + gyzx*x91 + gyzx*x97 + gyzxx*x43 - gyzxy*x52 + gyzxz*x2 + gyzy*x139 - gyzy*x156 - gyzy*x81 + gyzy*x86 + gyzz*x43*x60 - gzzxy*x2 + ixx*ixy*x35
               + ixx*x63 + ixy*x125 + (1.0/4.0)*ixy*x145 + ixy*x180 + ixz*x184 + ixz*x187 - iyz*x170 + izz*x182 - x12*x73 + x12*x99 + x120*x5 - x133*x44 + x135*x44 + x140*x44 + x146*x67 + x147*x52 + x150*x178 + x150*x38 + x150*x78 + x151*x78 + x152*x153 - x152*x80
               + x155*x62 + x155*x73 + x157*x160 + x157*x30 + x158*x87 + x158*x95 + x159*x41 + x161*x6 + x163*x87 + x163*x95 + x164*x6 + x165*x178 - x165*x38 - x165*x78 - x166*x78 - x167*x87 - x167*x95 - x169*x71 - x171*x92 - x172*x90 - x174*x175 + x174*x185 + x176*x38
               + x177*x70 + x179*x70 + x181*x192 + x181*x76 + x183*x92 + x188*x96 + x189*x5 - x19*x62 + x19*x72 + x191*x5 + x193*x5 + x194*x75 + x31*x65 + x5*x77 + x50*x70 + x54*x72 + x54*x73 + x58*x98 + x64*x79;
         x196 := gxyx*x113;
         x197 := gzzy*x121;
         x198 := x119*x90;
         x199 := (gzzx ** 2);
         x200 := izz*x199;
         x201 := gxzz*ixx;
         x202 := gzzx*x31;
         x203 := x1*x54;
         x204 := gxxz*x84;
         x205 := ixz*x6;
         x206 := ixy*x3;
         x207 := ixz*x1;
         x208 := gzzx*ixx;
         x209 := gxyz*ixy;
         x210 := ixz*x4;
         x211 := gxzz*x2;
         x212 := gxyz*gyzx;
         x213 := gyzy*x104;
         x214 := gyzz*x101;
         x215 := gxyz*x44;
         x216 := gyzx*x62;
         x217 := gzzx*x82;
         x218 := ixx*x0;
         x219 := ixx*x70;
         x220 := ixx*x7;
         x221 := izz*x96;
         x222 := gxzy*x174;
         x223 := gyzx*ixy;
         x224 := iyy*x90;
         x225 := x150*x92;
         x226 := x141*x23;
         x227 := x36*x7;
         x228 := gyyz*x7;
         x229 := izz*x223;
         x230 := gxxyz*x31 + gxxz*x109 + gxxz*x114 + gxxz*x143 + gxxz*x34*x5 + gxxz*x42 + gxxzz*x43 - gxyxz*x31 - gxyy*x1*x64 + gxyyz*x1 + gxyz*x110 + gxyz*x148 + gxyz*x196 + gxyz*x197 + gxyz*x198 + gxyz*x213 + gxyz*x214 + gxyz*x226 + gxyz*x59 + gxyzz*x52 + gxzx*x117
               - gxzxy*x31 - gxzxz*ixz - gxzy*x110 - gxzy*x118 - gxzy*x149 - gxzy*x197 - gxzy*x198 - gxzy*x213 - gxzy*x214 - gxzyy*x1 - gxzyz*x52 + gxzz*x204 - gxzz*x31*x54 - gxzz*x85 - gyyxz*x1 + gyzx*x110 + gyzx*x118 + gyzx*x149 + gyzx*x198 + gyzx*x213 + gyzx*x226
               + gyzx*x59 + gyzx*x75*x96 + gyzxx*x31 + gyzxy*x1 - gyzxz*x52 + gyzz*x139 + gyzz*x156 - gyzz*x168 - gyzz*x203 + gyzz*x61 + gyzz*x81 + gzzx*x111 + gzzx*x126 + gzzx*x127 + gzzx*x186 + gzzx*x215 + gzzx*x85 + gzzxx*x43 + gzzxy*x52 - gzzy*x1*x12 - gzzy*x156
               + gzzy*x168 + gzzy*x203 + gzzy*x23*x7 - gzzy*x61 + gzzy*x86 + ixy*x184 + ixy*x187 + ixz*x131 + (1.0/4.0)*ixz*x200 + iyy*x176*x70 + iyz*x182 + iyz*x216 + iyz*x225 + iyz*x227 + x100*x87 + x100*x95 + x106*x206 + x129*x7 + x13*x219 + x133*x41 - x135*x41
               + x140*x41 - x141*x188 - x146*x16 + x146*x55 + x16*x205 + x165*x37 + x169*x93 + x172*x174 + x175*x96 + x177*x4 - x179*x4 + x179*x92 - x185*x96 + x19*x202 - x190*x219 + x194*x69 + x201*x53 + x202*x54 + x205*x55 + x206*x87 - x206*x95 + x207*x212 + x207*x28
               + x208*x38 + x208*x78 + x209*x210 + x209*x211 + x209*x40 + x209*x74 - x211*x64 + x217*x6 + x218*x223 - x218*x64 + x219*x55 - x220*x60 - x221*x55 - x222*x69 + x223*x68 + x224*x55 + x228*x75 + x229*x7 + x50*x92 - x64*x68;
         x231 := (1.0/2.0)*gxyy;
         x232 := (1.0/2.0)*gxzy;
         x233 := (1.0/2.0)*gyzy;
         x234 := (1.0/2.0)*ixx;
         x235 := gyzy*ixz;
         x236 := ixx*x10;
         x237 := gyzz*ixz;
         x238 := iyz*x17;
         x239 := gxyy*gxzz;
         x240 := ixx*x37;
         x241 := gxxy*x31;
         x242 := ixx*x5;
         x243 := ixx*iyz;
         x244 := iyz*x51;
         x245 := ixy*x146;
         x246 := gxyy*x56;
         x247 := x19*x31;
         x248 := gzzx*x84;
         x249 := gzzy*x2;
         x250 := gyyx*x32;
         x251 := gxyz*x41;
         x252 := gyyz*x2;
         x253 := ixz*x62;
         x254 := ixz*x153;
         x255 := gxzy*x207;
         x256 := gxzy*gyyz;
         x257 := ixz*x223;
         x258 := gyyy*ixy;
         x259 := gxxx*x136;
         x260 := (gxyy ** 2);
         x261 := (gyyz ** 2);
         x262 := (1.0/4.0)*x261;
         x263 := (gyzy ** 2);
         x264 := gzzx*x116;
         x265 := gyyx*x113;
         x266 := (gzzy ** 2);
         x267 := (1.0/4.0)*x266;
         x268 := gzzy*ixy;
         x269 := ixx*x52;
         x270 := gxzz*x84;
         x271 := gyzz*x2;
         x272 := x1*x235;
         x273 := gxzy*x44;
         x274 := gxzz*x157;
         x275 := x208*x52;
         x276 := gyyx*x41;
         x277 := gzzx*x44;
         x278 := ixz*x141;
         x279 := izz*x278;
         x280 := -gxxyz*x234 + gxxz*x46*x96 + gxyxz*x234 - gxyy*x275 - gxyyz*x31 + gxyz*x109 + gxyz*x114 + gxyz*x143 + gxyz*x264 + gxyz*x270 + gxyz*x279 + gxyz*x42 + gxyzz*x43 + gxzxy*x234 + gxzy*x109 + gxzy*x117 + gxzy*x143 + gxzy*x246 + gxzy*x265 + gxzy*x272 + gxzy*x45
               + gxzyy*x31 - gxzyz*x43 + gxzz*x250 - gyyx*x201*x52 + gyyx*x275 + gyyxz*x31 + gyyz*x110 + gyyz*x196 + gyyz*x198 + gyyz*x240 + gyyz*x251 + gyyz*x254 + gyyz*x274 + gyyz*x277 + gyyzz*x52 - gyzx*x109 - gyzx*x143 - gyzx*x246 - gyzx*x264 - gyzx*x265 - gyzx*x270
               - gyzx*x272 - gyzx*x279 - gyzxx*x234 - gyzxy*x31 - gyzxz*x43 + gyzy*iyz*x99 + gyzy*x102 + gyzy*x122 - gyzy*x277 + gyzy*x66 - gyzyz*iyz + gyzz*x105 - gyzz*x154 - gyzz*x162 + gyzz*x247 - gyzz*x276 - gzzx*x250 + gzzxy*x43 + gzzy*ixx*x38 + gzzy*x111
               + gzzy*x137 + gzzy*x138 + gzzy*x154 + gzzy*x173 - gzzy*x186 + gzzy*x215 + gzzy*x273 + gzzy*x276 + gzzy*x91 + gzzyy*x52 + ixx*x184 + ixy*x13*x90 + ixy*x165*x70 + ixy*x21*x96 + ixz*x13*x174 + ixz*x176*x92 + ixz*x182 - ixz*x216 - ixz*x225 + x100*x151
               + x106*x219 + x129*x96 + x133*x32 + x135*x32 - x140*x32 - x146*x223 + x146*x64 - x151*x211 + x152*x202 + x165*x206 + x165*x210 - x165*x74 + x166*x68 - x166*x74 + x169*x222 + x171*x7 + x174*x191 + x174*x77 + x176*x210 - x176*x68 + x183*x7 + x190*x242
               + x209*x221 + x209*x271 + x21*x278 - x219*x223 + x219*x64 - x223*x224 - x223*x271 + x224*x64 + x229*x96 + x23*x267 + x237*x98 + x239*x269 - x239*x32 + x243*x93 + x249*x87 + x249*x95 + x256*x41 + x262*x69 + x268*x79 + x269*x29 + x29*x32 + x34*x71 + x58*x99;
         x281 := (1.0/2.0)*gxyz;
         x282 := (1.0/2.0)*gxzz;
         x283 := (1.0/2.0)*gyzz;
         x284 := gyzz*ixy;
         x285 := gzzx*x32;
         x286 := (gxzz ** 2);
         x287 := (gyzz ** 2);
         Rab (xx) := (1.0/4.0)*(gxxx ** 2)*x136 + gxxx*x109 + gxxx*x114 + gxxx*x117 + gxxx*x42 + gxxx*x45 - gxxy*x110 - gxxy*x122 - gxxy*x59 - gxxy*x66 - gxxyy*x1 - gxxyz*iyz - gxxz*x111 - gxxz*x127 - gxxz*x83 - gxxz*x85 - gxxz*x91 - gxxz*x94 - gxxz*x97 - gxxzz*x2
                     + gxyx*x102 + gxyx*x118 + gxyx*x25 - gxyx*x57 + gxyxy*iyy + gxyxz*iyz + gxyz*x20 + gxyz*x33 + gxyz*x61 + gxyz*x81 + gxzx*gyyy*x104 + gxzx*ixx*x14 + gxzx*x105 + gxzx*x107 + gxzx*x126 + gxzx*x24 + gxzx*x26 + gxzxy*iyz + gxzxz*izz + gxzy*x15 + gxzy*x33
                     + gxzy*x86 - gyyxx*x1 + gyzx*x15 + gyzx*x20 - gyzxx*iyz - gzzxx*x2 - ixxx*x0 + ixy*x63 - ixyx*x3 - ixzx*x4 + iyy*x11 + iyy*x125 + iyy*x39 - iyyx*x5 - iyzx*x6 + izz*x131 - izzx*x7 + x0*x47 + x0*x50 + x1*x65 + x10*x100 + x10*x22 + x10*x8 + x10*x9 + x10*x98
                     + x100*x9 + x101*x133 + x101*x135 + x103*x14 - x113*x123 + x113*x132 - x116*x130 + x116*x134 + x12*x13 + x12*x16 + x120*x3 - x120*x70 - x121*x133 - x121*x135 + x129*x4 - x129*x92 - x13*x54 + x132*x56 + x134*x84 + x14*x17 + x14*x99 + x18*x19 + x18*x21
                     - x19*x80 - x21*x55 + x27*x28 + x27*x29 + x27*x30 - x29*x75 + x34*x35 + x36*x37 + x37*x48 + x40*x8 + x51*x53 - x51*x78 - x54*x55 + x55*x79 - x60*x62 - x60*x64 - x60*x72 - x60*x73 - x67*x68 - x67*x74 - x68*x89 - x69*x71 - x70*x77 - x74*x89 - x75*x76
                     - x75*x93 + x8*x9 - x87*x88 - x88*x95 + x9*x98;
         Rab (xy) := -ixxy*x0 - ixyy*x3 - ixzy*x4 - iyyy*x5 - iyzy*x6 - izzy*x7 + x195;
         Rab (xz) := -ixxz*x0 - ixyz*x3 - ixzz*x4 - iyyz*x5 - iyzz*x6 - izzz*x7 + x230;
         Rab (xy) := -ixxx*x70 - ixyx*x231 - ixzx*x232 - iyyx*x90 - iyzx*x233 - izzx*x96 + x195;
         Rab (yy) := -gxxyy*x234 + gxyxy*ixx + gxyy*ixx*x241 + gxyy*x108 + gxyy*x236 + gxyy*x248 + gxyyz*ixz + gxyz*x250 - gxzyy*ixz - gyyx*x109 + gyyx*x114 - gyyx*x117 - gyyx*x143 - gyyx*x245 - gyyx*x246 - gyyx*x253 - gyyx*x257 + gyyx*x42 - gyyx*x45 - gyyxx*x234 - gyyxz*ixz
                     + (1.0/4.0)*(gyyy ** 2)*x119 + gyyy*x110 + gyyy*x122 + gyyy*x196 + gyyy*x238 + gyyy*x240 + gyyy*x251 + gyyy*x254 + gyyy*x255 + gyyy*x66 + gyyz*ixz*x98 - gyyz*x111 + gyyz*x127 - gyyz*x137 - gyyz*x164 - gyyz*x173 - gyyz*x229 + gyyz*x247 - gyyz*x83
                     + gyyz*x94 - gyyz*x97 - gyyzz*x2 + gyzxy*ixz + gyzy*ixx*x53 - gyzy*x105 + gyzy*x107 + gyzy*x161 + gyzy*x204 + gyzy*x229 + gyzy*x24 + gyzy*x26 + gyzy*x43*x8 + gyzy*x83 + gyzyz*izz - gzzyy*x2 + ixx*x11 + ixx*x155*x6 + ixx*x180 - ixx*x39 - ixxy*x70
                     - ixy*x48*x5 - ixyy*x231 - ixz*x170 - ixzy*x232 - iyyy*x90 - iyzy*x233 - izzy*x96 + x100*x235 + x101*x263 + x103*x152 - x113*x144 + x113*x260 + x115*x212 + x115*x239 + x115*x30 - x116*x133 - x116*x140 - x121*x261 + x121*x263 + x124*x145 - x129*x174
                     + x129*x233 + x133*x84 + x140*x84 + x147*x43 + x151*x19 + x152*x17 + x152*x223 + x152*x64 - x152*x99 - x155*x165 - x159*x32 + x165*x19 - x166*x79 + x19*x237 + x193*x90 + x209*x252 - x212*x46 + x22*x235 + x231*x259 + x235*x9 - x235*x98 - x242*x67
                     - x243*x71 + x244*x90 + x249*x58 - x252*x87 - x252*x95 - x256*x44 + x258*x40 + x258*x74 + x260*x56 + x262*x75 + x36*x58 - x46*x93 - x47*x5 - x5*x50 + x77*x90;
         Rab (yz) := -ixxz*x70 - ixyz*x231 - ixzz*x232 - iyyz*x90 - iyzz*x233 - izzz*x96 + x280;
         Rab (xz) := -ixxx*x92 - ixyx*x281 - ixzx*x282 - iyyx*x174 - iyzx*x283 - izzx*x141 + x230;
         Rab (yz) := -ixxy*x92 - ixyy*x281 - ixzy*x282 - iyyy*x174 - iyzy*x283 - izzy*x141 + x280;
         Rab (zz) := gxxz*x201*x43 - gxxzz*x234 - gxyz*x285 - gxyzz*ixy + gxzxz*ixx + gxzy*x285 + gxzyz*ixy + gxzz*gyyx*x56 + gxzz*x1*x258 + gxzz*x108 + gxzz*x236 - gxzz*x248 + gxzz*x257 - gyyzz*x1 + gyzxz*ixy + gyzyz*iyy + gyzz*iyy*x58 - gyzz*x102 + gyzz*x238 + gyzz*x25
                     + gyzz*x31*x8 + gyzz*x51*x52 + gyzz*x57 - gzzx*x109 - gzzx*x114 + gzzx*x117 - gzzx*x143 - gzzx*x179 - gzzx*x245 - gzzx*x253 - gzzx*x257 - gzzx*x42 + gzzx*x45 - gzzxx*x234 - gzzxy*ixy - gzzy*x110 + gzzy*x122 - gzzy*x191 - gzzy*x196 - gzzy*x198 - gzzy*x217
                     - gzzy*x240 - gzzy*x251 - gzzy*x254 + gzzy*x255 - gzzy*x274 - gzzy*x59 + gzzy*x66 - gzzyy*x1 + (1.0/4.0)*(gzzz ** 2)*x128 + gzzz*ixx*x178 + gzzz*x111 + gzzz*x127 + gzzz*x137 + gzzz*x174*x75 + gzzz*x24 + gzzz*x247 + gzzz*x273 + gzzz*x46*x92 + gzzz*x91
                     + gzzz*x94 - ixxz*x92 - ixy*x227 - ixyz*x281 - ixzz*x282 - iyyz*x174 - iyzz*x283 - izzz*x141 + x10*x284 + x100*x268 - x100*x284 + x101*x287 + x103*x209 - x106*x64 + x106*x72 - x106*x73 + x112*x160 + x112*x239 + x112*x28 - x113*x135 - x113*x140 - x116*x199
                     + x116*x286 + x120*x283 - x121*x266 + x121*x287 + x124*x200 + x135*x56 + x140*x56 + x151*x16 + x16*x176 + x16*x237 - x160*x34 - x169*x228 - x177*x7 - x189*x96 + x190*x237 - x192*x34 + x201*x241 + x201*x9 + x206*x49 + x208*x40 + x209*x87 + x209*x95
                     + x218*x49 + x22*x284 + x22*x49 - x220*x67 - x244*x96 + x259*x282 + x267*x75 - x268*x68 + x284*x98 + x286*x84 - x34*x76 + x49*x68 + x64*x95 + x87*x99 + x95*x99;
      end set_3d_ricci;


      Procedure set_3d_dot_N is
      begin
         dot_N := 0.0;
      end set_3d_dot_N;

      Procedure set_3d_dot_gab is
         x0, x1, x2, x3 : Real;
      begin
         x0 := 2.0*N;
         x1 := -Kxy*x0;
         x2 := -Kxz*x0;
         x3 := -Kyz*x0;
         dot_gab (xx) := -Kxx*x0;
         dot_gab (xy) := x1;
         dot_gab (xz) := x2;
         dot_gab (xy) := x1;
         dot_gab (yy) := -Kyy*x0;
         dot_gab (yz) := x3;
         dot_gab (xz) := x2;
         dot_gab (yz) := x3;
         dot_gab (zz) := -Kzz*x0;
      end set_3d_dot_gab;

      Procedure set_3d_dot_Kab is
         x0,  x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,  x10, x11, x12, x13, x14, x15, x16,
         x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33,
         x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50,
         x51, x52, x53 : Real;
      begin
         x0  := Nx*ixy;
         x1  := (1.0/2.0)*gxxy;
         x2  := Nx*ixz;
         x3  := (1.0/2.0)*gxxz;
         x4  := Ny*iyy;
         x5  := Ny*iyz;
         x6  := Nz*iyz;
         x7  := Nz*izz;
         x8  := Kxx*Kxy;
         x9  := 4.0*N;
         x10 := ixy*x9;
         x11 := Kxz*x9;
         x12 := Kxx*ixz;
         x13 := Kxy*iyz;
         x14 := 2.0*N;
         x15 := ixx*x14;
         x16 := (Kxy ** 2);
         x17 := x14*x16;
         x18 := (Kxz ** 2);
         x19 := x14*x18;
         x20 := N*trK;
         x21 := Nx*ixx;
         x22 := (1.0/2.0)*x2;
         x23 := (1.0/2.0)*gyyx;
         x24 := Ny*ixy;
         x25 := (1.0/2.0)*x5;
         x26 := Nz*ixz;
         x27 := (1.0/2.0)*x7;
         x28 := Kyy*x14;
         x29 := Kxx*ixy;
         x30 := Kyz*x14;
         x31 := Kxy*ixz;
         x32 := Kxz*x14;
         x33 := iyy*x14;
         x34 := Kxy*x33;
         x35 := Kxz*iyz;
         x36 := izz*x14;
         x37 := Kxz*x36;
         x38 := Kxy*x20 - Kyy*x34 - Kyz*x37 - Nxy + Rxy*N - gxyz*x22 - gxyz*x25 - gxyz*x27 + gxzy*x22 + gxzy*x25 + gxzy*x27 + gyzx*x22 + gyzx*x25 + gyzx*x27 - ixy*x17 + x0*x23 + x1*x21 + x1*x24 + x1*x26 - x12*x30 - x13*x30 - x15*x8 + x23*x4 + x23*x6 - x28*x29 - x28*x35
              - x31*x32;
         x39 := (1.0/2.0)*x0;
         x40 := (1.0/2.0)*x4;
         x41 := (1.0/2.0)*x6;
         x42 := Kxz*x15;
         x43 := Kzz*x14;
         x44 := Kxy*ixy;
         x45 := -Kxx*x42 + Kxz*x20 - Kyz*x34 - Kzz*x37 - Nxz + Rxz*N + gxyz*x39 + gxyz*x40 + gxyz*x41 - gxzy*x39 - gxzy*x40 - gxzy*x41 + gyzx*x39 + gyzx*x40 + gyzx*x41 + gzzx*x22 + gzzx*x25 + gzzx*x27 - ixz*x19 - x12*x43 - x13*x43 + x21*x3 + x24*x3 + x26*x3 - x29*x30
              - x30*x35 - x32*x44;
         x46 := Kyz*x9;
         x47 := iyz*x46;
         x48 := (Kyz ** 2);
         x49 := x14*x48;
         x50 := (1.0/2.0)*x21;
         x51 := (1.0/2.0)*x24;
         x52 := (1.0/2.0)*x26;
         x53 := -Kxy*x42 - Kxz*ixy*x28 - Kxz*ixz*x30 - Kyy*Kyz*x33 - Kyz*Kzz*x36 + Kyz*x20 - Kzz*iyz*x28 - Nyz + Ryz*N + gxyz*x50 + gxyz*x51 + gxyz*x52 + gxzy*x50 + gxzy*x51 + gxzy*x52 + gyyz*x39 + gyyz*x40 + gyyz*x41 - gyzx*x50 - gyzx*x51 - gyzx*x52 + gzzy*x22
              + gzzy*x25 + gzzy*x27 - iyz*x49 - x30*x44 - x31*x43;
         dot_Kab (xx) := -(Kxx ** 2)*x15 + Kxx*N*trK + (1.0/2.0)*Nx*gxxx*ixx + Nx*gxyx*ixy + Nx*gxzx*ixz - Nxx + (1.0/2.0)*Ny*gxxx*ixy + Ny*gxyx*iyy + Ny*gxzx*iyz + (1.0/2.0)*Nz*gxxx*ixz + Nz*gxyx*iyz + Nz*gxzx*izz + Rxx*N - iyy*x17 - izz*x19 - x0*x1 - x1*x4 - x1*x6 - x10*x8
                         - x11*x12 - x11*x13 - x2*x3 - x3*x5 - x3*x7;
         dot_Kab (xy) := x38;
         dot_Kab (xz) := x45;
         dot_Kab (xy) := x38;
         dot_Kab (yy) := -Kxy*Kyy*x10 - (Kyy ** 2)*x33 + Kyy*N*trK - Kyy*x47 + Nx*gxyy*ixx + (1.0/2.0)*Nx*gyyy*ixy + Nx*gyzy*ixz + Ny*gxyy*ixy + (1.0/2.0)*Ny*gyyy*iyy + Ny*gyzy*iyz - Nyy + Nz*gxyy*ixz + (1.0/2.0)*Nz*gyyy*iyz + Nz*gyzy*izz + Ryy*N - gyyz*x22 - gyyz*x25
                         - gyyz*x27 - izz*x49 - x15*x16 - x21*x23 - x23*x24 - x23*x26 - x31*x46;
         dot_Kab (yz) := x53;
         dot_Kab (xz) := x45;
         dot_Kab (yz) := x53;
         dot_Kab (zz) := -Kxz*Kyz*x10 - (Kzz ** 2)*x36 - Kzz*ixz*x11 + Kzz*N*trK - Kzz*x47 + Nx*gxzz*ixx + Nx*gyzz*ixy + (1.0/2.0)*Nx*gzzz*ixz + Ny*gxzz*ixy + Ny*gyzz*iyy + (1.0/2.0)*Ny*gzzz*iyz + Nz*gxzz*ixz + Nz*gyzz*iyz + (1.0/2.0)*Nz*gzzz*izz - Nzz + Rzz*N - gzzx*x50
                         - gzzx*x51 - gzzx*x52 - gzzy*x39 - gzzy*x40 - gzzy*x41 - x15*x18 - x33*x48;
      end set_3d_dot_Kab;

      procedure set_data is
         det : Real;
         i, j, k : Integer;
         x, y, z : Real;
      begin

         i := point.i;
         j := point.j;
         k := point.k;

         x := point.x;
         y := point.y;
         z := point.z;

         N   := ADMBase.N (i,j,k);     -- lapse
         gab := ADMBase.gab (i,j,k);   -- 3-metric
         Kab := ADMBase.Kab (i,j,k);   -- extrinsic curvature

         iab := symm_inverse (gab);    -- 3-metric inverse

         -- second order centred finite differences for the lapse

         declare
            N : LapseGridArray renames ADMBase.N;
         begin

            d1N (1) := (N (i+1,j,k) - N (i-1,j,k)) / (two_dx);
            d1N (2) := (N (i,j+1,k) - N (i,j-1,k)) / (two_dy);
            d1N (3) := (N (i,j,k+1) - N (i,j,k-1)) / (two_dz);

            d2N (xx) := (N (i+1,j,k) - 2.0*N (i,j,k) + N (i-1,j,k)) / (dxdx);
            d2N (yy) := (N (i,j+1,k) - 2.0*N (i,j,k) + N (i,j-1,k)) / (dydy);
            d2N (zz) := (N (i,j,k+1) - 2.0*N (i,j,k) + N (i,j,k-1)) / (dzdz);

            d2N (xy) := (N (i+1,j+1,k) + N (i-1,j-1,k) - N (i+1,j-1,k) - N (i-1,j+1,k)) / (four_dxdy);
            d2N (xz) := (N (i+1,j,k+1) + N (i-1,j,k-1) - N (i+1,j,k-1) - N (i-1,j,k+1)) / (four_dxdz);
            d2N (yz) := (N (i,j+1,k+1) + N (i,j-1,k-1) - N (i,j+1,k-1) - N (i,j-1,k+1)) / (four_dydz);

         end;

         -- second order centred finite differences for the 3-metric

         declare
            gab : MetricGridArray renames ADMBase.gab;
         begin

            d1gab (1) := (gab (i+1,j,k) - gab (i-1,j,k)) / (two_dx);
            d1gab (2) := (gab (i,j+1,k) - gab (i,j-1,k)) / (two_dy);
            d1gab (3) := (gab (i,j,k+1) - gab (i,j,k-1)) / (two_dz);

            d2gab (xx) := (gab (i+1,j,k) - 2.0*gab (i,j,k) + gab (i-1,j,k)) / (dxdx);
            d2gab (yy) := (gab (i,j+1,k) - 2.0*gab (i,j,k) + gab (i,j-1,k)) / (dydy);
            d2gab (zz) := (gab (i,j,k+1) - 2.0*gab (i,j,k) + gab (i,j,k-1)) / (dzdz);

            d2gab (xy) := (gab (i+1,j+1,k) + gab (i-1,j-1,k) - gab (i+1,j-1,k) - gab (i-1,j+1,k)) / (four_dxdy);
            d2gab (xz) := (gab (i+1,j,k+1) + gab (i-1,j,k-1) - gab (i+1,j,k-1) - gab (i-1,j,k+1)) / (four_dxdz);
            d2gab (yz) := (gab (i,j+1,k+1) + gab (i,j-1,k-1) - gab (i,j+1,k-1) - gab (i,j-1,k+1)) / (four_dydz);

         end;

         d1iab (1) := - symm_raise_indices (d1gab (1), iab);
         d1iab (2) := - symm_raise_indices (d1gab (2), iab);
         d1iab (3) := - symm_raise_indices (d1gab (3), iab);

         trK := ixx*Kxx + iyy*Kyy + izz*Kzz + 2.0*( ixy*Kxy + ixz*Kxz + iyz*Kyz );

         set_3d_ricci;

      end set_data;

   begin

      set_data;

      set_3d_dot_N;
      set_3d_dot_gab;
      set_3d_dot_Kab;

      declare
         i, j, k : Integer;
      begin

         i := point.i;
         j := point.j;
         k := point.k;

         ADMBase.dot_N   (i,j,k) := dot_N;     -- lapse
         ADMBase.dot_gab (i,j,k) := dot_gab;   -- 3-metric
         ADMBase.dot_Kab (i,j,k) := dot_Kab;   -- extrinsic curvature

      end;

   end set_time_derivatives;

   procedure set_finite_diff_factors is
   begin

      dxdx := dx * dx;
      dydy := dy * dy;
      dzdz := dz * dz;

      two_dx := 2.0 * dx;
      two_dy := 2.0 * dy;
      two_dz := 2.0 * dz;

      four_dxdy := 4.0 * dx * dy;
      four_dxdz := 4.0 * dx * dz;
      four_dydz := 4.0 * dy * dz;

   end set_finite_diff_factors;

   procedure set_time_derivatives_intr (params : slave_params_record) is

      b : Integer;

      the_task  : Integer := params.slave_id;

      beg_point : Integer := params.params (3);
      end_point : Integer := params.params (4);

   begin

      for a in beg_point .. end_point loop

         b := interior (a);

         set_time_derivatives (grid_point_list(b));

      end loop;

   end set_time_derivatives_intr;

   procedure set_time_derivatives_bndry_ns (params : slave_params_record) is

      b, i, j, k : Integer;

      the_task  : Integer := params.slave_id;

      beg_north : Integer := params.params (7);
      end_north : Integer := params.params (8);
      beg_south : Integer := params.params (9);
      end_south : Integer := params.params (10);

   begin

      -- apply periodic boundary conditions

      -- north/south boundaries

      for a in beg_north .. end_north loop

         b := north_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;  -- equals num_z

           dot_N (i,j,k) :=   dot_N (i,j,2);
         dot_gab (i,j,k) := dot_gab (i,j,2);
         dot_Kab (i,j,k) := dot_Kab (i,j,2);

      end loop;

      for a in beg_south .. end_south loop

         b := south_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;  -- equals 1

           dot_N (i,j,k) :=   dot_N (i,j,num_z-1);
         dot_gab (i,j,k) := dot_gab (i,j,num_z-1);
         dot_Kab (i,j,k) := dot_Kab (i,j,num_z-1);

      end loop;

   end set_time_derivatives_bndry_ns;

   procedure set_time_derivatives_bndry_ew (params : slave_params_record) is

      b, i, j, k : Integer;

      the_task  : Integer := params.slave_id;

      beg_east  : Integer := params.params (11);
      end_east  : Integer := params.params (12);
      beg_west  : Integer := params.params (13);
      end_west  : Integer := params.params (14);

   begin

      -- apply periodic boundary conditions

      -- east/west boundaries

      for a in beg_east .. end_east loop

         b := east_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;  -- equals num_y
         k := grid_point_list (b).k;

           dot_N (i,j,k) :=   dot_N (i,2,k);
         dot_gab (i,j,k) := dot_gab (i,2,k);
         dot_Kab (i,j,k) := dot_Kab (i,2,k);

      end loop;

      for a in beg_west .. end_west loop

         b := west_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;  -- equals 1
         k := grid_point_list (b).k;

           dot_N (i,j,k) :=   dot_N (i,num_y-1,k);
         dot_gab (i,j,k) := dot_gab (i,num_y-1,k);
         dot_Kab (i,j,k) := dot_Kab (i,num_y-1,k);

      end loop;

   end set_time_derivatives_bndry_ew;

   procedure set_time_derivatives_bndry_fb (params : slave_params_record) is

      b, i, j, k : Integer;

      the_task  : Integer := params.slave_id;

      beg_front : Integer := params.params (15);
      end_front : Integer := params.params (16);
      beg_back  : Integer := params.params (17);
      end_back  : Integer := params.params (18);

   begin

      -- apply periodic boundary conditions

      -- front/back boundaries

      for a in beg_front .. end_front loop

         b := front_bndry (a);

         i := grid_point_list (b).i;  -- equals num_x
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;

           dot_N (i,j,k) :=   dot_N (2,j,k);
         dot_gab (i,j,k) := dot_gab (2,j,k);
         dot_Kab (i,j,k) := dot_Kab (2,j,k);

      end loop;

      for a in beg_back .. end_back loop

         b := back_bndry (a);

         i := grid_point_list (b).i;  -- equals 1
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;

           dot_N (i,j,k) :=   dot_N (num_x-1,j,k);
         dot_gab (i,j,k) := dot_gab (num_x-1,j,k);
         dot_Kab (i,j,k) := dot_Kab (num_x-1,j,k);

      end loop;

   end set_time_derivatives_bndry_fb;

   procedure set_time_derivatives
   is
      params : slave_params_record := (1,
                                       (1,grid_point_num,
                                        1,interior_num,
                                        1,boundary_num,
                                        1,north_bndry_num,
                                        1,south_bndry_num,
                                        1,east_bndry_num,
                                        1,west_bndry_num,
                                        1,front_bndry_num,
                                        1,back_bndry_num));
   begin

      set_time_derivatives_intr (params);
      set_time_derivatives_bndry_fb (params);
      set_time_derivatives_bndry_ew (params);
      set_time_derivatives_bndry_ns (params);

   end set_time_derivatives;

end ADMBase.Time_Derivs;
