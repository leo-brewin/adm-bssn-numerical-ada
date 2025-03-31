package body BSSNBase.Time_Derivs is

   dxdx, dydy, dzdz : Real;
   two_dx, two_dy, two_dz : Real;
   four_dxdy, four_dxdz, four_dydz : Real;

   procedure set_time_derivatives (point : GridPoint) is

      -- these are local variables used by the basic codes such as dot-gab.ad, dot-kab.ad etc.
      -- renames are good as they don't involve a memory copy

      ----------------------------------------------------------------------
      -- basic local data

      N         : Real;                                     -- lapse
      phi       : Real;                                     -- conformal factor
      trK       : Real;                                     -- trace of K
      Gi        : GammaPointArray;                          -- trace of conformal connection, GammaBar^{i}
      gBar      : MetricPointArray;                         -- gBar_{ab}
      iBar      : MetricPointArray;                         -- gBar^{ab}
      ABar      : ExtcurvPointArray;                        -- ABar_{ab}
      BBar      : ExtcurvPointArray;                        -- ABar^{ab}

      ----------------------------------------------------------------------
      -- time derivatives of local data

      dot_N     : Real;
      dot_phi   : Real;
      dot_trK   : Real;
      dot_Gi    : GammaPointArray;
      dot_gBar  : MetricPointArray;
      dot_ABar  : ExtcurvPointArray;

      ----------------------------------------------------------------------
      -- partial derivatives of local data

      d1N       : Array (1..3) of Real;                     -- 1st & 2nd derivs of N
      d2N       : Array (symmetric) of Real;

      d1phi     : Array (1..3) of Real;                     -- 1st & 2nd derivs of phi
      d2phi     : Array (symmetric) of Real;

      d1trK     : Array (1..3) of Real;                     -- 1st derivs of trK

      d1Gi      : Array (1..3) of GammaPointArray;          -- 1st derivs of GammaBar^i

      d1gBar    : Array (1..3) of MetricPointArray;         -- 1st & 2nd derivs of gBar_{ab}
      d2gBar    : Array (symmetric) of MetricPointArray;

      d1iBar    : Array (1..3) of MetricPointArray;         -- 1st derivs of gBar^{ab}

      d1ABar    : Array (1..3) of ExtcurvPointArray;        -- 1st derivs of ABar_{ab}

      ----------------------------------------------------------------------
      -- miscellaneous data

      R         : Real;                                     -- Ricci scalar
      Rab       : Array (symmetric) of Real;                -- Ricci tensor

      ----------------------------------------------------------------------
      -- constraints

      Ham       : Real;                                     -- Hamiltonian constraint
      Mom       : MomConstraintPointArray;                  -- Momentum constraints

      ----------------------------------------------------------------------
      -- mappings between the Cadabra and Ada objects
      -- LCB: the following renames could be eliminated by using a set of sed-scripts
      --      to substitute the symbols in the foo.ad files

      -- \partial_{a}{\phi} --> dphi_{a} --> phia --> phi_a --> d1phi (a)

      phi_x : Real renames d1phi (1);
      phi_y : Real renames d1phi (2);
      phi_z : Real renames d1phi (3);

      -- \partial_{a b}{\phi} --> dphi_{a b} --> phiab --> phi_ab --> d2phi (ab)

      phi_xx : Real renames d2phi (xx);
      phi_xy : Real renames d2phi (xy);
      phi_xz : Real renames d2phi (xz);
      phi_yy : Real renames d2phi (yy);
      phi_yz : Real renames d2phi (yz);
      phi_zz : Real renames d2phi (zz);

      -- \partial_{a}{trK} --> d1trK_{a} --> trKa --> trK_a --> d1trK (a)

      trK_x : Real renames d1trK (1);
      trK_y : Real renames d1trK (2);
      trK_z : Real renames d1trK (3);

      -- GammaBar^{a} --> Gia --> Gi_a --> Gi (a)

      Gi_x : Real renames Gi (1);
      Gi_y : Real renames Gi (2);
      Gi_z : Real renames Gi (3);

      -- \partial_{a}{GammaBar^{b}} --> dGi^{b}_{a} --> Giba --> Gi_ba --> d1Gi (a)(b)

      Gi_xx : Real renames d1Gi (1)(1);
      Gi_xy : Real renames d1Gi (2)(1);
      Gi_xz : Real renames d1Gi (3)(1);

      Gi_yx : Real renames d1Gi (1)(2);
      Gi_yy : Real renames d1Gi (2)(2);
      Gi_yz : Real renames d1Gi (3)(2);

      Gi_zx : Real renames d1Gi (1)(3);
      Gi_zy : Real renames d1Gi (2)(3);
      Gi_zz : Real renames d1Gi (3)(3);

      -- \partial_{a}{N} --> dN_{a} --> Na --> N_a --> d1N (a)

      N_x : Real renames d1N (1);
      N_y : Real renames d1N (2);
      N_z : Real renames d1N (3);

      -- \partial_{a b}{N} --> dN_{a b} --> Nab --> N_ab --> d2N (ab)

      N_xx : Real renames d2N (xx);
      N_xy : Real renames d2N (xy);
      N_xz : Real renames d2N (xz);
      N_yy : Real renames d2N (yy);
      N_yz : Real renames d2N (yz);
      N_zz : Real renames d2N (zz);

      -- R_{a b} --> Rab --> Rab (ab)

      Rxx : Real renames Rab (xx);
      Rxy : Real renames Rab (xy);
      Rxz : Real renames Rab (xz);
      Ryy : Real renames Rab (yy);
      Ryz : Real renames Rab (yz);
      Rzz : Real renames Rab (zz);

      -- ABar_{a b} --> ABarab --> ABar_ab --> ABar (ab)

      ABar_xx : Real Renames ABar (xx);
      ABar_xy : Real Renames ABar (xy);
      ABar_xz : Real Renames ABar (xz);
      ABar_yy : Real Renames ABar (yy);
      ABar_yz : Real Renames ABar (yz);
      ABar_zz : Real Renames ABar (zz);

      -- ABar^{a b} --> BBarab --> BBar_ab --> BBar (ab)

      BBar_xx : Real Renames BBar (xx);
      BBar_xy : Real Renames BBar (xy);
      BBar_xz : Real Renames BBar (xz);
      BBar_yy : Real Renames BBar (yy);
      BBar_yz : Real Renames BBar (yz);
      BBar_zz : Real Renames BBar (zz);

      -- gBar_{a b} --> gBarab --> gBar_ab --> gBar (ab)

      gBar_xx : Real Renames gBar (xx);
      gBar_xy : Real Renames gBar (xy);
      gBar_xz : Real Renames gBar (xz);
      gBar_yy : Real Renames gBar (yy);
      gBar_yz : Real Renames gBar (yz);
      gBar_zz : Real Renames gBar (zz);

      -- gBar^{a b} --> iBarab --> iBar_ab --> iBar (ab)

      iBar_xx : Real renames iBar (xx);
      iBar_xy : Real renames iBar (xy);
      iBar_xz : Real renames iBar (xz);
      iBar_yy : Real renames iBar (yy);
      iBar_yz : Real renames iBar (yz);
      iBar_zz : Real renames iBar (zz);

      -- \partial_{a}{gBar_{b c}} --> dgBar_{b c a} --> gBarbca --> gBar_bca --> d1gBar(a)(bc)

      gBar_xxx : Real renames d1gBar (1)(xx);
      gBar_xxy : Real renames d1gBar (2)(xx);
      gBar_xxz : Real renames d1gBar (3)(xx);

      gBar_xyx : Real renames d1gBar (1)(xy);
      gBar_xyy : Real renames d1gBar (2)(xy);
      gBar_xyz : Real renames d1gBar (3)(xy);

      gBar_xzx : Real renames d1gBar (1)(xz);
      gBar_xzy : Real renames d1gBar (2)(xz);
      gBar_xzz : Real renames d1gBar (3)(xz);

      gBar_yyx : Real renames d1gBar (1)(yy);
      gBar_yyy : Real renames d1gBar (2)(yy);
      gBar_yyz : Real renames d1gBar (3)(yy);

      gBar_yzx : Real renames d1gBar (1)(yz);
      gBar_yzy : Real renames d1gBar (2)(yz);
      gBar_yzz : Real renames d1gBar (3)(yz);

      gBar_zzx : Real renames d1gBar (1)(zz);
      gBar_zzy : Real renames d1gBar (2)(zz);
      gBar_zzz : Real renames d1gBar (3)(zz);

      -- \partial_{a}{gBar^{b c}} --> dgBar^{b c}_{a} --> iBarbca --> iBar_bca --> d1iBar(a)(bc)

      iBar_xxx : Real renames d1iBar (1)(xx);
      iBar_xxy : Real renames d1iBar (2)(xx);
      iBar_xxz : Real renames d1iBar (3)(xx);

      iBar_xyx : Real renames d1iBar (1)(xy);
      iBar_xyy : Real renames d1iBar (2)(xy);
      iBar_xyz : Real renames d1iBar (3)(xy);

      iBar_xzx : Real renames d1iBar (1)(xz);
      iBar_xzy : Real renames d1iBar (2)(xz);
      iBar_xzz : Real renames d1iBar (3)(xz);

      iBar_yyx : Real renames d1iBar (1)(yy);
      iBar_yyy : Real renames d1iBar (2)(yy);
      iBar_yyz : Real renames d1iBar (3)(yy);

      iBar_yzx : Real renames d1iBar (1)(yz);
      iBar_yzy : Real renames d1iBar (2)(yz);
      iBar_yzz : Real renames d1iBar (3)(yz);

      iBar_zzx : Real renames d1iBar (1)(zz);
      iBar_zzy : Real renames d1iBar (2)(zz);
      iBar_zzz : Real renames d1iBar (3)(zz);

      -- \partial_{a}{ABar_{b c}} --> dABar_{b c a} --> ABarbca --> ABar_bca --> d1ABar(a)(bc)

      ABar_xxx : Real renames d1ABar (1)(xx);
      ABar_xxy : Real renames d1ABar (2)(xx);
      ABar_xxz : Real renames d1ABar (3)(xx);

      ABar_xyx : Real renames d1ABar (1)(xy);
      ABar_xyy : Real renames d1ABar (2)(xy);
      ABar_xyz : Real renames d1ABar (3)(xy);

      ABar_xzx : Real renames d1ABar (1)(xz);
      ABar_xzy : Real renames d1ABar (2)(xz);
      ABar_xzz : Real renames d1ABar (3)(xz);

      ABar_yyx : Real renames d1ABar (1)(yy);
      ABar_yyy : Real renames d1ABar (2)(yy);
      ABar_yyz : Real renames d1ABar (3)(yy);

      ABar_yzx : Real renames d1ABar (1)(yz);
      ABar_yzy : Real renames d1ABar (2)(yz);
      ABar_yzz : Real renames d1ABar (3)(yz);

      ABar_zzx : Real renames d1ABar (1)(zz);
      ABar_zzy : Real renames d1ABar (2)(zz);
      ABar_zzz : Real renames d1ABar (3)(zz);

      -- \partial_{a b}{gBar_{c d}} --> dgBar_{c d a b} --> gBarcdab --> gBar_cdab --> d2gBar(ab)(cd)

      gBar_xxxx : Real renames d2gBar (xx)(xx);
      gBar_xyxx : Real renames d2gBar (xx)(xy);
      gBar_xzxx : Real renames d2gBar (xx)(xz);
      gBar_yyxx : Real renames d2gBar (xx)(yy);
      gBar_yzxx : Real renames d2gBar (xx)(yz);
      gBar_zzxx : Real renames d2gBar (xx)(zz);

      gBar_xxxy : Real renames d2gBar (xy)(xx);
      gBar_xyxy : Real renames d2gBar (xy)(xy);
      gBar_xzxy : Real renames d2gBar (xy)(xz);
      gBar_yyxy : Real renames d2gBar (xy)(yy);
      gBar_yzxy : Real renames d2gBar (xy)(yz);
      gBar_zzxy : Real renames d2gBar (xy)(zz);

      gBar_xxxz : Real renames d2gBar (xz)(xx);
      gBar_xyxz : Real renames d2gBar (xz)(xy);
      gBar_xzxz : Real renames d2gBar (xz)(xz);
      gBar_yyxz : Real renames d2gBar (xz)(yy);
      gBar_yzxz : Real renames d2gBar (xz)(yz);
      gBar_zzxz : Real renames d2gBar (xz)(zz);

      gBar_xxyy : Real renames d2gBar (yy)(xx);
      gBar_xyyy : Real renames d2gBar (yy)(xy);
      gBar_xzyy : Real renames d2gBar (yy)(xz);
      gBar_yyyy : Real renames d2gBar (yy)(yy);
      gBar_yzyy : Real renames d2gBar (yy)(yz);
      gBar_zzyy : Real renames d2gBar (yy)(zz);

      gBar_xxyz : Real renames d2gBar (yz)(xx);
      gBar_xyyz : Real renames d2gBar (yz)(xy);
      gBar_xzyz : Real renames d2gBar (yz)(xz);
      gBar_yyyz : Real renames d2gBar (yz)(yy);
      gBar_yzyz : Real renames d2gBar (yz)(yz);
      gBar_zzyz : Real renames d2gBar (yz)(zz);

      gBar_xxzz : Real renames d2gBar (zz)(xx);
      gBar_xyzz : Real renames d2gBar (zz)(xy);
      gBar_xzzz : Real renames d2gBar (zz)(xz);
      gBar_yyzz : Real renames d2gBar (zz)(yy);
      gBar_yzzz : Real renames d2gBar (zz)(yz);
      gBar_zzzz : Real renames d2gBar (zz)(zz);

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
         x182, x183, x184, x185 : Real;
      begin
         x0   := 2.0*phi_xx;
         x1   := (1.0/2.0)*iBar_xx;
         x2   := (1.0/2.0)*iBar_yy;
         x3   := (1.0/2.0)*iBar_zz;
         x4   := (1.0/2.0)*gBar_xyx;
         x5   := (1.0/2.0)*iBar_xzx;
         x6   := iBar_xy*phi_x;
         x7   := iBar_yy*phi_y;
         x8   := iBar_yz*phi_z;
         x9   := iBar_xz*phi_x;
         x10  := iBar_yz*phi_y;
         x11  := iBar_zz*phi_z;
         x12  := (phi_x ** 2);
         x13  := gBar_xx*iBar_xx;
         x14  := 4.0*gBar_xx;
         x15  := iBar_xy*phi_xy;
         x16  := iBar_xz*phi_xz;
         x17  := 2.0*phi_yy;
         x18  := gBar_xx*iBar_yy;
         x19  := iBar_yz*phi_yz;
         x20  := 2.0*phi_zz;
         x21  := gBar_xx*iBar_zz;
         x22  := gBar_xxy*iBar_xz;
         x23  := iBar_xx*iBar_yy;
         x24  := gBar_xxy*iBar_yz;
         x25  := gBar_xzy*iBar_xy;
         x26  := gBar_xxz*iBar_yz;
         x27  := gBar_xyz*iBar_yz;
         x28  := iBar_xx*iBar_zz;
         x29  := gBar_xyz*iBar_xy;
         x30  := gBar_xyx*iBar_yz;
         x31  := gBar_xzx*iBar_xy;
         x32  := gBar_xyx*iBar_xz;
         x33  := gBar_xzy*iBar_xz;
         x34  := iBar_zz*x29;
         x35  := gBar_xyz*gBar_xzy;
         x36  := iBar_yy*iBar_zz;
         x37  := gBar_xzy*iBar_yz;
         x38  := gBar_xzx*x37;
         x39  := iBar_xy*phi_y;
         x40  := 24.0*gBar_xx;
         x41  := phi_x*x40;
         x42  := iBar_xz*phi_z;
         x43  := phi_y*x8;
         x44  := iBar_xxx + iBar_xyy + iBar_xzz;
         x45  := (1.0/2.0)*gBar_xxx;
         x46  := iBar_xyx + iBar_yyy + iBar_yzz;
         x47  := (1.0/2.0)*gBar_xxy;
         x48  := iBar_xzx + iBar_yzy + iBar_zzz;
         x49  := (1.0/2.0)*gBar_xxz;
         x50  := 12.0*x12;
         x51  := (phi_y ** 2);
         x52  := 12.0*x51;
         x53  := (phi_z ** 2);
         x54  := 12.0*x53;
         x55  := (iBar_xy ** 2);
         x56  := (gBar_xxy ** 2);
         x57  := (iBar_xz ** 2);
         x58  := (gBar_xxz ** 2);
         x59  := (gBar_xyx ** 2);
         x60  := (iBar_yz ** 2);
         x61  := (gBar_xyz ** 2);
         x62  := (gBar_xzx ** 2);
         x63  := (gBar_xzy ** 2);
         x64  := (1.0/2.0)*x55;
         x65  := (1.0/2.0)*x57;
         x66  := (1.0/2.0)*x60;
         x67  := 2.0*gBar_xx;
         x68  := phi_x*x44;
         x69  := phi_y*x46;
         x70  := phi_z*x48;
         x71  := (1.0/2.0)*gBar_xzx;
         x72  := (1.0/2.0)*gBar_xyy;
         x73  := (1.0/2.0)*gBar_xyz;
         x74  := (1.0/2.0)*gBar_xzy;
         x75  := (1.0/2.0)*gBar_xzz;
         x76  := (1.0/2.0)*iBar_xyx;
         x77  := (1.0/2.0)*iBar_yzx;
         x78  := (1.0/2.0)*gBar_yzz;
         x79  := x49*x57;
         x80  := gBar_xyz*x66;
         x81  := gBar_xzx*x65;
         x82  := gBar_xzy*x66;
         x83  := (1.0/2.0)*iBar_yz;
         x84  := gBar_xy*iBar_xx;
         x85  := gBar_xy*iBar_yy;
         x86  := gBar_xy*iBar_zz;
         x87  := 2.0*gBar_xy;
         x88  := 4.0*gBar_xy;
         x89  := iBar_xz*x29;
         x90  := iBar_yy*x1;
         x91  := gBar_xxy*x90;
         x92  := iBar_xz*iBar_yy;
         x93  := x47*x92;
         x94  := x1*x24;
         x95  := gBar_yzy*iBar_yz;
         x96  := iBar_xy*x47;
         x97  := gBar_xyy*iBar_xy;
         x98  := iBar_xz*x49;
         x99  := x1*x26;
         x100 := gBar_yyz*iBar_yz;
         x101 := iBar_zz*x1;
         x102 := gBar_xxz*x101;
         x103 := iBar_xy*iBar_zz;
         x104 := gBar_xyx*x90;
         x105 := x1*x27;
         x106 := iBar_xy*x4;
         x107 := gBar_yzx*iBar_xz;
         x108 := iBar_xy*x107;
         x109 := x2*x32;
         x110 := gBar_xyz*iBar_xz;
         x111 := iBar_xy*x37;
         x112 := gBar_xzx*x101;
         x113 := iBar_xy*x3;
         x114 := gBar_yyx*iBar_xy;
         x115 := iBar_yz*x73;
         x116 := iBar_zz*x2;
         x117 := gBar_xyz*x116;
         x118 := iBar_xz*x71;
         x119 := x3*x31;
         x120 := gBar_xzy*x116;
         x121 := iBar_yz*x74;
         x122 := 24.0*gBar_xy;
         x123 := phi_x*x122;
         x124 := -1.0/2.0*Gi_xx*gBar_xy - 1.0/2.0*Gi_xy*gBar_xx - 1.0/2.0*Gi_yx*gBar_yy - 1.0/2.0*Gi_yy*gBar_xy - 1.0/2.0*Gi_zx*gBar_yz - 1.0/2.0*Gi_zy*gBar_xz - 1.0/2.0*gBar_xxy*gBar_xyy*iBar_xx*iBar_yy - 1.0/2.0*gBar_xxy*gBar_xyz*iBar_xx*iBar_yz - 1.0/2.0*gBar_xxy*gBar_yyx*x55
               - 1.0/2.0*gBar_xxy*gBar_yyz*iBar_xy*iBar_yz - 1.0/2.0*gBar_xxy*gBar_yzx*iBar_xy*iBar_xz - 1.0/2.0*gBar_xxy*gBar_yzy*iBar_xz*iBar_yy - gBar_xxy*iBar_xx*phi_x - gBar_xxy*iBar_xy*phi_y - gBar_xxy*iBar_xz*phi_z - 1.0/2.0*gBar_xxz*gBar_xyy*iBar_xx*iBar_yz
               - 1.0/2.0*gBar_xxz*gBar_xyz*iBar_xx*iBar_zz - 1.0/2.0*gBar_xxz*gBar_yyx*iBar_xy*iBar_xz - 1.0/2.0*gBar_xxz*gBar_yyz*iBar_xy*iBar_zz - 1.0/2.0*gBar_xxz*gBar_yzx*x57 - 1.0/2.0*gBar_xxz*gBar_yzy*iBar_xz*iBar_yz - 1.0/2.0*gBar_xyx*gBar_xyy*x55
               - 1.0/2.0*gBar_xyx*gBar_xyz*iBar_xy*iBar_xz - 1.0/2.0*gBar_xyx*gBar_yyx*iBar_xx*iBar_yy - 1.0/2.0*gBar_xyx*gBar_yyz*iBar_xz*iBar_yy - 1.0/2.0*gBar_xyx*gBar_yzx*iBar_xx*iBar_yz - 1.0/2.0*gBar_xyx*gBar_yzy*iBar_xy*iBar_yz + gBar_xyx*x105 + gBar_xyxx*x1 + gBar_xyxy*iBar_xy
               + gBar_xyxz*iBar_xz - 1.0/2.0*gBar_xyy*gBar_xyz*iBar_xy*iBar_yz - 1.0/2.0*gBar_xyy*gBar_xzx*iBar_xy*iBar_xz + gBar_xyy*gBar_xzx*iBar_yz*x1 - 1.0/2.0*gBar_xyy*gBar_xzy*iBar_xz*iBar_yy + gBar_xyy*x104 + gBar_xyy*x110*x2 + gBar_xyy*x47*x55 + gBar_xyyy*x2 + gBar_xyyz*iBar_yz
               - 1.0/2.0*gBar_xyz*gBar_xzx*x57 - 1.0/2.0*gBar_xyz*gBar_xzy*iBar_xz*iBar_yz - 1.0/2.0*gBar_xyz*gBar_yyx*iBar_xz*iBar_yy - 1.0/2.0*gBar_xyz*gBar_yyz*iBar_yy*iBar_zz - 1.0/2.0*gBar_xyz*gBar_yzx*iBar_xz*iBar_yz - 1.0/2.0*gBar_xyz*gBar_yzy*x60 + gBar_xyz*x10 + gBar_xyz*x11
               + gBar_xyz*x112 + gBar_xyz*x5 + gBar_xyz*x79 + gBar_xyz*x9 + gBar_xyzz*x3 - 1.0/2.0*gBar_xzx*gBar_yyx*iBar_xx*iBar_yz - 1.0/2.0*gBar_xzx*gBar_yyz*iBar_xz*iBar_yz - 1.0/2.0*gBar_xzx*gBar_yzx*iBar_xx*iBar_zz - 1.0/2.0*gBar_xzx*gBar_yzy*iBar_xy*iBar_zz
               - 1.0/2.0*gBar_xzy*gBar_yyx*iBar_xy*iBar_yz - 1.0/2.0*gBar_xzy*gBar_yyz*x60 - 1.0/2.0*gBar_xzy*gBar_yzx*iBar_xy*iBar_zz - 1.0/2.0*gBar_xzy*gBar_yzy*iBar_yy*iBar_zz - gBar_xzy*iBar_xz*phi_x - gBar_xzy*iBar_yz*phi_y - gBar_xzy*iBar_zz*phi_z - gBar_yyx*iBar_xy*phi_x
               - gBar_yyx*iBar_yy*phi_y - gBar_yyx*iBar_yz*phi_z + gBar_yyx*x2*x33 + gBar_yyx*x4*x55 + gBar_yyx*x76 + gBar_yyx*x91 + gBar_yyx*x99 + gBar_yyz*x119 + gBar_yyz*x120 + gBar_yyz*x77 + gBar_yyz*x80 + gBar_yyz*x93 - gBar_yzx*iBar_xz*phi_x - gBar_yzx*iBar_yz*phi_y
               - gBar_yzx*iBar_zz*phi_z + gBar_yzx*x102 + gBar_yzx*x29*x3 + gBar_yzx*x5 + gBar_yzx*x81 + gBar_yzx*x94 + gBar_yzy*x103*x49 + gBar_yzy*x109 + gBar_yzy*x117 + gBar_yzy*x82 + iBar_xxx*x4 - 1.0/2.0*iBar_xy*iBar_zz*x61 + iBar_xyy*x47 + iBar_xz*x61*x83 + iBar_xzy*x49
               + iBar_yyy*x72 + iBar_yzy*x73 + iBar_yzy*x74 + iBar_zzx*x78 + iBar_zzy*x75 - 12.0*phi_x*phi_y + 2.0*phi_xy + x0*x84 + x100*x106 + x100*x98 + x107*x121 + x108*x4 + x111*x72 + x113*x35 + x114*x115 + x114*x118 + x118*x95 + x122*x43 + x123*x39 + x123*x42 + x15*x88
               + x16*x88 + x17*x85 + x19*x88 + x20*x86 + x4*x44 + x46*x72 + x47*x89 + x48*x73 + x50*x84 + x52*x85 + x54*x86 + x68*x87 + x69*x87 + x70*x87 + x95*x96 + x97*x98;
         x125 := (1.0/2.0)*iBar_yyx;
         x126 := gBar_xz*iBar_xx;
         x127 := gBar_xz*iBar_yy;
         x128 := gBar_xz*iBar_zz;
         x129 := 2.0*gBar_xz;
         x130 := 4.0*gBar_xz;
         x131 := iBar_xy*iBar_xz;
         x132 := gBar_zzy*iBar_yz;
         x133 := gBar_yzz*iBar_yz;
         x134 := gBar_zzy*x103;
         x135 := iBar_xz*x2;
         x136 := gBar_xzz*iBar_xz;
         x137 := gBar_yzx*iBar_yz;
         x138 := iBar_xy*x137;
         x139 := x29*x3;
         x140 := gBar_zzx*iBar_xz;
         x141 := 24.0*gBar_xz;
         x142 := phi_x*x141;
         x143 := -1.0/2.0*Gi_xx*gBar_xz - 1.0/2.0*Gi_xz*gBar_xx - 1.0/2.0*Gi_yx*gBar_yz - 1.0/2.0*Gi_yz*gBar_xy - 1.0/2.0*Gi_zx*gBar_zz - 1.0/2.0*Gi_zz*gBar_xz - 1.0/2.0*gBar_xxy*gBar_xzy*iBar_xx*iBar_yy - 1.0/2.0*gBar_xxy*gBar_xzz*iBar_xx*iBar_yz - 1.0/2.0*gBar_xxy*gBar_yzx*x55
               - 1.0/2.0*gBar_xxy*gBar_yzz*iBar_xy*iBar_yz - 1.0/2.0*gBar_xxy*gBar_zzx*iBar_xy*iBar_xz - 1.0/2.0*gBar_xxy*gBar_zzy*iBar_xz*iBar_yy - 1.0/2.0*gBar_xxz*gBar_xzy*iBar_xx*iBar_yz - 1.0/2.0*gBar_xxz*gBar_xzz*iBar_xx*iBar_zz - 1.0/2.0*gBar_xxz*gBar_yzx*iBar_xy*iBar_xz
               - 1.0/2.0*gBar_xxz*gBar_yzz*iBar_xy*iBar_zz - 1.0/2.0*gBar_xxz*gBar_zzx*x57 - 1.0/2.0*gBar_xxz*gBar_zzy*iBar_xz*iBar_yz - gBar_xxz*iBar_xx*phi_x - gBar_xxz*iBar_xy*phi_y - gBar_xxz*iBar_xz*phi_z - 1.0/2.0*gBar_xyx*gBar_xzy*x55 - 1.0/2.0*gBar_xyx*gBar_xzz*iBar_xy*iBar_xz
               - 1.0/2.0*gBar_xyx*gBar_yzx*iBar_xx*iBar_yy - 1.0/2.0*gBar_xyx*gBar_yzz*iBar_xz*iBar_yy - 1.0/2.0*gBar_xyx*gBar_zzx*iBar_xx*iBar_yz - 1.0/2.0*gBar_xyx*gBar_zzy*iBar_xy*iBar_yz - 1.0/2.0*gBar_xyz*gBar_xzy*iBar_xy*iBar_yz - 1.0/2.0*gBar_xyz*gBar_xzz*iBar_xy*iBar_zz
               - 1.0/2.0*gBar_xyz*gBar_yzx*iBar_xz*iBar_yy - 1.0/2.0*gBar_xyz*gBar_yzz*iBar_yy*iBar_zz - 1.0/2.0*gBar_xyz*gBar_zzx*iBar_xz*iBar_yz - 1.0/2.0*gBar_xyz*gBar_zzy*x60 - gBar_xyz*iBar_xy*phi_x - gBar_xyz*iBar_yy*phi_y - gBar_xyz*iBar_yz*phi_z - 1.0/2.0*gBar_xzx*gBar_xzy*iBar_xy*iBar_xz
               - 1.0/2.0*gBar_xzx*gBar_xzz*x57 - 1.0/2.0*gBar_xzx*gBar_yzx*iBar_xx*iBar_yz - 1.0/2.0*gBar_xzx*gBar_yzz*iBar_xz*iBar_yz - 1.0/2.0*gBar_xzx*gBar_zzx*iBar_xx*iBar_zz - 1.0/2.0*gBar_xzx*gBar_zzy*iBar_xy*iBar_zz + gBar_xzxx*x1 + gBar_xzxy*iBar_xy + gBar_xzxz*iBar_xz
               - 1.0/2.0*gBar_xzy*gBar_xzz*iBar_xz*iBar_yz - 1.0/2.0*gBar_xzy*gBar_yzx*iBar_xy*iBar_yz - 1.0/2.0*gBar_xzy*gBar_yzz*x60 - 1.0/2.0*gBar_xzy*gBar_zzx*iBar_xy*iBar_zz - 1.0/2.0*gBar_xzy*gBar_zzy*iBar_yy*iBar_zz + gBar_xzy*x104 + gBar_xzy*x47*x55 + gBar_xzy*x6 + gBar_xzy*x7
               + gBar_xzy*x8 + gBar_xzyy*x2 + gBar_xzyz*iBar_yz + gBar_xzz*x1*x30 + gBar_xzz*x112 + gBar_xzz*x131*x47 + gBar_xzz*x25*x3 + gBar_xzz*x79 + gBar_xzzz*x3 - gBar_yzx*iBar_xy*phi_x - gBar_yzx*iBar_yy*phi_y - gBar_yzx*iBar_yz*phi_z + gBar_yzx*x2*x33 + gBar_yzx*x4*x55 + gBar_yzx*x76
               + gBar_yzx*x91 + gBar_yzx*x99 + gBar_yzy*x125 + gBar_yzz*x119 + gBar_yzz*x120 + gBar_yzz*x80 + gBar_yzz*x93 - gBar_zzx*iBar_xz*phi_x - gBar_zzx*iBar_yz*phi_y - gBar_zzx*iBar_zz*phi_z + gBar_zzx*x102 + gBar_zzx*x131*x4 + gBar_zzx*x139 + gBar_zzx*x5 + gBar_zzx*x81 + gBar_zzx*x94
               + gBar_zzy*x109 + gBar_zzy*x117 + gBar_zzy*x77 + gBar_zzy*x82 + iBar_xxx*x71 + iBar_xy*x63*x83 + iBar_xyx*x74 + iBar_xyz*x47 - 1.0/2.0*iBar_xz*iBar_yy*x63 + iBar_xzz*x49 + iBar_yyz*x72 + iBar_yzz*x73 + iBar_yzz*x74 + iBar_zzz*x75 - 12.0*phi_x*phi_z + 2.0*phi_xz + x0*x126
               + x1*x38 + x106*x133 + x108*x71 + x115*x136 + x118*x132 + x121*x140 + x126*x50 + x127*x17 + x127*x52 + x128*x20 + x128*x54 + x129*x68 + x129*x69 + x129*x70 + x130*x15 + x130*x16 + x130*x19 + x132*x96 + x133*x98 + x134*x49 + x135*x35 + x138*x73 + x141*x43
               + x142*x39 + x142*x42 + x25*x98 + x44*x71 + x46*x74 + x48*x75;
         x144 := (1.0/2.0)*iBar_yzy;
         x145 := iBar_xx*phi_x;
         x146 := gBar_yy*iBar_xx;
         x147 := 4.0*gBar_yy;
         x148 := gBar_yy*iBar_yy;
         x149 := gBar_yy*iBar_zz;
         x150 := 24.0*gBar_yy;
         x151 := phi_x*x150;
         x152 := (1.0/2.0)*gBar_yyx;
         x153 := (1.0/2.0)*gBar_yyy;
         x154 := (1.0/2.0)*gBar_yyz;
         x155 := (gBar_xyy ** 2);
         x156 := (gBar_yyx ** 2);
         x157 := (gBar_yyz ** 2);
         x158 := (gBar_yzx ** 2);
         x159 := (gBar_yzy ** 2);
         x160 := 2.0*gBar_yy;
         x161 := (1.0/2.0)*gBar_yzy;
         x162 := (1.0/2.0)*gBar_yzx;
         x163 := (1.0/2.0)*gBar_zzx;
         x164 := (1.0/2.0)*x131;
         x165 := gBar_yz*iBar_xx;
         x166 := gBar_yz*iBar_yy;
         x167 := gBar_yz*iBar_zz;
         x168 := 2.0*gBar_yz;
         x169 := 4.0*gBar_yz;
         x170 := iBar_yz*x1;
         x171 := 24.0*gBar_yz;
         x172 := phi_x*x171;
         x173 := -1.0/2.0*Gi_xy*gBar_xz - 1.0/2.0*Gi_xz*gBar_xy - 1.0/2.0*Gi_yy*gBar_yz - 1.0/2.0*Gi_yz*gBar_yy - 1.0/2.0*Gi_zy*gBar_zz - 1.0/2.0*Gi_zz*gBar_yz - 1.0/2.0*gBar_xyy*gBar_xzy*iBar_xx*iBar_yy + gBar_xyy*gBar_xzy*x64 - 1.0/2.0*gBar_xyy*gBar_xzz*iBar_xx*iBar_yz
               - 1.0/2.0*gBar_xyy*gBar_yzx*x55 + gBar_xyy*gBar_yzx*x90 - 1.0/2.0*gBar_xyy*gBar_yzz*iBar_xy*iBar_yz + gBar_xyy*gBar_yzz*x135 - 1.0/2.0*gBar_xyy*gBar_zzx*iBar_xy*iBar_xz + gBar_xyy*gBar_zzx*x170 - 1.0/2.0*gBar_xyy*gBar_zzy*iBar_xz*iBar_yy
               - 1.0/2.0*gBar_xyz*gBar_xzy*iBar_xx*iBar_yz - 1.0/2.0*gBar_xyz*gBar_xzz*iBar_xx*iBar_zz + gBar_xyz*gBar_xzz*x65 - 1.0/2.0*gBar_xyz*gBar_yzx*iBar_xy*iBar_xz - 1.0/2.0*gBar_xyz*gBar_yzz*iBar_xy*iBar_zz + gBar_xyz*gBar_zzx*x101 - 1.0/2.0*gBar_xyz*gBar_zzx*x57
               - 1.0/2.0*gBar_xyz*gBar_zzy*iBar_xz*iBar_yz - gBar_xyz*iBar_xx*phi_x - gBar_xyz*iBar_xy*phi_y - gBar_xyz*iBar_xz*phi_z - 1.0/2.0*gBar_xzy*gBar_yyx*x55 + gBar_xzy*gBar_yyx*x90 - 1.0/2.0*gBar_xzy*gBar_yyz*iBar_xy*iBar_yz - 1.0/2.0*gBar_xzy*gBar_yzx*iBar_xy*iBar_xz
               - 1.0/2.0*gBar_xzy*gBar_yzy*iBar_xz*iBar_yy - gBar_xzy*iBar_xx*phi_x - gBar_xzy*iBar_xy*phi_y - gBar_xzy*iBar_xz*phi_z - 1.0/2.0*gBar_xzz*gBar_yyx*iBar_xy*iBar_xz + gBar_xzz*gBar_yyx*x170 - 1.0/2.0*gBar_xzz*gBar_yyz*iBar_xy*iBar_zz + gBar_xzz*gBar_yzx*x101
               - 1.0/2.0*gBar_xzz*gBar_yzx*x57 - 1.0/2.0*gBar_xzz*gBar_yzy*iBar_xz*iBar_yz + gBar_xzz*gBar_yzy*x113 + gBar_xzz*x131*x72 - 1.0/2.0*gBar_yyx*gBar_yzx*iBar_xx*iBar_yy + gBar_yyx*gBar_yzx*x64 - 1.0/2.0*gBar_yyx*gBar_yzz*iBar_xz*iBar_yy - 1.0/2.0*gBar_yyx*gBar_zzx*iBar_xx*iBar_yz
               - 1.0/2.0*gBar_yyx*gBar_zzy*iBar_xy*iBar_yz + gBar_yyx*gBar_zzy*x135 - 1.0/2.0*gBar_yyz*gBar_yzx*iBar_xz*iBar_yy - 1.0/2.0*gBar_yyz*gBar_yzz*iBar_yy*iBar_zz + gBar_yyz*gBar_yzz*x66 - 1.0/2.0*gBar_yyz*gBar_zzx*iBar_xz*iBar_yz + gBar_yyz*gBar_zzx*x113 + gBar_yyz*gBar_zzy*x116
               - 1.0/2.0*gBar_yyz*gBar_zzy*x60 - gBar_yyz*iBar_xy*phi_x - gBar_yyz*iBar_yy*phi_y - gBar_yyz*iBar_yz*phi_z + gBar_yyz*x2*x33 - 1.0/2.0*gBar_yzx*gBar_yzy*iBar_xy*iBar_yz + gBar_yzx*gBar_yzy*x135 - 1.0/2.0*gBar_yzx*gBar_yzz*iBar_xz*iBar_yz + gBar_yzx*gBar_yzz*x113
               - 1.0/2.0*gBar_yzx*gBar_zzx*iBar_xx*iBar_zz + gBar_yzx*gBar_zzx*x65 - 1.0/2.0*gBar_yzx*gBar_zzy*iBar_xy*iBar_zz + gBar_yzx*x1*x37 + gBar_yzx*x105 + gBar_yzx*x145 + gBar_yzx*x39 + gBar_yzx*x42 + gBar_yzxx*x1 + gBar_yzxy*iBar_xy + gBar_yzxz*iBar_xz + gBar_yzy*gBar_yzz*x116
               - 1.0/2.0*gBar_yzy*gBar_yzz*x60 - 1.0/2.0*gBar_yzy*gBar_zzx*iBar_xy*iBar_zz - 1.0/2.0*gBar_yzy*gBar_zzy*iBar_yy*iBar_zz + gBar_yzy*gBar_zzy*x66 + gBar_yzyy*x2 + gBar_yzyz*iBar_yz + gBar_yzzz*x3 + gBar_zzx*x131*x152 - gBar_zzy*iBar_xz*phi_x - gBar_zzy*iBar_yz*phi_y
               - gBar_zzy*iBar_zz*phi_z + gBar_zzy*x139 + gBar_zzy*x144 - 1.0/2.0*iBar_xx*iBar_yz*x158 + iBar_xxy*x71 + iBar_xxz*x4 + iBar_xy*x132*x72 + iBar_xy*x74*x95 + iBar_xyy*x162 + iBar_xyy*x74 + iBar_xyz*x152 + iBar_xz*x100*x75 + iBar_xz*x132*x162 + iBar_xz*x133*x73
               + iBar_xzy*x163 + iBar_xzz*x162 + iBar_xzz*x73 + iBar_yyy*x161 + iBar_yz*x114*x78 + iBar_yz*x140*x161 + iBar_yzz*x154 + iBar_zzz*x78 - 12.0*phi_y*phi_z + 2.0*phi_yz + x0*x165 + x138*x154 + x15*x169 + x158*x164 + x16*x169 + x161*x46 + x162*x44 + x164*x35 + x165*x50
               + x166*x17 + x166*x52 + x167*x20 + x167*x54 + x168*x68 + x168*x69 + x168*x70 + x169*x19 + x171*x43 + x172*x39 + x172*x42 + x48*x78;
         x174 := (1.0/2.0)*gBar_zzz;
         x175 := gBar_zz*iBar_xx;
         x176 := 4.0*gBar_zz;
         x177 := gBar_zz*iBar_yy;
         x178 := gBar_zz*iBar_zz;
         x179 := 24.0*gBar_zz;
         x180 := phi_x*x179;
         x181 := (gBar_xzz ** 2);
         x182 := (gBar_yzz ** 2);
         x183 := (gBar_zzx ** 2);
         x184 := (gBar_zzy ** 2);
         x185 := 2.0*gBar_zz;
         Rab (xx) := Gi_xx*gBar_xx + Gi_yx*gBar_xy + Gi_zx*gBar_xz + gBar_xxx*iBar_xx*phi_x - 3.0/4.0*gBar_xxx*iBar_xxx + gBar_xxx*iBar_xy*phi_y + gBar_xxx*iBar_xz*phi_z - gBar_xxxx*x1 - gBar_xxxy*iBar_xy - gBar_xxxz*iBar_xz + gBar_xxy*gBar_xxz*iBar_xx*iBar_yz - gBar_xxy*gBar_xyx*x23
                     + gBar_xxy*gBar_xyx*x55 + gBar_xxy*gBar_xyz*iBar_xy*iBar_yz + gBar_xxy*gBar_xzx*iBar_xy*iBar_xz + gBar_xxy*gBar_xzy*iBar_xz*iBar_yy - gBar_xxy*iBar_xyx - gBar_xxy*x6 - gBar_xxy*x7 - gBar_xxy*x8 - gBar_xxyy*x2 - gBar_xxyz*iBar_yz + gBar_xxz*gBar_xyx*iBar_xy*iBar_xz
                     + gBar_xxz*gBar_xyz*iBar_xy*iBar_zz - gBar_xxz*gBar_xzx*x28 + gBar_xxz*gBar_xzx*x57 + gBar_xxz*gBar_xzy*iBar_xz*iBar_yz - gBar_xxz*iBar_xy*x22 - gBar_xxz*iBar_xz*x27 - gBar_xxz*iBar_xzx - gBar_xxz*iBar_zz*x25 - gBar_xxz*x10 - gBar_xxz*x11 - gBar_xxz*x9 - gBar_xxzz*x3
                     + gBar_xyx*gBar_xyz*iBar_xz*iBar_yy + gBar_xyx*gBar_xzx*iBar_xx*iBar_yz + gBar_xyx*gBar_xzy*iBar_xy*iBar_yz - gBar_xyx*iBar_xx*x26 + 2.0*gBar_xyx*iBar_xy*phi_x + 2.0*gBar_xyx*iBar_yy*phi_y - gBar_xyx*iBar_yy*x33 + 2.0*gBar_xyx*iBar_yz*phi_z - gBar_xyy*iBar_yyx
                     + gBar_xyz*gBar_xzx*iBar_xz*iBar_yz + gBar_xyz*gBar_xzy*x60 - gBar_xyz*iBar_yy*x22 - gBar_xyz*iBar_yzx + gBar_xzx*gBar_xzy*iBar_xy*iBar_zz - gBar_xzx*iBar_xx*x24 + 2.0*gBar_xzx*iBar_xz*phi_x + 2.0*gBar_xzx*iBar_yz*phi_y + 2.0*gBar_xzx*iBar_zz*phi_z - gBar_xzx*x34 - gBar_xzx*x5
                     - gBar_xzy*iBar_yzx - gBar_xzz*iBar_zzx + (1.0/4.0)*gBar_yyx*iBar_yyx + (1.0/2.0)*gBar_yzx*iBar_yzx + (1.0/4.0)*gBar_zzx*iBar_zzx + (1.0/2.0)*iBar_xx*iBar_yy*x56 + (1.0/2.0)*iBar_xx*iBar_yy*x59 + (1.0/2.0)*iBar_xx*iBar_zz*x58 + (1.0/2.0)*iBar_xx*iBar_zz*x62 - iBar_xyx*x4
                     - iBar_xz*x38 + (1.0/2.0)*iBar_yy*iBar_zz*x61 + (1.0/2.0)*iBar_yy*iBar_zz*x63 - x0*x13 - x0 + 12.0*x12 - x13*x50 - x14*x15 - x14*x16 - x14*x19 - x17*x18 - x18*x52 - x20*x21 - x21*x54 - x24*x25 - x29*x30 - x31*x32 - x35*x36 - x39*x41 - x40*x43 - x41*x42 - x44*x45
                     - x46*x47 - x48*x49 - x56*x64 - x58*x65 - x59*x64 - x61*x66 - x62*x65 - x63*x66 - x67*x68 - x67*x69 - x67*x70;
         Rab (xy) := (1.0/4.0)*gBar_xxy*iBar_xxx + (1.0/2.0)*gBar_xzy*iBar_xzx - 1.0/4.0*gBar_yyy*iBar_yyx + (1.0/4.0)*gBar_zzy*iBar_zzx - iBar_xxy*x45 - iBar_xyy*x4 - iBar_xzy*x71 - x124;
         Rab (xz) := (1.0/4.0)*gBar_xxz*iBar_xxx + (1.0/2.0)*gBar_xyz*iBar_xyx + (1.0/4.0)*gBar_yyz*iBar_yyx - 1.0/4.0*gBar_zzz*iBar_zzx - iBar_xxz*x45 - iBar_xyz*x4 - iBar_xzz*x71 - x143;
         Rab (xy) := -1.0/4.0*gBar_xxx*iBar_xxy + (1.0/4.0)*gBar_yyx*iBar_yyy - gBar_yyy*x125 + (1.0/2.0)*gBar_yzx*iBar_yzy - gBar_yzy*x77 + (1.0/4.0)*gBar_zzx*iBar_zzy - iBar_xyx*x72 - x124;
         Rab (yy) := Gi_xy*gBar_xy + Gi_yy*gBar_yy + Gi_zy*gBar_yz + (1.0/4.0)*gBar_xxy*iBar_xxy - gBar_xyx*iBar_xxy + gBar_xyy*gBar_xyz*iBar_xx*iBar_yz - gBar_xyy*gBar_yyx*x23 + gBar_xyy*gBar_yyx*x55 + gBar_xyy*gBar_yyz*iBar_xy*iBar_yz - gBar_xyy*gBar_yyz*x92 + gBar_xyy*gBar_yzx*iBar_xy*iBar_xz
                     + gBar_xyy*gBar_yzy*iBar_xz*iBar_yy + 2.0*gBar_xyy*iBar_xx*phi_x - gBar_xyy*iBar_xx*x137 + 2.0*gBar_xyy*iBar_xy*phi_y + 2.0*gBar_xyy*iBar_xz*phi_z - gBar_xyy*x89 + gBar_xyz*gBar_yyx*iBar_xy*iBar_xz + gBar_xyz*gBar_yyz*iBar_xy*iBar_zz - gBar_xyz*gBar_yzx*x28 + gBar_xyz*gBar_yzx*x57
                     + gBar_xyz*gBar_yzy*iBar_xz*iBar_yz - gBar_xyz*iBar_xzy + (1.0/2.0)*gBar_xzy*iBar_xzy + gBar_yyx*gBar_yyz*iBar_xz*iBar_yy + gBar_yyx*gBar_yzx*iBar_xx*iBar_yz + gBar_yyx*gBar_yzy*iBar_xy*iBar_yz - gBar_yyx*gBar_yzy*x92 - gBar_yyx*iBar_xx*x27 - gBar_yyx*iBar_xyy - gBar_yyx*x145
                     - gBar_yyx*x39 - gBar_yyx*x42 - gBar_yyxx*x1 - gBar_yyxy*iBar_xy - gBar_yyxz*iBar_xz + gBar_yyy*iBar_xy*phi_x + gBar_yyy*iBar_yy*phi_y - 3.0/4.0*gBar_yyy*iBar_yyy + gBar_yyy*iBar_yz*phi_z - gBar_yyyy*x2 - gBar_yyyz*iBar_yz + gBar_yyz*gBar_yzx*iBar_xz*iBar_yz
                     - gBar_yyz*gBar_yzx*x103 - gBar_yyz*gBar_yzy*x36 + gBar_yyz*gBar_yzy*x60 - gBar_yyz*iBar_yzy - gBar_yyz*x10 - gBar_yyz*x11 - gBar_yyz*x9 - gBar_yyzz*x3 + gBar_yzx*gBar_yzy*iBar_xy*iBar_zz - gBar_yzx*iBar_xzy + 2.0*gBar_yzy*iBar_xz*phi_x + 2.0*gBar_yzy*iBar_yz*phi_y
                     + 2.0*gBar_yzy*iBar_zz*phi_z - gBar_yzy*x144 - gBar_yzy*x34 - gBar_yzz*iBar_zzy + (1.0/4.0)*gBar_zzy*iBar_zzy + (1.0/2.0)*iBar_xx*iBar_yy*x155 + (1.0/2.0)*iBar_xx*iBar_yy*x156 + (1.0/2.0)*iBar_xx*iBar_zz*x158 + (1.0/2.0)*iBar_xx*iBar_zz*x61 - iBar_xyy*x72
                     + (1.0/2.0)*iBar_yy*iBar_zz*x157 + (1.0/2.0)*iBar_yy*iBar_zz*x159 - x0*x146 - x100*x110 - x100*x114 - x107*x114 - x107*x95 - x146*x50 - x147*x15 - x147*x16 - x147*x19 - x148*x17 - x148*x52 - x149*x20 - x149*x54 - x150*x43 - x151*x39 - x151*x42 - x152*x44
                     - x153*x46 - x154*x48 - x155*x64 - x156*x64 - x157*x66 - x158*x65 - x159*x66 - x160*x68 - x160*x69 - x160*x70 - x17 + 12.0*x51 - x61*x65 - x95*x97;
         Rab (yz) := (1.0/4.0)*gBar_xxz*iBar_xxy + (1.0/2.0)*gBar_xyz*iBar_xyy + (1.0/4.0)*gBar_yyz*iBar_yyy - 1.0/4.0*gBar_zzz*iBar_zzy - iBar_xyz*x72 - iBar_yyz*x153 - iBar_yzz*x161 - x173;
         Rab (xz) := -1.0/4.0*gBar_xxx*iBar_xxz - gBar_xzz*x5 + (1.0/4.0)*gBar_yyx*iBar_yyz + (1.0/2.0)*gBar_yzx*iBar_yzz - gBar_yzz*x77 + (1.0/4.0)*gBar_zzx*iBar_zzz - iBar_zzx*x174 - x143;
         Rab (yz) := (1.0/4.0)*gBar_xxy*iBar_xxz + (1.0/2.0)*gBar_xzy*iBar_xzz - 1.0/4.0*gBar_yyy*iBar_yyz + (1.0/4.0)*gBar_zzy*iBar_zzz - iBar_xzy*x75 - iBar_yzy*x78 - iBar_zzy*x174 - x173;
         Rab (zz) := Gi_xz*gBar_xz + Gi_yz*gBar_yz + Gi_zz*gBar_zz + (1.0/4.0)*gBar_xxz*iBar_xxz + (1.0/2.0)*gBar_xyz*iBar_xyz - gBar_xzx*iBar_xxz + gBar_xzy*gBar_xzz*iBar_xx*iBar_yz - gBar_xzy*gBar_yzx*x23 + gBar_xzy*gBar_yzx*x55 + gBar_xzy*gBar_yzz*iBar_xy*iBar_yz + gBar_xzy*gBar_zzx*iBar_xy*iBar_xz
                     + gBar_xzy*gBar_zzy*iBar_xz*iBar_yy - gBar_xzy*iBar_xyz + gBar_xzz*gBar_yzx*iBar_xy*iBar_xz + gBar_xzz*gBar_yzz*iBar_xy*iBar_zz - gBar_xzz*gBar_yzz*iBar_xz*iBar_yz - gBar_xzz*gBar_zzx*x28 + gBar_xzz*gBar_zzx*x57 + gBar_xzz*gBar_zzy*iBar_xz*iBar_yz + 2.0*gBar_xzz*iBar_xx*phi_x
                     - gBar_xzz*iBar_xx*x137 + 2.0*gBar_xzz*iBar_xy*phi_y + 2.0*gBar_xzz*iBar_xz*phi_z - gBar_xzz*x134 + (1.0/4.0)*gBar_yyz*iBar_yyz + gBar_yzx*gBar_yzz*iBar_xz*iBar_yy + gBar_yzx*gBar_zzx*iBar_xx*iBar_yz + gBar_yzx*gBar_zzy*iBar_xy*iBar_yz - gBar_yzx*iBar_xyz - gBar_yzy*iBar_yyz
                     + gBar_yzz*gBar_zzx*iBar_xz*iBar_yz - gBar_yzz*gBar_zzx*x103 - gBar_yzz*gBar_zzy*x36 + gBar_yzz*gBar_zzy*x60 + 2.0*gBar_yzz*iBar_xy*phi_x + 2.0*gBar_yzz*iBar_yy*phi_y - gBar_yzz*iBar_yy*x33 + 2.0*gBar_yzz*iBar_yz*phi_z - gBar_yzz*x138 + gBar_zzx*gBar_zzy*iBar_xy*iBar_zz
                     - gBar_zzx*iBar_xx*x37 - gBar_zzx*iBar_xzz - gBar_zzx*x108 - gBar_zzx*x145 - gBar_zzx*x39 - gBar_zzx*x42 - gBar_zzxx*x1 - gBar_zzxy*iBar_xy - gBar_zzxz*iBar_xz - gBar_zzy*iBar_yy*x107 - gBar_zzy*iBar_yzz - gBar_zzy*x111 - 1.0/2.0*gBar_zzy*x46 - gBar_zzy*x6 - gBar_zzy*x7
                     - gBar_zzy*x8 - gBar_zzyy*x2 - gBar_zzyz*iBar_yz + gBar_zzz*iBar_xz*phi_x + gBar_zzz*iBar_yz*phi_y + gBar_zzz*iBar_zz*phi_z - 3.0/4.0*gBar_zzz*iBar_zzz - gBar_zzzz*x3 + (1.0/2.0)*iBar_xx*iBar_yy*x158 + (1.0/2.0)*iBar_xx*iBar_yy*x63 + (1.0/2.0)*iBar_xx*iBar_zz*x181
                     + (1.0/2.0)*iBar_xx*iBar_zz*x183 - iBar_xzz*x75 + (1.0/2.0)*iBar_yy*iBar_zz*x182 + (1.0/2.0)*iBar_yy*iBar_zz*x184 - iBar_yzz*x78 - x0*x175 - x132*x140 - x136*x25 - x15*x176 - x158*x64 - x16*x176 - x163*x44 - x17*x177 - x174*x48 - x175*x50 - x176*x19 - x177*x52
                     - x178*x20 - x178*x54 - x179*x43 - x180*x39 - x180*x42 - x181*x65 - x182*x66 - x183*x65 - x184*x66 - x185*x68 - x185*x69 - x185*x70 - x20 + 12.0*x53 - x63*x64;
      end set_3d_ricci;


      Procedure set_3d_dot_N is
      begin
         dot_N := 0.0;
      end set_3d_dot_N;

      Procedure set_3d_dot_Gi is
         x0,  x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,  x10, x11, x12, x13, x14, x15, x16,
         x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33,
         x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50 : Real;
      begin
         x0  := 2.0*N_x;
         x1  := 2.0*N_y;
         x2  := 2.0*N_z;
         x3  := 12.0*N;
         x4  := phi_x*x3;
         x5  := phi_y*x3;
         x6  := phi_z*x3;
         x7  := (4.0/3.0)*N;
         x8  := trK_x*x7;
         x9  := trK_y*x7;
         x10 := trK_z*x7;
         x11 := BBar_xx*N;
         x12 := gBar_xxx*x11;
         x13 := BBar_yy*N;
         x14 := gBar_yyy*x13;
         x15 := BBar_zz*N;
         x16 := gBar_zzz*x15;
         x17 := gBar_xxy*x11;
         x18 := gBar_xxz*x11;
         x19 := 2.0*x11;
         x20 := gBar_xyx*x19;
         x21 := gBar_xzx*x19;
         x22 := 2.0*N;
         x23 := BBar_xy*x22;
         x24 := gBar_xxy*x23;
         x25 := iBar_xz*x23;
         x26 := gBar_yyx*x23;
         x27 := BBar_xz*x22;
         x28 := gBar_xxz*x27;
         x29 := iBar_xy*x27;
         x30 := gBar_zzx*x27;
         x31 := 2.0*x13;
         x32 := gBar_xyy*x31;
         x33 := gBar_yyx*x13;
         x34 := gBar_yyz*x13;
         x35 := gBar_yzy*x31;
         x36 := BBar_yz*x22;
         x37 := iBar_xx*x36;
         x38 := iBar_xy*x36;
         x39 := gBar_zzy*x36;
         x40 := 2.0*x15;
         x41 := gBar_xzz*x40;
         x42 := gBar_yzz*x40;
         x43 := gBar_zzx*x15;
         x44 := gBar_zzy*x15;
         x45 := iBar_yz*x23;
         x46 := iBar_yy*x27;
         x47 := gBar_yyz*x36;
         x48 := iBar_zz*x23;
         x49 := iBar_yz*x27;
         x50 := iBar_xz*x36;
         dot_Gi (1) := -BBar_xx*x0 + BBar_xx*x4 - BBar_xy*x1 + BBar_xy*x5 - BBar_xz*x2 + BBar_xz*x6 - gBar_xyz*x25 + gBar_xyz*x29 + gBar_xyz*x37 + gBar_xzy*x25 - gBar_xzy*x29 + gBar_xzy*x37 + gBar_yyz*x38 + gBar_yzx*x25 + gBar_yzx*x29 - gBar_yzx*x37 + iBar_xx*x12 + iBar_xx*x24 + iBar_xx*x28
                     + iBar_xx*x32 - iBar_xx*x33 + iBar_xx*x41 - iBar_xx*x43 - iBar_xx*x8 + iBar_xy*x14 - iBar_xy*x17 + iBar_xy*x20 + iBar_xy*x26 + iBar_xy*x42 - iBar_xy*x44 - iBar_xy*x9 - iBar_xz*x10 + iBar_xz*x16 - iBar_xz*x18 + iBar_xz*x21 + iBar_xz*x30 - iBar_xz*x34 + iBar_xz*x35
                     + iBar_xz*x39;
         dot_Gi (2) := -BBar_xy*x0 + BBar_xy*x4 - BBar_yy*x1 + BBar_yy*x5 - BBar_yz*x2 + BBar_yz*x6 + gBar_xxz*x29 + gBar_xyz*x38 - gBar_xyz*x45 + gBar_xyz*x46 + gBar_xzy*x38 + gBar_xzy*x45 - gBar_xzy*x46 - gBar_yzx*x38 + gBar_yzx*x45 + gBar_yzx*x46 + iBar_xy*x12 + iBar_xy*x24 + iBar_xy*x32
                     - iBar_xy*x33 + iBar_xy*x41 - iBar_xy*x43 - iBar_xy*x8 + iBar_yy*x14 - iBar_yy*x17 + iBar_yy*x20 + iBar_yy*x26 + iBar_yy*x42 - iBar_yy*x44 + iBar_yy*x47 - iBar_yy*x9 - iBar_yz*x10 + iBar_yz*x16 - iBar_yz*x18 + iBar_yz*x21 + iBar_yz*x30 - iBar_yz*x34 + iBar_yz*x35
                     + iBar_yz*x39;
         dot_Gi (3) := -BBar_xz*x0 + BBar_xz*x4 - BBar_yz*x1 + BBar_yz*x5 - BBar_zz*x2 + BBar_zz*x6 + gBar_xxy*x25 - gBar_xyz*x48 + gBar_xyz*x49 + gBar_xyz*x50 + gBar_xzy*x48 - gBar_xzy*x49 + gBar_xzy*x50 + gBar_yzx*x48 + gBar_yzx*x49 - gBar_yzx*x50 + iBar_xz*x12 + iBar_xz*x28 + iBar_xz*x32
                     - iBar_xz*x33 + iBar_xz*x41 - iBar_xz*x43 - iBar_xz*x8 + iBar_yz*x14 - iBar_yz*x17 + iBar_yz*x20 + iBar_yz*x26 + iBar_yz*x42 - iBar_yz*x44 + iBar_yz*x47 - iBar_yz*x9 - iBar_zz*x10 + iBar_zz*x16 - iBar_zz*x18 + iBar_zz*x21 + iBar_zz*x30 - iBar_zz*x34 + iBar_zz*x35
                     + iBar_zz*x39;
      end set_3d_dot_Gi;

      Procedure set_3d_dot_phi is
      begin
         dot_phi := -1.0/6.0*N*trK;
      end set_3d_dot_phi;

      Procedure set_3d_dot_trK is
         x0, x1, x2, x3 : Real;
      begin
         x0 := 2.0*N;
         x1 := 2.0*iBar_xy;
         x2 := 2.0*iBar_xz;
         x3 := 2.0*iBar_yz;
         dot_trK := ABar_xx*BBar_xx*N + ABar_xy*BBar_xy*x0 + ABar_xz*BBar_xz*x0 + ABar_yy*BBar_yy*N + ABar_yz*BBar_yz*x0 + ABar_zz*BBar_zz*N + (1.0/3.0)*N*(trK ** 2) + (-2.0*N_x*iBar_xx*phi_x - N_x*phi_y*x1 - N_x*phi_z*x2 - N_x*(iBar_xxx + iBar_xyy + iBar_xzz) - N_xx*iBar_xx - N_xy*x1 - N_xz*x2
                  - 2.0*N_y*iBar_yy*phi_y - N_y*phi_x*x1 - N_y*phi_z*x3 - N_y*(iBar_xyx + iBar_yyy + iBar_yzz) - N_yy*iBar_yy - N_yz*x3 - 2.0*N_z*iBar_zz*phi_z - N_z*phi_x*x2 - N_z*phi_y*x3 - N_z*(iBar_xzx + iBar_yzy + iBar_zzz) - N_zz*iBar_zz)*exp(-4.0*phi);
      end set_3d_dot_trK;

      Procedure set_3d_dot_gBar is
         x0, x1, x2, x3 : Real;
      begin
         x0 := 2.0*N;
         x1 := -ABar_xy*x0;
         x2 := -ABar_xz*x0;
         x3 := -ABar_yz*x0;
         dot_gBar (xx) := -ABar_xx*x0;
         dot_gBar (xy) := x1;
         dot_gBar (xz) := x2;
         dot_gBar (xy) := x1;
         dot_gBar (yy) := -ABar_yy*x0;
         dot_gBar (yz) := x3;
         dot_gBar (xz) := x2;
         dot_gBar (yz) := x3;
         dot_gBar (zz) := -ABar_zz*x0;
      end set_3d_dot_gBar;

      Procedure set_3d_dot_ABar is
         x0,   x1,   x2,   x3,   x4,   x5,   x6,   x7,   x8,   x9,   x10,  x11,  x12,  x13,
         x14,  x15,  x16,  x17,  x18,  x19,  x20,  x21,  x22,  x23,  x24,  x25,  x26,  x27,
         x28,  x29,  x30,  x31,  x32,  x33,  x34,  x35,  x36,  x37,  x38,  x39,  x40,  x41,
         x42,  x43,  x44,  x45,  x46,  x47,  x48,  x49,  x50,  x51,  x52,  x53,  x54,  x55,
         x56,  x57,  x58,  x59,  x60,  x61,  x62,  x63,  x64,  x65,  x66,  x67,  x68,  x69,
         x70,  x71,  x72,  x73,  x74,  x75,  x76,  x77,  x78,  x79,  x80,  x81,  x82,  x83,
         x84,  x85,  x86,  x87,  x88,  x89,  x90,  x91,  x92,  x93,  x94,  x95,  x96,  x97,
         x98,  x99,  x100, x101, x102, x103, x104, x105, x106, x107, x108, x109, x110, x111,
         x112, x113 : Real;
      begin
         x0   := ABar_xx*ABar_xy;
         x1   := 4.0*N;
         x2   := iBar_xy*x1;
         x3   := ABar_xz*x1;
         x4   := ABar_xx*iBar_xz;
         x5   := ABar_xy*iBar_yz;
         x6   := 2.0*N;
         x7   := iBar_xx*x6;
         x8   := (ABar_xy ** 2);
         x9   := x6*x8;
         x10  := (ABar_xz ** 2);
         x11  := x10*x6;
         x12  := exp(-4.0*phi);
         x13  := N_x*iBar_xy;
         x14  := (1.0/2.0)*gBar_xxy;
         x15  := N_x*iBar_xz;
         x16  := (1.0/2.0)*gBar_xxz;
         x17  := N_y*iBar_yy;
         x18  := N_y*iBar_yz;
         x19  := N_z*iBar_yz;
         x20  := N_z*iBar_zz;
         x21  := (4.0/3.0)*gBar_xx;
         x22  := N_x*iBar_xx*phi_x;
         x23  := phi_y*x21;
         x24  := phi_z*x21;
         x25  := N_y*iBar_xy;
         x26  := phi_x*x21;
         x27  := N_z*iBar_xz;
         x28  := (1.0/3.0)*gBar_xx;
         x29  := Rxx*iBar_xx*N;
         x30  := (2.0/3.0)*gBar_xx;
         x31  := Rxy*N;
         x32  := iBar_xy*x31;
         x33  := Rxz*N;
         x34  := iBar_xz*x33;
         x35  := Ryy*iBar_yy*N;
         x36  := Ryz*N;
         x37  := iBar_yz*x36;
         x38  := Rzz*iBar_zz*N;
         x39  := iBar_xxx + iBar_xyy + iBar_xzz;
         x40  := iBar_xyx + iBar_yyy + iBar_yzz;
         x41  := iBar_xzx + iBar_yzy + iBar_zzz;
         x42  := ABar_yy*x6;
         x43  := ABar_xx*iBar_xy;
         x44  := ABar_yz*x6;
         x45  := ABar_xy*iBar_xz;
         x46  := ABar_xz*x6;
         x47  := iBar_yy*x6;
         x48  := ABar_xy*x47;
         x49  := ABar_xz*iBar_yz;
         x50  := iBar_zz*x6;
         x51  := ABar_xz*x50;
         x52  := 2.0*N_x;
         x53  := 2.0*phi_x;
         x54  := N_x*iBar_xx;
         x55  := (1.0/2.0)*x15;
         x56  := (1.0/2.0)*gBar_yyx;
         x57  := N_xx*iBar_xx;
         x58  := (1.0/3.0)*gBar_xy;
         x59  := N_xy*iBar_xy;
         x60  := (2.0/3.0)*gBar_xy;
         x61  := N_xz*iBar_xz;
         x62  := (1.0/2.0)*x18;
         x63  := N_yy*iBar_yy;
         x64  := N_yz*iBar_yz;
         x65  := (1.0/2.0)*x20;
         x66  := N_zz*iBar_zz;
         x67  := (4.0/3.0)*gBar_xy;
         x68  := phi_y*x67;
         x69  := phi_z*x67;
         x70  := phi_x*x67;
         x71  := N_x*x39;
         x72  := N_y*x40;
         x73  := N_z*x41;
         x74  := ABar_xy*N*trK - ABar_yy*x48 - ABar_yz*x51 - iBar_xy*x9 - x0*x7 + x12*(-N_xy + N_y*x53 - gBar_xyz*x55 - gBar_xyz*x62 - gBar_xyz*x65 + gBar_xzy*x55 + gBar_xzy*x62 + gBar_xzy*x65 + gBar_yzx*x55 + gBar_yzx*x62 + gBar_yzx*x65 + phi_y*x52 + x13*x56 - x13*x68 + x14*x25
               + x14*x27 + x14*x54 - x15*x69 + x17*x56 - x17*x68 - x18*x69 + x19*x56 - x19*x68 - x20*x69 - x22*x67 - x25*x70 - x27*x70 - x29*x58 + x31 - x32*x60 - x34*x60 - x35*x58 - x37*x60 - x38*x58 + x57*x58 + x58*x63 + x58*x66 + x58*x71 + x58*x72 + x58*x73 + x59*x60
               + x60*x61 + x60*x64) - x4*x44 - x42*x43 - x42*x49 - x44*x5 - x45*x46;
         x75  := ABar_xz*x7;
         x76  := ABar_zz*x6;
         x77  := ABar_xy*iBar_xy;
         x78  := (1.0/2.0)*x13;
         x79  := (1.0/3.0)*gBar_xz;
         x80  := (2.0/3.0)*gBar_xz;
         x81  := (1.0/2.0)*x17;
         x82  := (1.0/2.0)*x19;
         x83  := (4.0/3.0)*gBar_xz;
         x84  := phi_y*x83;
         x85  := phi_z*x83;
         x86  := phi_x*x83;
         x87  := -ABar_xx*x75 + ABar_xz*N*trK - ABar_yz*x48 - ABar_zz*x51 - iBar_xz*x11 + x12*(-N_xz + N_z*x53 + gBar_xyz*x78 + gBar_xyz*x81 + gBar_xyz*x82 - gBar_xzy*x78 - gBar_xzy*x81 - gBar_xzy*x82 + gBar_yzx*x78 + gBar_yzx*x81 + gBar_yzx*x82 + gBar_zzx*x55 + gBar_zzx*x62 + gBar_zzx*x65
               + phi_z*x52 - x13*x84 - x15*x85 + x16*x25 + x16*x27 + x16*x54 - x17*x84 - x18*x85 - x19*x84 - x20*x85 - x22*x83 - x25*x86 - x27*x86 - x29*x79 - x32*x80 + x33 - x34*x80 - x35*x79 - x37*x80 - x38*x79 + x57*x79 + x59*x80 + x61*x80 + x63*x79 + x64*x80
               + x66*x79 + x71*x79 + x72*x79 + x73*x79) - x4*x76 - x43*x44 - x44*x49 - x46*x77 - x5*x76;
         x88  := ABar_yz*x1;
         x89  := iBar_yz*x88;
         x90  := (ABar_yz ** 2);
         x91  := x6*x90;
         x92  := (4.0/3.0)*gBar_yy;
         x93  := phi_y*x92;
         x94  := phi_z*x92;
         x95  := phi_x*x92;
         x96  := (1.0/3.0)*gBar_yy;
         x97  := (2.0/3.0)*gBar_yy;
         x98  := (1.0/2.0)*x54;
         x99  := (1.0/3.0)*gBar_yz;
         x100 := (2.0/3.0)*gBar_yz;
         x101 := (1.0/2.0)*x25;
         x102 := (1.0/2.0)*x27;
         x103 := (4.0/3.0)*gBar_yz;
         x104 := phi_y*x103;
         x105 := phi_z*x103;
         x106 := phi_x*x103;
         x107 := -ABar_xy*x75 - ABar_xz*iBar_xy*x42 - ABar_xz*iBar_xz*x44 - ABar_yy*ABar_yz*x47 - ABar_yz*ABar_zz*x50 + ABar_yz*N*trK - ABar_zz*iBar_yz*x42 - iBar_yz*x91 + x12*(2.0*N_y*phi_z - N_yz + 2.0*N_z*phi_y + gBar_xyz*x101 + gBar_xyz*x102 + gBar_xyz*x98 + gBar_xzy*x101 + gBar_xzy*x102
               + gBar_xzy*x98 + gBar_yyz*x78 + gBar_yyz*x81 + gBar_yyz*x82 - gBar_yzx*x101 - gBar_yzx*x102 - gBar_yzx*x98 + gBar_zzy*x55 + gBar_zzy*x62 + gBar_zzy*x65 - x100*x32 - x100*x34 - x100*x37 + x100*x59 + x100*x61 + x100*x64 - x103*x22 - x104*x13 - x104*x17 - x104*x19
               - x105*x15 - x105*x18 - x105*x20 - x106*x25 - x106*x27 - x29*x99 - x35*x99 + x36 - x38*x99 + x57*x99 + x63*x99 + x66*x99 + x71*x99 + x72*x99 + x73*x99) - x44*x77 - x45*x76;
         x108 := (4.0/3.0)*gBar_zz;
         x109 := phi_y*x108;
         x110 := phi_z*x108;
         x111 := phi_x*x108;
         x112 := (1.0/3.0)*gBar_zz;
         x113 := (2.0/3.0)*gBar_zz;
         dot_ABar (xx) := -(ABar_xx ** 2)*x7 + ABar_xx*N*trK - iBar_yy*x9 - iBar_zz*x11 - x0*x2 + x12*((1.0/3.0)*N_x*gBar_xx*x39 + (1.0/2.0)*N_x*gBar_xxx*iBar_xx + N_x*gBar_xyx*iBar_xy + N_x*gBar_xzx*iBar_xz + 4.0*N_x*phi_x + (1.0/3.0)*N_xx*gBar_xx*iBar_xx - N_xx + (2.0/3.0)*N_xy*gBar_xx*iBar_xy
                          + (2.0/3.0)*N_xz*gBar_xx*iBar_xz + (1.0/3.0)*N_y*gBar_xx*x40 + (1.0/2.0)*N_y*gBar_xxx*iBar_xy + N_y*gBar_xyx*iBar_yy + N_y*gBar_xzx*iBar_yz + (1.0/3.0)*N_yy*gBar_xx*iBar_yy + (2.0/3.0)*N_yz*gBar_xx*iBar_yz + (1.0/3.0)*N_z*gBar_xx*x41 + (1.0/2.0)*N_z*gBar_xxx*iBar_xz
                          + N_z*gBar_xyx*iBar_yz + N_z*gBar_xzx*iBar_zz + (1.0/3.0)*N_zz*gBar_xx*iBar_zz + Rxx*N - x13*x14 - x13*x23 - x14*x17 - x14*x19 - x15*x16 - x15*x24 - x16*x18 - x16*x20 - x17*x23 - x18*x24 - x19*x23 - x20*x24 - x21*x22 - x25*x26 - x26*x27 - x28*x29 - x28*x35
                          - x28*x38 - x30*x32 - x30*x34 - x30*x37) - x3*x4 - x3*x5;
         dot_ABar (xy) := x74;
         dot_ABar (xz) := x87;
         dot_ABar (xy) := x74;
         dot_ABar (yy) := -ABar_xy*ABar_yy*x2 - (ABar_yy ** 2)*x47 + ABar_yy*N*trK - ABar_yy*x89 - iBar_zz*x91 + x12*(N_x*gBar_xyy*iBar_xx + (1.0/3.0)*N_x*gBar_yy*x39 + (1.0/2.0)*N_x*gBar_yyy*iBar_xy + N_x*gBar_yzy*iBar_xz + (1.0/3.0)*N_xx*gBar_yy*iBar_xx + (2.0/3.0)*N_xy*gBar_yy*iBar_xy
                          + (2.0/3.0)*N_xz*gBar_yy*iBar_xz + N_y*gBar_xyy*iBar_xy + (1.0/3.0)*N_y*gBar_yy*x40 + (1.0/2.0)*N_y*gBar_yyy*iBar_yy + N_y*gBar_yzy*iBar_yz + 4.0*N_y*phi_y + (1.0/3.0)*N_yy*gBar_yy*iBar_yy - N_yy + (2.0/3.0)*N_yz*gBar_yy*iBar_yz + N_z*gBar_xyy*iBar_xz + (1.0/3.0)*N_z*gBar_yy*x41
                          + (1.0/2.0)*N_z*gBar_yyy*iBar_yz + N_z*gBar_yzy*iBar_zz + (1.0/3.0)*N_zz*gBar_yy*iBar_zz + Ryy*N - gBar_yyz*x55 - gBar_yyz*x62 - gBar_yyz*x65 - x13*x93 - x15*x94 - x17*x93 - x18*x94 - x19*x93 - x20*x94 - x22*x92 - x25*x56 - x25*x95 - x27*x56 - x27*x95 - x29*x96
                          - x32*x97 - x34*x97 - x35*x96 - x37*x97 - x38*x96 - x54*x56) - x45*x88 - x7*x8;
         dot_ABar (yz) := x107;
         dot_ABar (xz) := x87;
         dot_ABar (yz) := x107;
         dot_ABar (zz) := -ABar_xz*ABar_yz*x2 - (ABar_zz ** 2)*x50 - ABar_zz*iBar_xz*x3 + ABar_zz*N*trK - ABar_zz*x89 - x10*x7 + x12*(N_x*gBar_xzz*iBar_xx + N_x*gBar_yzz*iBar_xy + (1.0/3.0)*N_x*gBar_zz*x39 + (1.0/2.0)*N_x*gBar_zzz*iBar_xz + (1.0/3.0)*N_xx*gBar_zz*iBar_xx
                          + (2.0/3.0)*N_xy*gBar_zz*iBar_xy + (2.0/3.0)*N_xz*gBar_zz*iBar_xz + N_y*gBar_xzz*iBar_xy + N_y*gBar_yzz*iBar_yy + (1.0/3.0)*N_y*gBar_zz*x40 + (1.0/2.0)*N_y*gBar_zzz*iBar_yz + (1.0/3.0)*N_yy*gBar_zz*iBar_yy + (2.0/3.0)*N_yz*gBar_zz*iBar_yz + N_z*gBar_xzz*iBar_xz + N_z*gBar_yzz*iBar_yz
                          + (1.0/3.0)*N_z*gBar_zz*x41 + (1.0/2.0)*N_z*gBar_zzz*iBar_zz + 4.0*N_z*phi_z + (1.0/3.0)*N_zz*gBar_zz*iBar_zz - N_zz + Rzz*N - gBar_zzx*x101 - gBar_zzx*x102 - gBar_zzx*x98 - gBar_zzy*x78 - gBar_zzy*x81 - gBar_zzy*x82 - x108*x22 - x109*x13 - x109*x17 - x109*x19 - x110*x15
                          - x110*x18 - x110*x20 - x111*x25 - x111*x27 - x112*x29 - x112*x35 - x112*x38 - x113*x32 - x113*x34 - x113*x37) - x47*x90;
      end set_3d_dot_ABar;

      procedure set_data is
         i, j, k : Integer;
         x, y, z : Real;
      begin

         i := point.i;
         j := point.j;
         k := point.k;

         x := point.x;
         y := point.y;
         z := point.z;

         N    := BSSNBase.N    (i,j,k);
         Gi   := BSSNBase.Gi   (i,j,k);
         phi  := BSSNBase.phi  (i,j,k);
         trK  := BSSNBase.trK  (i,j,k);
         gBar := BSSNBase.gBar (i,j,k);
         ABar := BSSNBase.ABar (i,j,k);

         iBar := symm_inverse (gBar);
         BBar := symm_raise_indices (ABar, iBar);

         -- second order centred finite differences for the lapse

         declare
            N : LapseGridArray renames BSSNBase.N;
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

         -- second order centred finite differences for the conformal factor phi

         declare
            phi : ConFactGridArray renames BSSNBase.phi;
         begin

            d1phi (1) := (phi (i+1,j,k) - phi (i-1,j,k)) / (two_dx);
            d1phi (2) := (phi (i,j+1,k) - phi (i,j-1,k)) / (two_dy);
            d1phi (3) := (phi (i,j,k+1) - phi (i,j,k-1)) / (two_dz);

            d2phi (xx) := (phi (i+1,j,k) - 2.0*phi (i,j,k) + phi (i-1,j,k)) / (dxdx);
            d2phi (yy) := (phi (i,j+1,k) - 2.0*phi (i,j,k) + phi (i,j-1,k)) / (dydy);
            d2phi (zz) := (phi (i,j,k+1) - 2.0*phi (i,j,k) + phi (i,j,k-1)) / (dzdz);

            d2phi (xy) := (phi (i+1,j+1,k) + phi (i-1,j-1,k) - phi (i+1,j-1,k) - phi (i-1,j+1,k)) / (four_dxdy);
            d2phi (xz) := (phi (i+1,j,k+1) + phi (i-1,j,k-1) - phi (i+1,j,k-1) - phi (i-1,j,k+1)) / (four_dxdz);
            d2phi (yz) := (phi (i,j+1,k+1) + phi (i,j-1,k-1) - phi (i,j+1,k-1) - phi (i,j-1,k+1)) / (four_dydz);

         end;

         -- second order centred finite differences for trK

         declare
            trK : TraceKGridArray renames BSSNBase.trK;
         begin

            d1trK (1) := (trK (i+1,j,k) - trK (i-1,j,k)) / (two_dx);
            d1trK (2) := (trK (i,j+1,k) - trK (i,j-1,k)) / (two_dy);
            d1trK (3) := (trK (i,j,k+1) - trK (i,j,k-1)) / (two_dz);

         end;

         -- second order centred finite differences for trK

         declare
            Gi : GammaGridArray renames BSSNBase.Gi;
         begin

            d1Gi (1) := (Gi (i+1,j,k) - Gi (i-1,j,k)) / (two_dx);
            d1Gi (2) := (Gi (i,j+1,k) - Gi (i,j-1,k)) / (two_dy);
            d1Gi (3) := (Gi (i,j,k+1) - Gi (i,j,k-1)) / (two_dz);

         end;

         -- second order centred finite differences for the 3-metric

         declare
            gBar : MetricGridArray renames BSSNBase.gBar;
         begin

            d1gBar (1) := (gBar (i+1,j,k) - gBar (i-1,j,k)) / (two_dx);
            d1gBar (2) := (gBar (i,j+1,k) - gBar (i,j-1,k)) / (two_dy);
            d1gBar (3) := (gBar (i,j,k+1) - gBar (i,j,k-1)) / (two_dz);

            d2gBar (xx) := (gBar (i+1,j,k) - 2.0*gBar (i,j,k) + gBar (i-1,j,k)) / (dxdx);
            d2gBar (yy) := (gBar (i,j+1,k) - 2.0*gBar (i,j,k) + gBar (i,j-1,k)) / (dydy);
            d2gBar (zz) := (gBar (i,j,k+1) - 2.0*gBar (i,j,k) + gBar (i,j,k-1)) / (dzdz);

            d2gBar (xy) := (gBar (i+1,j+1,k) + gBar (i-1,j-1,k) - gBar (i+1,j-1,k) - gBar (i-1,j+1,k)) / (four_dxdy);
            d2gBar (xz) := (gBar (i+1,j,k+1) + gBar (i-1,j,k-1) - gBar (i+1,j,k-1) - gBar (i-1,j,k+1)) / (four_dxdz);
            d2gBar (yz) := (gBar (i,j+1,k+1) + gBar (i,j-1,k-1) - gBar (i,j+1,k-1) - gBar (i,j-1,k+1)) / (four_dydz);

         end;

         d1iBar (1) := - symm_raise_indices (d1gBar (1), iBar);
         d1iBar (2) := - symm_raise_indices (d1gBar (2), iBar);
         d1iBar (3) := - symm_raise_indices (d1gBar (3), iBar);

         set_3d_ricci;

      end set_data;

   begin

      set_data;

      set_3d_dot_N;
      set_3d_dot_Gi;
      set_3d_dot_phi;
      set_3d_dot_trK;
      set_3d_dot_gBar;
      set_3d_dot_ABar;

      declare

         i, j, k : Integer;

      begin

         i := point.i;
         j := point.j;
         k := point.k;

         BSSNBase.dot_N    (i,j,k) := dot_N;
         BSSNBase.dot_Gi   (i,j,k) := dot_Gi;
         BSSNBase.dot_phi  (i,j,k) := dot_phi;
         BSSNBase.dot_trK  (i,j,k) := dot_trK;
         BSSNBase.dot_gBar (i,j,k) := dot_gBar;
         BSSNBase.dot_ABar (i,j,k) := dot_ABar;

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

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (i,j,2);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (i,j,2);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (i,j,2);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (i,j,2);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (i,j,2);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (i,j,2);

      end loop;

      for a in beg_south .. end_south loop

         b := south_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;  -- equals 1

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (i,j,num_z-1);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (i,j,num_z-1);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (i,j,num_z-1);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (i,j,num_z-1);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (i,j,num_z-1);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (i,j,num_z-1);

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

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (i,2,k);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (i,2,k);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (i,2,k);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (i,2,k);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (i,2,k);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (i,2,k);

      end loop;

      for a in beg_west .. end_west loop

         b := west_bndry (a);

         i := grid_point_list (b).i;
         j := grid_point_list (b).j;  -- equals 1
         k := grid_point_list (b).k;

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (i,num_y-1,k);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (i,num_y-1,k);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (i,num_y-1,k);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (i,num_y-1,k);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (i,num_y-1,k);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (i,num_y-1,k);

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

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (2,j,k);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (2,j,k);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (2,j,k);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (2,j,k);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (2,j,k);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (2,j,k);

      end loop;

      for a in beg_back .. end_back loop

         b := back_bndry (a);

         i := grid_point_list (b).i;  -- equals 1
         j := grid_point_list (b).j;
         k := grid_point_list (b).k;

            BSSNBase.dot_N (i,j,k) :=     BSSNBase.dot_N (num_x-1,j,k);
           BSSNBase.dot_Gi (i,j,k) :=    BSSNBase.dot_Gi (num_x-1,j,k);
          BSSNBase.dot_phi (i,j,k) :=   BSSNBase.dot_phi (num_x-1,j,k);
          BSSNBase.dot_trK (i,j,k) :=   BSSNBase.dot_trK (num_x-1,j,k);
         BSSNBase.dot_gBar (i,j,k) :=  BSSNBase.dot_gBar (num_x-1,j,k);
         BSSNBase.dot_ABar (i,j,k) :=  BSSNBase.dot_ABar (num_x-1,j,k);

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

end BSSNBase.Time_Derivs;
