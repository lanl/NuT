// T. M. Kelley (c) 2011 LANS LLC

#include "Mesh.hh"
#include "expect.hh"
#include "gtest/gtest.h"
#include "types.hh"
#include <algorithm>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>  // accumulate
#include <vector>

using nut::soft_equiv;
using test_aux::expect;
using test_aux::soft_expect;

using nut::cell_t;
using nut::geom_t;
using Sp1D = nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor>;

/* Make a mesh that is used in several places */
Sp1D
make_mesh()
{
  size_t constexpr n_cells(3);
  size_t constexpr n_bdys(n_cells + 1);
  geom_t const bdys_in[n_bdys] = {0.0, 1.0, 2.0, 35.0};
  nut::bdy_types::descriptor const bdy_ts[n_bdys] = {
      nut::bdy_types::REFLECTIVE, nut::bdy_types::CELL, nut::bdy_types::CELL,
      nut::bdy_types::VACUUM};
  Sp1D::vb bdys(&bdys_in[0], &bdys_in[n_bdys]);
  Sp1D::vbd bdy_types(&bdy_ts[0], &bdy_ts[n_bdys]);

  Sp1D mesh(bdys, bdy_types);
  return mesh;
};

TEST(nut_mesh, Sphere_1D_inst_init)
{
  Sp1D mesh{make_mesh()};
  EXPECT_TRUE(true);
  return;
}

TEST(nut_test, Sphere_1D_num_cells)
{
  Sp1D mesh{make_mesh()};
  EXPECT_EQ(mesh.num_cells(), 3u);
  return;
}  // test_2

TEST(nut_test, Sphere_1D_volume)
{
  Sp1D mesh{make_mesh()};

  geom_t const vols_exp[] = {4.1887902047863905, 29.321531433504735,
                             179560.86970857819};
  Sp1D::vb vols(mesh.num_cells());
  vols[0] = mesh.volume(1);
  vols[1] = mesh.volume(2);
  vols[2] = mesh.volume(3);

  bool passed = std::equal(vols.begin(), vols.end(), &vols_exp[0]);
  EXPECT_TRUE(passed);
  if(!passed) {
    std::cout << std::setprecision(16) << std::scientific << "volumes: ";
    std::copy(vols.begin(), vols.end(),
              std::ostream_iterator<geom_t>(std::cout, ","));
    std::cout << std::endl << "expected: ";
    std::copy(&vols_exp[0], &vols_exp[mesh.num_cells()],
              std::ostream_iterator<geom_t>(std::cout, ","));
    std::cout << std::endl;
  }
  return;
}  // test_3

TEST(nut_test, Sphere_1D_cell_across)
{
  Sp1D mesh{make_mesh()};
  Sp1D::Vector unused{0.0};

  cell_t const null{mesh.null_cell()};
  cell_t const cells_a_exp[] = {null, 2, 1, 3, 2, null};
  std::vector<cell_t> cells_across(mesh.num_cells() * 2);
  cells_across[0] = mesh.cell_across(1, 0, unused);
  cells_across[1] = mesh.cell_across(1, 1, unused);
  cells_across[2] = mesh.cell_across(2, 0, unused);
  cells_across[3] = mesh.cell_across(2, 1, unused);
  cells_across[4] = mesh.cell_across(3, 0, unused);
  cells_across[5] = mesh.cell_across(3, 1, unused);

  bool passed =
      std::equal(cells_across.begin(), cells_across.end(), &cells_a_exp[0]);
  EXPECT_TRUE(passed);
  if(!passed) {
    std::cout << std::setprecision(16) << std::scientific << "cells across:";
    std::copy(cells_across.begin(), cells_across.end(),
              std::ostream_iterator<cell_t>(std::cout, ","));
    std::cout << std::endl << "expected: ";
    std::copy(&cells_a_exp[0], &cells_a_exp[mesh.num_cells() * 2],
              std::ostream_iterator<cell_t>(std::cout, ","));
    std::cout << std::endl;
  }
  return;
}  // test_4

TEST(nut_test, Sphere_1D_distance_to_boundary_theta_0_cell_1)
{
  Sp1D mesh{make_mesh()};

  Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(0.5, 1.0, 1);

  bool passed1 = soft_expect(d_n_face.d, 0.5, "distance");
  bool passed2 = expect(d_n_face.face, 1u, "face");
  EXPECT_TRUE(passed1);
  EXPECT_TRUE(passed2);
  return;
}  // test_5

TEST(nut_test, Sphere_1D_distance_to_boundary_theta_pi_cell_1)
{
  Sp1D mesh{make_mesh()};

  Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(0.5, -1.0, 1);

  bool passed1 = soft_expect(d_n_face.d, 1.5, "distance");
  bool passed2 = expect(d_n_face.face, 1u, "face");
  EXPECT_TRUE(passed1);
  EXPECT_TRUE(passed2);
  return;
}  // test_6

TEST(nut_test, Sphere_1D_distance_to_boundary_theta_pi_cell_2)
{
  Sp1D mesh{make_mesh()};

  Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(1.5, -1.0, 2);

  bool passed1 = soft_expect(d_n_face.d, 0.5, "distance");
  bool passed2 = expect(d_n_face.face, 0u, "face");
  EXPECT_TRUE(passed1);
  EXPECT_TRUE(passed2);
  return;
}  // test_7

TEST(nut_test, Sphere_1D_distance_to_boundary_cell_2_grazing_inner_sphere)
{
  Sp1D mesh{make_mesh()};

  cell_t const cell = 2;

  geom_t const x = 2;
  geom_t const omega = -std::sqrt(3) / 2;
  geom_t const d_exp = std::sqrt(3);
  // This test has limited tolerance because of cancellation in
  // the code when computing the determinant.  The determinant is
  // 1.0*1.33333333333 - 4.0*0.33333333333. The exact value is
  // 1*4/3-4*1/3 = 0.  The floating point version fails pretty badly.
  geom_t const tol = 3.0e-8;

  Sp1D::d_to_b_t d_n_face = mesh.distance_to_bdy(x, omega, cell);

  bool passed1 = soft_expect(d_n_face.d, d_exp, "distance", tol);
  bool passed2 = expect(d_n_face.face, 0u, "face");
  EXPECT_TRUE(passed1);
  EXPECT_TRUE(passed2);
  return;
}  // test_8

namespace {

extern std::string const spherical_d_to_b_tests_cxx;

/* Read problems that were generated in Mathematica. Currently hard-coded below,
 * used to be read from file, but the configuration was always annoying. */
template <typename vb_t, typename vc_t>
bool
read_input_t9(vb_t & xs,
              vb_t & omegas,
              vb_t & r_los,
              vb_t & r_his,
              vb_t & d_exps,
              vc_t & f_exps)
{
  typedef typename vb_t::value_type geom_t;
  typedef typename vc_t::value_type cell_t;

  // std::string const fname("spherical-d-to-b-tests-cxx.txt");
  // std::ifstream inf(fname.c_str());
  // if(!inf.good()) {
  //   std::cerr << "mesh test_9 cannot read file " << fname << std::endl;
  //   std::cerr << "Are you running in the directory in which "
  //             << "the tests were built?" << std::endl;
  //   return false;
  // }

  std::stringstream instr(spherical_d_to_b_tests_cxx);

  uint32_t count(0);
  while(1) {
    geom_t x, omega, r_lo, r_hi, d_exp;
    cell_t f_exp;
    instr >> x >> omega >> r_lo >> r_hi >> d_exp >> f_exp;
    if(instr.eof()) break;
    // inf >> x >> omega >> r_lo >> r_hi >> d_exp >> f_exp;
    // if(inf.eof()) break;
    xs.push_back(x);
    omegas.push_back(omega);
    r_los.push_back(r_lo);
    r_his.push_back(r_hi);
    d_exps.push_back(d_exp);
    f_exps.push_back(f_exp - 1);
    count++;
  }
  std::cout << "read " << count << " input lines, will run that "
            << " many problems." << std::endl;
  // inf.close();
  return true;
}  // read_input
}  // namespace

TEST(
    nut_test,
    Sphere_1D_distance_to_boundary_impl_100_random_tests_created_with_Mathematica)
{
  Sp1D::vb xs, omegas, r_los, r_his, d_exps;
  std::vector<cell_t> f_exps;

  if(!read_input_t9(xs, omegas, r_los, r_his, d_exps, f_exps)) {
    std::cerr << "Failed to read input, this test SKIPPED" << std::endl;
    EXPECT_TRUE(false);
    return;
  }
  // This test has limited tolerance because of cancellation in
  // the code when computing the determinant.  The determinant is
  // 1.0*1.33333333333 - 4.0*0.33333333333. The exact value is
  // 1*4/3-4*1/3 = 0.  The floating point version fails pretty badly.
  geom_t const tol = 3.0e-8;
  test_aux::soft_eq_bound_tol<geom_t> s_eq(tol);
  Sp1D::vb ds(xs.size());
  std::vector<cell_t> fs(xs.size());
  for(size_t i = 0; i < xs.size(); ++i) {
    geom_t x = xs[i], o = omegas[i], rl = r_los[i], rh = r_his[i];
    Sp1D::d_to_b_t dnf = Sp1D::dist_to_bdy_impl(x, o, rl, rh, 0);
    ds[i] = (dnf.d);
    fs[i] = (dnf.face);
  }

  bool ds_passed = std::equal(ds.begin(), ds.end(), d_exps.begin(), s_eq);
  bool fs_passed = std::equal(fs.begin(), fs.end(), f_exps.begin());
  EXPECT_TRUE(ds_passed);
  EXPECT_TRUE(fs_passed);
  return;
}  // test_9

/* This test uses a different mesh. */
TEST(nut_test, Sphere_1D_sample_position_in_cell)
{
  size_t const n_cells(3);
  size_t const n_bdys(n_cells + 1);

  geom_t const bdys_in[n_bdys] = {0.0, 1.0, 2.0, 3.0};
  nut::bdy_types::descriptor const bdy_ts[n_bdys] = {
      nut::bdy_types::REFLECTIVE, nut::bdy_types::CELL, nut::bdy_types::CELL,
      nut::bdy_types::VACUUM};

  Sp1D::vb bdys(&bdys_in[0], &bdys_in[n_bdys]);
  Sp1D::vbd bdy_types(&bdy_ts[0], &bdy_ts[n_bdys]);

  Sp1D mesh(bdys, bdy_types);

  geom_t rns[] = {0.1};
  bool const silent{true};
  nut::Buffer_RNG<geom_t> rng(rns, 1, silent);

  geom_t position = mesh.sample_position(rng, 2);

  EXPECT_TRUE(soft_expect(position, 1.193483191927337, "sampled position"));
  return;
}  // test_10

namespace {
std::string const spherical_d_to_b_tests_cxx =
    "874.2153502626182 -0.03033663375173612 777.8031650087146 "
    "909.9127641144195 280.2767555431319 2\n"
    "739.5305386622146 -0.14307028584740067 709.965695760133 "
    "822.5307441896388 "
    "481.09949362525634 2\n"
    "359.6473848610964 0.17305031373270552 352.7618028981133 "
    "420.773585284999 "
    "164.86998406430524 2\n"
    "886.8311587844871 0.46121737024899545 839.1333934311679 "
    "1009.6636244204961 223.63126671365598 2\n"
    "568.3620502786575 -0.23930435259185412 527.8754766920204 "
    "704.1833666255529 573.4334199811067 2\n"
    "267.8464032914738 -0.547598222447438 224.80545606495684 "
    "301.60979199467147 129.10702781779267 1\n"
    "541.7900662787413 -0.952109897533806 487.72544704666416 "
    "598.1365295343455 "
    "57.11245315000766 1\n"
    "873.2379473868693 -0.6399762815968693 850.9067631131463 "
    "920.3726005914009 "
    "35.5802404125457 1\n"
    "297.33210065342246 -0.10290097025164746 283.11763768977244 "
    "331.1071262191506 179.46366142355302 2\n"
    "993.2972853077104 0.9596901736558299 948.8030613477008 "
    "994.3159170044139 "
    "1.061370613355976 2\n"
    "369.92912404107807 0.2667894503828223 356.5580084999324 "
    "478.6283247198687 "
    "220.64654004462236 2\n"
    "560.6447247186366 -0.5794229416159142 553.0717938206678 "
    "621.6786761089912 "
    "13.251804754872701 1\n"
    "899.7449706899623 -0.6744542760147199 833.6139841099259 "
    "925.2463822723373 "
    "103.22765304850122 1\n"
    "576.1085675798107 -0.2303599132833969 557.7034432814494 "
    "616.2125998943784 "
    "388.5041225518305 2\n"
    "728.697162194306 0.5946977255548012 712.1618249308217 "
    "730.2210779648419 "
    "2.557636736625827 2\n"
    "639.7883067592547 -0.7245540535212744 611.7384629372441 "
    "793.3535709388448 "
    "39.55193525321076 1\n"
    "419.4392097469225 0.43555160581358754 396.49028893339346 "
    "534.9313487866261 196.25451504984014 2\n"
    "593.9811855239203 -0.13704059445256567 465.86805398000547 "
    "598.0618489334039 188.59193491770094 2\n"
    "929.140163151597 0.457691478997134 887.3484743460783 "
    "942.1411018672693 "
    "27.701925975165096 2\n"
    "577.0762789223644 -0.13867542775000397 511.16229872520694 "
    "611.444383828537 297.40019771439046 2\n"
    "496.86431016881613 0.6486929958418699 487.43817216985803 "
    "625.0035367297136 175.32256876589003 2\n"
    "563.1405854998501 0.4441524176205376 477.96619403371324 "
    "592.9261281919214 "
    "61.318690416562546 2\n"
    "914.9358139255389 -0.9934030852025857 911.0413480167033 "
    "918.7646723739188 "
    "3.920439642179691 1\n"
    "345.25157808113113 0.9185813389608444 330.7245638314844 "
    "449.5085656549189 "
    "111.15511985311633 2\n"
    "633.272442356884 0.09213104265044159 598.3387203712086 "
    "755.6482330949289 "
    "358.0383389378323 2\n"
    "203.23283925176855 -0.3555650052542947 125.90986474959595 "
    "223.77135081892595 190.54983617251912 2\n"
    "827.127754638394 -0.23315253751447873 758.8398590082065 "
    "955.881734131246 "
    "709.3345991975053 2\n"
    "108.77273048851522 -0.257372751482126 91.70854596723325 "
    "123.6738393623219 "
    "93.16751863269833 2\n"
    "774.0890575925043 0.2962073160361278 734.5673846399011 "
    "859.107958201022 "
    "208.23276300202303 2\n"
    "626.4386861202061 -0.6192173635100073 550.9436898545746 "
    "743.1042988579782 "
    "139.74599478015597 1\n"
    "885.1471244849737 0.37031385644993353 848.0032383690232 "
    "890.4394860876546 "
    "14.03385545730026 2\n"
    "604.9366652501874 -0.00021110434061943906 604.1473698285702 "
    "621.0612958884399 140.7292004631697 2\n"
    "504.3198447299529 0.9655044814648552 498.0303651994668 "
    "523.7019107946758 "
    "20.04760141267793 2\n"
    "1055.4925624005605 -0.8274715268509594 935.5896522529015 "
    "1063.7471318467726 149.46075426452734 1\n"
    "377.1158268242932 0.8486297212850258 369.37373722650705 "
    "491.51944249340477 129.18462290947994 2\n"
    "611.9059629779958 0.6748723148214997 603.5809587375695 "
    "615.8190720491253 "
    "5.7764351581209725 2\n"
    "72.76913949302826 -0.22806008317363746 17.953994018444632 "
    "73.30774486940996 35.41318133767846 2\n"
    "552.7568349472464 0.22571874896267285 542.3283754489084 "
    "724.8908550099186 "
    "360.5092288454081 2\n"
    "775.3712648666592 0.46160471739372655 633.3888889502459 "
    "780.1982349509362 "
    "10.340120765348685 2\n"
    "310.30994589185786 0.23230160366269725 309.9057115618839 "
    "315.83414176833617 20.94816045969624 2\n"
    "202.709163829038 0.0993359878656781 102.93505186381708 "
    "291.7577591163991 "
    "190.66467188454362 2\n"
    "864.6031253007001 -0.10426192050766225 756.0270058111371 "
    "915.1802063102407 403.42205608895677 2\n"
    "516.5463987479951 -0.20851168119912789 438.8455799108792 "
    "602.9551813467388 436.8493037625859 2\n"
    "257.89291809134477 0.3637434774780983 111.86205994808051 "
    "274.53157786245765 39.078641429502184 2\n"
    "174.64208086107584 0.2936765526069718 165.94523248655833 "
    "313.3590019790171 213.89937725827036 2\n"
    "74.8815033821894 -0.15557970242348418 52.61469056758074 "
    "82.84783142113457 "
    "48.9628695838838 2\n"
    "621.8748372819153 0.31878560445123094 620.9458756418117 "
    "637.4785223065023 "
    "44.55462072806939 2\n"
    "902.7793943920236 -0.8247853234251146 817.3804632306337 "
    "999.5788339094678 "
    "106.21988151834995 1\n"
    "727.9653562971693 0.24853878471549518 723.6934931437747 "
    "757.8107764310984 "
    "96.70237770208604 2\n"
    "922.4580783054486 0.5978121517089487 873.8942336988812 "
    "1072.1573884405948 "
    "224.8789712307825 2\n"
    "803.5858151366704 0.9378139928277904 652.0201989373647 "
    "843.7829993823375 "
    "42.72364628126807 2\n"
    "377.6657537401995 -0.9542461169947254 364.8264657611521 "
    "419.7146714982823 "
    "13.478231712882788 1\n"
    "325.15915324751995 0.8671212870597 308.13800331767584 "
    "331.9311391527441 "
    "7.7836207138497775 2\n"
    "941.7235070767672 -0.21478705626605965 931.8918107910886 "
    "977.9232481585695 52.29554144103488 1\n"
    "441.9936243419802 0.010102069933997093 424.7077787330984 "
    "453.52176477543946 97.2382623337492 2\n"
    "436.4029597513562 -0.24951128210630102 427.16516873894716 "
    "461.518351256053 46.60570315306243 1\n"
    "971.4301226812836 0.10965105908230433 933.1074622942356 "
    "1052.8066597022294 313.0904481643627 2\n"
    "807.4987522221504 0.7993639668868235 760.8199321402294 "
    "817.4546241917888 "
    "12.412182637739814 2\n"
    "995.0797164156887 -0.5347063653870192 993.3205808603152 "
    "1040.6379023272104 3.2972177635803073 1\n"
    "312.47670342542983 0.13159076244904 175.6833039696328 "
    "352.6404422132688 "
    "127.41690513175446 2\n"
    "490.48725533089384 -0.15971771528341838 419.46464413066496 "
    "494.3498194479481 178.0444394891509 2\n"
    "281.96081077201745 -0.16657132545481712 270.20842352596674 "
    "293.2407112514598 140.20861578263052 2\n"
    "1073.4570124292115 -0.8751104126794789 906.3099501086626 "
    "1096.4359534893813 196.73080386865024 1\n"
    "195.73408828234932 -0.5800228373994445 190.46862046399406 "
    "309.3494515888048 9.340137381018337 1\n"
    "683.6898936686363 0.9701658212050903 673.3607708520717 "
    "690.2120043492506 "
    "6.720693823050835 2\n"
    "454.0137504005387 0.7067020956744936 424.5566524091232 "
    "547.1616557608023 "
    "122.09639811281332 2\n"
    "1015.70688847206 0.2105265569534276 919.1794286461802 "
    "1024.6432740640532 "
    "39.06596454482041 2\n"
    "1032.4372786091321 0.4086399003628842 942.912556883831 "
    "1133.3253617235673 "
    "207.7834930316393 2\n"
    "1019.7240308373591 -0.42276558284896826 986.7299288674196 "
    "1064.9038266270627 85.20003864312977 1\n"
    "223.6454060113848 -0.6358438769509664 199.37758920817942 "
    "255.81709055030308 42.42376196161549 1\n"
    "698.6563424412027 0.4727835183931348 606.6355661919615 "
    "713.6895042867209 "
    "30.711499378063273 2\n"
    "693.5208088428545 0.44991764719545646 645.7802861041371 "
    "701.7236363782354 "
    "17.83022569998685 2\n"
    "256.76753956290645 -0.907422159137401 248.88763953676903 "
    "269.05891768901193 8.713514211316491 1\n"
    "348.2718617837195 0.5108196692722697 347.5531616390008 "
    "365.54058775246546 "
    "31.801644517334598 2\n"
    "1018.8123312034563 0.24869876875646968 972.8882535836624 "
    "1134.1490575405187 305.65358503798717 2\n"
    "45.42101916463845 0.7220054570556762 37.319646309114205 "
    "45.52992360882166 "
    "0.1506707444469689 2\n"
    "676.4200698235035 0.2944646316821622 672.3024611709598 "
    "736.397692488805 "
    "153.53741434914176 2\n"
    "41.666002779726156 -0.30424783657626575 30.16005697477135 "
    "154.75749424163791 162.25796405575474 2\n"
    "725.8128195546789 0.788265596136204 700.6884933666099 "
    "747.0382131689896 "
    "26.69752520396909 2\n"
    "790.328071094304 -0.3471150484721406 768.1389947940991 "
    "794.0314629081256 "
    "72.64535303389815 1\n"
    "391.6086862497455 0.738283783419083 358.3778116903602 "
    "489.7384121047945 "
    "123.28408852287366 2\n"
    "365.75025700959446 0.7986633599060999 294.6709148764496 "
    "470.13381176563405 123.3180620592988 2\n"
    "793.0387078895052 0.36080038051320207 766.4641697208792 "
    "910.237468957348 "
    "244.42903154393287 2\n"
    "394.12585228406397 -0.4011307898771741 273.4182364763574 "
    "407.5029779676682 347.0876528637904 2\n"
    "568.7269499547582 -0.6629062707924609 558.7519484714205 "
    "581.1449813684758 "
    "15.222745600056692 1\n"
    "404.9029219219859 0.7866001206165816 240.07274550712054 "
    "417.90219194890665 16.370459010849615 2\n"
    "229.68576488684235 -0.26346561983331673 165.06160454933092 "
    "257.8816618610371 192.45891855557704 2\n"
    "339.81210636681783 0.19290211975504912 334.5667513807675 "
    "377.75761160862726 112.00348897976187 2\n"
    "594.4411408107612 -0.8889291467903377 587.404575471219 "
    "628.2433691558797 "
    "7.9284061594012325 1\n"
    "774.5154063117095 0.7359191640394425 734.4708451662807 "
    "844.4690434638976 "
    "91.93458426752541 2\n"
    "695.3121729097788 -0.4389622786691243 640.9038371014651 "
    "793.7668780732115 "
    "162.19342303737344 1\n"
    "834.0940337765742 0.3732645797244545 797.3852903915113 "
    "842.999448462057 "
    "23.126608468740027 2\n"
    "858.6321060648309 0.9053658745011735 848.2123573384874 "
    "932.0381039030967 "
    "80.3881549771462 2\n"
    "157.16305537399 -0.6858412569543297 114.00008534425137 "
    "267.8939941626765 "
    "350.03963438254135 2\n"
    "701.0675040572786 0.28061976982468906 695.24422907243 "
    "716.4229786931545 "
    "49.17363461002282 2\n"
    "487.46103100919476 0.09817565639601566 429.6818565371516 "
    "536.056367724158 "
    "180.24064303880704 2\n"
    "599.3093717409276 -0.32374204584576605 567.5243911917976 "
    "608.7272134495987 170.4290766405198 1\n"
    "558.803071176873 0.6168205244880611 543.3553712028813 "
    "703.2359860829624 "
    "204.0310269381856 2\n"
    "764.9246301685787 -0.30742003551143293 757.094767418394 "
    "842.1271101396636 "
    "26.874967409181487 1\n"
    "934.1071854494862 0.9215436527178302 915.7339252093132 "
    "936.202146795913 "
    "2.2728661894828077 2\n";
}

// End of file
