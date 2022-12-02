#ifndef PAU_CONSTS__H
#define PAU_CONSTS__H

#include <array>
#include "ioClass.h"
//glue in the binning for UEbbc
const ioBinVec bins_EAbbc30 {{ 0., 0., 30., 64000. }};
const ioBinVec bins_EAbbc50 {{ 0., 0., 50., 64000. }};


// charged track pT bins
const ioBinVec bins_trpt {{ // track pT bins
     0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,
     1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,
     3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4.0,  4.1,  4.2,  4.3,  4.4,
     4.6,  4.8,  5.0,  5.2,  5.4,  5.6,  5.8,  6.0,  6.2,  6.4,  6.6,  6.8,  7.1,  7.4,  7.7,
     8.0,  8.3,  8.6,  9.0,  9.4,  9.8, 10.3, 10.8, 11.3, 11.9, 12.5, 13.2, 14.0, 15.0
}};

const ioBinVec bins_etaEMW {{ -1., -.3, 0.3, 1. }}; // EMW for EastMidWest


// important: update: bin decile values selected for MB 500004 events, with good run list, with |Vz|<10, |Vz-VzVPD|<6
//
// 2022.06.15 Bins update for ZDCx set to range 4 to 22kHz
 const ioBinVec bins_EAbbc10 {{     0.0000,  2751.8790,  5376.4117,  8315.4606, 11603.2222, 15289.4399, 19479.7043, 24292.8207, 30059.4932, 37647.3152, 64000.0000 }};
/* const ioBinVec bins_EAbbc10 {{ 0, 2752.52, 5378.27, 8315.51, 11603.5, 15290.3, 19479.4, 24293.3, 30060.2, 37651.3, 64000 }}; */
const ioBinVec bins_EAbbc10_arch0 {{ 0, 2767.15, 5397.98, 8333.35, 11610.5, 15280.9, 19440.2, 24219.7, 29959, 37534.5, 64000 }};
/* const ioBinVec bins_EAbbc10 {{ 0, 2768.06, 5400.77, 8334.48, 11611.2, 15280.8, 19440.2, 24222.4, 29960.6, 37539.6, 64000 */
const ioBinVec bins_EAbbc3 {{     0.0000,  8315.4606, 24292.8207, 64000.0000 }};
/* const ioBinVec bins_EAbbc3  {{ 0, 8333.35, 24219.7, 64000 }}; */

const ioBinVec bins_EAbbc10_LEGACY {{ 0, 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000 }};
const ioBinVec bins_EAbbc3_LEGACY  {{ 0, 10126.1, 26718.1, 64000 }};

const std::array<double,3>  NchZDCx_pfit__DCA3 {-0.5518790, 4.96012e-4, -9.56337e-09};
const std::array<double,3>  NchZDCx_pfit       { 0.8394191, 8.241e-06,  -2.42942e-10};

 const std::array<double,3>  PU_pol2_fit      {  1.011172473, 4.540577178e-05, -7.810329033e-10 };
 const std::array<double,3>  east_PU_pol2_fit {  1.118514079, 4.53705939e-05, -8.667241617e-10  };
 const std::array<double,3>  west_PU_pol2_fit { 0.9173822694, 3.562268455e-05, -2.433180247e-10 };
 const std::array<double,3>  mid_PU_pol2_fit  {  1.006795055, 5.512093342e-05, -1.231373389e-09 };


// OK, track density binning should center bins on integer number of tracks. 
//  - rho(whole TPC)  \equiv n/(4pi)
//  - rho(transverse) \equiv n/(4/3pi)
// The min bias event top out at about 81 tracks, which is equivalent to 27 tracks transverse
const ioBinVec bins_rhoTPC      {{ 0., 0.5/4/M_PI, 0.5/4/M_PI, 81, (81+0.5)/4/M_PI }};
const ioBinVec bins_rhoTPCtrans {{ 0., 0.5/(4.0/3.0*M_PI), 0.5/(4.0/3.0*M_PI), 26, (26+0.5)/(4.0/3.0*M_PI), (81+0.5)/4/M_PI}};

const ioBinVec bins_UEtpc10 {{ 0.0000, 0.5472, 0.7240, 0.8761, 1.0230, 1.1744, 1.3392, 1.5305, 1.7728, 2.1406, 6.4856 }};
const ioBinVec bins_UEtpc3  {{ 0.0000, 0.8761, 1.5305, 6.4856 }};

/* // note that the inherent binning in UEtpc is in increments of 1/4pi (which is 0.0795774), */
/* // or, when using transverse, is in increments of 3/4pi (1/3 is transverse) */
/* // -- ideally put in the center of each bin */
/* // So, going out to 80 particles is like, (using x=1/4/M_PI) */
/* // for for 1/(4pi):  0, x/2, 3x/2, 5x/2 ... 161x/x */
/* const ioBinVec bins_tr_rho_per4pi   {{ 0., 1./8/M_PI, 1./8/M_PI, 81, (1.+81.*2.)/8./M_PI }};   //  | 0        1        2     ....    n-1      n */
/*                                                                                                //  |      |        |                       | */
/*                                                                                                //  rho = n/(4 pi) */
/* const ioBinVec bins_tr_rho_per4o3pi {{ 0., 3./8/M_PI, 3./8/M_PI, 27, (1.+27.*2.)*3/8./M_PI }}; // "4 over 3 pi   n/(4/3 pi) */

//glue in the bad tower list
const ioIntSet BAD_TOWERS_v0 {{
     35, 95, 96, 104, 214, 246, 308, 315, 493, 555, 562, 606, 740, 750, 758, 784,
     822, 882, 887, 897, 923, 924, 972, 991, 1001, 1028, 1125, 1130, 1132,
     1183, 1197, 1217, 1233, 1250, 1257, 1273, 1294, 1312, 1375, 1378, 1382,
     1392, 1537, 1709, 1720, 1909, 1943, 1984, 2023, 2027, 2043, 2063, 2065,
     2162, 2303, 2391, 2439, 2459, 2487, 2633, 2749, 2792, 2865, 3420, 3473,
     3513, 3588, 3692, 3720, 3737, 3821, 3831, 3834, 3838, 3840, 3859, 4013,
     4017, 4053, 4057, 4217, 4319, 4355, 4369, 4437, 4498, 4563, 4657, 4766
}};

//glue in the good runs list
const ioIntSet GOOD_RUNS_v0 {{
    16125035, 16125039, 16125043, 16125050, 16125053, 16125054, 16125056,
    16125057, 16126006, 16126007, 16126008, 16126009, 16126010, 16126011,
    16126012, 16126013, 16126015, 16126016, 16126017, 16126018, 16126030,
    16126031, 16126032, 16126033, 16126034, 16126035, 16126036, 16126037,
    16126038, 16126039, 16126040, 16126041, 16126042, 16126043, 16126044,
    16126050, 16126051, 16126052, 16126053, 16127007, 16127008, 16127010,
    16127011, 16127012, 16127013, 16127014, 16127015, 16127016, 16127017,
    16127018, 16127019, 16127020, 16127021, 16127023, 16127025, 16127026,
    16127027, 16127028, 16127029, 16127030, 16127031, 16127032, 16127033,
    16127047, 16127052, 16127053, 16127055, 16128001, 16128002, 16128003,
    16128007, 16128008, 16128012, 16128015, 16128016, 16128021, 16128022,
    16128028, 16128029, 16128030, 16128034, 16128035, 16128036, 16128037,
    16128039, 16128040, 16128041, 16128043, 16128048, 16128049, 16128050,
    16128051, 16128052, 16128055, 16128057, 16128058, 16128059, 16128064,
    16128065, 16129001, 16129002, 16129003, 16129004, 16129005, 16129006,
    16129007, 16129008, 16129009, 16129010, 16129011, 16129013, 16129014,
    16129022, 16129023, 16129025, 16129029, 16129030, 16129031, 16129032,
    16129033, 16129034, 16129035, 16129036, 16129037, 16129038, 16129039,
    16129043, 16129044, 16129045, 16129046, 16129047, 16129048, 16129049,
    16129050, 16129051, 16130001, 16130002, 16130003, 16130004, 16130005,
    16130033, 16130034, 16130035, 16130038, 16130040, 16130044, 16130045,
    16130046, 16130048, 16130049, 16130051, 16130052, 16131004, 16131005,
    16131006, 16131007, 16131008, 16131009, 16131010, 16131011, 16131012,
    16131013, 16131014, 16131028, 16131030, 16131039, 16131041, 16131042,
    16131043, 16131044, 16131047, 16131049, 16131050, 16132001, 16132003,
    16132021, 16132022, 16132023, 16132024, 16132025, 16132026, 16132027,
    16132030, 16132032, 16132048, 16132049, 16132050, 16132051, 16132053,
    16132054, 16133001, 16133002, 16133003, 16133004, 16133006, 16133007,
    16133008, 16133009, 16133010, 16133011, 16133013, 16133086, 16133087,
    16133088, 16133089, 16133091, 16133092, 16133093, 16134001, 16134002,
    16134003, 16134004, 16134005, 16134013, 16134014, 16134015, 16134016,
    16134017, 16134018, 16134019, 16134028, 16134029, 16134030, 16134043,
    16134044, 16134045, 16134046, 16134047, 16134048, 16134054, 16134055,
    16135001, 16135002, 16135003, 16135004, 16135005, 16135006, 16135007,
    16135008, 16135009, 16135010, 16135011, 16135013, 16135014, 16135025,
    16135032, 16135033, 16135034, 16135035, 16135036, 16135037, 16135038,
    16135039, 16135040, 16135041, 16135043, 16135048, 16135049, 16135053,
    16135056, 16135061, 16136001, 16136002, 16136003, 16136004, 16136011,
    16136013, 16136014, 16136016, 16136017, 16136035, 16136036, 16136043,
    16136044, 16136045, 16136046, 16136047, 16136048, 16136049, 16136050,
    16136051, 16136052, 16136053, 16136057, 16137001, 16137002, 16137003,
    16137004, 16137005, 16137006, 16137007, 16137008, 16137010, 16137011,
    16137012, 16137013, 16137017, 16137018, 16137019, 16137020, 16137021,
    16137022, 16137028, 16137029, 16137030, 16137031, 16137032, 16137033,
    16137038, 16137039, 16137040, 16137041, 16137042, 16137043, 16137045,
    16137046, 16137047, 16137048, 16137049, 16137052, 16138014, 16138015,
    16138016, 16138017, 16138021, 16138022, 16138023, 16138024, 16138025,
    16138026, 16138027, 16138028, 16138029, 16138040, 16138042, 16138044,
    16138045, 16138046, 16138047, 16138048, 16138049, 16138050, 16138051,
    16138052, 16138054, 16138059, 16138061, 16138062, 16138063, 16138064,
    16138065, 16138066, 16139002, 16139003, 16139004, 16139005, 16139007,
    16139008, 16139009, 16139010, 16139011, 16139012, 16139013, 16139021,
    16139022, 16139023, 16139024, 16139025, 16139026, 16139027, 16139028,
    16139029, 16139030, 16139031, 16139032, 16139033, 16139034, 16139054,
    16139055, 16139056, 16139057, 16140001, 16140002, 16140003, 16140004,
    16140008, 16140009, 16140010, 16140012, 16140013, 16141007, 16141008,
    16141009, 16141014, 16141016, 16141017, 16141018, 16141021, 16141022,
    16141029, 16141032, 16141033, 16141034, 16141035, 16141038, 16141040,
    16141049, 16141050, 16142001, 16142003, 16142004, 16142006, 16142007,
    16142008, 16142010, 16142012, 16142013, 16142014, 16142015, 16142022,
    16142024, 16142025, 16142027, 16142028, 16142032, 16142035, 16142036,
    16142037, 16142038, 16142039, 16142044, 16142061, 16142065, 16142066,
    16142067, 16142068, 16142069, 16142070, 16142071, 16142072, 16142076,
    16142077, 16142079, 16142080, 16143002, 16143003, 16143004, 16143005,
    16143006, 16143007, 16143009, 16143010, 16143011, 16143013, 16143014,
    16143015, 16143016, 16143018, 16143019, 16143020, 16143021, 16143022,
    16143031, 16143032, 16143033, 16143034, 16143035, 16143036, 16143038,
    16143039, 16143040, 16143041, 16143042, 16143043, 16143044, 16143045,
    16143052, 16143053, 16143054, 16143055, 16143056, 16143057, 16143058,
    16143059, 16144001, 16144002, 16144003, 16144004, 16144005, 16144007,
    16144013, 16144014, 16144015, 16144016, 16144019, 16144020, 16144021,
    16144023, 16144024, 16144026, 16144038, 16144039, 16144040, 16144042,
    16144056, 16144061, 16144065, 16144070, 16144071, 16144073, 16144075,
    16144076, 16145003, 16145010, 16145011, 16145012, 16145013, 16145014,
    16145015, 16145016, 16145017, 16145018, 16145019, 16145020, 16145021,
    16145022, 16145023, 16145035, 16145036, 16145039, 16145040, 16145041,
    16145042, 16145043, 16145044, 16145045, 16145046, 16145048, 16145049,
    16145050, 16145055, 16145056, 16145057, 16145058, 16145059, 16145060,
    16145061, 16146003, 16146004, 16146005, 16146007, 16146008, 16146015,
    16146016, 16146051, 16146052, 16146097, 16146098, 16146100, 16146106,
    16146107, 16146109, 16146111, 16146112, 16146114, 16146115, 16146116,
    16146117, 16146118, 16146119, 16147006, 16147007, 16147008, 16147009,
    16147010, 16147012, 16147013, 16147014, 16147056, 16147057, 16147058,
    16147059, 16147061, 16147062, 16147063, 16147064, 16147065, 16147066,
    16147067, 16147068, 16147069, 16147070, 16148002, 16148003, 16148004,
    16148005, 16148006, 16148007, 16148008, 16148009, 16148010, 16148019,
    16148020, 16148021, 16148022, 16149044, 16149049, 16149050, 16149052,
    16149053, 16149054, 16149055, 16149056, 16150004, 16150010, 16150011,
    16150012, 16150013, 16150014, 16150015, 16150018, 16150019, 16150020,
    16150021, 16150022, 16150023, 16150024, 16150044, 16150045, 16150046,
    16150048, 16150049, 16150052, 16150054, 16150055, 16150056, 16150057,
    16150058, 16150059, 16150061, 16150065, 16150066, 16151001, 16151002,
    16151022, 16151023, 16151024, 16151025, 16151026, 16151027, 16151028,
    16151029, 16151030, 16151031, 16151032, 16151034, 16151035, 16151041,
    16151042, 16151043, 16151044, 16151045, 16151046, 16151048, 16151049,
    16151050, 16151053, 16151055, 16151056, 16151062, 16152007, 16152008,
    16152009, 16152019, 16152021, 16152022, 16152024, 16152025, 16152028,
    16152029, 16152030, 16152031, 16152032, 16152034, 16152035, 16152045,
    16152046, 16152047, 16152049, 16152051, 16152052, 16152054, 16152057,
    16152058, 16152059, 16153001, 16153007, 16153009, 16153010, 16153011,
    16153012, 16153013, 16153014, 16153015, 16153016, 16153019, 16153021,
    16153022, 16153032, 16153033, 16153034, 16153036, 16153037, 16153039,
    16153040, 16153041, 16153042, 16153043, 16153044, 16153047, 16153048,
    16154003, 16154004, 16154005, 16154006, 16154008, 16154009, 16154011,
    16154022, 16155001, 16155002, 16155003, 16155004, 16155006, 16155007,
    16155008, 16155011, 16155012, 16155013, 16155014, 16155015, 16155016,
    16155019, 16155028, 16155032, 16155034, 16155038, 16155041, 16155042,
    16155047, 16155048, 16156005, 16156007, 16156008, 16156009, 16156011,
    16156013, 16156014, 16156015, 16156017, 16156019, 16156020, 16156021,
    16156022, 16156043, 16156048, 16156049, 16156050, 16156051, 16156052,
    16156053, 16156054, 16156056, 16156057, 16157004, 16157005, 16157006,
    16157007, 16157013, 16157015, 16157016, 16157018, 16157019, 16157020,
    16157021, 16157025, 16157026, 16157027, 16157035, 16157036, 16157037,
    16157038, 16157039, 16157040, 16157042, 16157043, 16157044, 16157045,
    16157046, 16157062, 16157063, 16157064, 16157065, 16157066, 16157068,
    16157069, 16157072, 16157073, 16158001, 16158002, 16158005, 16158006,
    16158016, 16158017, 16158018, 16158019, 16158020, 16158023, 16158024,
    16158027, 16158030, 16158031, 16158033, 16158045, 16158046, 16158047,
    16158048, 16158049, 16158051, 16158062, 16158064, 16158066, 16158069,
    16158070, 16159002, 16159003, 16159004, 16159011, 16159012, 16159013,
    16159014, 16159017, 16159018, 16159020, 16159023, 16159024
}};


#endif