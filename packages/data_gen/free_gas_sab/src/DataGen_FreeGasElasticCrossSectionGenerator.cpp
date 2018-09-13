//---------------------------------------------------------------------------//
//!
//! \file   DataGen_FreeGasElasticCrossSectionGenerator.cpp
//! \author Eli Moll
//! \brief  Free gas elastic cross section generator
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "DataGen_FreeGasElasticCrossSectionGenerator.hpp"

#include <cmath>

namespace DataGen{

//  Constructor
FreeGasElasticCrossSectionGenerator::FreeGasElasticCrossSectionGenerator(
	   double kT,
	   std::vector<double> E,
     double A,
     int beta_num,
     int alpha_num,
     double beta_max_multiplier,
     double zero_tolerance )
  : d_kT( kT ),
    d_E( E ),
    d_A( A ),
    d_beta_num( beta_num ),
    d_alpha_num( alpha_num ),
    d_beta_max_multiplier( beta_max_multiplier ),
    d_zero_tolerance( zero_tolerance )
{
  // Make sure the values are valid
  testPrecondition( A > 0.0 );
  testPrecondition( kT > 0.0 );

  setBaseCrossSection();
  setBaseAngularDistribution();
  d_sab_function.reset( new FreeGasElasticSAlphaBetaFunction(
						                    d_cross_section, 
						                    d_angular_distribution,
						                    d_A,
						                    d_kT ) );
}

// Set base cross section
void FreeGasElasticCrossSectionGenerator::setBaseCrossSection()
{
  //std::vector<double> energy_vector = {1E-11,1.03125E-11,1.0625E-11,1.09375E-11,1.125E-11,1.15625E-11,1.1875E-11,1.21875E-11,1.25E-11,1.28125E-11,1.3125E-11,1.34375E-11,1.375E-11,1.4375E-11,1.5E-11,1.5625E-11,1.625E-11,1.6875E-11,1.75E-11,1.8125E-11,1.875E-11,1.9375E-11,2E-11,2.09375E-11,2.1875E-11,2.28125E-11,2.375E-11,2.46875E-11,2.5625E-11,2.65625E-11,2.75E-11,2.84375E-11,2.9375E-11,3.03125E-11,3.125E-11,3.21875E-11,3.3125E-11,3.40625E-11,3.5E-11,3.59375E-11,3.6875E-11,3.78125E-11,3.875E-11,3.96875E-11,4.0625E-11,4.25E-11,4.4375E-11,4.625E-11,4.8125E-11,5E-11,5.15625E-11,5.3125E-11,5.46875E-11,5.625E-11,5.78125E-11,5.9375E-11,6.09375E-11,6.25E-11,6.40625E-11,6.5625E-11,6.71875E-11,6.875E-11,7.1875E-11,7.5E-11,7.8125E-11,8.125E-11,8.4375E-11,8.75E-11,9.0625E-11,9.375E-11,9.6875E-11,1E-10,1.03125E-10,1.0625E-10,1.09375E-10,1.125E-10,1.15625E-10,1.1875E-10,1.21875E-10,1.25E-10,1.28125E-10,1.3125E-10,1.34375E-10,1.375E-10,1.4375E-10,1.5E-10,1.5625E-10,1.625E-10,1.6875E-10,1.75E-10,1.8125E-10,1.875E-10,1.9375E-10,2E-10,2.09375E-10,2.1875E-10,2.28125E-10,2.375E-10,2.46875E-10,2.5625E-10,2.65625E-10,2.75E-10,2.84375E-10,2.9375E-10,3.03125E-10,3.125E-10,3.21875E-10,3.3125E-10,3.40625E-10,3.5E-10,3.59375E-10,3.6875E-10,3.78125E-10,3.875E-10,3.96875E-10,4.0625E-10,4.25E-10,4.4375E-10,4.625E-10,4.8125E-10,5E-10,5.15625E-10,5.3125E-10,5.46875E-10,5.625E-10,5.78125E-10,5.9375E-10,6.09375E-10,6.25E-10,6.40625E-10,6.5625E-10,6.71875E-10,6.875E-10,7.1875E-10,7.5E-10,7.8125E-10,8.125E-10,8.4375E-10,8.75E-10,9.0625E-10,9.375E-10,9.6875E-10,1E-09,1.03125E-09,1.0625E-09,1.09375E-09,1.125E-09,1.15625E-09,1.1875E-09,1.21875E-09,1.25E-09,1.28125E-09,1.3125E-09,1.34375E-09,1.375E-09,1.4375E-09,1.5E-09,1.5625E-09,1.625E-09,1.6875E-09,1.75E-09,1.8125E-09,1.875E-09,1.9375E-09,2E-09,2.09375E-09,2.1875E-09,2.28125E-09,2.375E-09,2.46875E-09,2.5625E-09,2.65625E-09,2.75E-09,2.84375E-09,2.9375E-09,3.03125E-09,3.125E-09,3.21875E-09,3.3125E-09,3.40625E-09,3.5E-09,3.59375E-09,3.6875E-09,3.78125E-09,3.875E-09,3.96875E-09,4.0625E-09,4.25E-09,4.4375E-09,4.625E-09,4.8125E-09,5E-09,5.15625E-09,5.3125E-09,5.46875E-09,5.625E-09,5.78125E-09,5.9375E-09,6.09375E-09,6.25E-09,6.40625E-09,6.5625E-09,6.71875E-09,6.875E-09,7.1875E-09,7.5E-09,7.8125E-09,8.125E-09,8.4375E-09,8.75E-09,9.0625E-09,9.375E-09,9.6875E-09,1E-08,1.03125E-08,1.0625E-08,1.09375E-08,1.125E-08,1.15625E-08,1.1875E-08,1.21875E-08,1.25E-08,1.28125E-08,1.3125E-08,1.34375E-08,1.375E-08,1.4375E-08,1.5E-08,1.5625E-08,1.625E-08,1.6875E-08,1.75E-08,1.8125E-08,1.875E-08,1.9375E-08,2E-08,2.06625E-08,2.1325E-08,2.19875E-08,2.265E-08,2.33125E-08,2.3975E-08,2.46375E-08,2.53E-08,2.645782E-08,2.761563E-08,2.877344E-08,2.993125E-08,3.108907E-08,3.224688E-08,3.340469E-08,3.45625E-08,3.552735E-08,3.649219E-08,3.745704E-08,3.842188E-08,3.938673E-08,4.035157E-08,4.131641E-08,4.228125E-08,4.421094E-08,4.614063E-08,4.807032E-08,5E-08,5.15625E-08,5.3125E-08,5.46875E-08,5.625E-08,5.78125E-08,5.9375E-08,6.09375E-08,6.25E-08,6.40625E-08,6.5625E-08,6.71875E-08,6.875E-08,7.1875E-08,7.5E-08,7.8125E-08,8.125E-08,8.4375E-08,8.75E-08,9.0625E-08,9.375E-08,9.6875E-08,1E-07,1.03125E-07,1.0625E-07,1.09375E-07,1.125E-07,1.15625E-07,1.1875E-07,1.21875E-07,1.25E-07,1.28125E-07,1.3125E-07,1.34375E-07,1.375E-07,1.4375E-07,1.5E-07,1.5625E-07,1.625E-07,1.6875E-07,1.75E-07,1.8125E-07,1.875E-07,1.9375E-07,2E-07,2.09375E-07,2.1875E-07,2.28125E-07,2.375E-07,2.46875E-07,2.5625E-07,2.65625E-07,2.75E-07,2.84375E-07,2.9375E-07,3.03125E-07,3.125E-07,3.21875E-07,3.3125E-07,3.40625E-07,3.5E-07,3.59375E-07,3.6875E-07,3.78125E-07,3.875E-07,3.96875E-07,4.0625E-07,4.25E-07,4.4375E-07,4.625E-07,5E-07,5.3125E-07,5.625E-07,5.9375E-07,6.25E-07,6.875E-07,7.5E-07,8.125E-07,8.75E-07,9.375E-07,1E-06,1.0625E-06,1.125E-06,1.1875E-06,1.25E-06,1.375E-06,1.5E-06,1.625E-06,1.75E-06,1.875E-06,2E-06,2.1875E-06,2.375E-06,2.5625E-06,2.75E-06,2.9375E-06,3.125E-06,3.3125E-06,3.5E-06,3.875E-06,4.25E-06,4.625E-06,5E-06,5.3125E-06,5.625E-06,5.9375E-06,6.25E-06,6.875E-06,7.5E-06,8.125E-06,8.75E-06,9.375E-06,1E-05,1.0625E-05,1.125E-05,1.1875E-05,1.25E-05,1.375E-05,1.5E-05,1.625E-05,1.75E-05,1.875E-05,2E-05,2.1875E-05,2.375E-05,2.5625E-05,2.75E-05,2.9375E-05,3.125E-05,3.3125E-05,3.5E-05,3.875E-05,4.25E-05,4.625E-05,5E-05,5.3125E-05,5.625E-05,5.9375E-05,6.25E-05,6.875E-05,7.5E-05,8.125E-05,8.75E-05,9.375E-05,0.0001,0.00010625,0.0001125,0.00011875,0.000125,0.0001375,0.00015,0.0001625,0.000175,0.0001875,0.0002,0.00021875,0.0002375,0.00025625,0.000275,0.00029375,0.0003125,0.00033125,0.00035,0.0003875,0.000425,0.0004625,0.0005,0.00053125,0.0005625,0.00059375,0.000625,0.0006875,0.00075,0.0008125,0.000875,0.0009375,0.001,0.0010625,0.001125,0.0011875,0.00125,0.001375,0.0015,0.001625,0.00175,0.001875,0.002,0.002125,0.00225,0.002375,0.0025,0.00275,0.003,0.00325,0.0035,0.00375,0.004,0.00425,0.0045,0.00475,0.005,0.00525,0.0055,0.00575,0.006,0.00625,0.0065,0.00675,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.01,0.01125,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.1125,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.75,0.8,0.85,0.9,0.95,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20};
  //std::vector<double> xs_vector     = {37.134595704,36.879616558,36.63594252,36.402811346,36.179448673,35.965217456,35.759488032,35.561742637,35.371448198,35.188170529,35.011462607,34.840959521,34.676313854,34.363268937,34.070005477,33.794525904,33.535101197,33.290225603,33.058593384,32.839054577,32.630596472,32.432302651,32.24338955,31.97604298,31.726079903,31.491702316,31.271346956,31.063676393,30.867518529,30.681842891,30.505759893,30.338455416,30.179239307,30.027474852,29.882595261,29.744101522,29.611536314,29.484488513,29.362589277,29.245494999,29.132902928,29.02453051,28.920120829,28.819438104,28.722239132,28.53715164,28.363944535,28.20136169,28.048382242,27.90408784,27.789907244,27.680804143,27.576414658,27.476407646,27.380484851,27.288377526,27.199837962,27.114642264,27.032584281,26.953478261,26.877156032,26.803448697,26.66332579,26.532067224,26.408776437,26.292682891,26.183115273,26.079482198,25.981276807,25.88803188,25.799352658,25.714874245,25.63428294,25.557286746,25.483627853,25.41307724,25.345412803,25.280459804,25.21802779,25.157966567,25.100128147,25.04438441,24.990607201,24.938687502,24.840014254,24.747621138,24.660880051,24.57917809,24.50133626,24.427719172,24.357961039,24.291738997,24.228765265,24.168792461,24.083946754,24.004664255,23.930356112,23.860554713,23.794797707,23.732729363,23.674023027,23.618390467,23.565567914,23.515349485,23.46751678,23.421906521,23.378336714,23.336679942,23.296796441,23.258574494,23.221897001,23.186297568,23.151608072,23.118207612,23.086028183,23.054981651,22.996065543,22.94102264,22.889440366,22.841011637,22.795407245,22.759395009,22.725044165,22.692238155,22.660871749,22.630853648,22.602085677,22.574499726,22.548011297,22.522554281,22.497557149,22.47298377,22.449292813,22.404368909,22.362417265,22.323170873,22.286366585,22.251767732,22.219187798,22.188471296,22.159443809,22.13199142,22.10597054,22.080343903,22.055317493,22.031464339,22.008723272,21.987003812,21.966259475,21.946409811,21.927396277,21.90920187,21.891761638,21.875027396,21.858971403,21.828724645,21.799027891,21.770569807,21.744050444,21.719297604,21.696164112,21.674519254,21.654244857,21.635241526,21.616879287,21.588854362,21.563111484,21.539437231,21.517643026,21.497562557,21.479051467,21.460056811,21.441515203,21.424255694,21.408184948,21.393220087,21.379286228,21.366313043,21.353435365,21.339758314,21.326913207,21.314856621,21.303543063,21.292928829,21.282978239,21.2736563,21.264675029,21.244608391,21.226705922,21.21078101,21.196664037,21.182214403,21.16986815,21.158567388,21.148250216,21.13885841,21.130340678,21.120769815,21.111057286,21.102125602,21.093925566,21.08642775,21.079595589,21.072940737,21.057326301,21.044042683,21.03288083,21.021020363,21.009184732,20.99910986,20.990655878,20.980367336,20.970867947,20.962763837,20.955948834,20.947408511,20.939383762,20.932495229,20.926667945,20.920107522,20.913003427,20.90685499,20.901605136,20.897204705,20.890684925,20.8849631,20.875799091,20.866254293,20.857248764,20.850362873,20.84130336,20.834427962,20.827765548,20.820476833,20.815000437,20.808735604,20.802315036,20.797598298,20.791316697,20.785971297,20.782071108,20.776243674,20.771616058,20.768223788,20.759795204,20.753627845,20.746609672,20.740902704,20.734925717,20.729703782,20.724435312,20.719737103,20.715602134,20.712678216,20.70832814,20.705148628,20.701903637,20.698397079,20.695892109,20.692476647,20.687342841,20.68161837,20.676513384,20.672181712,20.668327609,20.665024358,20.661681095,20.658541754,20.655562079,20.652618074,20.649886009,20.647163609,20.644612974,20.642122737,20.639685974,20.637444314,20.633116691,20.629090005,20.625121065,20.621341367,20.617819758,20.614545226,20.611501303,20.608674718,20.605881521,20.603095827,20.59918384,20.597006051,20.59489342,20.592839106,20.590848256,20.588914848,20.587034712,20.585217486,20.583453456,20.581751801,20.580097563,20.578489206,20.57542801,20.572539293,20.569810305,20.567244236,20.564807331,20.562500588,20.560317971,20.558237736,20.55625931,20.554372434,20.551696261,20.549198469,20.546853456,20.544650923,20.542574914,20.540607061,20.538756821,20.536992993,20.535317509,20.533722922,20.532206485,20.530754495,20.529364753,20.528039064,20.526766215,20.525543213,20.524376341,20.523244222,20.522156548,20.521116517,20.520101932,20.519124288,20.517275163,20.515549626,20.513924891,20.510945137,20.508715835,20.506672805,20.504785008,20.503048613,20.499935511,20.497207319,20.494812889,20.492668593,20.490748083,20.489005589,20.487418711,20.485963899,20.4846304,20.483397305,20.48117619,20.47924003,20.477537268,20.476013783,20.474636886,20.473396614,20.471738769,20.470278024,20.468972353,20.467802828,20.466745716,20.465782024,20.464900623,20.464086784,20.462629448,20.461369573,20.46026537,20.45927285,20.458521921,20.457832158,20.457198959,20.456608765,20.455538536,20.454596941,20.453764426,20.453006012,20.452312811,20.45167988,20.451102804,20.450558724,20.450058823,20.449579883,20.4675322,20.4651954,20.4631604,20.4613748,20.4597973,20.4583793,20.4565007,20.4548555,20.4534013,20.4521001,20.4509239,20.4498467,20.4488451,20.4479257,20.4462623,20.4447891,20.4434535,20.442227,20.4412842,20.4403898,20.4395417,20.4387279,20.4371967,20.4357752,20.4344258,20.4331531,20.431922,20.4307385,20.4295905,20.4284658,20.4273727,20.4262998,20.4242096,20.4221776,20.4201982,20.4182576,20.4163424,20.4144496,20.4116652,20.4089209,20.4062028,20.4035082,20.4008248,20.3981706,20.3955246,20.3928956,20.3876652,20.3824639,20.3772777,20.3721237,20.3678628,20.3636078,20.3593681,20.3551229,20.346674,20.3382278,20.3298018,20.3213942,20.3129935,20.3045985,20.296298,20.2880218,20.2797393,20.2714601,20.2549198,20.2383786,20.2218647,20.2053467,20.1888436,20.1723446,20.1561885,20.1400555,20.1239251,20.1077971,20.0755469,20.0433133,20.0110649,19.9788307,19.94661,19.9144123,19.8830965,19.8518029,19.8205112,19.7892312,19.7579427,19.7266655,19.6953794,19.6641244,19.6337,19.6033164,19.5729236,19.5425416,19.4817593,19.4210093,19.3619107,19.3028337,19.2437681,19.1847438,19.044008,18.9033222,18.6220408,18.3589895,18.0961984,17.8709633,17.6459383,17.4209033,17.1958783,16.9708433,16.7458083,16.5207833,16.2960783,16.1172116,15.938675,15.7601283,15.5815917,15.403055,15.2245184,15.0459717,14.8677251,14.7222512,14.5770673,14.4318835,14.2866996,13.9963219,13.7062041,13.465069,13.2241739,12.9832889,12.7427238,12.2876665,11.8329391,10.9238445,10.2833569,9.64358528,9.22042561,8.79761494,8.37480327,7.95229561,7.68310802,7.41422244,7.14533586,6.87663428,6.68870753,6.50096578,6.31322304,6.12560229,5.98583894,5.84619659,5.70655525,5.5669989,5.45819576,5.34947862,5.24076048,5.13210634,4.95682329,4.78165124,4.63655432,4.4915424,4.36882155,4.2461807,3.85048846,3.54178322,3.29134987,3.08222339,2.90368177,2.74858002,2.61195516,2.49023421,2.3807732,2.28155813,2.19103001,2.10795386,2.03133769,1.9603705,1.89438531,1.8328231,1.77521289,1.72115368,1.67029947,1.62235327,1.51358777,1.41819131,1.33374287,1.25840048,1.19073011,1.12959677,1.07408446,1.02344716,0.977066588,0.934429025,0.895099875,0.858710735,0.824946302,0.793535076,0.764241754,0.736861535,0.711214819,0.687143905,0.664509492,0.643187979,0.623069267,0.604055154,0.58605774,0.568997726,0.552804011,0.537412095,0.522763779,0.508805961,0.495490443,0.482773424};
  //std::vector<double> energy_vector = {1e-11, 1e-1, 20};
  //std::vector<double> xs_vector     = {20.43, 20.43, 0.5};
  
  std::vector<double> energy_vector = {1e-11, 1.154357e-11, 20};
  std::vector<double> xs_vector = {4.790229528, 4.790246677, 0.1};
  
  const Teuchos::Array<double> energy_array( energy_vector );
  const Teuchos::Array<double> xs_array( xs_vector );
  d_cross_section.reset( new Utility::TabularDistribution<Utility::LinLin>( energy_array, xs_array ) );
}

// Set base angular distribution
void FreeGasElasticCrossSectionGenerator::setBaseAngularDistribution()
{
  // Initialize the scattering probability distribution
  Teuchos::RCP<Utility::TabularOneDDistribution> isotropic_distribution(
			  new Utility::UniformDistribution( -1.0, 1.0, 0.5 ) );

  // Initialize the scattering distribution
  MonteCarlo::NuclearScatteringAngularDistribution::AngularDistribution
    distribution( 2 );

  distribution[0].first = 0.0;
  distribution[0].second = isotropic_distribution;
  
  distribution[1].first = 20.0;
  distribution[1].second = isotropic_distribution;

  d_angular_distribution.reset( 
			 new MonteCarlo::NuclearScatteringAngularDistribution(
							      distribution ) );
}

// Calculate Analytical Cross Section for Isotropic Scattering and Unity Sigma
double FreeGasElasticCrossSectionGenerator::analyticCrossSectionValue( 
        double alpha, 
        double beta,
        double E )
{
  double pi3 = Utility::PhysicalConstants::pi*
    Utility::PhysicalConstants::pi*
    Utility::PhysicalConstants::pi;

  if( alpha > 0.0 || beta > 0.0 )
  {
    return (d_kT*(d_A+1)*(d_A+1)/(16*sqrt(pi3)*d_A*E*sqrt(alpha)))*
      exp( -(alpha + beta)*(alpha + beta)/(4*alpha) );
  }
  else
    return std::numeric_limits<double>::infinity();
}

// Calculate cross section
double FreeGasElasticCrossSectionGenerator::crossSectionValue( 
        double beta_int,
        double E )
{
  double pi3 = Utility::PhysicalConstants::pi*
    Utility::PhysicalConstants::pi*
    Utility::PhysicalConstants::pi;

  return (d_A+1)*(d_A+1)*(d_A+1)*(d_A+1)*(d_kT/E)/(4*d_A*sqrt(pi3))*beta_int;
}

//  Constructruct full double differential cross section at a given energy
void FreeGasElasticCrossSectionGenerator::doubleDifferentialCrossSectionValue( 
        double E,
        DoubleDifferentialCrossSection& double_differential_sigma )
{
  double beta_min = Utility::calculateBetaMin( E, d_kT );
  double beta_max = d_beta_max_multiplier*beta_min;
  double beta_spread = (beta_max - beta_min)/(d_beta_num - 1.0);

  // Loop over all beta for given energy
  for ( int j = 0; j < d_beta_num; ++j )
  {
    double beta = beta_min + j*beta_spread;

    // Correct for zero value beta
    if (beta > -1.0*d_zero_tolerance && beta < 0.0 )
    {
      beta = -1*d_zero_tolerance;
    }
    else if (beta < 1.0*d_zero_tolerance && beta > 0.0 )
    {
      beta = d_zero_tolerance;
    }
    else if (beta <= beta_min)
    {
      beta = beta_min - beta_min*1e-3;
    }

    double alpha_min = Utility::calculateAlphaMin( E, 
                                                   beta, 
                                                   d_A, 
                                                   d_kT);

    double alpha_max = Utility::calculateAlphaMax( E, 
                                                   beta, 
                                                   d_A, 
                                                   d_kT);

    double alpha_spread = (alpha_max - alpha_min)/(d_alpha_num - 1.0);

    for (int k = 0; k < d_alpha_num; ++k )
    {
      double alpha = alpha_min + k*alpha_spread;

      double sab = (*d_sab_function)( alpha, beta, E );

      double value = crossSectionValue( sab, E );

      std::pair<double,double> beta_alpha( beta, alpha );
      double_differential_sigma[beta_alpha] = value;    
    }
  }
}

// Integrate over energy and angle for a total cross section value at a given energy
void FreeGasElasticCrossSectionGenerator::totalCrossSectionValue( 
           double E )
{
  if (E > 500*d_kT/d_A)
  {
    d_total_cross_section[E] = d_cross_section->evaluate( E );
  }
  else
  {
    d_beta_function.reset( new DataGen::FreeGasElasticMarginalBetaFunction(
                  d_cross_section, 
                  d_angular_distribution,
                  d_A,
                  d_kT,
                  E ) );

    d_total_cross_section[E] = this->crossSectionValue( d_beta_function->getNormalizationConstant(), E );
  }
}

void FreeGasElasticCrossSectionGenerator::energyCrossSectionValue(
           double E )
{
  d_beta_function.reset( new DataGen::FreeGasElasticMarginalBetaFunction(
						    d_cross_section, 
						    d_angular_distribution,
						    d_A,
						    d_kT,
						    E ) );

  double beta_min = Utility::calculateBetaMin( E, d_kT );
  double beta_max = d_beta_max_multiplier*beta_min;
  double beta_spacing = (beta_max - beta_min)/(d_beta_num - 1.0);

  DifferentialEnergyCrossSection beta_pdf;

  // Loop over beta
  for (int i = 0; i < d_beta_num; ++i)
  {
    if (true)
    {
      double beta = beta_min + i*beta_spacing;
      double pdf  = d_beta_function->operator()( beta );
      beta_pdf.push_back( std::make_pair(beta, pdf) ); 
    }
    else
    {
      double beta = beta_min + i*beta_spacing;

      d_alpha_function.reset( new DataGen::FreeGasElasticMarginalAlphaFunction(
						    d_cross_section, 
						    d_angular_distribution,
						    d_A,
						    d_kT,
						    beta,
                E ) );

      double pdf = d_alpha_function->getNormalizationConstant();

      beta_pdf.push_back( std::make_pair(beta, pdf) );
    }
  }

  d_beta_pdf_map[E] = beta_pdf;
}

// Get total cross section
void FreeGasElasticCrossSectionGenerator::getTotalCrossSection( 
    boost::unordered_map< double, double >& total_cross_section )
{
  for (int i = 0; i < d_E.size(); ++i)
  {
    this->totalCrossSectionValue( d_E[i] );
  }

  total_cross_section = d_total_cross_section;
}

// Calculate cross sections for all energies
void FreeGasElasticCrossSectionGenerator::getDifferentialEnergyCrossSectionMap(
  DifferentialEnergyCrossSectionMap& energy_cross_section_map )
{
  // Loop over all energies
  for ( int i = 0; i < d_E.size(); ++i) 
  {
    energyCrossSectionValue( d_E[i] );
  }

  energy_cross_section_map = d_beta_pdf_map;
}

} // end DataGen namespace

//---------------------------------------------------------------------------//
// end DataGen_FreeGasElasticCrossSectionGenerator.cpp
//---------------------------------------------------------------------------//
