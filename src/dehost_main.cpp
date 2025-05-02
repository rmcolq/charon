#include <unordered_set>
#include <iostream>
#include <algorithm>

#include "dehost_main.hpp"
#include "classify_stats.hpp"
#include "index.hpp"
#include "load_index.hpp"
#include "utils.hpp"
#include "version.h"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>
#include <seqan3/alphabet/quality/phred_base.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

const std::vector<float> default_pos_data = {0.603448, 0.480363, 0.411268, 0.444444, 0.392523, 0.31383 ,
                                             0.281447, 0.53539 , 0.553792, 0.537906, 0.273994, 0.693023,
                                             0.680751, 0.428571, 0.254005, 0.696133, 0.251462, 0.594203,
                                             0.361582, 0.515789, 0.724771, 0.312821, 0.457627, 0.728614,
                                             0.573222, 0.440476, 0.381633, 0.234127, 0.673988, 0.368421,
                                             0.548387, 0.304688, 0.365348, 0.280899, 0.482315, 0.445283,
                                             0.449541, 0.59589 , 0.671362, 0.59854 , 0.445415, 0.394265,
                                             0.525597, 0.201365, 0.702703, 0.47505 , 0.509174, 0.692308,
                                             0.230573, 0.553616, 0.653846, 0.690476, 0.432203, 0.674377,
                                             0.298387, 0.626712, 0.589347, 0.545045, 0.624679, 0.240658,
                                             0.254003, 0.746667, 0.357143, 0.606909, 0.344558, 0.409836,
                                             0.266998, 0.588235, 0.537815, 0.50641 , 0.371951, 0.415929,
                                             0.689266, 0.309735, 0.705645, 0.549296, 0.80754 , 0.398979,
                                             0.459941, 0.491647, 0.292453, 0.43871 , 0.649254, 0.309154,
                                             0.617021, 0.601671, 0.383721, 0.628571, 0.487936, 0.440367,
                                             0.382465, 0.660297, 0.742489, 0.737143, 0.282759, 0.490196,
                                             0.247191, 0.690994, 0.731602, 0.666667, 0.535714, 0.322368,
                                             0.694389, 0.273333, 0.783972, 0.307692, 0.284091, 0.661333,
                                             0.366492, 0.733728, 0.352113, 0.803922, 0.269911, 0.366599,
                                             0.642197, 0.662921, 0.294118, 0.234146, 0.345309, 0.364407,
                                             0.352941, 0.288462, 0.547619, 0.712121, 0.53211 , 0.510823,
                                             0.631579, 0.490132, 0.651007, 0.653846, 0.773931, 0.550265,
                                             0.56044 , 0.45977 , 0.445743, 0.643761, 0.232422, 0.564171,
                                             0.443966, 0.638356, 0.314516, 0.276398, 0.51746 , 0.745342,
                                             0.366093, 0.471049, 0.205198, 0.643423, 0.581197, 0.453686,
                                             0.587719, 0.305322, 0.473118, 0.582593, 0.595486, 0.573822,
                                             0.370192, 0.681661, 0.513514, 0.289655, 0.349364, 0.35873 ,
                                             0.571429, 0.28125 , 0.785425, 0.689055, 0.675   , 0.255814,
                                             0.542553, 0.333333, 0.544643, 0.668831, 0.534483, 0.237838,
                                             0.372308, 0.548387, 0.227818, 0.33101 , 0.640931, 0.55814 ,
                                             0.355993, 0.643443, 0.333333, 0.323288, 0.634146, 0.352   ,
                                             0.243902, 0.333333, 0.663968, 0.246212, 0.540323, 0.565217,
                                             0.715847, 0.655431, 0.707736, 0.407792, 0.52924 , 0.511292,
                                             0.281195, 0.2875  , 0.323296, 0.230392, 0.262222, 0.371245,
                                             0.206025, 0.438525, 0.660377, 0.30426 , 0.261044, 0.360902,
                                             0.34    , 0.695   , 0.455462, 0.621891, 0.41003 , 0.349481,
                                             0.656863, 0.482143, 0.379032, 0.389474, 0.658451, 0.295181,
                                             0.745614, 0.23622 , 0.520408, 0.570213, 0.476923, 0.336864,
                                             0.431227, 0.389381, 0.368601, 0.489879, 0.228682, 0.413636,
                                             0.511905, 0.34375 , 0.688235, 0.717131, 0.340996, 0.521739,
                                             0.302752, 0.724551, 0.247059, 0.634551, 0.375204, 0.585271,
                                             0.488372, 0.707071, 0.698305, 0.384615, 0.620743, 0.43738 ,
                                             0.8     , 0.275   , 0.434783, 0.335079, 0.28789 , 0.242915,
                                             0.657588, 0.664134, 0.255172, 0.636364, 0.349593, 0.849057,
                                             0.371681, 0.535032, 0.383721, 0.597403, 0.694268, 0.588235,
                                             0.481013, 0.555556, 0.5     , 0.37931 , 0.259843, 0.681259,
                                             0.545455, 0.550725, 0.338583, 0.3085  , 0.437318, 0.540123,
                                             0.464481, 0.597598, 0.481928, 0.488889, 0.565543, 0.54549 ,
                                             0.428571, 0.34891 , 0.62069 , 0.656428, 0.480392, 0.738295,
                                             0.475   , 0.409574, 0.5625  , 0.582677, 0.640572, 0.307263,
                                             0.414634, 0.641566, 0.486188, 0.679245, 0.566845, 0.53528 ,
                                             0.48318 , 0.422096, 0.362694, 0.261866, 0.207268, 0.756184,
                                             0.431373, 0.781818, 0.479079, 0.505464, 0.318792, 0.209524,
                                             0.55268 , 0.666667, 0.365482, 0.745455, 0.662791, 0.670455,
                                             0.611511, 0.562197, 0.737557, 0.315315, 0.520661, 0.35514 ,
                                             0.476489, 0.310502, 0.511416, 0.642857, 0.378277, 0.44898 ,
                                             0.665635, 0.494518, 0.487179, 0.614568, 0.528736, 0.405263,
                                             0.275556, 0.504808, 0.264591, 0.518817, 0.587156, 0.652278,
                                             0.551429, 0.497797, 0.54646 , 0.546392, 0.409938, 0.447531,
                                             0.438247, 0.363636, 0.637631, 0.57037 , 0.276042, 0.554479,
                                             0.56    , 0.367521, 0.348039, 0.259259, 0.210818, 0.485477,
                                             0.285714, 0.39801 , 0.546921, 0.603306, 0.521739, 0.409812,
                                             0.509729, 0.349593, 0.54023 , 0.487805, 0.516468, 0.606383,
                                             0.746606, 0.549133, 0.394265, 0.776722, 0.380762, 0.661348,
                                             0.493056, 0.293147, 0.566265, 0.241573, 0.257956, 0.577519,
                                             0.563506, 0.694981, 0.550725, 0.48169 , 0.654391, 0.505208,
                                             0.241158, 0.530488, 0.223785, 0.689655, 0.507246, 0.421731,
                                             0.449541, 0.651007, 0.272727, 0.642412, 0.398098, 0.281588,
                                             0.429379, 0.647059, 0.716667, 0.668966, 0.727811, 0.608355,
                                             0.63785 , 0.428571, 0.647376, 0.549223, 0.409091, 0.479365,
                                             0.317333, 0.49373 , 0.706186, 0.357287, 0.534884, 0.603306,
                                             0.630435, 0.506944, 0.653846, 0.385714, 0.362126, 0.439252,
                                             0.61    , 0.367742, 0.474638, 0.37247 , 0.819549, 0.65762 ,
                                             0.66416 , 0.479532, 0.532338, 0.379562, 0.205524, 0.706564,
                                             0.466667, 0.555556, 0.70801 , 0.719397, 0.333333, 0.791667,
                                             0.831707, 0.537081, 0.434343, 0.349398, 0.367537, 0.371939,
                                             0.78392 , 0.343612, 0.541667, 0.631016, 0.390625, 0.7     ,
                                             0.566524, 0.305   , 0.472868, 0.693182, 0.413534, 0.653696,
                                             0.220264, 0.397959, 0.371069, 0.214674, 0.495114, 0.719298,
                                             0.351955, 0.213483, 0.289796, 0.515464, 0.30599 , 0.574257,
                                             0.617883, 0.376226, 0.408582, 0.339426, 0.671815, 0.333333,
                                             0.364179, 0.541126, 0.633188, 0.506045, 0.702512, 0.355014,
                                             0.494845, 0.66    , 0.55814 , 0.375   , 0.238889, 0.594444,
                                             0.597701, 0.37234 };
const std::vector<float> default_neg_data = {0.029661  , 0.00869565, 0.0830769 , 0.0914513 , 0.00950119,
                                             0.144928  , 0.0286494 , 0.0164141 , 0.0425532 , 0.032345  ,
                                             0.0186514 , 0.0281899 , 0.        , 0.0379651 , 0.0427928 ,
                                             0.0571696 , 0.0196319 , 0.00830484, 0.0254545 , 0.0104167 ,
                                             0.0192781 , 0.0434917 , 0.0165426 , 0.0256917 , 0.0288875 ,
                                             0.0231214 , 0.0139555 , 0.0181997 , 0.0180437 , 0.0222222 ,
                                             0.0177083 , 0.0252525 , 0.0381125 , 0.0156128 , 0.048603  ,
                                             0.0373519 , 0.0379747 , 0.0123748 , 0.0380165 , 0.042953  ,
                                             0.0284091 , 0.0450205 , 0.014876  , 0.0477674 , 0.0108873 ,
                                             0.0485075 , 0.0283688 , 0.0397799 , 0.037058  , 0.0195178 ,
                                             0.0285344 , 0.032835  , 0.0380864 , 0.0271024 , 0.        ,
                                             0.0509554 , 0.0440744 , 0.0226349 , 0.0388787 , 0.0307692 ,
                                             0.0197628 , 0.0183486 , 0.0220588 , 0.0194175 , 0.00157233,
                                             0.0150637 , 0.0372737 , 0.0365034 , 0.0361446 , 0.0120919 ,
                                             0.0185449 , 0.0505495 , 0.0140728 , 0.039557  , 0.0632184 ,
                                             0.0312891 , 0.0134626 , 0.0209111 , 0.0341297 , 0.0286807 ,
                                             0.0166852 , 0.0230701 , 0.0101138 , 0.0255102 , 0.0425532 ,
                                             0.0276442 , 0.0402415 , 0.0368196 , 0.0164688 , 0.0488798 ,
                                             0.0283768 , 0.0268714 , 0.0496157 , 0.0496806 , 0.0158311 ,
                                             0.0191657 , 0.0465511 , 0.0460564 , 0.0228311 , 0.0272793 ,
                                             0.0475241 , 0.0167959 , 0.0221843 , 0.0322767 , 0.0318091 ,
                                             0.0171265 , 0.0280222 , 0.0237781 , 0.0108481 , 0.0288136 ,
                                             0.0408627 , 0.0212766 , 0.0322148 , 0.0272989 , 0.0347826 ,
                                             0.00402901, 0.0259009 , 0.0351459 , 0.00761905, 0.0261438 ,
                                             0.0463215 , 0.0265372 , 0.0582691 , 0.0325203 , 0.0304878 ,
                                             0.0283688 , 0.0220264 , 0.00813638, 0.0259434 , 0.0219928 ,
                                             0.0142292 , 0.0324519 , 0.0206349 , 0.0231729 , 0.0179541 ,
                                             0.0447154 , 0.0357436 , 0.0221402 , 0.0304124 , 0.0559384 ,
                                             0.00854701, 0.0185567 , 0.0107759 , 0.016129  , 0.0297824 ,
                                             0.0211216 , 0.0352645 , 0.046285  , 0.0254181 , 0.0440367 ,
                                             0.0274964 , 0.0229277 , 0.0542232 , 0.00228833, 0.0342029 ,
                                             0.0397839 , 0.0209068 , 0.0293501 , 0.0262069 , 0.0157853 ,
                                             0.0365854 , 0.026362  , 0.0125953 , 0.0194384 , 0.0250284 ,
                                             0.0417973 , 0.0236829 , 0.00610583, 0.0164331 , 0.026455  ,
                                             0.047619  , 0.0228571 , 0.0458221 , 0.0965909 , 0.0348162 ,
                                             0.0821918 , 0.0198265 , 0.026455  , 0.0153846 , 0.0615977 ,
                                             0.0296296 , 0.034555  , 0.0374805 , 0.0253921 , 0.0186992 ,
                                             0.0233236 , 0.0119363 , 0.0366881 , 0.00961538, 0.0453846 ,
                                             0.0459094 , 0.0391798 , 0.0248075 , 0.0330113 , 0.00591716,
                                             0.0397306 , 0.0323607 , 0.0278035 , 0.0193237 , 0.0242152 ,
                                             0.00440529, 0.0294826 , 0.0188679 , 0.0543933 , 0.0327363 ,
                                             0.0365233 , 0.0515856 , 0.042086  , 0.012907  , 0.0109664 ,
                                             0.0351605 , 0.0251852 , 0.0494148 , 0.0401606 , 0.        ,
                                             0.0369686 , 0.0411622 , 0.0257143 , 0.0416667 , 0.0400943 ,
                                             0.0409836 , 0.00296296, 0.0714286 , 0.0270856 , 0.0117521 ,
                                             0.0174472 , 0.0354167 , 0.0237056 , 0.0292683 , 0.0253723 ,
                                             0.025     , 0.0163043 , 0.0229983 , 0.0418502 , 0.0259179 ,
                                             0.0386047 , 0.0355295 , 0.0451613 , 0.0353982 , 0.0310711 ,
                                             0.0250344 , 0.0125604 , 0.0237691 , 0.0406835 , 0.0524476 ,
                                             0.0307692 , 0.0248619 , 0.0264085 , 0.03003   , 0.0219273 ,
                                             0.0214827 , 0.039629  , 0.0393462 , 0.0207164 , 0.0193896 ,
                                             0.0312947 , 0.0485252 , 0.0265004 , 0.0298322 , 0.0178571 ,
                                             0.0333333 , 0.00632911, 0.0328767 , 0.0355556 , 0.0331825 ,
                                             0.0197568 , 0.0413793 , 0.0215983 , 0.0182584 , 0.0211161 ,
                                             0.0290276 , 0.00200803, 0.035284  , 0.0368509 , 0.0279971 ,
                                             0.02334   , 0.0245637 , 0.00180505, 0.0298013 , 0.0173069 ,
                                             0.0167254 , 0.0321569 , 0.0357616 , 0.0438791 , 0.024147  ,
                                             0.0399344 , 0.0285962 , 0.00518623, 0.0311218 , 0.0288498 ,
                                             0.0236686 , 0.012583  , 0.0243531 , 0.0218015 , 0.0121639 ,
                                             0.041868  , 0.0243223 , 0.0260317 , 0.0256233 , 0.0325203 ,
                                             0.024186  , 0.0125116 , 0.0138889 , 0.0221271 , 0.0445293 ,
                                             0.0166667 , 0.0292373 , 0.00485437, 0.0443599 , 0.0174766 ,
                                             0.023511  , 0.0411594 , 0.0189873 , 0.0137931 , 0.0388889 ,
                                             0.0490798 , 0.0394737 , 0.0539568 , 0.0342466 , 0.0679612 ,
                                             0.0250696 , 0.03125   , 0.0365034 , 0.0265252 , 0.060241  ,
                                             0.0166008 , 0.0295602 , 0.0162665 , 0.0100925 , 0.0506757 ,
                                             0.0150754 , 0.0737705 , 0.0371367 , 0.026362  , 0.0289958 ,
                                             0.0297619 , 0.0301003 , 0.0219895 , 0.0134804 , 0.0175439 ,
                                             0.0443319 , 0.0137457 , 0.0286561 , 0.0546875 , 0.0655738 ,
                                             0.0326741 , 0.0251232 , 0.0363743 , 0.0313901 , 0.0442404 ,
                                             0.0262935 , 0.019601  , 0.0291751 , 0.0316978 , 0.0287356 ,
                                             0.0226586 , 0.0341207 , 0.00504356, 0.0359589 , 0.0232751 ,
                                             0.0324374 , 0.0125962 , 0.0216401 , 0.0282589 , 0.0180505 ,
                                             0.0345912 , 0.0228188 , 0.0260374 , 0.0195787 , 0.00486618,
                                             0.0508732 , 0.0280469 , 0.0224719 , 0.0320856 , 0.0444444 ,
                                             0.0365535 , 0.018797  , 0.0759494 , 0.072     , 0.0235199 ,
                                             0.0193548 , 0.0453202 , 0.0414972 , 0.019802  , 0.019802  ,
                                             0.0317003 , 0.0341047 , 0.025544  , 0.0518672 , 0.0271859 ,
                                             0.0528355 , 0.02398   , 0.030303  , 0.0219146 , 0.0397456 ,
                                             0.010582  , 0.0218182 , 0.010846  , 0.0348174 , 0.032967  ,
                                             0.011117  , 0.0404494 , 0.0176829 , 0.052     , 0.0307692 ,
                                             0.0404355 , 0.0392216 , 0.00652667, 0.0440465 , 0.0325444 ,
                                             0.0322581 , 0.0220096 , 0.030146  , 0.0356802 , 0.0349938 ,
                                             0.0285036 , 0.0176194 , 0.0571429 , 0.0108509 , 0.0279416 ,
                                             0.0334225 , 0.0136482 , 0.0535294 , 0.0213592 , 0.027933  ,
                                             0.0296296 , 0.0430925 , 0.0200236 , 0.0278884 , 0.0463054 ,
                                             0.0368098 , 0.0205479 , 0.027833  , 0.0247748 , 0.0292398 ,
                                             0.0149489 , 0.0326087 , 0.024183  , 0.0265487 , 0.0101597 ,
                                             0.028125  , 0.0333333 , 0.024296  , 0.0395074 , 0.0255272 ,
                                             0.0282575 , 0.0175867 , 0.0225443 , 0.0253317 , 0.0226131 ,
                                             0.0478723 , 0.0238837 , 0.0229253 , 0.0208153 , 0.0377277 ,
                                             0.0350765 , 0.0526316 , 0.0409556 , 0.0325581 , 0.0384906 ,
                                             0.0438909 , 0.0250765 , 0.0383352 , 0.021057  , 0.0140392 ,
                                             0.0248447 , 0.0277778 , 0.0128023 , 0.0438523 , 0.0338753 ,
                                             0.0164701 , 0.0372493 , 0.0367422 , 0.0209832 , 0.0272031 ,
                                             0.0305114 , 0.0208228 , 0.0246575 , 0.0269278 , 0.034039  ,
                                             0.0348558 , 0.0205266 , 0.0447761 , 0.0470446 , 0.028169  ,
                                             0.0107239 , 0.021978  , 0.00425532, 0.0169492 , 0.0470163 ,
                                             0.0326693 , 0.0257069 , 0.0371163 , 0.0258519 , 0.0670444 ,
                                             0.0401786 , 0.00708215, 0.0335475 , 0.030137  , 0.0210084};

void setup_dehost_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<DehostArguments>();
    auto* dehost_subcommand = app.add_subcommand(
        "dehost", "Dehost read file into host and other using index.");

    dehost_subcommand->add_option("<fastaq>", opt->read_file, "Fasta/q file")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    dehost_subcommand->add_option("<fastaq>", opt->read_file2, "Paired Fasta/q file")
            ->transform(make_absolute)
            ->check(CLI::ExistingFile.description(""))
            ->type_name("FILE");

    dehost_subcommand->add_option("--db", opt->db, "Prefix for the index.")
            ->type_name("FILE")
            ->required()
            ->check(CLI::ExistingPath.description(""));

    dehost_subcommand->add_option("-e,--extract", opt->category_to_extract, "Reads from this category in the index will be extracted to file.")
            ->type_name("STRING");

    dehost_subcommand->add_option("-p,--prefix", opt->prefix, "Prefix for the output files.")
            ->type_name("FILE")
            ->check(CLI::NonexistentPath.description(""))
            ->default_str("<prefix>");

    dehost_subcommand
            ->add_option("--chunk_size", opt->chunk_size, "Read file is read in chunks of this size, to be processed in parallel within a chunk.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--lo_hi_threshold", opt->lo_hi_threshold, "Threshold used during model fitting stage to decide if read should be used to train lo or hi distribution.")
            ->type_name("FLOAT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--num_reads_to_fit", opt->num_reads_to_fit, "Number of reads to use to train each distribution in the model.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand->add_option("-d,--dist", opt->dist, "Probability distribution to use for modelling.")
            ->type_name("STRING");

    dehost_subcommand
            ->add_option("--min_length", opt->min_length, "Minimum read length to classify.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--min_quality", opt->min_quality, "Minimum read quality to classify.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--min_compression", opt->min_compression, "Minimum read gzip compression ratio to classify (a measure of how much information is in the read.")
            ->type_name("FLOAT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--confidence", opt->confidence_threshold, "Minimum difference between the top 2 unique hit counts.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--host_unique_prop_lo_threshold", opt->host_unique_prop_lo_threshold, "Require non-host reads to have unique host proportion below this threshold for classification.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--min_proportion_diff", opt->min_proportion_difference, "Minimum difference between the proportion of (non-unique) kmers found in each category.")
            ->type_name("FLOAT")
            ->capture_default_str();
    dehost_subcommand
            ->add_option("--min_probability_diff", opt->min_prob_difference, "Minimum difference between the probability found in each category.")
            ->type_name("FLOAT")
            ->capture_default_str();

    dehost_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    dehost_subcommand
            ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    dehost_subcommand->callback([opt]() { dehost_main(*opt); });
}

void dehost_reads(const DehostArguments& opt, const Index& index){
    PLOG_INFO << "Dehosting file " << opt.read_file;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<my_traits> fin{opt.read_file};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    using outfile_field_ids = decltype(fin)::field_ids;
    using outfile_format = decltype(fin)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    for (auto && chunk : fin | seqan3::views::chunk(opt.chunk_size))
    {
        // You can use a for loop:
        for (auto & record : chunk)
        {
            records.push_back(std::move(record));
        }

#pragma omp parallel for firstprivate(agent, hash_adaptor) num_threads(opt.threads) shared(result)
        for (auto i=0; i<records.size(); ++i){

            const record_type & record = records[i];
            const auto read_id = split(record.id(), " ")[0];
            const uint32_t read_length = std::ranges::size(record.sequence());
            if (read_length > std::numeric_limits<uint32_t>::max()){
                PLOG_WARNING << "Ignoring read " << record.id() << " as too long!";
                continue;
            }
            if (read_length == 0){
                PLOG_WARNING << "Ignoring read " << record.id() << " as has zero length!";
                continue;
            }
            auto qualities = record.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
            float mean_quality = 0;
            if (std::ranges::size(qualities) > 0)
                mean_quality = static_cast< float >( sum )/ static_cast< float >(std::ranges::size(qualities));
            PLOG_VERBOSE << "Mean quality of read  " << record.id() << " is " << mean_quality;

            float compression_ratio = get_compression_ratio(sequence_to_string(record.sequence()));
            PLOG_VERBOSE << "Found compression ratio of read  " << record.id() << " is " << compression_ratio;

            auto read = ReadEntry(read_id, read_length, mean_quality, compression_ratio, result.input_summary());
            for (auto && value : record.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            PLOG_VERBOSE << "Finished adding raw hash counts for read " << read_id;

            read.post_process(result.input_summary());
#pragma omp critical(add_read_to_results)
            result.add_read(read, record, true);
        }
        records.clear();
    }
    result.complete(true);
    result.print_summary();
}


void dehost_paired_reads(const DehostArguments& opt, const Index& index){
    PLOG_INFO << "Dehosting files " << opt.read_file << " and " << opt.read_file2;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<my_traits> fin1{opt.read_file};
    seqan3::sequence_file_input<my_traits> fin2{opt.read_file2};
    using record_type = decltype(fin1)::record_type;
    std::vector<record_type> records1{};
    std::vector<record_type> records2{};

    using outfile_field_ids = decltype(fin1)::field_ids;
    using outfile_format = decltype(fin1)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    for (auto && chunk : fin1 | seqan3::views::chunk(opt.chunk_size))
    {
        for (auto & record : chunk)
        {
            records1.push_back(std::move(record));
        }

        // loop in the second file and get same amount of reads
        for ( auto& record2 : fin2 | std::views::take( opt.chunk_size ) )
        {
            records2.push_back( std::move( record2 ));
        }

#pragma omp parallel for firstprivate(agent, hash_adaptor) num_threads(opt.threads) shared(result)
        for (auto i=0; i<records1.size(); ++i){

            const auto & record1 = records1[i];
            const auto & record2 = records2[i];

            auto id1 = record1.id();
            id1.erase(id1.size() - 1);
            auto id2 = record2.id();
            id2.erase(id2.size() - 1);
            if (id1 != id2){
                std::cout << id1 << " " << id2;
                throw std::runtime_error("Your pairs don't match for read ids.");
            }
            const auto read_id = split(record1.id(), " ")[0];
            const uint32_t read_length = std::ranges::size(record1.sequence()) + std::ranges::size(record2.sequence());
            if (read_length > std::numeric_limits<uint32_t>::max()){
                PLOG_WARNING << "Ignoring read " << record1.id() << " as too long!";
                continue;
            }
            if (read_length == 0){
                PLOG_WARNING << "Ignoring read " << record1.id() << " as has zero length!";
                continue;
            }
            auto qualities1 = record1.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto qualities2 = record2.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto sum = std::accumulate(qualities1.begin(), qualities1.end(), 0);
            sum = std::accumulate(qualities2.begin(), qualities2.end(), sum);
            float mean_quality = 0;
            if (std::ranges::size(qualities1) + std::ranges::size(qualities2) > 0)
                mean_quality = static_cast< float >( sum )/ static_cast< float >(std::ranges::size(qualities1) + std::ranges::size(qualities2));
            PLOG_VERBOSE << "Mean quality of read  " << record1.id() << " is " << mean_quality;

            auto combined_record = sequence_to_string(record1.sequence()) + sequence_to_string(record2.sequence());
            float compression_ratio = get_compression_ratio(combined_record);
            PLOG_VERBOSE << "Found compression ratio of read  " << record1.id() << " is " << compression_ratio;

            auto read = ReadEntry(read_id, read_length, mean_quality, compression_ratio, result.input_summary());
            for (auto && value : record1.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            for (auto && value : record2.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            PLOG_VERBOSE << "Finished adding raw hash counts for read " << read_id;

            read.post_process(result.input_summary());
#pragma omp critical(add_read_to_results)
            result.add_paired_read(read, record1, record2);
        }
        records1.clear();
        records2.clear();
    }
    result.complete();
    result.print_summary();
}


int dehost_main(DehostArguments & opt)
{
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, opt.log_file.c_str(), 10000000, 5);

    if (!ends_with(opt.db, ".idx")) {
        opt.db += ".idx";
    }

    if (opt.read_file2 != "") {
        opt.is_paired = true;
        opt.min_length = 80;
    }

    auto args = opt.to_string();
    LOG_INFO << "Running charon dehost\n\nCharon version: " << SOFTWARE_VERSION << "\n" << args;

    auto index = Index();
    load_index(index, opt.db);
    auto host_index = index.get_host_index();
    LOG_INFO << "Found host at index " << +host_index << " in the index categories";

    opt.run_extract = (opt.category_to_extract != "");
    const auto categories = index.categories();
    if (opt.run_extract and opt.category_to_extract != "all" and std::find(categories.begin(), categories.end(), opt.category_to_extract) == categories.end())
    {
        std::string options = "";
        for (auto i: categories)
            options += i + " ";
        PLOG_ERROR << "Cannot extract " << opt.category_to_extract << ", please chose one of [ all " << options << "]";
        return 1;
    } else if (opt.run_extract){
        if (opt.prefix == "")
            opt.prefix = "charon";
        std::vector<std::string> to_extract;
        if (opt.category_to_extract == "all")
            to_extract = categories;
        else
            to_extract.push_back(opt.category_to_extract);

        const auto extension = get_extension(opt.read_file);
        for (const auto & category : to_extract){
            const auto category_index = index.get_category_index(category);
            if (opt.is_paired) {
                opt.extract_category_to_file[category_index].push_back(opt.prefix + "_" + category + "_1" + extension + ".gz");
                opt.extract_category_to_file[category_index].push_back(opt.prefix + "_" + category + "_2" + extension + ".gz");
            } else {
                opt.extract_category_to_file[category_index].push_back(opt.prefix + "_" + category + extension + ".gz");
            }
        }
    }

    if (opt.dist != "gamma" and opt.dist != "beta" and opt.dist != "kde")
    {
        PLOG_ERROR << "Supported distributions are [gamma , beta, kde]";
        return 1;
    }


    if (opt.is_paired)
        dehost_paired_reads(opt, index);
    else
        dehost_reads(opt, index);

    return 0;
}