#ifndef BICYCLE_3D_TRACK_H
# define BICYCLE_3D_TRACK_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include <slap/slap.h>

float Xref_data[300] = {0.0,
                          0.0,
                          0.0,
                          0.09999315158682012,
                          0.0010135154712834142,
                          0.02027100355086725,
                          0.19994521603621784,
                          0.00405364543211875,
                          0.0405420071017345,
                          0.2998151230934555,
                          0.00911914069458163,
                          0.06081301065260175,
                          0.39956183626222797,
                          0.01620791984912891,
                          0.081084014203469,
                          0.4991443696665392,
                          0.025317070119848864,
                          0.10135501775433625,
                          0.5985218048917784,
                          0.03644284856131937,
                          0.1216260213052035,
                          0.6976533077980773,
                          0.0495806835965824,
                          0.14189702485607075,
                          0.7964981452990376,
                          0.06472517689560273,
                          0.162168028406938,
                          0.8950157020989371,
                          0.08187010559343916,
                          0.18243903195780525,
                          0.9931654973815347,
                          0.10100842484721667,
                          0.2027100355086725,
                          1.090907201443618,
                          0.12213227073084898,
                          0.22298103905953975,
                          1.1882006522664592,
                          0.14523296346632192,
                          0.243252042610407,
                          1.2850058720183684,
                          0.17030101099021,
                          0.26352304616127425,
                          1.3812830834815657,
                          0.19732611285396068,
                          0.2837940497121415,
                          1.476992726396621,
                          0.2262971644563437,
                          0.30406505326300876,
                          1.5720954737177444,
                          0.25720226160632625,
                          0.324336056813876,
                          1.6665522477722514,
                          0.2900287054144994,
                          0.34460706036474326,
                          1.7603242363175584,
                          0.3247630075110454,
                          0.3648780639156105,
                          1.8533729084891155,
                          0.36139089558810267,
                          0.38514906746647776,
                          1.945660030632719,
                          0.3998973192642496,
                          0.405420071017345,
                          2.037147682014702,
                          0.4402664562686995,
                          0.42569107456821226,
                          2.1277982704035456,
                          0.4824817189426637,
                          0.4459620781190795,
                          2.2175745475165076,
                          0.5265257610552124,
                          0.46623308166994676,
                          2.3064396243249248,
                          0.572380484930833,
                          0.486504085220814,
                          2.3943569862118967,
                          0.6200270488857558,
                          0.5067750887716813,
                          2.4812905079761234,
                          0.6694458749699924,
                          0.5270460923225486,
                          2.567204468675734,
                          0.7206166570119066,
                          0.5473170958734159,
                          2.6520635663060035,
                          0.77351836896201,
                          0.5675880994242832,
                          2.73583293230493,
                          0.8281292735325554,
                          0.5878591029751505,
                          2.8184781458807096,
                          0.8844269311293783,
                          0.6081301065260178,
                          2.899965248155223,
                          0.9423882090723148,
                          0.6284011100768851,
                          2.980260756117722,
                          1.0019892911004087,
                          0.6486721136277525,
                          3.059331676382983,
                          1.0632056871580016,
                          0.6689431171786198,
                          3.137145518748272,
                          1.126012243457685,
                          0.6892141207294871,
                          3.213670309543554,
                          1.190383152815979,
                          0.7094851242803544,
                          3.2888746047694566,
                          1.2562919652574922,
                          0.7297561278312217,
                          3.362727503017598,
                          1.323711598883203,
                          0.750027131382089,
                          3.435198658167956,
                          1.3926143509983984,
                          0.7702981349329563,
                          3.506258291858077,
                          1.4629719094956992,
                          0.7905691384838236,
                          3.5758772057189905,
                          1.5347553644884901,
                          0.8108401420346909,
                          3.644026793372803,
                          1.6079352201899784,
                          0.8311111455855582,
                          3.710679052187047,
                          1.682481407032998,
                          0.8513821491364255,
                          3.7758065947809465,
                          1.758363294025581,
                          0.8716531526872928,
                          3.8393826602788814,
                          1.835549701337216,
                          0.8919241562381601,
                          3.901381125306416,
                          1.9140089131106266,
                          0.9121951597890274,
                          3.9617765147243817,
                          1.9937086904938013,
                          0.9324661633398947,
                          4.020544012096599,
                          2.074616284886922,
                          0.952737166890762,
                          4.0776594698869415,
                          2.156698451398747,
                          0.9730081704416293,
                          4.133099419381542,
                          2.239921462506923,
                          0.9932791739924967,
                          4.186841080332078,
                          2.3242511219166047,
                          1.013550177543364,
                          4.238862370316167,
                          2.409652778611698,
                          1.0338211810942313,
                          4.2891419138110205,
                          2.4960913410929453,
                          1.0540921846450986,
                          4.337659050976637,
                          2.5835312917970037,
                          1.0743631881959659,
                          4.384393846144917,
                          2.671936701690597,
                          1.0946341917468332,
                          4.429327096011225,
                          2.761271245033736,
                          1.1149051952977005,
                          4.472440337525015,
                          2.85149821430595,
                          1.1351761988485678,
                          4.513715855476286,
                          2.9425805352893857,
                          1.155447202399435,
                          4.553136689774759,
                          3.034480782302588,
                          1.1757182059503024,
                          4.590686642418769,
                          3.127161193578693,
                          1.1959892095011697,
                          4.626350284151009,
                          3.220583686781721,
                          1.216260213052037,
                          4.660112960798413,
                          3.314709874654587,
                          1.2365312166029043,
                          4.691960799293546,
                          3.40950108079241,
                          1.2568022201537716,
                          4.721880713375036,
                          3.504918355534626,
                          1.277073223704639,
                          4.74986040896472,
                          3.600922491969386,
                          1.2973442272555062,
                          4.7758883892192685,
                          3.6974740420436554,
                          1.3176152308063735,
                          4.799953959254236,
                          3.794533332772402,
                          1.3378862343572409,
                          4.822047230538584,
                          3.8920604825402054,
                          1.3581572379081082,
                          4.842159124957872,
                          3.9900154174885927,
                          1.3784282414589755,
                          4.8602813785444505,
                          4.088357887982368,
                          1.3986992450098428,
                          4.876406544873123,
                          4.1870474851481685,
                          1.41897024856071,
                          4.890527998120876,
                          4.286043657478451,
                          1.4392412521115774,
                          4.902639935789429,
                          4.385305727494083,
                          1.4595122556624447,
                          4.912737381089476,
                          4.484792908458708,
                          1.479783259213312,
                          4.920816184985648,
                          4.584464321137983,
                          1.5000542627641793,
                          4.9268730279013475,
                          4.68427901059685,
                          1.5203252663150466,
                          4.930905421082762,
                          4.784195963027885,
                          1.540596269865914,
                          4.932911707621491,
                          4.884174122603864,
                          1.5608672734167812,
                          4.932891063135371,
                          4.984172408347573,
                          1.5811382769676485,
                          4.930843496107212,
                          5.084149731011967,
                          1.6014092805185158,
                          4.92676984788131,
                          5.184065009963722,
                          1.6216802840693831,
                          4.920671792317743,
                          5.283877190063242,
                          1.6419512876202504,
                          4.912551835104582,
                          5.383545258534201,
                          1.6622222911711177,
                          4.902413312728299,
                          5.4830282618156705,
                          1.682493294721985,
                          4.890260391102817,
                          5.58228532238992,
                          1.7027642982728524,
                          4.876098063857729,
                          5.681275655578966,
                          1.7230353018237197,
                          4.85993215028642,
                          5.779958586302978,
                          1.743306305374587,
                          4.841769292954927,
                          5.878293565793651,
                          1.7635773089254543,
                          4.821616954972505,
                          5.976240188255671,
                          1.7838483124763216,
                          4.799483416925044,
                          6.07375820746944,
                          1.804119316027189,
                          4.775377773472578,
                          6.170807553328224,
                          1.8243903195780562,
                          4.749309929612293,
                          6.267348348302945,
                          1.8446613231289235,
                          4.721290596608572,
                          6.363340923827836,
                          1.8649323266797908,
                          4.691331287591737,
                          6.458745836600235,
                          1.8852033302306581,
                          4.659444312827312,
                          6.55352388478782,
                          1.9054743337815254,
                          4.6256427746577415,
                          6.647636124136623,
                          1.9257453373323927,
                          4.589940562118643,
                          6.741043883973203,
                          1.94601634088326,
                          4.5523523452318075,
                          6.833708783094411,
                          1.9662873444341273,
                          4.5128935689773,
                          6.925592745538203,
                          1.9865583479849946,
                          4.471580446947118,
                          7.0166580162290355,
                          2.006829351535862};

float Uref_data[198] = {
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2,
    1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0, 0.2, 1.0,
    0.2, 1.0, 0.2};


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef BICYCLE_3D_TRACK_H