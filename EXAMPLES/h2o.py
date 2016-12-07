"""
Example for the water molecule. Here H is the hamiltonian of the final SCF
interation while A is the canonical transformation matrix. This
calculation was run using pyQuante 1.6.5

The resuling eigenvalues refer to the molecular orbital energies while the
eigenvectors are the molecular orbitals.

USEAGE:
    python h2o.py
"""
import sys
sys.path.append('../')
import symmetricEigenSolver as sym
import numpy as np

np.set_printoptions(linewidth=500)

# Hamiltonian for H2O
H = np.array(
    [
        [-20.5511687037, -0.334947304214, 2.41069498212e-11, 0.0, 0.0236675537507, 0.0736280581338, -9.28768041396e-12, 0.0, -0.00788755257093, 0.0276330589412, 0.0, -1.29976925453e-11, 0.0354873979148, 0.0, 0.0619049772543, -0.00337416551197, 0.00315137627033, -0.00604297784234, 0.0, 0.00551623166563, -0.00345285040238, 0.00165107532799, 0.00598014176699, 0.0, 0.005307940673],
        [-0.334947304214, -0.901527501167, -4.20879733988e-11, 0.0, 0.0574467280185, -0.890012748528, 1.88356908791e-11, 0.0, 0.0450283264511, -0.208433543391, 0.0, -1.87571048793e-11, -0.117571015329, 0.0, -0.30650831233, -0.0766313549084, 0.0626188136642, 0.0546879764716, 0.0, -0.0594076130945, -0.0697178939822, 0.0629120873273, -0.0602691730481, 0.0, -0.0599222424355],
        [2.41069498212e-11, -4.20879733988e-11, -0.117564528177, 0.0, -3.00378928054e-10, 1.10224603795e-10, -0.745528716822, 0.0, -1.48209283692e-10, 1.06551446888e-10, 0.0, 0.257067607745, -5.59899422178e-11, 0.0, -1.05478673706e-10, -0.145232737763, 0.144871667957, 0.117237413951, 0.0, -0.176754979936, 0.167206766317, -0.140106920653, 0.11140397162, 0.0, 0.165667936213],
        [0.0, 0.0, 0.0, -0.165361804483, 0.0, 0.0, 0.0, -0.660848309349, 0.0, 0.0, 4.21156890864e-12, 0.0, 0.0, 0.0735875824957, 0.0, 0.0, 0.0, 0.0, -0.0541151304205, 0.0, 0.0, 0.0, 0.0, -0.0557623078024, 0.0],
        [0.0236675537507, 0.0574467280185, -3.00378928054e-10, 0.0, -0.140023349165, 0.0853300137284, 2.22029241764e-10, 0.0, -0.712602423802, 0.0490707083454, 0.0, -1.06306490247e-10, -0.0925523592099, 0.0, 0.084322930421, 0.133475310714, -0.0987809873098, -0.135836678209, 0.0, 0.104780480072, 0.123646644412, -0.101470856582, 0.149698846569, 0.0, 0.100036212966],
        [0.0736280581338, -0.890012748528, 1.10224603795e-10, 0.0, 0.0853300137284, 0.791734717987, -5.71276786963e-13, 0.0, 0.181873482133, 0.163940663844, 0.0, 1.74616341492e-11, 0.44944505451, 0.0, 0.616000409239, -0.224130669636, -0.100620489614, 0.0424033087971, 0.0, -0.0601134927613, -0.223305412676, -0.0932661808094, -0.0534520362169, 0.0, -0.0201907425537],
        [-9.28768041396e-12, 1.88356908791e-11, -0.745528716822, 0.0, 2.22029241764e-10, -5.71276786963e-13, 0.63378081335, 0.0, 7.87494976615e-11, 7.75746336078e-11, 0.0, 0.234939967916, 2.96218023625e-11, 0.0, -5.86427354111e-11, -0.159958876948, -0.279920239487, -0.137852943435, 0.0, 0.163583041463, 0.140369796325, 0.282404154457, -0.134081481646, 0.0, -0.107907449131],
        [0.0, 0.0, 0.0, -0.660848309349, 0.0, 0.0, 0.0, 0.842517888492, 0.0, 0.0, -1.39215946743e-11, 0.0, 0.0, 0.0500797113091, 0.0, 0.0, 0.0, 0.0, -0.00147271709509, 0.0, 0.0, 0.0, 0.0, -0.00151754418805, 0.0],
        [-0.00788755257093, 0.0450283264511, -1.48209283692e-10, 0.0, -0.712602423802, 0.181873482133, 7.87494976615e-11, 0.0, 0.723045774809, 0.0470953814175, 0.0, -6.24668509542e-11, -0.0900982680982, 0.0, 0.0468770828047, 0.12341738794, 0.197127538825, 0.109679900323, 0.0, -0.0953972287448, 0.13620830877, 0.202845501976, -0.107852424806, 0.0, -0.13698593361],
        [0.0276330589412, -0.208433543391, 1.06551446888e-10, 0.0, 0.0490707083454, 0.163940663844, 7.75746336078e-11, 0.0, 0.0470953814175, 2.04941014589, 0.0, 6.44442312624e-11, 0.218043023011, 0.0, 0.313257713678, -0.446927889443, 0.224953086641, -0.20789526734, 0.0, -0.2275506523, -0.40099422093, 0.208614777312, 0.222700091732, 0.0, -0.252781866854],
        [0.0, 0.0, 0.0, 4.21156890864e-12, 0.0, 0.0, 0.0, -1.39215946743e-11, 0.0, 0.0, 2.03241784147, 0.0, 0.0, 4.08266339264e-13, 0.0, 0.0, 0.0, 0.0, -0.339557881919, 0.0, 0.0, 0.0, 0.0, 0.329527592782, 0.0],
        [-1.29976925453e-11, -1.87571048793e-11, 0.257067607745, 0.0, -1.06306490247e-10, 1.74616341492e-11, 0.234939967916, 0.0, -6.24668509542e-11, 6.44442312624e-11, 0.0, 1.77307286887, 2.89033641712e-12, 0.0, 6.71168097454e-12, 0.736284740583, 0.108704045358, -0.140880711696, 0.0, -0.0507863330845, -0.735937052012, -0.0884976400444, -0.102732359263, 0.0, -0.0395672915019],
        [0.0354873979148, -0.117571015329, -5.59899422178e-11, 0.0, -0.0925523592099, 0.44944505451, 2.96218023625e-11, 0.0, -0.0900982680982, 0.218043023011, 0.0, 2.89033641712e-12, 2.3458284115, 0.0, 0.493861154689, 0.265305386353, 0.461627121142, -0.117884333014, 0.0, -0.027070796847, 0.293483674637, 0.433876727116, 0.143638358954, 0.0, -0.142425353246],
        [0.0, 0.0, 0.0, 0.0735875824957, 0.0, 0.0, 0.0, 0.0500797113091, 0.0, 0.0, 4.08266339264e-13, 0.0, 0.0, 2.06666563547, 0.0, 0.0, 0.0, 0.0, 0.251257238041, 0.0, 0.0, 0.0, 0.0, 0.258905103586, 0.0],
        [0.0619049772543, -0.30650831233, -1.05478673706e-10, 0.0, 0.084322930421, 0.616000409239, -5.86427354111e-11, 0.0, 0.0468770828047, 0.313257713678, 0.0, 6.71168097454e-12, 0.493861154689, 0.0, 2.8225229499, -0.115911152074, 0.703325788764, 0.118775024115, 0.0, 0.168454868877, -0.102551304011, 0.606217471791, -0.159810965354, 0.0, 0.102677675006],
        [-0.00337416551197, -0.0766313549084, -0.145232737763, 0.0, 0.133475310714, -0.224130669636, -0.159958876948, 0.0, 0.12341738794, -0.446927889443, 0.0, 0.736284740583, 0.265305386353, 0.0, -0.115911152074, 1.95254026734, -0.0669654481899, -0.257465831449, 0.0, 0.283968432137, -0.129916823063, 0.0106297006192, 0.0385119280047, 0.0, -0.0266532361383],
        [0.00315137627033, 0.0626188136642, 0.144871667957, 0.0, -0.0987809873098, -0.100620489614, -0.279920239487, 0.0, 0.197127538825, 0.224953086641, 0.0, 0.108704045358, 0.461627121142, 0.0, 0.703325788764, -0.0669654481899, 0.941163591108, 0.17979308133, 0.0, -0.181022891294, 0.014958123733, 0.297571137997, 0.0113930882051, 0.0, 0.0163937483419],
        [-0.00604297784234, 0.0546879764716, 0.117237413951, 0.0, -0.135836678209, 0.0424033087971, -0.137852943435, 0.0, 0.109679900323, -0.20789526734, 0.0, -0.140880711696, -0.117884333014, 0.0, 0.118775024115, -0.257465831449, 0.17979308133, 3.21219219412, 0.0, -0.440294061894, -0.0390784998201, 0.0791577728344, -0.159052818881, 0.0, 0.159064466745],
        [0.0, 0.0, 0.0, -0.0541151304205, 0.0, 0.0, 0.0, -0.00147271709509, 0.0, 0.0, -0.339557881919, 0.0, 0.0, 0.251257238041, 0.0, 0.0, 0.0, 0.0, 2.81532970395, 0.0, 0.0, 0.0, 0.0, 0.0400189336766, 0.0],
        [0.00551623166563, -0.0594076130945, -0.176754979936, 0.0, 0.104780480072, -0.0601134927613, 0.163583041463, 0.0, -0.0953972287448, -0.2275506523, 0.0, -0.0507863330845, -0.027070796847, 0.0, 0.168454868877, 0.283968432137, -0.181022891294, -0.440294061894, 0.0, 3.2665208867, -0.154990910505, -0.332862221867, -0.13686564435, 0.0, 0.137214753331],
        [-0.00345285040238, -0.0697178939822, 0.167206766317, 0.0, 0.123646644412, -0.223305412676, 0.140369796325, 0.0, 0.13620830877, -0.40099422093, 0.0, -0.735937052012, 0.293483674637, 0.0, -0.102551304011, -0.129916823063, 0.014958123733, -0.0390784998201, 0.0, -0.154990910505, 1.98129639503, -0.0330485995634, 0.281325706101, 0.0, 0.265821930149],
        [0.00165107532798, 0.0629120873273, -0.140106920653, 0.0, -0.101470856582, -0.0932661808094, 0.282404154457, 0.0, 0.202845501976, 0.208614777312, 0.0, -0.0884976400444, 0.433876727116, 0.0, 0.606217471791, 0.0106297006192, 0.297571137997, 0.0791577728344, 0.0, -0.332862221867, -0.0330485995634, 0.928143446102, -0.172026720806, 0.0, -0.246731159112],
        [0.00598014176699, -0.0602691730481, 0.11140397162, 0.0, 0.149698846569, -0.0534520362169, -0.134081481646, 0.0, -0.107852424806, 0.222700091732, 0.0, -0.102732359263, 0.143638358954, 0.0, -0.159810965354, 0.0385119280047, 0.0113930882051, -0.159052818881, 0.0, -0.13686564435, 0.281325706101, -0.172026720806, 3.28056404505, 0.0, 0.388523579568],
        [0.0, 0.0, 0.0, -0.0557623078024, 0.0, 0.0, 0.0, -0.00151754418805, 0.0, 0.0, 0.329527592782, 0.0, 0.0, 0.258905103586, 0.0, 0.0, 0.0, 0.0, 0.0400189336766, 0.0, 0.0, 0.0, 0.0, 2.81772994584, 0.0],
        [0.005307940673, -0.0599222424355, 0.165667936213, 0.0, 0.100036212966, -0.0201907425537, -0.107907449131, 0.0, -0.13698593361, -0.252781866854, 0.0, -0.0395672915019, -0.142425353246, 0.0, 0.102677675006, -0.0266532361383, 0.0163937483419, 0.159064466745, 0.0, 0.137214753331, 0.265821930149, -0.246731159112, 0.388523579568, 0.0, 3.36666336117]
    ],
    dtype=np.double
)

val, vec = sym.getEig(H)

print 'eigenvalues:'
print val

print 'eigenvectors:'
print vec

# The canonical transformation matrix
A = np.array(
    [
        [1.0, -0.240344769484, 0.0, 0.0, 9.45056202824e-18, 0.0183117251501, 0.0, 0.0, -4.76144022074e-18, 0.129457719476, 0.0, 0.0, 0.18146086361, 0.0, 0.307084507763, 0.0312912591008, 0.205273235627, 0.0147122963636, 0.0, -0.0157012043583, 0.0418332077401, 0.190246053252, -0.0150864799073, 0.0, -0.0504211709448],
        [0.0, 1.02847732509, 0.0, 0.0, 8.60875693935e-18, -1.18722981282, 0.0, 0.0, -9.0753596195e-18, -0.0777877030292, 0.0, 0.0, -0.109035010249, 0.0, -0.184518919314, 0.316357595077, 0.644185137657, 0.0486301297194, 0.0, -0.0518988732844, 0.347813187937, 0.609832273153, -0.0373878808983, 0.0, -0.185933315396],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -0.579695019023, 0.0, 0.0, 0.0, 0.0, -1.17002210696e-17, 0.0, 0.0, 0.0, -0.00463718761047, 0.230575200313, 0.343618293339, 0.0, -0.379649803743, 0.0369980137594, -0.194428968102, 0.337241345167, 0.0, 0.350503647031],
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -0.579695019023, 0.0, 0.0, 0.0, 0.0, 0.0, -1.17002210696e-17, 0.0, 0.0, 0.0, 0.0, -0.0146904610367, 0.0, 0.0, 0.0, 0.0, -0.0151376150064, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, -2.55228334971e-18, 0.0, 0.0, -0.579695019023, -1.69795185818e-17, 0.0, 0.0, -2.38001883395e-17, 0.0, -6.32755744785e-17, 0.00349636635094, -0.173850065914, -0.272508221334, 0.0, 0.273670095228, -0.0190304483647, -0.206414194923, 0.279383678732, 0.0, 0.296689629555],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.54909296201, 0.0, 0.0, 3.95456536376e-18, -0.947868171155, 0.0, 0.0, -1.32862665604, 0.0, -2.24842235705, -1.03177958868, -3.46731921042, 0.546196171417, 0.0, -0.582909526509, -1.15988756213, -3.0901734121, -0.578758279386, 0.0, 0.0655750987194],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15587469696, 0.0, 0.0, 0.0, 0.0, 2.3329490576e-17, 0.0, 0.0, 0.0, -0.855033432896, -0.916636336495, 0.381674488079, 0.0, -0.689041405203, 0.878694624371, 1.12252540318, 0.417829817809, 0.0, 0.902977243701],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15587469696, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3329490576e-17, 0.0, 0.0, 0.0, 0.0, -0.319951423596, 0.0, 0.0, 0.0, 0.0, -0.329690229535, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15587469696, 2.28383567836e-17, 0.0, 0.0, 3.20125208613e-17, 0.0, 1.0003250044e-16, 0.644681728415, 0.691129346532, -0.58017651805, 0.0, 0.245542594692, 0.675539578971, 0.609064718128, 0.65084120872, 0.0, 0.0289506909274],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.41000515404, 0.0, 0.0, 0.485239404145, 0.0, 0.82116606636, -0.335023108689, 0.634339112454, 0.208581968665, 0.0, -0.799926194515, -0.233762494849, 0.684721622527, -0.185323799683, 0.0, -0.909843662052],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.378562477752, 0.0, 0.0, 0.0, 0.0, 0.367380021687, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.599680220404, 0.0718080760009, -0.614069960996, 0.0, 0.464586384818, -0.626431288582, -0.144237290431, -0.583337415693, 0.0, -0.548769822519],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.49116458306, 0.0, 0.82116606636, 0.353768886184, 0.716817783976, -0.0993642822453, 0.0, 0.106043194271, 0.379896816447, 0.649379057008, 0.113937302376, 0.0, -0.0441049483653],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.285430140017, 0.0, 0.0, 0.0, 0.0, 0.294118173693, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.70231769136, -0.0378040081965, 0.669929300409, 0.527508580335, 0.0, 0.0143582688506, -0.0303499755317, 0.614469827224, -0.593984326397, 0.0, -0.0272327396586],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.83131320896, -0.194357108597, -0.856776771862, 0.0, 0.914366245216, -0.0864690848896, -0.0514634887782, 0.0920039404776, 0.0, -0.162171810458],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.33501688423, -0.527004956377, 0.0, 0.562428346564, 0.0827175513442, -0.370664829215, -0.0588649831391, 0.0, -0.373621818105],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.39893015463, 0.0, -0.476744124573, -0.0174011574562, 0.108608980069, -0.0941552328531, 0.0, 0.0907816749163],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1541573403, 0.0, 0.0, 0.0, 0.0, 0.034611822936, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.34779624047, -0.130189513292, -0.210292926473, -0.119377967971, 0.0, 0.093178972796],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.84405379976, -0.184409226565, 0.882095170651, 0.0, 0.890994706664],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3820108036, 0.570188449488, 0.0, 0.578417099014],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.4126840305, 0.0, 0.466711472998],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.15467620762, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.37864446118]
    ],
    dtype=np.double
)

vec_trans = A.dot(vec)

print 'eigenvectors (transformed)'
print vec_trans