from ase.data import s22
from ase.embedding.multiase.lammps import COMPASS, CHARMM, OPLSAA, ReaxFF
from ase.embedding.multiase.lammps.amber import AMBER
from ase.embedding.multiase.lammps.pcff import PCFF
from ase.embedding.multiase.utils import get_datafile
import unittest
import numpy as np

def floatEqual(val1, val2):
	return np.max(np.abs((val1-val2)/val1)) < 2e-5

class PhenolDimerForces(unittest.TestCase):
	def test_compass(self):
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = COMPASS(get_datafile('compass.frc'), debug=True)
		
		correct_potential_energy = -0.8709613598457432
		correct_forces = np.array([[ 0.35112694,  0.25463791,  0.05644749],
			[ 0.72228116,  0.0895421 ,  0.40237507],
			[-0.57104885, -0.14180625, -0.37549366],
			[-0.58755757, -0.07901503, -0.30178597],
			[-0.0170888 ,  0.45440809, -0.28923813],
			[ 0.21728754,  0.47817162, -0.15168806],
			[ 0.35739913,  0.14291464,  0.12856372],
			[ 0.08421959, -0.55165643,  0.39717225],
			[ 0.50162293,  0.03809237,  0.28184412],
			[ 0.21723637, -0.42526915,  0.39161687],
			[-0.30955074, -0.52814441,  0.1298048 ],
			[-0.54157861, -0.1383445 , -0.25536383],
			[-0.33960857,  0.35724171, -0.43102226],
			[ 0.34620208, -0.02740724,  0.45917381],
			[-0.49434643,  1.08521701,  0.49173591],
			[-0.16693748, -0.33332728, -0.54045982],
			[ 0.62644216, -0.50698273,  0.1855112 ],
			[ 0.02029852,  0.28945408,  0.2759518 ],
			[-0.24916276,  0.3939607 ,  0.08951348],
			[-0.43254608,  0.37100895, -0.1199759 ],
			[-0.07833031, -0.40428916, -0.45046629],
			[-0.49351384,  0.10938595, -0.42777082],
			[-0.1822294 , -0.3577326 , -0.50842675],
			[ 0.29216867, -0.54694709, -0.17945497],
			[ 0.49935065, -0.24628122,  0.30982784],
			[ 0.22785798,  0.22316858,  0.43161028]])

		self.assertTrue(floatEqual(atoms.get_potential_energy(), correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(correct_potential_energy))
		self.assertTrue(floatEqual(atoms.get_forces(), correct_forces),
			'Force test failed. Output is %s' % repr(correct_forces))

	def test_pcff(self):
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = PCFF(get_datafile('pcff.frc'), debug=True)
		
		correct_potential_energy = -0.6984590910167954
		correct_forces = np.array([[ 0.38203557,  0.31769929,  0.04481593],
			[ 0.06645072,  0.103403  ,  0.08574124],
			[ 0.00229426, -0.10892152, -0.09467121],
			[ 0.22153722,  0.10880227,  0.09930943],
			[ 0.3400227 , -0.16116659,  0.31053077],
			[-0.22498033, -0.3186932 ,  0.06136411],
			[-0.46089536, -0.00924397, -0.28932139],
			[-0.28527205,  0.0557619 , -0.20593482],
			[-0.23553429, -0.14936895, -0.05274159],
			[-0.11395132,  0.11342922, -0.14521467],
			[ 0.09841353,  0.17775813, -0.04823693],
			[ 0.18713908,  0.01462723,  0.10630103],
			[-0.06184675, -0.2062761 ,  0.08826763],
			[ 0.16000573, -0.27071082, -0.02609756],
			[-0.50141911,  1.13370675,  0.51524359],
			[ 0.1432776 , -0.02164376,  0.01262255],
			[-0.07181789, -0.21984863, -0.30902257],
			[-0.29430696, -0.16109547, -0.45144199],
			[ 0.19343122, -0.37887476, -0.12902685],
			[ 0.27851158,  0.06016552,  0.35355663],
			[ 0.1822294 ,  0.08455133,  0.26004845],
			[ 0.16078021, -0.20149043, -0.02997253],
			[ 0.08616924,  0.0722537 ,  0.1520211 ],
			[-0.09563302,  0.1601445 ,  0.03820269],
			[-0.14750256,  0.01951328, -0.13697983],
			[-0.00913885, -0.21448015, -0.20936275]])
			
		potential_energy = atoms.get_potential_energy()
		forces = atoms.get_forces()
		
		self.assertTrue(floatEqual(potential_energy, correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(potential_energy))
		self.assertTrue(floatEqual(forces, correct_forces),
			'Force test failed. Output is\n%s' % repr(forces))
		
	def test_charmm(self):
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = CHARMM(get_datafile('par_all36_cgenff.prm'), auto_charges=True, debug=True)
			
		correct_potential_energy = 0.13048314753274043
		correct_forces = np.array([[ -1.14463884e+00,  -1.68207618e+00,   3.29438386e-01],
			[  7.32632172e-01,   1.23258558e+00,  -2.10086066e-01],
			[  6.14478002e-02,  -3.51428322e-02,  -1.11381564e-01],
			[ -3.83892857e-01,   8.85109028e-02,  -2.69391246e-01],
			[ -2.85321917e-02,   2.78799956e-01,  -1.93727826e-01],
			[  2.50792383e-01,   5.02099931e-01,  -1.46366854e-01],
			[  2.67473252e-01,   5.31778323e-02,   1.19425171e-01],
			[  4.38307000e-01,  -2.35313133e-01,   4.25593946e-01],
			[ -5.97756804e-02,  -1.48735402e-01,   5.85194224e-02],
			[ -4.37864686e-02,   5.72293404e-03,  -3.22041337e-02],
			[ -1.94106828e-02,   1.09826094e-03,  -5.47744985e-03],
			[  2.54398108e-02,  -1.46024278e-02,   2.67202660e-02],
			[ -1.74992799e-01,  -8.94137434e-02,  -4.01153071e-02],
			[ -4.07283454e-01,   9.06656651e-01,   5.03734758e-01],
			[  4.46953802e-01,  -1.07517389e+00,  -5.36400936e-01],
			[  1.89104779e-02,   1.54977664e-02,  -6.15510068e-02],
			[  5.78260304e-01,  -1.36903073e-01,   4.56580633e-01],
			[  4.72586323e-02,   1.69949554e-01,   2.11128106e-01],
			[ -2.75001694e-01,   4.50292838e-01,   1.19754739e-01],
			[ -3.09711622e-01,   2.23511159e-01,  -1.03908195e-01],
			[ -1.82374237e-01,  -1.64251511e-01,  -3.31384567e-01],
			[ -1.18465957e-03,  -1.46414989e-01,  -1.39178821e-01],
			[  2.86580776e-02,  -1.30676421e-02,   1.20760786e-02],
			[  2.21930104e-02,  -1.87663790e-02,  -1.36623808e-02],
			[  2.58436173e-04,  -3.01109919e-02,  -3.58351835e-02], 
			[  1.11995166e-01,  -1.37929935e-01,  -3.22990144e-02]])

		potential_energy = atoms.get_potential_energy()
		forces = atoms.get_forces()
		
		self.assertTrue(floatEqual(potential_energy, correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(potential_energy))
		self.assertTrue(floatEqual(forces, correct_forces),
			'Force test failed. Output is\n%s' % repr(forces))

	def test_oplsaa(self):
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = OPLSAA(get_datafile('gromacs_top'), debug=True)
		
		correct_potential_energy = 0.08010372652676295
		correct_forces = np.array([[ -5.44488341e-01,   1.61736393e-01,  -4.22837290e-01],
       [  1.07106731e+00,  -2.81245691e-01,   8.18176536e-01],
       [ -4.61021117e-01,   2.37460090e-01,  -5.56734364e-01],
       [  3.23237319e-01,   3.13988190e-01,   3.61013523e-02],
       [  2.73089337e-01,  -2.31379142e-01,   3.05502702e-01],
       [ -1.23124997e-01,  -1.06339186e-01,  -6.92923663e-03],
       [ -4.34486623e-01,  -1.55414774e-01,  -1.83176906e-01],
       [  1.96855245e-01,   2.46456406e-01,  -1.87242725e-02],
       [ -5.08778003e-02,  -1.42329223e-01,   4.85747328e-02],
       [ -5.24731656e-02,   2.74903691e-02,  -4.94159964e-02],
       [ -7.73945152e-04,   5.27723779e-03,  -3.43772390e-03],
       [  6.03931853e-02,  -1.52109996e-03,   3.70595689e-02],
       [ -1.43192602e-01,  -1.03796749e-01,  -2.51707365e-02],
       [  9.20173241e-01,  -2.34231199e-01,   9.53381471e-01],
       [ -9.38121643e-01,   7.61985333e-01,  -3.16286054e-01],
       [ -5.03288108e-01,  -5.84227205e-02,  -6.35266753e-01],
       [  5.44670470e-02,   6.22665145e-02,   7.99590682e-02],
       [ -1.74277292e-01,  -2.47938591e-01,  -3.90135987e-01],
       [  6.65903491e-02,  -1.76886943e-01,  -8.28223999e-02],
       [  2.89269784e-01,  -5.22077773e-03,   3.06853928e-01],
       [  8.80763945e-02,   2.33088121e-01,   3.00527105e-01],
       [ -8.05648647e-03,  -1.38594273e-01,  -1.27830435e-01],
       [  3.24510056e-02,   9.33698505e-03,   3.42864346e-02],
       [ -1.20272506e-03,  -6.85846642e-03,  -1.08188232e-02],
       [ -3.25913752e-02,  -1.90343692e-02,  -4.83995419e-02],
       [  9.22991908e-02,  -1.49871108e-01,  -4.24363273e-02]])
		
		potential_energy = atoms.get_potential_energy()
		forces = atoms.get_forces()
		
		self.assertTrue(floatEqual(potential_energy, correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(potential_energy))
		self.assertTrue(floatEqual(forces, correct_forces),
			'Force test failed. Output is\n%s' % repr(forces))

	def test_amber(self):
		""" No Coulomb interactions! """
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = AMBER(get_datafile('gromacs_top'), debug=True)
		
		correct_potential_energy = 0.3324552179661166
		correct_forces = np.array([[-0.49542186,  0.37058441, -0.51020468],
       [ 0.72921508, -0.31332949,  0.63050104],
       [-0.35354753,  0.27442235, -0.31369939],
       [ 0.63123823,  0.40587369,  0.14572593],
       [ 0.4113909 , -0.46782494,  0.53476611],
       [-0.24513423, -0.31201122,  0.04369237],
       [-0.76459585, -0.21723681, -0.35159658],
       [ 0.15155233,  0.50109822, -0.20145834],
       [-0.0932519 , -0.08532494,  0.02256945],
       [-0.06985697, -0.00636169, -0.03790135],
       [ 0.00617748,  0.01371325, -0.00475019],
       [ 0.03777204, -0.0378543 ,  0.0452886 ],
       [-0.18549298, -0.14524676, -0.02796516],
       [ 1.18299873, -0.4691172 ,  0.53836099],
       [-1.15412257,  1.00876177, -0.28264548],
       [-0.31591095,  0.18178232, -0.12320002],
       [-0.18158674,  0.19132199, -0.05114579],
       [-0.30875588, -0.43867993, -0.7026459 ],
       [ 0.19808895, -0.41808979, -0.16460493],
       [ 0.56060678, -0.13405752,  0.48329726],
       [ 0.19004058,  0.39510942,  0.53258923],
       [ 0.04055402, -0.15973731, -0.11038506],
       [ 0.05862436, -0.02846233,  0.03145194],
       [-0.00785263,  0.00732359, -0.0017853 ],
       [-0.01925375, -0.04043039, -0.05409628],
       [-0.00347931, -0.07622802, -0.07016138]])

		potential_energy = atoms.get_potential_energy()
		forces = atoms.get_forces()
		
		self.assertTrue(floatEqual(potential_energy, correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(potential_energy))
		self.assertTrue(floatEqual(forces, correct_forces),
			'Force test failed. Output is\n%s' % repr(forces))
	
	def test_reaxff(self):
		atoms = s22.create_s22_system('Phenol_dimer')
		atoms.calc = ReaxFF(get_datafile('ffield.reax'), debug=True)
		
		correct_potential_energy = -142.4691415172731
		correct_forces = np.array([[ 2.07723156,  3.96541299, -1.15132559],
			[-3.2265494 , -0.20475922, -1.49226283],
			[ 2.60729701,  0.3025817 ,  1.33967323],
			[ 5.04814524,  2.42645569,  1.66324749],
			[ 2.1388216 , -3.76116807,  3.61281947],
			[-2.80884034, -4.91705556,  1.26078092],
			[-5.08652247, -1.16750474, -2.46529692],
			[-0.73403717,  4.58245815, -3.24360016],
			[ 0.57589696,  0.03272494,  0.34879092],
			[ 0.22167946, -0.63391813,  0.53482248],
			[-0.39215936, -0.67218695,  0.16863779],
			[-0.75861594, -0.27159978, -0.30523254],
			[-0.45115578,  0.36083747, -0.49195707],
			[-1.93148048, -0.10359034, -2.32416844],
			[-1.29829954,  2.23213681,  0.54021264],
			[ 1.73759958,  1.86858085,  3.20233054],
			[-4.27793374,  3.30169938, -1.55869664],
			[-1.68371535, -3.23798885, -4.58176432],
			[ 2.71178281, -4.94754052, -1.52588302],
			[ 4.49239091, -2.09399613,  2.87009647],
			[ 0.39408212,  4.17231613,  4.14201763],
			[-0.52669605,  0.0963932 , -0.53510001],
			[-0.16398092, -0.54798349, -0.68089446],
			[ 0.37304923, -0.6820306 , -0.22145353],
			[ 0.68719527, -0.39223828,  0.3756307 ],
			[ 0.27481566,  0.29199609,  0.51855661]])

		potential_energy = atoms.get_potential_energy()
		forces = atoms.get_forces()
		
		self.assertTrue(floatEqual(potential_energy, correct_potential_energy),
			'Potential energy test failed. Output is %s' % repr(potential_energy))
		self.assertTrue(floatEqual(forces, correct_forces),
			'Force test failed. Output is\n%s' % repr(forces))
			
if __name__ == '__main__':
	unittest.main()
