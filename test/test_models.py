import unittest
from pythermo import pythermo as pt

class Testing(unittest.TestCase):
    """
        Here the PBubble function is exposed to unit tests
    """
    def test_PBubble_output(self):
        Thermo = pt.Model()
        Thermo.ChooseAModel(1)
        Thermo.NoPureComp(1)
        Thermo.CritProps(1, 647.29, 220.64000, 0.3449)
        Thermo.CPAParams(1, 14.515, 1017.3, 0.6736)
        Thermo.AssocParams(1, 22, 69.2, 2003.2)
        Thermo.Setup_Thermo()
        P, LnK, ierr = Thermo.PBubble(373.15,1) # T = 371.15K, moles = 1
        Thermo.Finishup_Thermo()
        self.assertGreater(P,0)
        self.assertAlmostEqual(P,1.01325,1)
        self.assertNotEqual(ierr,1)

    def test_PBubble_values(self):
        self.assertRaises(ValueError, pt.Model().PBubble, -10,1,1)
        self.assertRaises(ValueError, pt.Model().PBubble, 0,1,1)
        self.assertRaises(ValueError, pt.Model().PBubble, 10,1,-10)
        self.assertRaises(ValueError, pt.Model().PBubble, 10,1,0)

    def test_PBubble_types(self):
        self.assertRaises(TypeError, pt.Model().PBubble, 3+5j, 1,1)
        self.assertRaises(TypeError, pt.Model().PBubble, "some string", 1,1)
        self.assertRaises(TypeError, pt.Model().PBubble, True, 1,1)


    """
        Here the TBubble function is exposed to unit tests
    """
    def test_TBubble_output(self):
        Thermo = pt.Model()
        Thermo.ChooseAModel(1)
        Thermo.NoPureComp(1)
        Thermo.CritProps(1, 647.29, 220.64000, 0.3449)
        Thermo.CPAParams(1, 14.515, 1017.3, 0.6736)
        Thermo.AssocParams(1, 22, 69.2, 2003.2)
        Thermo.Setup_Thermo()
        T, LnK, ierr = Thermo.TBubble(1.01325,1,500) # P = 1.01325bar, moles = 1, Tini = 500K
        Thermo.Finishup_Thermo()
        self.assertGreater(T,0)
        self.assertAlmostEqual(T,373.15,0)
        self.assertNotEqual(ierr,1)

    def test_TBubble_values(self):
        self.assertRaises(ValueError, pt.Model().TBubble, -10,1,1)
        self.assertRaises(ValueError, pt.Model().TBubble, 0,1,1)
        self.assertRaises(ValueError, pt.Model().TBubble, 10,1,-10)
        self.assertRaises(ValueError, pt.Model().TBubble, 10,1,0)

    def test_TBubble_types(self):
        self.assertRaises(TypeError, pt.Model().TBubble, 3+5j, 1,1)
        self.assertRaises(TypeError, pt.Model().TBubble, "some string", 1,1)
        self.assertRaises(TypeError, pt.Model().TBubble, True, 1,1)


if __name__ == '__main__':
    unittest.main()