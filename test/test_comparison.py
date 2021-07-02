import unittest
from pythermo import pythermo as pt

class Testing(unittest.TestCase):
    """
        Here the ComparisonFuncs initialization procedure is tested
    """
    def test_ComparisonFuncs_values(self):
        self.assertRaises(ValueError, pt.ComparisonFuncs, pt.Model, "bussemand")

    def test_ComparisonFuncs_types(self):
        self.assertRaises(TypeError, pt.ComparisonFuncs, pt.Model, 2)
        self.assertRaises(TypeError, pt.ComparisonFuncs, pt.Model, 2+3j)
        self.assertRaises(TypeError, pt.ComparisonFuncs, pt.Model, True)
        self.assertRaises(TypeError, pt.ComparisonFuncs, 1, "ARD")


    """
        Here the ComparisonFuncs initialization procedure is tested
    """
    def test_PBubble_comparison_values(self):
        Thermo = pt.Model()
        Comp = pt.ComparisonFuncs(Thermo,"ARD")
        self.assertRaises(ValueError, Comp.PBubble_comparison, -1, 1, 1, 1)
        self.assertRaises(ValueError, Comp.PBubble_comparison, 1, -1, 1, 1)
        self.assertRaises(ValueError, Comp.PBubble_comparison, 1, 1, -1, 1)
        self.assertRaises(ValueError, Comp.PBubble_comparison, 1, 1, 1, -1)


    def test_PBubble_comparison_values(self):
        Thermo = pt.Model()
        Comp = pt.ComparisonFuncs(Thermo,"ARD")
        self.assertRaises(TypeError, Comp.PBubble_comparison, "1", 1, 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, "1", 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, "1", 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, 1, "1")

        self.assertRaises(TypeError, Comp.PBubble_comparison, True, 1, 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, True, 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, True, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, 1, True)

        self.assertRaises(TypeError, Comp.PBubble_comparison, 1+1j, 1, 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1+1j, 1, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, 1+1j, 1)
        self.assertRaises(TypeError, Comp.PBubble_comparison, 1, 1, 1, 1+1j)


if __name__ == '__main__':
    unittest.main()