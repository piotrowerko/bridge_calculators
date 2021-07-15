import unittest
from char_geom import CharGeom


class TestCharGeom(unittest.TestCase):
    """Testy dla klasy CharGeom"""

    def test_area(self):
        """Sprawdzenie pola wierzchni"""
        moje_cw_test = CharGeom()
        moje_cw_test.znajdz_srodek_c((3, 1.5, 1.5, 2.5, 4, 2.5))
        print(f'pole: {sum(moje_cw_test.A)}; środek ciężkości: {round((moje_cw_test.ex), 2)}')
        self.assertEqual(sum(moje_cw_test.A), 18.25)


    def test_ex(self):
        moje_cw_test = CharGeom()
        moje_cw_test.znajdz_srodek_c((3, 1.5, 1.5, 2.5, 4, 2.5))
        self.assertEqual(round((moje_cw_test.ex), 2), 3.63)


if __name__ == '__main__':
    unittest.main()