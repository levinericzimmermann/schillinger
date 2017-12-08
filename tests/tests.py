from schillinger import schillinger
import unittest


class TestSchillinger(unittest.TestCase):
    def test_set_range(self):
        self.assertEqual(schillinger.setrange(1, 9),
                         set((0, 1, 2, 3, 4, 5, 6, 7, 8)))
        self.assertEqual(schillinger.setrange(2, 9),
                         set((0, 2, 4, 6, 8)))
        self.assertEqual(schillinger.setrange(3, 13),
                         set((0, 3, 6, 9, 12)))

    def test_synchronize(self):
        self.assertEqual(schillinger.synchronize(3, 2), [2, 1, 1, 2])
        self.assertEqual(schillinger.synchronize(5, 3, 2),
                         [2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1,
                          1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 2])

    def test_synchronize_complementary(self):
        self.assertEqual(schillinger.synchronize_complementary(2, 3, 5),
                         [6, 4, 2, 3, 3, 2, 4, 6])
