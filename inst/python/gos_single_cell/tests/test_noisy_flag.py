"""Synthetic tests for the data-driven noisy-patient flag: python3 -m pytest test_noisy_flag.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gos_sc_noisy as n  # noqa: E402


class Noisy(unittest.TestCase):
    def test_default_empty(self):
        self.assertEqual(n.noisy_patients({}), set())
        self.assertFalse(n.is_noisy("BWH69", {}))

    def test_env_set(self):
        e = {"NOISY_PATIENTS": "BWH69, MGH999,"}
        self.assertEqual(n.noisy_patients(e), {"BWH69", "MGH999"})
        self.assertTrue(n.is_noisy("MGH999", e))
        self.assertFalse(n.is_noisy("BWH70", e))

    def test_caveat_cleaned_set_when_not_noisy(self):
        c = n.noisy_caveat(env={})
        self.assertIn("cleaned site set (5,675 sites)", c)
        self.assertNotIn("noisy", c)

    def test_caveat_noisy(self):
        c = n.noisy_caveat(env={"NOISY_PATIENTS": "BWH69"})
        self.assertIn("BWH69 SNV calls are noisy", c)

    def test_caveat_none(self):
        self.assertIsNone(n.noisy_caveat(env={"CLEANED_SNV_PATIENTS": ""}))
        self.assertIsNone(n.noisy_caveat(["BWH70"], env={}))


if __name__ == "__main__":
    unittest.main()
