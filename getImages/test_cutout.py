from getpass import getpass, getuser
import hashlib
import md5
import os
from unittest import TestCase
import vos
import get_stamps

__author__ = 'jjk'

TEST_DIR = "asteroid_families/TEST/TEST_stamps/"
TEST_STAMP = "TEST_1667879p_21.123633_11.869728.fits"
TEST_NAME = 'TEST'
TEST_OBSERVATION = '1667879p'
TEST_MD5 = 'fc1d9caeae25ba087674305a91f464bf'


class TestCutout(TestCase):

    def test_cutout(self):
        ra = 21.1236333333
        dec = 11.8697277778
        radius = (20*0.185)/3600.0
        username = raw_input('Username: ')
        password = getpass('Password: ')
        output_name = os.path.join(TEST_DIR, TEST_STAMP)
        if os.path.exists(output_name):
            os.unlink(output_name)
        get_stamps.cutout(TEST_NAME, TEST_OBSERVATION, ra, dec, radius, username, password, TEST_NAME, test=True)
        out_md5 = hashlib.md5(open(output_name).read())
        self.assertEqual(out_md5.hexdigest(), TEST_MD5)

    def test_cutout_multi(self):
        test_observation = '1616690p'
        test_stamp = 'TEST_1616690p_216.966498_-13.832557.fits'
        test_md5 = '5eb99521f96054f3df6301b8ed3fe3df'
        ra = 216.966497722
        dec = -13.8325573937
        radius = (20*0.185)/3600.0
        username = raw_input('Username: ')
        password = getpass('Password: ')
        output_name = os.path.join(TEST_DIR, test_stamp)
        if os.path.exists(output_name):
            os.unlink(output_name)
        get_stamps.cutout(TEST_NAME, test_observation, ra, dec, radius, username, password, TEST_NAME, test=True)
        out_md5 = hashlib.md5(open(output_name).read())
        self.assertEqual(out_md5.hexdigest(), test_md5)
