#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import unittest
from analysis import list_matching_subdirectories
from analysis import find_measurements_files
from analysis import reduce_to_simulation_name
from analysis import write_file
from analysis import read_lines_from_file
from analysis import verify_bestSeros
from analysis import get_valid_header
from analysis import find_element_in_header
from analysis import get_data_array
from analysis import find_minimum_line
from analysis import collect_info
from analysis import convert_to_str
from analysis import write_info

cwd = os.getcwd()+'/'
filename = 'assets/example_measurements/measurementsSerosSimple'
filename_long = 'assets/example_measurements/measurementsSeros0'
filename_long_error = 'assets/example_measurements/error_measurementsSeros0'

file_content_long = [
'iT time seros rhobarInletAvg rhobarInletAvgSmooth rhobarInletAvgSmoothGrad rhobarGlobalAvg rhobarGlobalAvgSmooth rhobarGlobalAvgSmoothGrad kineticEnergyAvg kineticEnergyAvgSmooth kineticEnergyAvgSmoothGrad diameter fluidCells wallCells volumeToSurfaceRatio',
'0 0 0 0.000779342 0.000779342 0 0.00039115 0.00039115 0 7.37704e-07 7.37704e-07 0 1 11676 1108 10.5379',
'1500 15 0 0.000623192 0.000623192 0 0.000314544 0.000314544 0 1.10697e-06 1.10697e-06 0 1 11676 1108 10.5379',
'3000 30 1 0.000606805 0.000621881 -1.31093e-06 0.000304119 0.00031371 -8.34035e-07 1.13213e-06 1.10898e-06 2.01274e-09 1.00104 11224 1064 10.5489',
'4500 45 2 0.000617643 0.000621542 -3.39008e-07 0.00031024 0.000313433 -2.77651e-07 1.14868e-06 1.11216e-06 3.17605e-09 0.996222 11044 1052 10.4981',
'6000 60 3 0.000627091 0.000621986 4.43946e-07 0.000315152 0.00031357 1.37557e-07 1.17104e-06 1.11687e-06 4.71029e-09 1.00122 10804 1024 10.5508',
'7500 75 4 0.000641786 0.00062357 1.58402e-06 0.000321621 0.000314214 6.44074e-07 1.19306e-06 1.12296e-06 6.09573e-09 0.986394 10644 1024 10.3945',
'9000 90 5 0.000644422 0.000625238 1.66813e-06 0.00032264 0.000314888 6.74093e-07 1.20807e-06 1.12977e-06 6.80867e-09 1.00193 10516 996 10.5582',
'10500 105 6 0.000645274 0.000626841 1.60288e-06 0.000322606 0.000315506 6.17423e-07 1.21067e-06 1.13624e-06 6.47153e-09 0.99888 10484 996 10.5261',
'12000 120 7 0.00064861 0.000628582 1.7415e-06 0.000324583 0.000316232 7.2621e-07 1.21661e-06 1.14267e-06 6.42954e-09 0.997933 10432 992 10.5161',
'13500 135 8 0.00064008 0.000629502 9.19771e-07 0.000320509 0.000316574 3.42169e-07 1.21232e-06 1.14824e-06 5.57148e-09 1.00197 10432 988 10.5587',
'15000 150 9 0.000644888 0.000630733 1.23085e-06 0.000323263 0.000317109 5.35086e-07 1.21569e-06 1.15364e-06 5.39594e-09 1.0018 10388 984 10.5569',
'16500 165 10 0.000639051 0.000631399 6.65447e-07 0.000319987 0.000317339 2.30182e-07 1.21317e-06 1.1584e-06 4.76214e-09 0.997748 10388 988 10.5142',
'18000 180 11 0.00064195 0.000632243 8.44117e-07 0.000321139 0.000317643 3.03971e-07 1.21636e-06 1.16304e-06 4.6368e-09 0.993906 10348 988 10.4737',
'19500 195 12 0.000642659 0.000633076 8.33318e-07 0.000321763 0.000317973 3.29564e-07 1.22076e-06 1.16766e-06 4.61771e-09 1.00185 10304 976 10.5574',
'21000 210 13 0.000639331 0.000633576 5.00376e-07 0.000320407 0.000318168 1.94687e-07 1.21788e-06 1.17167e-06 4.01767e-09 0.997759 10304 980 10.5143',
]

file_header_long = 'iT time seros rhobarInletAvg rhobarInletAvgSmooth rhobarInletAvgSmoothGrad rhobarGlobalAvg rhobarGlobalAvgSmooth rhobarGlobalAvgSmoothGrad kineticEnergyAvg kineticEnergyAvgSmooth kineticEnergyAvgSmoothGrad diameter fluidCells wallCells volumeToSurfaceRatio'

file_content = [
    'iT time\n',
    '0 0\n',
    '1000 0.1\n',
    ]

file_header = 'iT time'

file_data = [
    [0, 0],
    [1000, 0.1],
    ]

file_data_min = [
    [1000, 0.1, 0, 1020],
    [2000, 0.2, 1, 1000],
    [3000, 0.3, 2, 1030],
    ]


class Test_file_writer(unittest.TestCase):

    def setUp(self):
        self.write_dir = cwd+'assets/tmp_test/'
        self.file_full_path = self.write_dir+'test.file'

    def tearDown(self):
        os.remove(self.file_full_path)

    def test_write_file(self):
        write_file(file_content, self.file_full_path)
        self.assertTrue(os.path.isfile(self.file_full_path))

        # get content of written file:
        content = []
        f = open(self.file_full_path, "r")
        content = [line for line in f]
        f.close()
        self.assertTrue(content, file_content)


class TestListDirectories(unittest.TestCase):

    def test_list_matching_subdirectories(self):
        parent_dir = cwd + 'assets/example_directories/'
        pattern = '[2-3]D_'
        directories = [
                      parent_dir + '2D_blub/',
                      parent_dir + '3D_foo/',
                      ]

        result = list_matching_subdirectories(parent_dir, pattern)
        self.assertEqual(result, directories)


class TestFilesFinder(unittest.TestCase):

    def test_find_measurements_files(self):
        parent_dir = cwd + 'assets/example_directories/'
        file_pattern = 'measurementsSeros0'
        result_files, result_dirs = find_measurements_files(parent_dir)

        files = [
                parent_dir + '2D_blub/' + file_pattern,
                parent_dir + '3D_foo/' + file_pattern,
                ]

        directories = [
                          parent_dir + '2D_blub/',
                          parent_dir + '3D_foo/',
                          ]

        self.assertEqual(result_files, files)
        self.assertEqual(result_dirs, directories)


class TestSimulationNames(unittest.TestCase):

    def test_reduce_to_simulation_name(self):

        parent_dir = cwd + 'assets/example_directories/'
        directories = [
                          parent_dir + '2D_blub/',
                          parent_dir + '3D_foo/',
                          ]
        simulation_names = [
                '2D_blub',
                '3D_foo',
                ]
        result = reduce_to_simulation_name(directories)

        self.assertEqual(result, simulation_names)

class TestFileReader(unittest.TestCase):

    def test_read_file_linewise(self):
        content = read_lines_from_file(filename)

        self.assertEqual(content, file_content)


class TestFileHeader(unittest.TestCase):

    def test_get_valid_header(self):
        header = get_valid_header(file_content)

        self.assertEqual(header, file_header)


class TestFileData(unittest.TestCase):

    def test_get_data_array(self):
        data = get_data_array(file_content)

        self.assertEqual(data, file_data)


class TestFindElement(unittest.TestCase):

    def test_find_element_in_header(self):
        position_iT = find_element_in_header(file_header, 'iT')
        position_time = find_element_in_header(file_header, 'time')
        position_rhobarGlobalAvg = find_element_in_header(file_header_long, 'rhobarGlobalAvg')

        self.assertEqual(position_iT, 0)
        self.assertEqual(position_time, 1)
        self.assertEqual(position_rhobarGlobalAvg, 6)


class TestFindMinimumline(unittest.TestCase):

    def test_find_minimum_line(self):
        minimum_line = find_minimum_line(file_data_min, 3)

        self.assertEqual(minimum_line, 1)


class TestSerosVerification(unittest.TestCase):

    def test_verify_bestSeros(self):
        simulation_dir = cwd + 'assets/example_measurements/'
        bestSeros_count = verify_bestSeros(simulation_dir, 42)

        self.assertEqual(bestSeros_count, 42)


class TestCollectInfo(unittest.TestCase):

    def test_collect_info_on_rhobar(self):
        info_inlet, info_width, best_seros = collect_info(filename_long)

        self.assertEqual(best_seros, 1)
        self.assertEqual(info_inlet, [0.000779342,0.000606805])


class TestConvertInfo(unittest.TestCase):

    def test_convert_to_str(self):
        info = (1,2,3)
        str_info = convert_to_str(info)

        self.assertEqual(str_info, '1 2 3')


class TestWriteInfo(unittest.TestCase):

    write_dir = cwd+'assets/tmp_test/'

    def test_write_info_rhobar(self):
        info_rhobar = collect_info(filename_long)
        str_info_rhobar = convert_to_str(info_rhobar)
        write_file(str_info_rhobar, self.write_dir+"info_rhobar")

    def tearDown(self):
        os.remove(self.write_dir+"info_rhobar")


if __name__ == "__main__":
    unittest.main(verbosity=2)



