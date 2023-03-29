
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
from analysis import get_init_value
from analysis import find_minimum_line
from analysis import get_seros_at_min
from analysis import get_seros_for_min_value
from analysis import collect_info
from analysis import convert_to_str
from analysis import write_info

from create_smoothing import extract_column
from create_smoothing import smooth_column
from create_smoothing import smooth_data_from_file

cwd = os.getcwd()+'/'
filename_gradSmooth = 'assets/example_measurements_with_gradSmooth/measurementsSeros0'

class TestExtractColumn(unittest.TestCase):

    def test_smooth_data(self):
        multiple_columns = [
            [1000, 0.1, 0, 1020],
            [2000, 0.2, 1, 1000],
            [3000, 0.3, 2, 1030],
            ]
        first_column = [1000,2000,3000]
        second_column = [.1,.2,.3]

        extracted_first_column = extract_column(multiple_columns, 0)
        extracted_second_column = extract_column(multiple_columns, 1)
        self.assertEqual(extracted_first_column, first_column )
        self.assertEqual(extracted_second_column, second_column )

class TestSmoothColumn(unittest.TestCase):

    def test_no_smooth_simple_column(self):
        column = [0,1,1]
        factor = 1
        smoothed_column = [0,1,1]

        smoothed = smooth_column(column, factor)
        self.assertEqual(smoothed_column, smoothed)

    def test_smooth_simple_column(self):
        column = [0,1,1]
        factor = .1
        smoothed_column = [0,.1,.1+.1*.9]

        smoothed = smooth_column(column, factor)
        self.assertEqual(smoothed_column, smoothed)

    def test_smooth_real_column(self):
        rhobarAvg = [0.00323776,0.00312972]
        rhobarAvgSmooth = [0.00323776,0.00323344]
        factor = 0.04

        smoothed = smooth_column(rhobarAvg, factor)
        rounded_smoothed = [round(x,8) for x in smoothed]
        self.assertEqual(rhobarAvgSmooth, rounded_smoothed)

class TestSmoothDataFromFile(unittest.TestCase):

    def test_smooth_measurements_file(self):

        file_content = read_lines_from_file(filename_gradSmooth)
        header = get_valid_header(file_content)
        data = get_data_array(file_content)
        position_in_header = find_element_in_header(header, 'rhobarGlobalAvgSmooth')
        rhobarGlobalAvgSmooth = extract_column(data, position_in_header)

        smoothed = smooth_data_from_file(filename_gradSmooth, 'rhobarGlobalAvg', 0.04)
        rounded_smoothed = [round(x,8) for x in smoothed]
        self.assertEqual(rhobarGlobalAvgSmooth[0], rounded_smoothed[0])
        self.assertEqual(rhobarGlobalAvgSmooth[1], rounded_smoothed[1])
        self.assertEqual(rhobarGlobalAvgSmooth[2], rounded_smoothed[2])




if __name__ == "__main__":
    unittest.main(verbosity=2)

