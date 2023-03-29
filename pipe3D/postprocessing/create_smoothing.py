
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse # https://docs.python.org/3.7/howto/argparse.html
import os

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



def main():
    parser = argparse.ArgumentParser(description='create smoothed data from simulation data.')
    parser.add_argument('dir',
            help="directory containing simulation(s)")
    parser.add_argument('factor',
            help='smoothing factor (smaller number means more aggressive smoothing')
    parser.add_argument('--key',
            help='what measurement should be smoothed')

    args = parser.parse_args()

    if args.key == None:
        key = 'rhobarGlobalAvgSmoothGrad'
    else:
        key = args.key

    if args.dir[-1] == '/':
        directory = args.dir
    else:
        directory = args.dir+'/'

    files, sub_directories = find_measurements_files(directory)
    if sub_directories == []:
        simulation_names = reduce_to_simulation_name([directory])
    else:
        simulation_names = reduce_to_simulation_name(sub_directories)

    smoothed_data = smooth_data_from_file(files[0], key, float(args.factor))
    output=[str(x) for x in smoothed_data]
    write_file(output, directory+key+'Smooth')

'''
    header = generate_header()
    output = [header]
    for i,f in enumerate(files):
        info = collect_info(f, key)
        if info == None:
            continue
        print(simulation_names[i])
        print_info(info, key)
        if key == 'rhobarGlobalAvg':
            bestSeros_count = verify_bestSeros(sub_directories[i], info[3])
        output.append(simulation_names[i]+ ' ' + convert_to_str(info))

'''

def extract_column(data, element):
    column_data = []

    for line in data:
        current_value = line[element]
        column_data.append(current_value)

    return column_data

def smooth_column(column, factor):
    smoothed_column = []
    smoothed_column.append(column[0])
    for i,element in enumerate(column[1:]):
        smoothed_element = element*factor + smoothed_column[i]*(1-factor)
        smoothed_column.append(smoothed_element)

    return smoothed_column

def smooth_data_from_file(file_name, key, factor):
    file_content = read_lines_from_file(file_name)
    header = get_valid_header(file_content)
    data = get_data_array(file_content)
    position_in_header = find_element_in_header(header, key)
    column = extract_column(data, position_in_header)
    smoothed = smooth_column(column, factor)
    return smoothed



if __name__ == "__main__":
    main()
