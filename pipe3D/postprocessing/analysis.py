
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import argparse # https://docs.python.org/3.7/howto/argparse.html
import os

from shutil import copy
import re

def main():
    parser = argparse.ArgumentParser(description='process simulation data.')
    parser.add_argument('dir',
            help="directory containing simulation(s)")
    args = parser.parse_args()

    if args.dir[-1] == '/':
        directory = args.dir
    else:
        directory = args.dir+'/'

    files, sub_directories = find_measurements_files(directory)
    if sub_directories == []:
        simulation_names = reduce_to_simulation_name([directory])
    else:
        simulation_names = reduce_to_simulation_name(sub_directories)

    header = generate_header()
    output = [header]
    for i,f in enumerate(files):
        print(simulation_names[i])
        info_inlet, info_width, best_seros = collect_info(f)
        if best_seros == None:
            continue
        verify_bestSeros(sub_directories[i], best_seros)
        stuff = (simulation_names[i]+ ' ' + str(best_seros)
                + ' ' + convert_to_str(info_inlet)
                + ' 0'
                + ' ' + convert_to_str(info_width))
        output.append(stuff)

    write_file(output, directory+'rhobarInletAvg_analysis')

def write_file(file_content, file_full_path):
    f = open(file_full_path, "w+")
    for line in file_content:
        f.write(line+'\n')
    f.close()
    return

def list_matching_subdirectories(parent_dir, pattern):
    result = glob.glob(parent_dir + pattern + '*/')
    return result

def find_measurements_files(parent_dir, number=0):
    file_pattern = 'measurementsSeros'+str(number)
    directory_pattern = '[2-3]D_'
    sub_directories = []
    sub_directories_with_files = []
    files = glob.glob(parent_dir + file_pattern)

    if len(files) == 0:
        sub_directories = list_matching_subdirectories(parent_dir, directory_pattern)
        for dir in sub_directories:
            found_files = glob.glob(dir + file_pattern)
            if found_files != []:
                files.append(found_files[0])
                sub_directories_with_files.append(dir)
    return files, sub_directories_with_files

def reduce_to_simulation_name(directories):
    simulation_names = [split(dir, '/')[-2] for dir in directories]
    return simulation_names


def read_lines_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    if not lines:
        print('no data in file %s'%filename)
    return lines

def split(line, delimiter=' '):
    elements = [item.strip() for item in line.split(delimiter)]
    return elements

def verify_bestSeros(directory, seros_at_min):
    files = glob.glob(directory+"bestSeros")
    bestSeros_count = 0
    if files != []:
        bestSeros_count = int(read_lines_from_file(files[0])[0].strip())
        if bestSeros_count != seros_at_min:
            print("   does not correspond with bestSeros file!")
    else:
        print("   could not verify with bestSeros file")
    return bestSeros_count

def get_valid_header(file_content):
    for line in file_content:
        if line.startswith('iT'):
            return line.strip()
        else:
            return ''

def get_data_array(file_content):
    if len(file_content) <= 1:
        return None

    data_array = []
    for line in file_content[1:]:
        elements = split(line)
        if elements:
            float_elements = [float(item) for item in elements if item]
        if float_elements:
            data_array.append(float_elements)

    return data_array

def find_element_in_header(header, element):
    elements = split(header)
    for column,item in enumerate(elements):
        if item == element:
            return column

def find_minimum_line(data, column):
    minimum_line = 0
    minimum = data[0][column]

    for number,line in enumerate(data):
        current_value = line[column]
        if current_value <= minimum:
            minimum = current_value
            minimum_line = number

    return minimum_line

def collect_info(file_name):
    file_content = read_lines_from_file(file_name)
    header = get_valid_header(file_content)
    data = get_data_array(file_content)
    if data is None:
        return None, None, None

    column_inlet  = find_element_in_header(header, "rhobarInletAvg")
    column_fluid  = find_element_in_header(header, "fluidCells")
    column_width  = find_element_in_header(header, "diameter")
    min_value_line = find_minimum_line(data, column_inlet)

    init_inlet  = data[0][column_inlet]
    init_fluid  = data[0][column_fluid]
    init_width  = data[0][column_width]
    end_inlet   = data[min_value_line][column_inlet]
    end_fluid   = data[min_value_line][column_fluid]
    end_width   = data[min_value_line][column_width]
    best_seros  = data[min_value_line][2]

    if init_inlet <0:
        return None

    info_inlet  = [init_inlet ,end_inlet ]
    info_width  = [init_width ,end_width ]
    return info_inlet, info_width, best_seros

def generate_header():
    header = 'simulation best_seros init_inlet end_inlet space init_width end_width'
    return header

def convert_to_str(info):
    line = ''
    for item in info:
        line+= str(item)+' '
    return line.strip()

def write_info(info, destination):
    str_info = convert_to_str(info)
    write_file([str_info], destination)
    return


if __name__ == "__main__":
    main()
