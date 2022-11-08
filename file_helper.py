import os.path
from subprocess import call

def generate_output_path(filename):
    return "output/" + filename

def save_file(file_name, content): 
	f = open(file_name, "w")
	f.write(content)
	f.close()

def create_bash_file(file_name, command):
    f = open(file_name, "w")
    f.write("#!/bin/sh\n")
    f.write(command)
    f.close()

def run_bash_file(file_name):
    call(file_name, shell=True)

def run_bash_file_with_arguments(file_name, args):
    call(['bash', file_name, *args])

def file_exists(file):
    return os.path.exists(file)

def delete_file(path):
    os.remove(path)
