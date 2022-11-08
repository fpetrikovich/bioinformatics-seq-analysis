import argparse
import time

from orf_calculator import run_exercise_1
from file_helper import generate_output_path
from constants import ORFS_FILE_SUFFIX, FASTA_EXTENSION

def main():
    # Get the current time for the output files
    str_time = time.strftime("%Y%m%d-%H%M%S");

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', '--exercise', required=True)   # Ejercicio para correr
    ######################### Exercise 1 params #########################
    parser.add_argument('-gb', '--genbank', help='identifier of genbank input file',
                        type=str, required=False)

    args = parser.parse_args()

    input_file = ""
    output_file = ""
    item = 0

    # Param parsing and setup
    try:
        item = int(args.exercise)
        if item == 1 and args.genbank != None: 
            input_file = args.genbank
            output_file = generate_output_path(str(args.genbank))

        else:
            print("[ERROR] Invalid combination of params. Check the manual.")
            exit(1)

    except Exception as e:
        print("[ERROR] " + str(e))
        exit(1)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    start_time = time.time()
    try:
        if item == 1:
            nucleotide_file = output_file + '_nucleotides' + FASTA_EXTENSION
            proteins_file = output_file + ORFS_FILE_SUFFIX + FASTA_EXTENSION
            run_exercise_1(input_file, nucleotide_file, proteins_file)

    except Exception as e:
        print("[ERROR] Unknown error when running the exercise.")
        print("[ERROR][MESSAGE] " + str(e))
        exit(1)
    
    execution_time = int(time.time() - start_time)
    print("[DONE] Execution took %i mins %i seconds" % (int(execution_time / 60), execution_time % 60))


if __name__ == '__main__':
    main()
