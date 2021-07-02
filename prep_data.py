"""Entrypoint script to process and stitch data."""
import glob
import os

SOLUTIONS_DIR = 'data/solutions'


def main():
    """Entry point."""
    pass


if __name__ == '__main__':
    for solution_dir in glob.glob(os.path.join(SOLUTIONS_DIR, '*')):
        print(solution_dir)
