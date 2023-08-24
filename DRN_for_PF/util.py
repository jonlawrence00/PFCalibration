import subprocess, os, sys

def check_git():
    #check if git working directory is clean (ie no changes to be committed)
    status = subprocess.check_output(['git', 'status', '--porcelain'], cwd=os.path.dirname(__file__)).decode('ascii').strip()
    if status == '':
        #working directory is clean
        pass
    else:
        print("ERROR: you have uncomitted changes in the git reposity. Commit your changes and try again")
        sys.exit(-1)

    git_hash = subprocess.check_output(['git', 'describe', '--always'], cwd=os.path.dirname(__file__)).decode('ascii').strip()
    print("Current git commit:",git_hash)
    print()

def print_params(args):
    print("Parameters:")
    for arg in vars(args):
        print(arg+":", getattr(args, arg))
    print()
