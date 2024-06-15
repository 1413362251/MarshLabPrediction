import subprocess
import sys

def install_packages():
    packages = [
        'selenium',
        'webdriver-manager',
        'pandas',
        'requests',
        'lxml'
    ]
    #pip install packages
    subprocess.check_call([sys.executable, "-m", "pip", "install"] + packages)

if __name__ == "__main__":
    install_packages()