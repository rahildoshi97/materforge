import sys


def error(msg: str):
    msg = "Required Python packages could not be found:\n" + msg
    print(msg, file=sys.stderr)
    sys.exit(1)


print("Checking Python environment - ", end="")

try:
    import lbmpy
except ImportError:
    error("lbmpy is not installed in the current Python environment.")

try:
    import materforge
except ImportError:
    error("The materforge Python package is not installed in the current Python environment.")

try:
    import pystencils
except ImportError:
    error("pystencils is not installed in the current Python environment.")

try:
    import pystencilssfg
except ImportError:
    error("pystencils-sfg is not installed in the current Python environment.")

print("found required packages pystencils, pystencils-sfg, materforge", end="")

sys.exit(0)
