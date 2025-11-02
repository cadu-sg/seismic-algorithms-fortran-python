This project demonstrates how to call, from within Python, Fortran subroutines related to velocity analysis.


#### How to run the test example

Install all packages specified in Pipfile.lock (exact versions)
```
pipenv sync --dev
```

Activate the virtual environment
```
pipenv shell
```

You need to generate an "entension module" file named `seismic_algorithms.cpython-312-x86_64-linux-gnu.so`. When this file is present, we can import it from Python as if it were a regular module (`import seismic_algorithms`) and then call the Fortran functions contained inside it.

Run the following command to generate the extension module file
```
f2py -c -m seismic_algorithms seismic_algorithms.f90
```

If the command was run successfully, the extension module should have been created. Now the test programs are ready to be used.

#### Test programs

**`example_bokeh.ipynb`**

Shows how to create a semblance plot and a very simple point drawer.

**`velocity_analysis.py`**

Shows how to implement a simple velocity analysis tool, which includes a semblance plot, a more polished picking tool, and NMO correction. It is implemented as a Bokeh application, Python code that should be run by a Bokeh server with the following command:
```
bokeh serve --show velocity_analysis.py
```
