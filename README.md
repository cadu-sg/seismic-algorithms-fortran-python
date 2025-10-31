This project allows us to call the `semblance` Fortran subroutine from Python.


#### How to run the test example

Install all packages specified in Pipfile.lock (exact versions)
```
pipenv sync --dev
```

Activate the virtual environment
```
pipenv shell
```

You need to generate an "entension module" file named `seismic_algorithms.cpython-312-x86_64-linux-gnu.so`. When this file is present, we can import it from Pyhton as if it were a regular module (`import seismic_algorithms`) and then call the Fortran functions contained inside it.

Run the following command to generate the extension module file
```
f2py -c -m seismic_algorithms seismic_algorithms.f90
```

The extension module should be present. You can now run the `example_bokeh.ipynb` Jupyter Notebook.
